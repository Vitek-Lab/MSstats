## ── Memory-monitoring helpers ─────────────────────────────────────────────────

# RSS of the current process in MB.
# On Linux reads /proc/self/status (VmRSS); elsewhere falls back to gc() counts.
.current_rss_mb <- function() {
    if (file.exists("/proc/self/status")) {
        ln <- readLines("/proc/self/status", warn = FALSE)
        m  <- grep("^VmRSS:", ln, value = TRUE)
        if (length(m))
            return(as.numeric(gsub("[^0-9]", "", m[1L])) / 1024)
    }
    g <- gc(reset = FALSE)
    (g["Ncells", "used"] * 8L + g["Vcells", "used"] * 8L) / 1024^2
}

# Cross-platform peak-RSS reader. Reflects the true lifetime peak of the
# calling process, regardless of when you call it — no polling required.
# All three branches report OS-level peak resident/working-set memory,
# so figures are comparable across platforms.
.peak_rss_mb <- function() {
  if (.Platform$OS.type == "windows") {
    if (!exists(".peakRSS_windows_impl", mode = "function")) {
      Rcpp::cppFunction(
        depends  = "Rcpp",
        includes = c(
          "#include <windows.h>",
          "#include <psapi.h>"),
        code = "
        double peakRSS_windows_impl() {
          PROCESS_MEMORY_COUNTERS pmc;
          if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
            return (double) pmc.PeakWorkingSetSize / (1024.0*1024.0); // bytes -> MB
          }
          return NA_REAL;
        }")
      assign(".peakRSS_windows_impl", peakRSS_windows_impl, envir = .GlobalEnv)
    }
    return(.peakRSS_windows_impl())
  }
  
  if (file.exists("/proc/self/status")) {
    # Linux: VmHWM = kernel-maintained peak resident set size
    ln <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
    if (length(ln)) return(as.numeric(sub("\\D+(\\d+).*", "\\1", ln)) / 1024)
  }
  
  # macOS (also works as a Linux fallback): POSIX getrusage() ru_maxrss
  # is likewise a lifetime high-water mark, just different units per OS.
  if (!exists(".rusage_maxrss_mb_impl", mode = "function")) {
    Rcpp::cppFunction(
      depends  = "Rcpp",
      includes = "#include <sys/resource.h>",
      code = "
      double rusage_maxrss_mb_impl() {
        struct rusage ru; getrusage(RUSAGE_SELF, &ru);
      #ifdef __APPLE__
        return (double) ru.ru_maxrss / (1024.0*1024.0); // bytes -> MB
      #else
        return (double) ru.ru_maxrss / 1024.0;           // KB -> MB
      #endif
      }")
    assign(".rusage_maxrss_mb_impl", rusage_maxrss_mb_impl, envir = .GlobalEnv)
  }
  .rusage_maxrss_mb_impl()
}

# Print a formatted memory report to stderr via message().
# checkpoints:    named numeric vector of RSS snapshots (MB).
# worker_peak_mb: numeric vector of per-worker peak RSS values (may be NA).
# elapsed:        total wall-clock seconds (NULL to omit).
.print_memory_report <- function(function_name, checkpoints, worker_peak_mb = NULL,
                                 elapsed = NULL) {
    rule_width <- 65L
    rule <- strrep("─", rule_width)
    format_mb    <- function(x) if (is.na(x)) "    n/a" else sprintf("%7.1f", x)
    format_delta <- function(now, previous_value) {
        if (is.na(now) || is.na(previous_value)) return("")
        sprintf("  (%+.1f MB)", now - previous_value)
    }
    lines <- c(
        rule,
        sprintf(" MSstats Memory Report — %s", function_name),
        rule,
        sprintf("  %-36s  %7s  %s", "Checkpoint", "RSS MB", "Delta")
    )
    previous_value <- NA_real_
    for (checkpoint_name in names(checkpoints)) {
        checkpoint_value <- checkpoints[[checkpoint_name]]
        lines <- c(lines,
                   sprintf("  %-36s  %s%s", checkpoint_name,
                           format_mb(checkpoint_value),
                           format_delta(checkpoint_value, previous_value)))
        previous_value <- checkpoint_value
    }
    if (!is.null(worker_peak_mb)) {
        observed_peaks <- worker_peak_mb[!is.na(worker_peak_mb)]
        if (length(observed_peaks)) {
            lines <- c(lines, "",
                       sprintf("  Worker RSS  min / mean / max :  %.1f / %.1f / %.1f MB",
                               min(observed_peaks), mean(observed_peaks), max(observed_peaks)))
        }
    }
    if (!is.null(elapsed))
        lines <- c(lines, sprintf("  Total elapsed: %.1f s", as.numeric(elapsed)))
    lines <- c(lines, rule)
    message(paste(lines, collapse = "\n"))
}

## ── Protein-slot pack/unpack helpers ───────────────────────────────────────────
##
## Packed double-vector layout for one protein slot (all column-major matrices).
## FL = n_feature_labels (unique FEATURE × LABEL pairs), R = n_runs:
##
##  pos 1                          : slot_index  (protein index, cast to double)
##  pos 2          .. FL*R+1       : newABUNDANCE  (FL × R)
##  pos FL*R+2     .. 2*FL*R+1    : ABUNDANCE      (FL × R; NA for unlabeled/TMP)
##  pos 2*FL*R+2   .. 3*FL*R+1   : censored        (FL × R; 0.0/1.0)
##  pos 3*FL*R+2   .. 4*FL*R+1   : cen             (FL × R; 1-censored event flag)
##  pos 4*FL*R+2   .. 5*FL*R+1   : ANOMALYSCORES   (FL × R; NA if unused)
##  pos 5*FL*R+2   .. 5*FL*R+FL+1 : n_obs          (FL; per feature-label)
##  pos 5*FL*R+FL+2.. 5*FL*R+FL+R+1: n_obs_run     (R;  per run)
##  pos 5*FL*R+FL+R+2..5*FL*R+FL+2R+1: prop_features (R; per run)
##
## Total length: 1 + 5*FL*R + FL + 2*R
## ─────────────────────────────────────────────────────────────────────────────

#' Build the packed double vector and lightweight metadata for one protein slot
#'
#' @param protein_dt data.table rows belonging to one protein (or protein × label) slot
#' @param slot_index integer position of this slot in the global protein list
#' @param all_runs character vector of all run names in global order
#' @return list with elements \code{packed} (double vector) and \code{meta} (list)
#' @keywords internal
.pack_protein_slot <- function(protein_dt, slot_index, all_runs) {

    n_runs <- length(all_runs)

    # ── Unique (FEATURE, LABEL) → "effective features", stable ordering ───────
    has_peptide <- "PEPTIDE" %in% colnames(protein_dt)
    feature_label_dt <- unique(protein_dt[, .(
        FEATURE = as.character(FEATURE),
        LABEL   = as.character(LABEL),
        PEPTIDE = if (has_peptide) as.character(PEPTIDE) else as.character(FEATURE)
    )])
    data.table::setorder(feature_label_dt, FEATURE, LABEL)
    n_feature_labels <- nrow(feature_label_dt)

    # ── Index maps: feature-label → row index, run name → column index ─────────
    feature_row_idx <- match(
        paste(as.character(protein_dt$FEATURE), as.character(protein_dt$LABEL), sep = "\t"),
        paste(feature_label_dt$FEATURE, feature_label_dt$LABEL, sep = "\t"))
    run_col_idx_lookup <- seq_len(n_runs); names(run_col_idx_lookup) <- all_runs
    run_col_idx  <- run_col_idx_lookup[as.character(protein_dt$RUN)]
    row_is_valid <- !is.na(feature_row_idx) & !is.na(run_col_idx)

    # ── Helper: allocate n_feature_labels × n_runs matrix, fill from a protein_dt column ──
    scatter_into_matrix <- function(col_vec, fill = NA_real_) {
        m <- matrix(fill, nrow = n_feature_labels, ncol = n_runs)
        m[cbind(feature_row_idx[row_is_valid], run_col_idx[row_is_valid])] <-
            as.double(col_vec[row_is_valid])
        m
    }

    # ── Build the five numeric matrices ───────────────────────────────────────
    new_abundance_mat <- scatter_into_matrix(protein_dt$newABUNDANCE)

    has_ABUNDANCE <- "ABUNDANCE" %in% colnames(protein_dt)
    abundance_mat <- if (has_ABUNDANCE) scatter_into_matrix(protein_dt$ABUNDANCE) else
        matrix(NA_real_, n_feature_labels, n_runs)

    has_censored <- "censored" %in% colnames(protein_dt)
    censored_mat <- if (has_censored) scatter_into_matrix(as.double(protein_dt$censored)) else
        matrix(0.0, n_feature_labels, n_runs)   # treat as non-censored when column absent

    has_cen <- "cen" %in% colnames(protein_dt)
    event_mat <- if (has_cen) scatter_into_matrix(protein_dt$cen) else
        matrix(NA_real_, n_feature_labels, n_runs)

    has_anom <- "ANOMALYSCORES" %in% colnames(protein_dt) &&
        !all(is.na(protein_dt$ANOMALYSCORES))
    anomaly_scores_mat <- if (has_anom) scatter_into_matrix(protein_dt$ANOMALYSCORES) else
        matrix(NA_real_, n_feature_labels, n_runs)

    # ── Per-feature scalar: n_obs (constant within FEATURE × LABEL) ───────────
    n_obs_by_feature_label <- protein_dt[, .(n_obs = as.double(n_obs[1L])),
                      by = .(FEATURE = as.character(FEATURE),
                             LABEL   = as.character(LABEL))]
    feature_label_with_nobs <- n_obs_by_feature_label[feature_label_dt, on = c("FEATURE", "LABEL")]
    data.table::setorder(feature_label_with_nobs, FEATURE, LABEL)
    n_obs_vec   <- feature_label_with_nobs$n_obs

    # ── Per-run scalars: n_obs_run, prop_features (constant within RUN) ───────
    run_level_scalars <- protein_dt[, .(n_obs_run    = as.double(n_obs_run[1L]),
                       prop_features = as.double(prop_features[1L])),
                   by = .(RUN = as.character(RUN))]
    run_scalars_all_runs <- run_level_scalars[data.table::data.table(RUN = all_runs), on = "RUN"]
    n_obs_run_vec     <- run_scalars_all_runs$n_obs_run
    prop_features_vec <- run_scalars_all_runs$prop_features

    # ── Pack ─────────────────────────────────────────────────────────────────
    packed <- c(
        as.double(slot_index),
        as.vector(new_abundance_mat),
        as.vector(abundance_mat),
        as.vector(censored_mat),
        as.vector(event_mat),
        as.vector(anomaly_scores_mat),
        n_obs_vec,
        n_obs_run_vec,
        prop_features_vec
    )

    # ── Metadata (string labels, flags; kept in main-process RAM) ────────────
    meta <- list(
        PROTEIN           = as.character(protein_dt$PROTEIN[1L]),
        feature_label_dt  = as.data.frame(feature_label_dt),   # n_feature_labels × 3: FEATURE, LABEL, PEPTIDE
        runs              = all_runs,
        n_feature_labels  = n_feature_labels,
        n_runs            = n_runs,
        is_labeled_ref    = "is_labeled_ref" %in% colnames(protein_dt) &&
            isTRUE(any(protein_dt$is_labeled_ref, na.rm = TRUE)),
        has_ABUNDANCE     = has_ABUNDANCE,
        has_censored      = has_censored,
        has_cen           = has_cen,
        has_anom          = has_anom,
        add_ref_covariate = "ref_covariate" %in% colnames(protein_dt)
    )

    list(packed = packed, meta = meta)
}


#' Reconstruct a per-protein data.table from a packed double vector
#'
#' @param packed double vector produced by \code{.pack_protein_slot}
#' @param meta metadata list from \code{.pack_protein_slot}
#' @return data.table compatible with \code{MSstatsSummarizeSingleTMP} /
#'   \code{MSstatsSummarizeSingleLinear}
#' @keywords internal
.unpack_protein_slot <- function(packed, meta) {

    n_feature_labels <- meta$n_feature_labels
    n_runs           <- meta$n_runs
    matrix_len       <- n_feature_labels * n_runs

    # ── Unpack sections (1-indexed; position 1 is slot_index header) ─────────
    cursor <- 2L
    read_next_matrix <- function() {
        m <- matrix(packed[cursor:(cursor + matrix_len - 1L)], nrow = n_feature_labels, ncol = n_runs)
        cursor <<- cursor + matrix_len
        m
    }
    new_abundance_mat  <- read_next_matrix()
    abundance_mat      <- read_next_matrix()
    censored_mat       <- read_next_matrix()
    event_mat          <- read_next_matrix()
    anomaly_scores_mat <- read_next_matrix()

    n_obs_vec         <- packed[cursor:(cursor + n_feature_labels - 1L)]; cursor <- cursor + n_feature_labels
    n_obs_run_vec     <- packed[cursor:(cursor + n_runs - 1L)]; cursor <- cursor + n_runs
    prop_features_vec <- packed[cursor:(cursor + n_runs - 1L)]

    # ── Column-major melt: rep/each match the matrix column-major ordering ────
    feature_label_dt <- meta$feature_label_dt   # data.frame: FEATURE, LABEL, PEPTIDE (n_feature_labels rows)
    runs   <- meta$runs
    n_rows <- n_feature_labels * n_runs

    protein_dt <- data.table::data.table(
        PROTEIN       = rep(meta$PROTEIN, n_rows),
        FEATURE       = rep(feature_label_dt$FEATURE,  times = n_runs),
        LABEL         = rep(feature_label_dt$LABEL,    times = n_runs),
        PEPTIDE       = rep(feature_label_dt$PEPTIDE,  times = n_runs),
        RUN           = rep(runs,         each  = n_feature_labels),
        newABUNDANCE  = as.vector(new_abundance_mat),
        n_obs         = as.integer(rep(n_obs_vec,     times = n_runs)),
        n_obs_run     = as.integer(rep(n_obs_run_vec, each  = n_feature_labels)),
        prop_features = rep(prop_features_vec, each = n_feature_labels)
    )

    # Optional: ABUNDANCE (labeled linear model uses raw log-intensities)
    if (meta$has_ABUNDANCE) {
        protein_dt[, ABUNDANCE := as.vector(abundance_mat)]
    }

    # censored: stored as 0.0/1.0 doubles; NA → treat as non-censored
    protein_dt[, censored := {
        v <- as.vector(censored_mat)
        if (meta$has_censored) as.logical(v > 0.5) else rep(FALSE, n_rows)
    }]

    # cen: survival event indicator (1 = observed, 0 = left-censored)
    if (meta$has_cen) {
        protein_dt[, cen := as.vector(event_mat)]
    }

    # ANOMALYSCORES: NA when unused (linear model skips anomaly weighting)
    protein_dt[, ANOMALYSCORES := as.vector(anomaly_scores_mat)]

    # SRM-specific columns derived from LABEL and RUN
    if (meta$is_labeled_ref) {
        protein_dt[, is_labeled_ref := (LABEL == "H")]
        if (meta$add_ref_covariate) {
            protein_dt[, ref_covariate := factor(
                data.table::fifelse(LABEL == "L", as.character(RUN), "0"))]
        }
    }

    protein_dt
}





#' Build the per-record worker closure for
#' \code{MSstatsSummarizeWithMultipleCores}
#'
#' Defined at package top level — not nested via \code{local()} inside
#' \code{MSstatsSummarizeWithMultipleCores} — so the returned closure's
#' enclosing environment chain is this factory's own (small) call frame plus
#' the package namespace. \code{BiocParallel}/\code{matter} serialize a
#' closure's entire enclosing environment chain to ship it to each socket
#' worker, not just the variables the closure body actually references. A
#' closure built via \code{local()} inside
#' \code{MSstatsSummarizeWithMultipleCores} would have that function's own
#' evaluation frame in its chain — which holds \code{input},
#' \code{protein_records}, and other run-scale objects — so every worker
#' would receive and retain a serialized copy of them even though
#' \code{.worker} never touches them. Building the closure here instead keeps
#' its captured state to only the scalar run parameters.
#'
#' @keywords internal
.build_summarize_worker <- function(
        use_TMP, impute, censored_symbol, remove50missing,
        aft_iterations, equal_variance
) {
    unpack_fn        <- .unpack_protein_slot
    use_TMP_         <- use_TMP
    impute_          <- impute
    censored_symbol_ <- censored_symbol
    remove50missing_ <- remove50missing
    aft_iterations_  <- aft_iterations
    equal_variance_  <- equal_variance

    function(record) {
        meta   <- record$meta
        protein_dt <- unpack_fn(record$packed, meta)
        result <- if (use_TMP_) {
            MSstatsSummarizeSingleTMP(
                protein_dt, impute_, censored_symbol_,
                remove50missing_, aft_iterations_)
        } else {
            MSstatsSummarizeSingleLinear(
                protein_dt, impute_, censored_symbol_,
                remove50missing_, aft_iterations_,
                equal_variances = equal_variance_)
        }
        # Normalize column types/levels at the worker boundary. Both
        # MSstatsSummarizeSingleTMP and MSstatsSummarizeSingleLinear call a
        # bare factor(RUN)/factor(FEATURE) on their input. For
        # MSstatsSummarizeWithSingleCore that's a no-op order-wise, since RUN
        # and FEATURE arrive already factors there and bare factor() keeps
        # existing levels. Here protein_dt comes from .unpack_protein_slot,
        # which reconstructs RUN/FEATURE as plain character — so the bare
        # factor() call inside those functions re-sorts alphabetically
        # instead ("1","10","11","2",... for >=10 runs). Re-leveling against
        # meta$runs restores the numeric run order that levels(input$RUN)
        # established upstream; as.character() first strips the alphabetical
        # level order factor() just introduced. droplevels() then matches
        # SingleCore's factor(RUN), which only ever sees — and so only ever
        # keeps — the runs present for that protein, since meta$runs carries
        # every run across the whole input.
        #
        # FEATURE is re-leveled with a fresh base R sort() rather than
        # inheriting meta$feature_label_dt's existing row order: that order
        # comes from .pack_protein_slot's data.table::setorder(feature_label_dt,
        # FEATURE, LABEL), and data.table sorts character columns in the
        # C-locale (byte order: all uppercase before any lowercase) for
        # platform-independence — whereas SingleCore's bare factor(FEATURE)
        # goes through base R's factor()/sort(), which use the session's
        # collation locale (e.g. en_US.UTF-8, where case is interleaved:
        # "a" < "A" < "b" < "B"). Any FEATURE strings with mixed case (e.g.
        # modification tags) sort differently between the two, so
        # re-deriving the level order with a base R sort() reproduces
        # SingleCore's order instead of inheriting data.table's.
        feature_levels <- sort(unique(meta$feature_label_dt$FEATURE))
        for (idx in 1:2) {
            if (!is.null(result[[idx]])) {
                if ("RUN" %in% colnames(result[[idx]]))
                    result[[idx]][, RUN := droplevels(factor(as.character(RUN), levels = meta$runs))]
                if ("FEATURE" %in% colnames(result[[idx]]))
                    result[[idx]][, FEATURE := droplevels(factor(as.character(FEATURE), levels = feature_levels))]
            }
        }
        if (!is.null(result[[2L]]) && "cen" %in% colnames(result[[2L]]))
            result[[2L]][, cen := as.integer(cen)]
        result
    }
}

#' Per-worker peak-RSS query task for \code{MSstatsSummarizeWithMultipleCores}
#'
#' Dispatched once per worker (via \code{seq_len(bpnworkers(BPPARAM))}) after
#' the main summarization \code{bplapply} call, while the persistent workers
#' are still alive, so each worker reports its own true lifetime-peak RSS
#' rather than a snapshot taken mid-run.
#'
#' @keywords internal
#' @noRd
.report_worker_peak <- function(i) {
    list(worker = i, pid = Sys.getpid(), peak_mb = .peak_rss_mb())
}

# Per-worker warm-up task for MSstatsSummarizeWithMultipleCores: loads
# MSstats once per persistent worker process and pins data.table to a single
# thread. Defined at top level rather than inline inside
# MSstatsSummarizeWithMultipleCores: even though its body never references
# `input`/`protein_records`, a closure defined inline there would still carry
# that function's evaluation frame in its enclosing environment chain, and
# BiocParallel would serialize that whole frame — including those run-scale
# objects — to every worker just to ship this no-op task.
.warmup_worker <- function(i) {
    library(MSstats, quietly = TRUE, warn.conflicts = FALSE)
    data.table::setDTthreads(1)
    NULL
}


#' Feature-level data summarization via socket-dispatched protein records
#'
#' Fixes a hidden RAM cost present in an earlier iteration of this function,
#' where the worker closure captured \code{meta_list} — the metadata for
#' \emph{every} protein — by reference, so each worker received and retained
#' the full metadata set for the whole run regardless of how many proteins
#' were actually assigned to it. Bounding the per-task \emph{packed-vector}
#' payload (e.g. via \code{SnowfastParam(tasks = ...)}) did nothing to reduce
#' that, since the metadata was baked into the closure once, not re-sliced
#' per task.
#'
#' This version pairs each protein's packed double vector with only its own
#' metadata into a single list element (\code{list(packed = ..., meta = ...)})
#' and dispatches that combined list through \code{bplapply}.
#' \code{BiocParallel} then only serializes and sends the records actually
#' assigned to a given task, so a worker never holds metadata for proteins
#' outside its own task(s) — unlike that earlier iteration, RAM scales with
#' the batch actually being processed.
#'
#' Progress is reported by turning on \code{SnowfastParam}'s built-in
#' \code{progressbar}, not by having workers talk back to the parent.
#' \code{BiocParallel} already ticks that progress bar from inside the
#' manager process itself, once per task result it collects — the same
#' receive that would happen regardless, over the same single
#' \code{bplapply} call. Workers never see the flag and send nothing extra
#' because of it, so this adds no serialization and no IPC beyond what an
#' unmonitored run already does. Reporting granularity therefore tracks
#' however many tasks the run is already split into (see
#' \code{max_proteins_per_worker} below): with the default \code{tasks = 0}
#' that's one step per worker.
#'
#' @param input feature-level data processed by dataProcess subfunctions
#' @param method summarization method: "linear" or "TMP" 
#' @param equal_variance only for summaryMethod = "linear". Default is TRUE. 
#' Logical variable for whether the model should account for heterogeneous variation 
#' among intensities from different features. Default is TRUE, which assume equal
#' variance among intensities from features. FALSE means that we cannot assume 
#' equal variance among intensities from features, then we will account for
#' heterogeneous variation from different features.
#' @param censored_symbol Indicates how censored missing values are encoded in
#' the 'Intensity' column. 'NA' (default) treats all NA intensities as
#' left-censored (i.e., below the limit of detection). '0' treats zero
#' intensities as left-censored; in this case NA intensities are assumed to be
#' missing at random and are not censored. Skyline output should use '0'. NULL
#' assumes that all missing values are missing at random — no values are treated
#' as censored, and imputation is disabled (impute is ignored).
#' @param remove50missing only for summaryMethod = "TMP". TRUE removes the proteins
#' where every run has at least 50\% missing values for each peptide. FALSE is default.
#' @param impute only for summaryMethod = "TMP" and censored_symbol = 'NA' or '0'.
#' TRUE (default) imputes censored missing values using an Accelerated Failure
#' Time model. FALSE excludes censored observations from summarization entirely,
#' treating them as missing at random; no imputed values are introduced.
#' Has no effect when censored_symbol = NULL, since no values are considered censored.
#' @param numberOfCores Number of cores for parallel processing. When > 1, 
#' a logfile named `MSstats_dataProcess_log_progress.log` is created to 
#' track progress. Only works for Linux & Mac OS. Default is 1.
#' @param max_proteins_per_worker integer; caps how many protein records
#'   (packed vector + its own metadata) are bundled into a single
#'   \code{bplapply} task sent to one worker. Translated into
#'   \code{SnowfastParam(tasks = ...)}: with \code{N} proteins this becomes
#'   \code{tasks = ceiling(N / max_proteins_per_worker)}. Only applied when
#'   \code{BPPARAM} is \code{NULL}; ignored if the caller supplies
#'   \code{BPPARAM} directly. Default \code{0} reproduces \code{tasks = 0}:
#'   \code{X} is divided as evenly as possible across \code{numberOfCores}
#'   workers.
#'
#' @return A named list with one element per protein slot, keyed by protein
#'   (or protein \eqn{\times} label) identifier.
#'
#' @importFrom matter SnowfastParam
#' @importFrom BiocParallel bplapply bpstart bpstop bpisup bpnworkers bpprogressbar
#' @importFrom data.table data.table fifelse setDTthreads
#' @importFrom stats median predict
#'
#' @export
MSstatsSummarizeWithMultipleCores <- function(
        input,
        method,
        impute,
        censored_symbol,
        remove50missing,
        equal_variance,
        numberOfCores  = 1L,
        aft_iterations = 90L,
        verbose        = FALSE,
        BPPARAM        = NULL,
        track_memory   = FALSE,
        max_proteins_per_worker = 0L
) {
    # ── 0. Single-core fallback ────────────────────────────────────────────────
    if (numberOfCores <= 1L && is.null(BPPARAM)) {
        return(MSstatsSummarizeWithSingleCore(
            input, method, impute, censored_symbol,
            remove50missing, equal_variance, aft_iterations))
    }
    
    start_time <- proc.time()[["elapsed"]]
    memory_checkpoints <- list()
    
    # ── 1. Split input by protein slot ────────────────────────────────────────
    is_labeled_reference <- "is_labeled_ref" %in% colnames(input) &&
        any(input$is_labeled_ref, na.rm = TRUE)
    split_keys <- if (is_labeled_reference) list(input$PROTEIN) else
        list(input$PROTEIN, input$LABEL)
    protein_indices <- split(seq_len(nrow(input)), split_keys)
    protein_ids     <- names(protein_indices)
    num_proteins    <- length(protein_indices)
    
    # Sort on input$RUN's own type (integer/numeric, typically) before
    # converting to character. Sorting the character form directly — as a
    # naive `sort(as.character(...))` would — collates lexicographically
    # ("1","10","11","2",...) whenever RUN isn't already a pre-existing
    # factor, e.g. when `input` comes straight from fread() and RUN reads in
    # as integer. This mirrors what factor() itself does internally for a
    # non-factor input: unique() + order() on the original values, then
    # as.character() only at the end.
    all_runs <- if (is.factor(input$RUN)) levels(input$RUN) else
        as.character(sort(unique(input$RUN)))
    
    getOption("MSstatsLog")("INFO",
                            paste0("Packing ", num_proteins, " proteins × ",
                                   length(all_runs), " runs into per-protein records"))

    # ── 2. Pack each protein into a (packed vector, own metadata) record ──────
    #
    # Unlike earlier iterations' parallel packed_list/meta_list arrays, each
    # protein's metadata travels bundled with its own packed vector.
    # bplapply/BPPARAM slice this combined list into tasks, so a worker's
    # closure never needs — and never receives — metadata for proteins
    # outside its assigned task(s).
    protein_records <- vector("list", num_proteins)
    for (slot_index in seq_len(num_proteins)) {
        protein_records[[slot_index]] <- .pack_protein_slot(
            input[protein_indices[[slot_index]], ], slot_index, all_runs)
    }
    
    payload_mb <- sum(vapply(protein_records,
                             function(r) length(r$packed), integer(1))) * 8 / 1024^2
    getOption("MSstatsLog")("INFO",
                            paste0("Dispatching via sockets (",
                                   format(round(payload_mb, 1)),
                                   " MB total packed payload; metadata sharded per task)"))
    
    # ── 3. Worker closure ─────────────────────────────────────────────────────
    #
    # No captured meta_list: each task's records already carry their own
    # metadata, so the closure only needs the scalar run parameters.
    #
    # Built via .build_summarize_worker() (defined at package top level)
    # rather than local() here, so the closure's enclosing environment chain
    # never includes this function's own frame — which holds `input`,
    # `protein_records`, etc. — and BiocParallel doesn't serialize those
    # run-scale objects to every worker.
    use_TMP <- identical(method, "TMP")

    worker_fn <- .build_summarize_worker(
        use_TMP, impute, censored_symbol, remove50missing,
        aft_iterations, equal_variance)
    
    # ── 4. Dispatch ───────────────────────────────────────────────────────────
    #
    # tasks controls how many protein records are bundled into one bplapply
    # task (i.e. one message sent to one worker). tasks == 0 (default) leaves
    # BiocParallel's own behavior in place: X divided as evenly as possible
    # over numberOfCores workers. When max_proteins_per_worker > 0, tasks is
    # sized so no task exceeds that many records — and because metadata now
    # travels with each record instead of being fully captured by the
    # closure, this actually bounds peak worker RAM.
    if (is.null(BPPARAM)) {
        tasks <- if (max_proteins_per_worker > 0L) {
            as.integer(ceiling(num_proteins / max_proteins_per_worker))
        } else {
            0L
        }
        getOption("MSstatsLog")("INFO",
                                paste0("Dispatching as ",
                                       if (tasks > 0L) tasks else numberOfCores,
                                       " task(s) (max_proteins_per_worker = ",
                                       max_proteins_per_worker, ")"))
        BPPARAM <- matter::SnowfastParam(
            workers       = numberOfCores,
            tasks         = tasks,
            progressbar   = TRUE,
            force.GC      = TRUE,
            stop.on.error = FALSE)
    }
    
    # ── Cluster setup ─────────────────────────────────────────────────────────
    # Load MSstats once per persistent worker process instead of once per
    # protein record, and pin data.table to a single thread per worker —
    # otherwise each of the numberOfCores workers would independently
    # auto-detect its own DT thread pool, oversubscribing the node's cores.
    started_here <- !BiocParallel::bpisup(BPPARAM)
    if (started_here) {
        BiocParallel::bpstart(BPPARAM)
        on.exit(BiocParallel::bpstop(BPPARAM), add = TRUE)
    }
    # Progress bar is on BPPARAM itself, so it would otherwise also print for
    # these two trivial one-task-per-worker calls. Toggle it off around them
    # and restore whatever it was (set above) so only the main summarization
    # bplapply — the one call whose progress is actually informative — shows
    # a bar.
    show_progress <- BiocParallel::bpprogressbar(BPPARAM)
    
    BiocParallel::bpprogressbar(BPPARAM) <- FALSE
    BiocParallel::bplapply(
        seq_len(BiocParallel::bpnworkers(BPPARAM)),
        .warmup_worker, BPPARAM = BPPARAM)
    BiocParallel::bpprogressbar(BPPARAM) <- show_progress
    
    results <- BiocParallel::bplapply(protein_records, worker_fn, BPPARAM = BPPARAM)
    names(results) <- protein_ids

    worker_peaks <- NULL
    if (track_memory) {
        BiocParallel::bpprogressbar(BPPARAM) <- FALSE
        worker_peaks <- BiocParallel::bplapply(
            seq_len(BiocParallel::bpnworkers(BPPARAM)),
            .report_worker_peak, BPPARAM = BPPARAM)   # must run BEFORE bpstop() while workers are alive
        BiocParallel::bpprogressbar(BPPARAM) <- show_progress
        memory_checkpoints[["parent peak (main)"]] <- .peak_rss_mb()
        worker_peak_mb <- vapply(worker_peaks, function(x) x$peak_mb, numeric(1L))
        .print_memory_report(
            "MSstatsSummarizeWithMultipleCores",
            memory_checkpoints, worker_peak_mb,
            elapsed = proc.time()[["elapsed"]] - start_time)
    }

    getOption("MSstatsLog")("INFO", "Summarization complete.")
    results
}
