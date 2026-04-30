#' Feature-level data summarization with multiple cores via matter out-of-core storage
#'
#' A drop-in replacement for \code{MSstatsSummarizeWithMultipleCores} that avoids
#' copying the full feature-level data to every worker.  Instead it:
#' \enumerate{
#'   \item Serializes each protein's data slice to a disk-backed
#'         \code{matter_list}, so workers memory-map only their assigned slice.
#'   \item Dispatches work with \code{chunkLapply} + \code{SnowfastParam}
#'         (GC between tasks, dynamic port selection, worker timeouts).
#'   \item Sizes chunks adaptively from \code{mem_realized} so that the
#'         number of concurrent proteins per worker scales with the actual
#'         serialized footprint rather than being fixed at one-per-task.
#'   \item Optionally streams results to disk via \code{outpath} so the main
#'         process never accumulates the full result list in RAM.
#' }
#'
#' @param input feature-level data prepared by \code{MSstatsPrepareForSummarization}.
#' @param method summarization method: \code{"TMP"} (Tukey median polish) or
#'   \code{"linear"}.
#' @param impute logical; passed to \code{MSstatsSummarizeSingleTMP} /
#'   \code{MSstatsSummarizeSingleLinear}.
#' @param censored_symbol character or NULL; how censored values are encoded.
#' @param remove50missing logical; remove proteins where every run has \eqn{\ge}50\%
#'   missing values.
#' @param equal_variance logical; assume equal feature variances (linear method).
#' @param numberOfCores integer number of parallel workers.  Values \eqn{\le}1
#'   fall back to \code{MSstatsSummarizeWithSingleCore}.
#' @param aft_iterations integer; maximum iterations for AFT survival model.
#' @param backing_path character path for the matter backing file that stores
#'   serialized per-protein input slices.  A \code{tempfile()} is used when
#'   \code{NULL} (default).
#' @param result_path character path for an optional matter backing file that
#'   streams serialized per-protein results to disk.  When \code{NULL} (default)
#'   results accumulate in the main process memory, matching the original
#'   behavior.
#' @param ram_overhead_factor numeric multiplier applied to the serialized
#'   per-protein size to estimate its in-memory footprint during model fitting.
#'   Default \code{8}; increase for experiments with many features or heavy AFT
#'   convergence.
#' @param verbose logical; passed to \code{chunkLapply} for chunk-progress
#'   messages.
#'
#' @return A named list with one element per protein (or protein × label group).
#'   Each element is \code{list(ProteinLevelData, SurvivalData)} matching the
#'   format of \code{MSstatsSummarizeWithMultipleCores}.
#'
#' @importFrom matter matter_list mem_realized SnowfastParam chunkLapply
#'
#' @export
MSstatsSummarizeWithMultipleCoresV2 <- function(
    input,
    method,
    impute,
    censored_symbol,
    remove50missing,
    equal_variance,
    numberOfCores    = 1L,
    aft_iterations   = 90L,
    backing_path     = NULL,
    result_path      = NULL,
    ram_overhead_factor = 8,
    verbose          = FALSE
) {
    # ── 0. Single-core fallback ────────────────────────────────────────────────
    if (numberOfCores <= 1L) {
        return(MSstatsSummarizeWithSingleCore(
            input, method, impute, censored_symbol,
            remove50missing, equal_variance, aft_iterations))
    }

    # ── 1. Split input by protein (replicating original logic) ────────────────
    is_labeled_reference <- "is_labeled_ref" %in% colnames(input) &&
        any(input$is_labeled_ref, na.rm = TRUE)

    split_keys <- if (is_labeled_reference) {
        list(input$PROTEIN)
    } else {
        list(input$PROTEIN, input$LABEL)
    }
    protein_indices <- split(seq_len(nrow(input)), split_keys)
    protein_ids     <- names(protein_indices)
    num_proteins    <- length(protein_indices)

    getOption("MSstatsLog")("INFO",
        paste0("V2: serializing ", num_proteins,
               " proteins to out-of-core backing store"))

    # ── 2. Build matter_list (one element = raw-serialized protein data.table) ─
    #
    # Each element is a raw vector produced by serialize().  Workers memory-map
    # the backing file and read only their assigned slice — O(protein_size) per
    # worker instead of O(full_dataset) via clusterExport.
    if (is.null(backing_path)) {
        backing_path <- tempfile(fileext = ".bin")
    }

    protein_raw_list <- lapply(protein_indices, function(idx) {
        serialize(input[idx, ], connection = NULL)
    })

    protein_matter <- matter::matter_list(
        protein_raw_list,
        type  = "raw",
        path  = backing_path,
        names = protein_ids
    )

    # Drop the in-memory copy; workers will read from disk.
    rm(protein_raw_list)
    gc(verbose = FALSE)

    getOption("MSstatsLog")("INFO",
        paste0("V2: backing store written to ", backing_path,
               " (", format(matter::mem_realized(protein_matter)), ")"))

    # ── 3. Adaptive chunk sizing ──────────────────────────────────────────────
    #
    # mem_realized() returns the on-disk footprint without reading any data.
    # Multiplying by ram_overhead_factor approximates the in-memory cost of
    # deserializing a protein slice plus the residuals of AFT / LME model fits.
    per_protein_bytes <- as.numeric(matter::mem_realized(protein_matter)) /
        num_proteins
    effective_bytes_per_protein <- per_protein_bytes * ram_overhead_factor

    # Estimate headroom from gc() — conservative: 50% of currently used Vcells.
    gc_stats     <- gc(reset = FALSE)
    used_bytes   <- sum(gc_stats[, "used"]) * 8   # Vcells are 8-byte doubles
    worker_budget_bytes <- max(256 * 1024^2,       # floor: 256 MB
                               used_bytes * 0.5)   # half of current session RAM

    n_concurrent <- max(1L,
        floor(worker_budget_bytes /
              (numberOfCores * effective_bytes_per_protein)))

    getOption("MSstatsLog")("INFO",
        paste0("V2: chunk size = ", n_concurrent,
               " proteins/worker (estimated ",
               format(round(effective_bytes_per_protein / 1024^2, 1)),
               " MB per protein)"))

    # ── 4. Parallel backend ───────────────────────────────────────────────────
    #
    # SnowfastParam over SnowParam / makeCluster because:
    #   • force.GC = TRUE: GC between tasks prevents residual AFT / LME objects
    #     from piling up across proteins in a single worker session.
    #   • Random port: no collision when multiple MSstats sessions run in
    #     parallel (e.g., Bioconductor build machines).
    #   • stop.on.error = FALSE: a single failed protein does not abort the run;
    #     callers already handle NULL results.
    BPPARAM <- matter::SnowfastParam(
        workers       = numberOfCores,
        force.GC      = TRUE,
        stop.on.error = FALSE
    )

    # ── 5. Per-protein worker function ────────────────────────────────────────
    #
    # Capture all scalar parameters in a local environment so BiocParallel
    # serializes only the values (not the full calling frame) to each worker.
    # Using MSstats:: ensures the worker only needs the package installed, not
    # a prior library() call.
    use_TMP <- identical(method, "TMP")

    .processProtein <- local({
        use_TMP_          <- use_TMP
        impute_           <- impute
        censored_symbol_  <- censored_symbol
        remove50missing_  <- remove50missing
        aft_iterations_   <- aft_iterations
        equal_variance_   <- equal_variance

        function(raw_bytes) {
            library(MSstats, quietly = TRUE, warn.conflicts = FALSE)
            single_protein <- unserialize(raw_bytes)
            if (use_TMP_) {
                MSstatsSummarizeSingleTMP(
                    single_protein,
                    impute_,
                    censored_symbol_,
                    remove50missing_,
                    aft_iterations_)
            } else {
                MSstatsSummarizeSingleLinear(
                    single_protein,
                    impute_,
                    censored_symbol_,
                    remove50missing_,
                    aft_iterations_,
                    equal_variances = equal_variance_)
            }
        }
    })

    # ── 6. Dispatch via chunkLapply ───────────────────────────────────────────
    #
    # chunkopts$chunksize controls how many proteins are batched into each
    # bplapply task (i.e., how many each worker processes before returning
    # results to the manager).  Adaptive sizing from Step 3 ensures that
    # peak worker-RAM stays within budget regardless of protein size.
    #
    # When result_path is supplied, chunk_writer streams serialized results to
    # a matter_list on disk — the main process never holds the full result list.
    if (!is.null(result_path)) {
        # Wrap FUN to return serialized bytes so chunk_writer can store them in
        # a matter_list (which requires atomic vector elements).
        .processProteinStreaming <- local({
            inner_ <- .processProtein
            function(raw_bytes) serialize(inner_(raw_bytes), connection = NULL)
        })

        result_matter <- matter::chunkLapply(
            protein_matter,
            FUN      = .processProteinStreaming,
            outpath  = result_path,
            chunkopts = list(chunksize = n_concurrent),
            verbose  = verbose,
            BPPARAM  = BPPARAM
        )

        getOption("MSstatsLog")("INFO",
            paste0("V2: results streamed to ", result_path,
                   "; deserializing into list"))

        # Materialize: result_matter is a matter_list of raw-serialized results.
        summarized_results <- lapply(seq_len(num_proteins), function(i) {
            unserialize(result_matter[[i]])
        })
        names(summarized_results) <- protein_ids

    } else {
        summarized_results <- matter::chunkLapply(
            protein_matter,
            FUN      = .processProtein,
            chunkopts = list(chunksize = n_concurrent),
            verbose  = verbose,
            BPPARAM  = BPPARAM
        )
        names(summarized_results) <- protein_ids
    }

    getOption("MSstatsLog")("INFO", "V2: summarization complete.")
    summarized_results
}


## ── Memory-monitoring helpers ─────────────────────────────────────────────────

# RSS of the current process in MB.
# On Linux reads /proc/self/status (VmRSS); elsewhere falls back to gc() counts.
.memMB <- function() {
    if (file.exists("/proc/self/status")) {
        ln <- readLines("/proc/self/status", warn = FALSE)
        m  <- grep("^VmRSS:", ln, value = TRUE)
        if (length(m))
            return(as.numeric(gsub("[^0-9]", "", m[1L])) / 1024)
    }
    g <- gc(reset = FALSE)
    (g["Ncells", "used"] * 8L + g["Vcells", "used"] * 8L) / 1024^2
}

# Print a formatted memory report to stderr via message().
# checkpoints: named numeric vector of RSS snapshots (MB).
# worker_mems: numeric vector of per-task RSS values from workers (may be NA).
# elapsed:     total wall-clock seconds (NULL to omit).
.printMemReport <- function(fn_name, checkpoints, worker_mems = NULL,
                            elapsed = NULL) {
    w  <- 65L
    hr <- strrep("─", w)
    fmt_mb    <- function(x) if (is.na(x)) "    n/a" else sprintf("%7.1f", x)
    fmt_delta <- function(now, prev) {
        if (is.na(now) || is.na(prev)) return("")
        sprintf("  (%+.1f MB)", now - prev)
    }
    lines <- c(
        hr,
        sprintf(" MSstats Memory Report — %s", fn_name),
        hr,
        sprintf("  %-36s  %7s  %s", "Checkpoint", "RSS MB", "Delta")
    )
    prev <- NA_real_
    for (nm in names(checkpoints)) {
        val <- checkpoints[[nm]]
        lines <- c(lines,
            sprintf("  %-36s  %s%s", nm, fmt_mb(val), fmt_delta(val, prev)))
        prev <- val
    }
    if (!is.null(worker_mems)) {
        wm <- worker_mems[!is.na(worker_mems)]
        if (length(wm)) {
            lines <- c(lines, "",
                sprintf("  Worker RSS  min / mean / max :  %.1f / %.1f / %.1f MB",
                        min(wm), mean(wm), max(wm)))
        }
    }
    if (!is.null(elapsed))
        lines <- c(lines, sprintf("  Total elapsed: %.1f s", as.numeric(elapsed)))
    lines <- c(lines, hr)
    message(paste(lines, collapse = "\n"))
}

## ── V3 internal helpers ───────────────────────────────────────────────────────
##
## Packed double-vector layout for one protein slot (all column-major matrices):
##
##  pos 1                          : slot_k  (protein index, cast to double)
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
#' @param dt data.table rows belonging to one protein (or protein × label) slot
#' @param slot_k integer position of this slot in the global protein list
#' @param all_runs character vector of all run names in global order
#' @return list with elements \code{packed} (double vector) and \code{meta} (list)
#' @keywords internal
.buildProteinSlotV3 <- function(dt, slot_k, all_runs) {

    R <- length(all_runs)

    # ── Unique (FEATURE, LABEL) → "effective features", stable ordering ───────
    has_peptide <- "PEPTIDE" %in% colnames(dt)
    flp <- unique(dt[, .(
        FEATURE = as.character(FEATURE),
        LABEL   = as.character(LABEL),
        PEPTIDE = if (has_peptide) as.character(PEPTIDE) else as.character(FEATURE)
    )])
    data.table::setorder(flp, FEATURE, LABEL)
    FL <- nrow(flp)

    # ── Index maps: feature-label → row index, run name → column index ─────────
    fi <- match(
        paste(as.character(dt$FEATURE), as.character(dt$LABEL), sep = "\t"),
        paste(flp$FEATURE, flp$LABEL, sep = "\t"))
    ri_map <- seq_len(R); names(ri_map) <- all_runs
    ri     <- ri_map[as.character(dt$RUN)]
    valid  <- !is.na(fi) & !is.na(ri)

    # ── Helper: allocate FL × R matrix and fill from dt column ───────────────
    fill_mat <- function(col_vec, fill = NA_real_) {
        m <- matrix(fill, nrow = FL, ncol = R)
        m[cbind(fi[valid], ri[valid])] <- as.double(col_vec[valid])
        m
    }

    # ── Build the five numeric matrices ───────────────────────────────────────
    mat_newABU <- fill_mat(dt$newABUNDANCE)

    has_ABUNDANCE <- "ABUNDANCE" %in% colnames(dt)
    mat_ABU <- if (has_ABUNDANCE) fill_mat(dt$ABUNDANCE) else
        matrix(NA_real_, FL, R)

    has_censored <- "censored" %in% colnames(dt)
    mat_cens <- if (has_censored) fill_mat(as.double(dt$censored)) else
        matrix(0.0, FL, R)   # treat as non-censored when column absent

    has_cen <- "cen" %in% colnames(dt)
    mat_cen <- if (has_cen) fill_mat(dt$cen) else
        matrix(NA_real_, FL, R)

    has_anom <- "ANOMALYSCORES" %in% colnames(dt) &&
        !all(is.na(dt$ANOMALYSCORES))
    mat_anom <- if (has_anom) fill_mat(dt$ANOMALYSCORES) else
        matrix(NA_real_, FL, R)

    # ── Per-feature scalar: n_obs (constant within FEATURE × LABEL) ───────────
    n_obs_by_fl <- dt[, .(n_obs = as.double(n_obs[1L])),
                       by = .(FEATURE = as.character(FEATURE),
                              LABEL   = as.character(LABEL))]
    flp_nobs    <- n_obs_by_fl[flp, on = c("FEATURE", "LABEL")]
    data.table::setorder(flp_nobs, FEATURE, LABEL)
    n_obs_vec   <- flp_nobs$n_obs

    # ── Per-run scalars: n_obs_run, prop_features (constant within RUN) ───────
    run_vals <- dt[, .(n_obs_run    = as.double(n_obs_run[1L]),
                        prop_features = as.double(prop_features[1L])),
                    by = .(RUN = as.character(RUN))]
    run_full <- run_vals[data.table::data.table(RUN = all_runs), on = "RUN"]
    n_obs_run_vec     <- run_full$n_obs_run
    prop_features_vec <- run_full$prop_features

    # ── Pack ─────────────────────────────────────────────────────────────────
    packed <- c(
        as.double(slot_k),
        as.vector(mat_newABU),
        as.vector(mat_ABU),
        as.vector(mat_cens),
        as.vector(mat_cen),
        as.vector(mat_anom),
        n_obs_vec,
        n_obs_run_vec,
        prop_features_vec
    )

    # ── Metadata (string labels, flags; kept in main-process RAM) ────────────
    meta <- list(
        PROTEIN           = as.character(dt$PROTEIN[1L]),
        feat_label_pep    = as.data.frame(flp),   # FL × 3: FEATURE, LABEL, PEPTIDE
        runs              = all_runs,
        FL                = FL,
        R                 = R,
        is_labeled_ref    = "is_labeled_ref" %in% colnames(dt) &&
            isTRUE(any(dt$is_labeled_ref, na.rm = TRUE)),
        has_ABUNDANCE     = has_ABUNDANCE,
        has_censored      = has_censored,
        has_cen           = has_cen,
        has_anom          = has_anom,
        add_ref_covariate = "ref_covariate" %in% colnames(dt)
    )

    list(packed = packed, meta = meta)
}


#' Reconstruct a per-protein data.table from a V3 packed double vector
#'
#' @param packed double vector produced by \code{.buildProteinSlotV3}
#' @param meta metadata list from \code{.buildProteinSlotV3}
#' @return data.table compatible with \code{MSstatsSummarizeSingleTMP} /
#'   \code{MSstatsSummarizeSingleLinear}
#' @keywords internal
.reconstructProteinDTV3 <- function(packed, meta) {

    FL    <- meta$FL
    R     <- meta$R
    n_mat <- FL * R

    # ── Unpack sections (1-indexed; position 1 is slot_k header) ─────────────
    i <- 2L
    extract_mat <- function() {
        m <- matrix(packed[i:(i + n_mat - 1L)], nrow = FL, ncol = R)
        i <<- i + n_mat
        m
    }
    mat_newABU <- extract_mat()
    mat_ABU    <- extract_mat()
    mat_cens   <- extract_mat()
    mat_cen    <- extract_mat()
    mat_anom   <- extract_mat()

    n_obs_vec       <- packed[i:(i + FL - 1L)]; i <- i + FL
    n_obs_run_vec   <- packed[i:(i + R  - 1L)]; i <- i + R
    prop_feat_vec   <- packed[i:(i + R  - 1L)]

    # ── Column-major melt: rep/each match the matrix column-major ordering ────
    flp    <- meta$feat_label_pep   # data.frame: FEATURE, LABEL, PEPTIDE (FL rows)
    runs   <- meta$runs
    n_rows <- FL * R

    dt <- data.table::data.table(
        PROTEIN       = rep(meta$PROTEIN, n_rows),
        FEATURE       = rep(flp$FEATURE,  times = R),
        LABEL         = rep(flp$LABEL,    times = R),
        PEPTIDE       = rep(flp$PEPTIDE,  times = R),
        RUN           = rep(runs,         each  = FL),
        newABUNDANCE  = as.vector(mat_newABU),
        n_obs         = as.integer(rep(n_obs_vec,     times = R)),
        n_obs_run     = as.integer(rep(n_obs_run_vec, each  = FL)),
        prop_features = rep(prop_feat_vec, each = FL)
    )

    # Optional: ABUNDANCE (labeled linear model uses raw log-intensities)
    if (meta$has_ABUNDANCE) {
        dt[, ABUNDANCE := as.vector(mat_ABU)]
    }

    # censored: stored as 0.0/1.0 doubles; NA → treat as non-censored
    dt[, censored := {
        v <- as.vector(mat_cens)
        if (meta$has_censored) as.logical(v > 0.5) else rep(FALSE, n_rows)
    }]

    # cen: survival event indicator (1 = observed, 0 = left-censored)
    if (meta$has_cen) {
        dt[, cen := as.vector(mat_cen)]
    }

    # ANOMALYSCORES: NA when unused (linear model skips anomaly weighting)
    dt[, ANOMALYSCORES := as.vector(mat_anom)]

    # SRM-specific columns derived from LABEL and RUN
    if (meta$is_labeled_ref) {
        dt[, is_labeled_ref := (LABEL == "H")]
        if (meta$add_ref_covariate) {
            dt[, ref_covariate := factor(
                data.table::fifelse(LABEL == "L", as.character(RUN), "0"))]
        }
    }

    dt
}


#' Feature-level data summarization via native matter numeric matrices (V3)
#'
#' Improves on \code{MSstatsSummarizeWithMultipleCoresV2} by eliminating
#' \code{serialize()}/\code{unserialize()} overhead.  Each protein's numeric
#' columns (\code{newABUNDANCE}, \code{ABUNDANCE}, \code{censored}, \code{cen},
#' \code{ANOMALYSCORES}) are stored as typed double vectors in a single
#' disk-backed \code{matter_list}, and the per-run/per-feature summary scalars
#' (\code{n_obs}, \code{n_obs_run}, \code{prop_features}) are appended to the
#' same vector.  Only lightweight string metadata (factor levels, protein names)
#' lives in the main-process R heap.
#'
#' Workers receive a matter descriptor (file path + offsets), memory-map the
#' backing file, read their slice as a raw double vector, and call
#' \code{.reconstructProteinDTV3} to rebuild the data.table — no R object
#' serialization anywhere in the hot path.
#'
#' @inheritParams MSstatsSummarizeWithMultipleCoresV2
#'
#' @return A named list with one element per protein slot, identical in
#'   structure to \code{MSstatsSummarizeWithMultipleCores}.
#'
#' @importFrom matter matter_list mem_realized SnowfastParam chunkLapply
#' @importFrom data.table data.table setorder fifelse
#'
#' @export
MSstatsSummarizeWithMultipleCoresV3 <- function(
    input,
    method,
    impute,
    censored_symbol,
    remove50missing,
    equal_variance,
    numberOfCores       = 1L,
    aft_iterations      = 90L,
    backing_path        = NULL,
    result_path         = NULL,
    ram_overhead_factor = 4,
    verbose             = FALSE
) {
    # ── 0. Single-core fallback ────────────────────────────────────────────────
    if (numberOfCores <= 1L) {
        return(MSstatsSummarizeWithSingleCore(
            input, method, impute, censored_symbol,
            remove50missing, equal_variance, aft_iterations))
    }

    # ── 1. Split input by protein slot ────────────────────────────────────────
    is_labeled_reference <- "is_labeled_ref" %in% colnames(input) &&
        any(input$is_labeled_ref, na.rm = TRUE)
    split_keys <- if (is_labeled_reference) {
        list(input$PROTEIN)
    } else {
        list(input$PROTEIN, input$LABEL)
    }
    protein_indices <- split(seq_len(nrow(input)), split_keys)
    protein_ids     <- names(protein_indices)
    num_proteins    <- length(protein_indices)

    # Global run ordering: columns of every per-protein FL × R matrix
    all_runs <- if (is.factor(input$RUN)) levels(input$RUN) else
        sort(unique(as.character(input$RUN)))
    R <- length(all_runs)

    getOption("MSstatsLog")("INFO",
        paste0("V3: building native matter matrices for ", num_proteins,
               " proteins × ", R, " runs"))

    # ── 2. Build per-protein packed vectors and collect metadata ──────────────
    #
    # Each protein's numeric columns are packed into a flat double vector and
    # stored in the matter_list backing file.  Only the string metadata (factor
    # levels, protein names) remains in main-process RAM.
    if (is.null(backing_path)) {
        backing_path <- tempfile(fileext = ".bin")
    }

    slots <- vector("list", num_proteins)
    for (k in seq_len(num_proteins)) {
        slots[[k]] <- .buildProteinSlotV3(
            dt       = input[protein_indices[[k]], ],
            slot_k   = k,
            all_runs = all_runs
        )
    }

    packed_list <- lapply(slots, `[[`, "packed")
    meta_list   <- lapply(slots, `[[`, "meta")
    rm(slots)

    protein_matter <- matter::matter_list(
        packed_list,
        type  = "double",
        path  = backing_path,
        names = protein_ids
    )
    rm(packed_list)
    gc(verbose = FALSE)

    getOption("MSstatsLog")("INFO",
        paste0("V3: backing store written to ", backing_path,
               " (", format(matter::mem_realized(protein_matter)), ")"))

    # ── 3. Adaptive chunk sizing ──────────────────────────────────────────────
    #
    # mem_realized() reads stored atom metadata (no I/O) for the total
    # on-disk footprint.  ram_overhead_factor accounts for the in-memory
    # expansion when reconstructing the data.table (column vectors, factor
    # levels, model fitting residuals).  V3 uses factor 4 (vs. 8 in V2)
    # because the double vector is already a compact representation — no
    # serialization framing overhead.
    per_protein_bytes <- as.numeric(matter::mem_realized(protein_matter)) /
        num_proteins
    effective_bytes   <- per_protein_bytes * ram_overhead_factor

    gc_stats <- gc(reset = FALSE)
    used_bytes <- sum(gc_stats[, "used"]) * 8
    worker_budget_bytes <- max(256 * 1024^2, used_bytes * 0.5)
    n_concurrent <- max(1L,
        floor(worker_budget_bytes / (numberOfCores * effective_bytes)))

    getOption("MSstatsLog")("INFO",
        paste0("V3: chunk size = ", n_concurrent,
               " proteins/worker (estimated ",
               format(round(effective_bytes / 1024^2, 2)), " MB per protein)"))

    # ── 4. Parallel backend ───────────────────────────────────────────────────
    BPPARAM <- matter::SnowfastParam(
        workers       = numberOfCores,
        force.GC      = TRUE,
        stop.on.error = FALSE
    )

    # ── 5. Per-protein worker function ────────────────────────────────────────
    #
    # The closure captures meta_list (string metadata, typically a few MB)
    # and the two reconstruction / summarization helpers.  Workers receive
    # the closure once; subsequent tasks in the same worker reuse it.
    #
    # .reconstructProteinDTV3 is captured so workers do not need MSstats
    # internals beyond what library(MSstats) provides.
    use_TMP <- identical(method, "TMP")

    .processProteinV3 <- local({
        meta_list_          <- meta_list
        reconstruct_        <- .reconstructProteinDTV3
        use_TMP_            <- use_TMP
        impute_             <- impute
        censored_symbol_    <- censored_symbol
        remove50missing_    <- remove50missing
        aft_iterations_     <- aft_iterations
        equal_variance_     <- equal_variance

        function(packed) {
            library(MSstats, quietly = TRUE, warn.conflicts = FALSE)
            k  <- as.integer(packed[1L])
            dt <- reconstruct_(packed, meta_list_[[k]])
            if (use_TMP_) {
                MSstatsSummarizeSingleTMP(
                    dt, impute_, censored_symbol_,
                    remove50missing_, aft_iterations_)
            } else {
                MSstatsSummarizeSingleLinear(
                    dt, impute_, censored_symbol_,
                    remove50missing_, aft_iterations_,
                    equal_variances = equal_variance_)
            }
        }
    })

    # ── 6. Dispatch ───────────────────────────────────────────────────────────
    if (!is.null(result_path)) {
        .processProteinV3Streaming <- local({
            inner_ <- .processProteinV3
            function(packed) serialize(inner_(packed), connection = NULL)
        })
        result_matter <- matter::chunkLapply(
            protein_matter,
            FUN       = .processProteinV3Streaming,
            outpath   = result_path,
            chunkopts = list(chunksize = n_concurrent),
            verbose   = verbose,
            BPPARAM   = BPPARAM
        )
        summarized_results <- lapply(seq_len(num_proteins), function(i) {
            unserialize(result_matter[[i]])
        })
        names(summarized_results) <- protein_ids
    } else {
        summarized_results <- matter::chunkLapply(
            protein_matter,
            FUN       = .processProteinV3,
            chunkopts = list(chunksize = n_concurrent),
            verbose   = verbose,
            BPPARAM   = BPPARAM
        )
        names(summarized_results) <- protein_ids
    }

    getOption("MSstatsLog")("INFO", "V3: summarization complete.")
    summarized_results
}


## ── V4: matrix-native TMP ─────────────────────────────────────────────────────
##
## MSstatsSummarizeSingleTMPV2 takes the V3 packed double vector + meta list
## directly and produces TMP results without ever building the full long-format
## data.table or calling dcast.
##
## Compared with MSstatsSummarizeSingleTMP the hot path is:
##   V3: packed → long-DT (melt) → filter → AFT fit → dcast → TMP
##   V4: packed → FL×R matrix → mask → (AFT fit if needed) → t(mat) → TMP
##
## The dcast savings are most significant for high-feature proteins (>100
## features) and large run counts (>50 runs), where the intermediate FL*R
## data.table allocation dominates.
## ─────────────────────────────────────────────────────────────────────────────


#' Summarize a single protein with TMP directly from a V3 packed double vector
#'
#' Bypasses the long-format \code{data.table} reconstruction and the
#' \code{dcast} inside \code{.fitTukey}.  The FL×R packed matrices are
#' operated on directly:
#' \enumerate{
#'   \item Row/column masking replaces the \code{n_obs}/\code{n_obs_run} filter.
#'   \item AFT survival fitting uses a lazily-built minimal \code{data.table}
#'     (constructed only when \code{impute=TRUE} and censored values are
#'     present).
#'   \item TMP is applied via a vectorized scatter into a
#'     \code{(LABEL×RUN) × FEATURE} wide matrix followed by
#'     \code{median_polish_summary} — no \code{dcast} round-trip.
#' }
#'
#' @param packed double vector produced by \code{.buildProteinSlotV3}
#' @param meta metadata list from \code{.buildProteinSlotV3}
#' @param impute logical; impute censored values with AFT survival model
#' @param censored_symbol \code{"0"}, \code{"NA"}, or \code{NULL}
#' @param remove50missing logical; skip proteins where all runs are >50\%
#'   missing
#' @param aft_iterations integer; max AFT iterations
#' @return \code{list(result_dt, survival_dt)} matching the format of
#'   \code{MSstatsSummarizeSingleTMP}
#' @keywords internal
MSstatsSummarizeSingleTMPV2 <- function(packed, meta,
    impute, censored_symbol, remove50missing, aft_iterations = 90L)
{
    FL      <- meta$FL
    R       <- meta$R
    n_mat   <- FL * R
    PROTEIN <- meta$PROTEIN
    flp     <- meta$feat_label_pep   # data.frame: FEATURE, LABEL, PEPTIDE
    runs    <- meta$runs

    # ── 1. Unpack matrices ─────────────────────────────────────────────────────
    # Layout (see top-of-file comment block):
    #   pos 1           : slot_k (skip)
    #   pos 2..FL*R+1   : newABUNDANCE (FL × R, column-major)
    #   pos FL*R+2..    : ABUNDANCE    (skip — not needed for TMP)
    #   pos 2*FL*R+2..  : censored     (FL × R)
    #   pos 3*FL*R+2..  : cen          (FL × R)
    #   pos 4*FL*R+2..  : ANOMALYSCORES(skip — not needed for TMP)
    #   then: n_obs(FL), n_obs_run(R), prop_features(R)
    i <- 2L
    mat_newABU <- matrix(packed[i:(i + n_mat - 1L)], nrow = FL, ncol = R); i <- i + n_mat
    i          <- i + n_mat                  # skip ABUNDANCE
    mat_cens   <- matrix(packed[i:(i + n_mat - 1L)], nrow = FL, ncol = R); i <- i + n_mat
    mat_cen    <- matrix(packed[i:(i + n_mat - 1L)], nrow = FL, ncol = R); i <- i + n_mat
    i          <- i + n_mat                  # skip ANOMALYSCORES
    n_obs_vec     <- packed[i:(i + FL - 1L)]; i <- i + FL
    n_obs_run_vec <- packed[i:(i + R  - 1L)]; i <- i + R
    prop_feat_vec <- packed[i:(i + R  - 1L)]

    # ── 2. Row/column masking (mirrors the n_obs / n_obs_run filter) ───────────
    keep_feat <- !is.na(n_obs_vec)    & n_obs_vec    > 1
    keep_runs <- !is.na(n_obs_run_vec) & n_obs_run_vec > 0

    mat_newABU <- mat_newABU[keep_feat, keep_runs, drop = FALSE]
    mat_cens   <- mat_cens  [keep_feat, keep_runs, drop = FALSE]
    mat_cen    <- mat_cen   [keep_feat, keep_runs, drop = FALSE]

    flp_filt       <- flp[keep_feat, , drop = FALSE]
    runs_filt      <- runs[keep_runs]
    prop_feat_filt <- prop_feat_vec[keep_runs]
    FL_filt        <- nrow(mat_newABU)
    R_filt         <- ncol(mat_newABU)

    # ── 3. Empty-matrix early exit ─────────────────────────────────────────────
    if (FL_filt == 0L || R_filt == 0L) {
        msg <- paste("Can't summarize for protein", PROTEIN,
                     "because all measurements are missing or censored.")
        try(getOption("MSstatsMsg")("INFO", msg), silent = TRUE)
        try(getOption("MSstatsLog")("INFO", msg), silent = TRUE)
        return(list(NULL, NULL))
    }

    # ── 4. remove50missing check (independent of imputation) ──────────────────
    if (remove50missing &&
        all(prop_feat_filt <= 0.5 | is.na(prop_feat_filt))) {
        msg <- paste("Can't summarize for protein", PROTEIN,
                     "because all runs have more than 50% missing values and",
                     "are removed with the option, remove50missing=TRUE.")
        try(getOption("MSstatsMsg")("INFO", msg), silent = TRUE)
        try(getOption("MSstatsLog")("INFO", msg), silent = TRUE)
        return(list(NULL, NULL))
    }

    # ── 5. AFT survival imputation (lazy: only if any censored present) ────────
    n_rows_filt <- FL_filt * R_filt
    orig_abu_vec <- as.vector(mat_newABU)

    surv_cols <- c("newABUNDANCE", "cen", "RUN", "FEATURE", "ref_covariate")
    any_censored <- meta$has_censored && any(mat_cens > 0.5, na.rm = TRUE)

    if (impute && any_censored && meta$has_cen) {
        # Build minimal long-format DT — only columns needed by .fitSurvival
        surv_dt <- data.table::data.table(
            newABUNDANCE = orig_abu_vec,
            cen          = as.vector(mat_cen),
            censored     = as.logical(as.vector(mat_cens) > 0.5),
            FEATURE      = factor(rep(flp_filt$FEATURE, times = R_filt)),
            RUN          = factor(rep(runs_filt,         each  = FL_filt)),
            LABEL        = rep(as.character(flp_filt$LABEL), times = R_filt)
        )
        if (meta$add_ref_covariate) {
            surv_dt[, ref_covariate := factor(
                data.table::fifelse(LABEL == "L", as.character(RUN), "0"))]
        }

        fit_cols <- intersect(surv_cols, colnames(surv_dt))
        fit_data <- if (meta$is_labeled_ref) {
            surv_dt[LABEL != "H", fit_cols, with = FALSE]
        } else {
            surv_dt[, fit_cols, with = FALSE]
        }

        converged <- TRUE
        survival_fit <- withCallingHandlers({
            .fitSurvival(fit_data, aft_iterations)
        }, warning = function(w) {
            if (grepl("converge", conditionMessage(w), ignore.case = TRUE)) {
                message("Convergence warning caught: ", conditionMessage(w))
                converged <<- FALSE
            }
        })

        predicted_all <- if (converged) {
            predict(survival_fit, newdata = surv_dt)
        } else {
            rep(NA_real_, n_rows_filt)
        }

        use_imputed <- if (meta$is_labeled_ref) {
            surv_dt$censored & surv_dt$LABEL != "H"
        } else {
            surv_dt$censored
        }

        imputed_abu <- ifelse(use_imputed, predicted_all, orig_abu_vec)
        surv_dt[, predicted    := ifelse(use_imputed, predicted_all, NA_real_)]
        surv_dt[, newABUNDANCE := imputed_abu]
        mat_newABU <- matrix(imputed_abu, nrow = FL_filt, ncol = R_filt)
        survival   <- surv_dt[, intersect(c(surv_cols, "LABEL", "predicted"),
                                          colnames(surv_dt)), with = FALSE]

    } else {
        # No imputation: build survival DT from original values
        surv_dt <- data.table::data.table(
            newABUNDANCE = orig_abu_vec,
            cen          = if (meta$has_cen) as.vector(mat_cen)
                           else rep(NA_real_, n_rows_filt),
            FEATURE      = rep(as.character(flp_filt$FEATURE), times = R_filt),
            RUN          = rep(runs_filt,                       each  = FL_filt),
            LABEL        = rep(as.character(flp_filt$LABEL),    times = R_filt),
            predicted    = NA_real_
        )
        survival <- surv_dt[, intersect(c(surv_cols, "LABEL", "predicted"),
                                        colnames(surv_dt)), with = FALSE]
    }

    # ── 6. Post-imputation .isSummarizable check ───────────────────────────────
    if (all(is.na(mat_newABU) | mat_newABU == 0)) {
        msg <- paste("Can't summarize for protein", PROTEIN,
                     "because all measurements are missing or censored.")
        try(getOption("MSstatsMsg")("INFO", msg), silent = TRUE)
        try(getOption("MSstatsLog")("INFO", msg), silent = TRUE)
        return(list(NULL, NULL))
    }

    # ── 7. TMP via vectorized scatter into (LABEL×RUN) × FEATURE wide matrix ───
    #
    # Mirrors .fitTukey: dcast(LABEL + RUN ~ FEATURE, value.var="newABUNDANCE")
    # Row order in wide_mat: (lbl1,run1),(lbl1,run2),...,(lbl2,run1),...
    # where labels are sorted — matches dcast alphabetical ordering for
    # character LABEL values.
    if (FL_filt == 1L) {
        # ── Single (FEATURE, LABEL) row: skip TMP, use values directly ────────
        if (meta$is_labeled_ref) {
            h_sel <- as.character(flp_filt$LABEL) == "H"
            l_sel <- as.character(flp_filt$LABEL) == "L"
            if (any(h_sel) && any(l_sel)) {
                h_vals   <- as.vector(mat_newABU[h_sel, , drop = FALSE])
                l_vals   <- as.vector(mat_newABU[l_sel, , drop = FALSE])
                h_median <- stats::median(h_vals, na.rm = TRUE)
                result <- data.table::data.table(
                    Protein        = PROTEIN,
                    LABEL          = "L",
                    RUN            = runs_filt,
                    LogIntensities = l_vals - h_vals + h_median
                )
            } else {
                result <- data.table::data.table(
                    Protein        = PROTEIN,
                    LABEL          = as.character(flp_filt$LABEL[1L]),
                    RUN            = runs_filt,
                    LogIntensities = as.vector(mat_newABU)
                )
            }
        } else {
            result <- data.table::data.table(
                Protein        = PROTEIN,
                LABEL          = rep(as.character(flp_filt$LABEL), times = R_filt),
                RUN            = rep(runs_filt, each = FL_filt),
                LogIntensities = as.vector(mat_newABU)
            )
        }
    } else {
        # ── Multi-feature: scatter into wide matrix, apply TMP ─────────────────
        unique_labels <- sort(unique(as.character(flp_filt$LABEL)))
        unique_feats  <- sort(unique(as.character(flp_filt$FEATURE)))
        n_lbl  <- length(unique_labels)
        n_feat <- length(unique_feats)

        lbl_idx_map  <- match(as.character(flp_filt$LABEL),   unique_labels)
        feat_idx_map <- match(as.character(flp_filt$FEATURE),  unique_feats)

        # Vectorized scatter: mat[fl, run] → wide[(lbl_idx-1)*R + run, feat_idx]
        fl_rep   <- rep(seq_len(FL_filt), times = R_filt)
        run_rep  <- rep(seq_len(R_filt),  each  = FL_filt)
        wide_row <- (lbl_idx_map[fl_rep] - 1L) * R_filt + run_rep
        wide_col <- feat_idx_map[fl_rep]

        wide_mat <- matrix(NA_real_, nrow = n_lbl * R_filt, ncol = n_feat)
        wide_mat[cbind(wide_row, wide_col)] <- as.vector(mat_newABU)

        tmp_fitted <- median_polish_summary(wide_mat)

        # Row → (LABEL, RUN) mapping: label index cycles every R_filt rows
        result_labels <- rep(unique_labels, each = R_filt)
        result_runs   <- rep(runs_filt,     times = n_lbl)

        if (meta$is_labeled_ref) {
            h_li <- match("H", unique_labels)
            l_li <- match("L", unique_labels)
            if (!is.na(h_li) && !is.na(l_li)) {
                h_rows  <- (h_li - 1L) * R_filt + seq_len(R_filt)
                l_rows  <- (l_li - 1L) * R_filt + seq_len(R_filt)
                h_med   <- stats::median(tmp_fitted[h_rows], na.rm = TRUE)
                result  <- data.table::data.table(
                    Protein        = PROTEIN,
                    LABEL          = "L",
                    RUN            = runs_filt,
                    LogIntensities = tmp_fitted[l_rows] - tmp_fitted[h_rows] + h_med
                )
            } else {
                result <- data.table::data.table(
                    Protein        = PROTEIN,
                    LABEL          = result_labels,
                    RUN            = result_runs,
                    LogIntensities = tmp_fitted
                )
            }
        } else {
            result <- data.table::data.table(
                Protein        = PROTEIN,
                LABEL          = result_labels,
                RUN            = result_runs,
                LogIntensities = tmp_fitted
            )
        }
    }

    list(result, survival)
}


#' Feature-level data summarization via matrix-native TMP (V4)
#'
#' Improves on \code{MSstatsSummarizeWithMultipleCoresV3} by eliminating
#' the long-format \code{data.table} reconstruction and the \code{dcast}
#' inside Tukey median polish fitting.  Workers extract the FL×R matrix
#' from the V3 packed double vector and call
#' \code{MSstatsSummarizeSingleTMPV2} directly.
#'
#' For the \code{method = "linear"} case this function delegates to
#' \code{MSstatsSummarizeWithMultipleCoresV3} transparently (V4 only
#' implements the TMP path).
#'
#' @inheritParams MSstatsSummarizeWithMultipleCoresV2
#'
#' @return A named list with one element per protein slot, identical in
#'   structure to \code{MSstatsSummarizeWithMultipleCores}.
#'
#' @importFrom matter matter_list mem_realized SnowfastParam chunkLapply
#' @importFrom data.table data.table fifelse
#' @importFrom stats median predict
#'
#' @export
MSstatsSummarizeWithMultipleCoresV4 <- function(
    input,
    method,
    impute,
    censored_symbol,
    remove50missing,
    equal_variance,
    numberOfCores       = 1L,
    aft_iterations      = 90L,
    backing_path        = NULL,
    result_path         = NULL,
    ram_overhead_factor = 4,
    verbose             = FALSE
) {
    # ── 0. Fallbacks ───────────────────────────────────────────────────────────
    if (numberOfCores <= 1L) {
        return(MSstatsSummarizeWithSingleCore(
            input, method, impute, censored_symbol,
            remove50missing, equal_variance, aft_iterations))
    }
    if (!identical(method, "TMP")) {
        # V4 implements TMP only; delegate linear to V3
        return(MSstatsSummarizeWithMultipleCoresV3(
            input, method, impute, censored_symbol,
            remove50missing, equal_variance, numberOfCores,
            aft_iterations, backing_path, result_path,
            ram_overhead_factor, verbose))
    }

    # ── 1. Split and build matter backing store (identical to V3) ─────────────
    is_labeled_reference <- "is_labeled_ref" %in% colnames(input) &&
        any(input$is_labeled_ref, na.rm = TRUE)
    split_keys <- if (is_labeled_reference) list(input$PROTEIN) else
        list(input$PROTEIN, input$LABEL)
    protein_indices <- split(seq_len(nrow(input)), split_keys)
    protein_ids     <- names(protein_indices)
    num_proteins    <- length(protein_indices)

    all_runs <- if (is.factor(input$RUN)) levels(input$RUN) else
        sort(unique(as.character(input$RUN)))
    R <- length(all_runs)

    getOption("MSstatsLog")("INFO",
        paste0("V4: building matter matrices for ", num_proteins,
               " proteins × ", R, " runs"))

    if (is.null(backing_path)) backing_path <- tempfile(fileext = ".bin")

    slots <- vector("list", num_proteins)
    for (k in seq_len(num_proteins)) {
        slots[[k]] <- .buildProteinSlotV3(
            dt       = input[protein_indices[[k]], ],
            slot_k   = k,
            all_runs = all_runs)
    }
    packed_list <- lapply(slots, `[[`, "packed")
    meta_list   <- lapply(slots, `[[`, "meta")
    rm(slots)

    protein_matter <- matter::matter_list(
        packed_list, type = "double",
        path = backing_path, names = protein_ids)
    rm(packed_list)
    gc(verbose = FALSE)

    getOption("MSstatsLog")("INFO",
        paste0("V4: backing store written to ", backing_path,
               " (", format(matter::mem_realized(protein_matter)), ")"))

    # ── 2. Adaptive chunk sizing ──────────────────────────────────────────────
    per_protein_bytes <- as.numeric(matter::mem_realized(protein_matter)) /
        num_proteins
    effective_bytes   <- per_protein_bytes * ram_overhead_factor
    gc_stats <- gc(reset = FALSE)
    used_bytes <- sum(gc_stats[, "used"]) * 8
    worker_budget_bytes <- max(256 * 1024^2, used_bytes * 0.5)
    n_concurrent <- max(1L,
        floor(worker_budget_bytes / (numberOfCores * effective_bytes)))

    getOption("MSstatsLog")("INFO",
        paste0("V4: chunk size = ", n_concurrent,
               " proteins/worker (estimated ",
               format(round(effective_bytes / 1024^2, 2)), " MB per protein)"))

    # ── 3. Parallel backend ───────────────────────────────────────────────────
    BPPARAM <- matter::SnowfastParam(
        workers       = numberOfCores,
        force.GC      = TRUE,
        stop.on.error = FALSE)

    # ── 4. Worker closure ─────────────────────────────────────────────────────
    #
    # Captures meta_list and MSstatsSummarizeSingleTMPV2.  Workers receive the
    # closure once; subsequent proteins in the same chunk reuse it.
    # MSstatsSummarizeSingleTMPV2 is defined in the MSstats namespace, so its
    # internal calls (.fitSurvival, median_polish_summary) resolve correctly
    # after library(MSstats) in the worker.
    .processProteinV4 <- local({
        meta_list_       <- meta_list
        summarize_v2_    <- MSstatsSummarizeSingleTMPV2
        impute_          <- impute
        censored_symbol_ <- censored_symbol
        remove50missing_ <- remove50missing
        aft_iterations_  <- aft_iterations

        function(packed) {
            library(MSstats, quietly = TRUE, warn.conflicts = FALSE)
            k <- as.integer(packed[1L])
            summarize_v2_(packed, meta_list_[[k]],
                          impute_, censored_symbol_,
                          remove50missing_, aft_iterations_)
        }
    })

    # ── 5. Dispatch ───────────────────────────────────────────────────────────
    if (!is.null(result_path)) {
        .streaming <- local({
            inner_ <- .processProteinV4
            function(packed) serialize(inner_(packed), connection = NULL)
        })
        result_matter <- matter::chunkLapply(
            protein_matter, FUN = .streaming,
            outpath   = result_path,
            chunkopts = list(chunksize = n_concurrent),
            verbose   = verbose,
            BPPARAM   = BPPARAM)
        summarized_results <- lapply(seq_len(num_proteins), function(i) {
            unserialize(result_matter[[i]])
        })
        names(summarized_results) <- protein_ids
    } else {
        summarized_results <- matter::chunkLapply(
            protein_matter, FUN = .processProteinV4,
            chunkopts = list(chunksize = n_concurrent),
            verbose   = verbose,
            BPPARAM   = BPPARAM)
        names(summarized_results) <- protein_ids
    }

    getOption("MSstatsLog")("INFO", "V4: summarization complete.")
    summarized_results
}


#' Feature-level data summarization via socket-dispatched packed double vectors (V5)
#'
#' Replaces the matter backing-file approach of V3/V4 with direct socket
#' transfer of compact per-protein double vectors.  Each protein's numeric
#' columns are packed into a flat \code{double} vector by
#' \code{.buildProteinSlotV3} in the main process; \code{bplapply} sends
#' exactly one protein's vector to each worker task through the socket channel.
#'
#' Compared with earlier versions:
#' \itemize{
#'   \item No temporary backing file — works identically on LUSTRE, NFS, or any
#'     filesystem because data never touches the disk during parallel dispatch.
#'   \item Smaller per-task payload than the original \code{parLapply} approach
#'     (packed doubles only — no factor levels, character columns, or R object
#'     framing).
#'   \item \code{SnowfastParam} provides GC-between-tasks and dynamic port
#'     selection over plain \code{makeCluster}.
#'   \item For \code{method = "TMP"}: workers call
#'     \code{MSstatsSummarizeSingleTMPV2} (no \code{dcast} round-trip).
#'   \item For \code{method = "linear"}: workers reconstruct the data.table
#'     via \code{.reconstructProteinDTV3} and call
#'     \code{MSstatsSummarizeSingleLinear}.
#' }
#'
#' @inheritParams MSstatsSummarizeWithMultipleCoresV2
#'
#' @return A named list with one element per protein slot, identical in
#'   structure to \code{MSstatsSummarizeWithMultipleCores}.
#'
#' @importFrom matter SnowfastParam
#' @importFrom BiocParallel bplapply
#' @importFrom data.table data.table fifelse
#' @importFrom stats median predict
#'
#' @export
MSstatsSummarizeWithMultipleCoresV5 <- function(
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
    track_memory   = FALSE
) {
    # ── 0. Single-core fallback ────────────────────────────────────────────────
    if (numberOfCores <= 1L && is.null(BPPARAM)) {
        return(MSstatsSummarizeWithSingleCore(
            input, method, impute, censored_symbol,
            remove50missing, equal_variance, aft_iterations))
    }

    t_start <- proc.time()[["elapsed"]]
    .snap   <- if (track_memory) .memMB else function() NA_real_
    mem_log <- list()
    mem_log[["baseline (main)"]] <- .snap()

    # ── 1. Split input by protein slot ────────────────────────────────────────
    is_labeled_reference <- "is_labeled_ref" %in% colnames(input) &&
        any(input$is_labeled_ref, na.rm = TRUE)
    split_keys <- if (is_labeled_reference) list(input$PROTEIN) else
        list(input$PROTEIN, input$LABEL)
    protein_indices <- split(seq_len(nrow(input)), split_keys)
    protein_ids     <- names(protein_indices)
    num_proteins    <- length(protein_indices)
    mem_log[["after protein split"]] <- .snap()

    all_runs <- if (is.factor(input$RUN)) levels(input$RUN) else
        sort(unique(as.character(input$RUN)))

    getOption("MSstatsLog")("INFO",
        paste0("V5: packing ", num_proteins, " proteins × ",
               length(all_runs), " runs into double vectors"))

    # ── 2. Pack each protein into a compact double vector ─────────────────────
    packed_list <- vector("list", num_proteins)
    meta_list   <- vector("list", num_proteins)
    for (k in seq_len(num_proteins)) {
        slot             <- .buildProteinSlotV3(
                                input[protein_indices[[k]], ],
                                k, all_runs)
        packed_list[[k]] <- slot$packed
        meta_list[[k]]   <- slot$meta
    }
    mem_log[["after packing (packed_list built)"]] <- .snap()

    payload_mb <- sum(lengths(packed_list)) * 8 / 1024^2
    getOption("MSstatsLog")("INFO",
        paste0("V5: dispatching via sockets (",
               format(round(payload_mb, 1)), " MB total payload)"))

    # ── 3. Worker closure ─────────────────────────────────────────────────────
    use_TMP <- identical(method, "TMP")

    .worker <- local({
        meta_list_       <- meta_list
        reconstruct_     <- .reconstructProteinDTV3
        use_TMP_         <- use_TMP
        impute_          <- impute
        censored_symbol_ <- censored_symbol
        remove50missing_ <- remove50missing
        aft_iterations_  <- aft_iterations
        equal_variance_  <- equal_variance
        track_memory_    <- track_memory

        function(packed) {
            library(MSstats, quietly = TRUE, warn.conflicts = FALSE)
            k <- as.integer(packed[1L])
            result <- if (use_TMP_) {
                MSstatsSummarizeSingleTMPV2(
                    packed, meta_list_[[k]],
                    impute_, censored_symbol_,
                    remove50missing_, aft_iterations_)
            } else {
                dt <- reconstruct_(packed, meta_list_[[k]])
                MSstatsSummarizeSingleLinear(
                    dt, impute_, censored_symbol_,
                    remove50missing_, aft_iterations_,
                    equal_variances = equal_variance_)
            }
            if (track_memory_)
                list(.r = result, .m = MSstats:::.memMB())
            else
                result
        }
    })

    # ── 4. Dispatch ───────────────────────────────────────────────────────────
    if (is.null(BPPARAM))
        BPPARAM <- matter::SnowfastParam(
            workers       = numberOfCores,
            force.GC      = TRUE,
            stop.on.error = FALSE)

    results <- BiocParallel::bplapply(packed_list, .worker, BPPARAM = BPPARAM)
    mem_log[["after bplapply"]] <- .snap()
    names(results) <- protein_ids

    worker_mems <- NULL
    if (track_memory) {
        worker_mems <- vapply(results, function(x)
            if (is.list(x) && !is.null(x$.m)) x$.m else NA_real_, numeric(1L))
        results <- lapply(results, function(x)
            if (is.list(x) && !is.null(x$.r)) x$.r else x)
        names(results) <- protein_ids
        .printMemReport(
            "MSstatsSummarizeWithMultipleCoresV5",
            mem_log, worker_mems,
            elapsed = proc.time()[["elapsed"]] - t_start)
    }

    getOption("MSstatsLog")("INFO", "V5: summarization complete.")
    results
}
