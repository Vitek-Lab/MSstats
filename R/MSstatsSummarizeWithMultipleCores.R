.peak_rss_mb <- function() {
  if (.Platform$OS.type != "windows" && file.exists("/proc/self/status")) {
    ln <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
    if (length(ln)) return(as.numeric(sub("\\D+(\\d+).*", "\\1", ln)) / 1024)
  }
  peak_rss_mb()
}

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

#' Pack one protein slot into a double vector plus metadata
#'
#' @param protein_dt data.table rows for one protein (or protein x label) slot
#' @param slot_index integer position of this slot in the global protein list
#' @param all_runs character vector of all run names in global order
#' @return list with elements \code{packed} (double vector) and \code{meta} (list)
#' @keywords internal
.pack_protein_slot <- function(protein_dt, slot_index, all_runs) {

    n_runs <- length(all_runs)

    has_peptide <- "PEPTIDE" %in% colnames(protein_dt)
    feature_label_dt <- unique(protein_dt[, .(
        FEATURE = as.character(FEATURE),
        LABEL   = as.character(LABEL),
        PEPTIDE = if (has_peptide) as.character(PEPTIDE) else as.character(FEATURE)
    )])
    data.table::setorder(feature_label_dt, FEATURE, LABEL)
    n_feature_labels <- nrow(feature_label_dt)

    feature_row_idx <- match(
        paste(as.character(protein_dt$FEATURE), as.character(protein_dt$LABEL), sep = "\t"),
        paste(feature_label_dt$FEATURE, feature_label_dt$LABEL, sep = "\t"))
    run_col_idx_lookup <- seq_len(n_runs); names(run_col_idx_lookup) <- all_runs
    run_col_idx  <- run_col_idx_lookup[as.character(protein_dt$RUN)]
    row_is_valid <- !is.na(feature_row_idx) & !is.na(run_col_idx)

    scatter_into_matrix <- function(col_vec, fill = NA_real_) {
        m <- matrix(fill, nrow = n_feature_labels, ncol = n_runs)
        m[cbind(feature_row_idx[row_is_valid], run_col_idx[row_is_valid])] <-
            as.double(col_vec[row_is_valid])
        m
    }

    new_abundance_mat <- scatter_into_matrix(protein_dt$newABUNDANCE)

    has_ABUNDANCE <- "ABUNDANCE" %in% colnames(protein_dt)
    abundance_mat <- if (has_ABUNDANCE) scatter_into_matrix(protein_dt$ABUNDANCE) else
        matrix(NA_real_, n_feature_labels, n_runs)

    has_censored <- "censored" %in% colnames(protein_dt)
    censored_mat <- if (has_censored) scatter_into_matrix(as.double(protein_dt$censored)) else
        matrix(0.0, n_feature_labels, n_runs)

    has_cen <- "cen" %in% colnames(protein_dt)
    event_mat <- if (has_cen) scatter_into_matrix(protein_dt$cen) else
        matrix(NA_real_, n_feature_labels, n_runs)

    has_anom <- "ANOMALYSCORES" %in% colnames(protein_dt) &&
        !all(is.na(protein_dt$ANOMALYSCORES))
    anomaly_scores_mat <- if (has_anom) scatter_into_matrix(protein_dt$ANOMALYSCORES) else
        matrix(NA_real_, n_feature_labels, n_runs)

    n_obs_by_feature_label <- protein_dt[, .(n_obs = as.double(n_obs[1L])),
                      by = .(FEATURE = as.character(FEATURE),
                             LABEL   = as.character(LABEL))]
    feature_label_with_nobs <- n_obs_by_feature_label[feature_label_dt, on = c("FEATURE", "LABEL")]
    data.table::setorder(feature_label_with_nobs, FEATURE, LABEL)
    n_obs_vec   <- feature_label_with_nobs$n_obs

    run_level_scalars <- protein_dt[, .(n_obs_run    = as.double(n_obs_run[1L]),
                       prop_features = as.double(prop_features[1L])),
                   by = .(RUN = as.character(RUN))]
    run_scalars_all_runs <- run_level_scalars[data.table::data.table(RUN = all_runs), on = "RUN"]
    n_obs_run_vec     <- run_scalars_all_runs$n_obs_run
    prop_features_vec <- run_scalars_all_runs$prop_features

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

    meta <- list(
        PROTEIN           = as.character(protein_dt$PROTEIN[1L]),
        feature_label_dt  = as.data.frame(feature_label_dt),
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

    feature_label_dt <- meta$feature_label_dt
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

    if (meta$has_ABUNDANCE) {
        protein_dt[, ABUNDANCE := as.vector(abundance_mat)]
    }

    protein_dt[, censored := {
        v <- as.vector(censored_mat)
        if (meta$has_censored) as.logical(v > 0.5) else rep(FALSE, n_rows)
    }]

    if (meta$has_cen) {
        protein_dt[, cen := as.vector(event_mat)]
    }

    protein_dt[, ANOMALYSCORES := as.vector(anomaly_scores_mat)]

    if (meta$is_labeled_ref) {
        protein_dt[, is_labeled_ref := (LABEL == "H")]
        if (meta$add_ref_covariate) {
            protein_dt[, ref_covariate := factor(
                data.table::fifelse(LABEL == "L", as.character(RUN), "0"))]
        }
    }

    protein_dt
}


#' Build the per-record worker closure for \code{MSstatsSummarizeWithMultipleCores}
#'
#' Defined at package top level so the closure only captures the scalar run
#' parameters, not the caller's run-scale objects.
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
        result
    }
}

#' Build the per-bundle worker closure for \code{MSstatsSummarizeWithMultipleCores}
#'
#' Defined at package top level, like \code{.build_summarize_worker}, so the
#' closure captures only \code{worker_fn} -- not the caller's run-scale
#' objects (e.g. \code{input}). A closure defined inline inside
#' \code{MSstatsSummarizeWithMultipleCores} would instead capture that
#' function's entire call frame, and \code{bpiterate} would re-serialize all
#' of it (including \code{input}) every time \code{FUN} is shipped to a
#' worker.
#'
#' @keywords internal
.build_bundle_worker <- function(worker_fn) {
    function(bundle) lapply(bundle, worker_fn)
}

#' Per-worker peak-RSS query task for \code{MSstatsSummarizeWithMultipleCores}
#'
#' @keywords internal
#' @noRd
.report_worker_peak <- function(i) {
    list(worker = i, pid = Sys.getpid(), peak_mb = .peak_rss_mb())
}

.warmup_worker <- function(i) {
    library(MSstats, quietly = TRUE, warn.conflicts = FALSE)
    data.table::setDTthreads(1)
    NULL
}


#' Feature-level data summarization via socket-dispatched protein records
#'
#' @param input feature-level data processed by dataProcess subfunctions
#' @param method summarization method: "linear" or "TMP"
#' @param impute only for method = "TMP"; imputes censored values via AFT model when TRUE
#' @param censored_symbol how censored values are encoded: 'NA', '0', or NULL for none
#' @param remove50missing only for method = "TMP"; drops proteins missing >=50\% per peptide in every run
#' @param equal_variance only for method = "linear"; assume equal variance among feature intensities
#' @param numberOfCores number of cores for parallel processing (Linux/Mac only)
#' @param aft_iterations number of AFT model iterations
#' @param verbose whether to print verbose output
#' @param BPPARAM optional \code{BiocParallelParam} instance
#' @param track_memory whether to report per-worker peak RSS memory usage
#' @param max_proteins_per_worker how many protein records are packed into a
#'   single dispatched \code{bpiterate} job; a throughput knob (fewer, larger
#'   jobs cut per-job socket overhead) rather than a memory cap -- parent
#'   memory is bounded by \code{bpnworkers(BPPARAM)} in-flight jobs regardless
#'   of this value, since records are packed lazily as workers ask for more
#'   work. 0 (default) packs one protein per job.
#'
#' @return A named list with one element per protein slot, keyed by protein
#'   (or protein \eqn{\times} label) identifier.
#'
#' @importFrom matter SnowfastParam
#' @importFrom BiocParallel bplapply bpiterate bpstart bpstop bpisup bpnworkers bpprogressbar
#' @importFrom data.table data.table fifelse setDTthreads
#' @importFrom stats median
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
        max_proteins_per_worker = 50L
) {
    if (numberOfCores <= 1L && is.null(BPPARAM)) {
        return(MSstatsSummarizeWithSingleCore(
            input, method, impute, censored_symbol,
            remove50missing, equal_variance, aft_iterations))
    }

    start_time <- proc.time()[["elapsed"]]
    memory_checkpoints <- list()

    is_labeled_reference <- "is_labeled_ref" %in% colnames(input) &&
        any(input$is_labeled_ref, na.rm = TRUE)
    split_keys <- if (is_labeled_reference) list(input$PROTEIN) else
        list(input$PROTEIN, input$LABEL)
    protein_indices <- split(seq_len(nrow(input)), split_keys)
    protein_ids     <- names(protein_indices)
    num_proteins    <- length(protein_indices)

    all_runs <- if (is.factor(input$RUN)) levels(input$RUN) else
        as.character(sort(unique(input$RUN)))

    ## bpiterate() uses the same continuous, capacity-bounded dispatch loop as
    ## bplapply() (BiocParallel:::.bploop_impl): at most bpnworkers(BPPARAM)
    ## jobs are ever in flight, and a free worker is refilled the instant it
    ## finishes -- a slow protein only stalls the worker that drew it, not its
    ## siblings. Calling bplapply() repeatedly over fixed-size batches (the
    ## previous approach here) does not have that property: it introduces a
    ## barrier at every batch boundary where idle workers wait for the whole
    ## batch, including any straggler, before the next batch is dispatched.
    ## Feeding bpiterate() packed records lazily -- one bundle at a time,
    ## built only when a worker asks for more work -- also means the packed
    ## payload for all num_proteins proteins and the collected result list
    ## are never both fully resident in the parent at once, unlike calling
    ## bplapply(protein_records, ...) on an eagerly-packed input list (see
    ## BiocParallel:::.lapplyReducer, which pre-allocates a length(X)-sized
    ## result list up front regardless of task granularity).
    bundle_size <- if (max_proteins_per_worker > 0L) max_proteins_per_worker else 1L
    n_bundles   <- ceiling(num_proteins / bundle_size)

    getOption("MSstatsLog")("INFO",
                            paste0("Summarizing ", num_proteins, " proteins × ",
                                   length(all_runs), " runs via ", n_bundles,
                                   " job(s) of up to ", bundle_size,
                                   " protein(s) each"))

    use_TMP <- identical(method, "TMP")

    worker_fn <- .build_summarize_worker(
        use_TMP, impute, censored_symbol, remove50missing,
        aft_iterations, equal_variance)
    bundle_worker_fn <- .build_bundle_worker(worker_fn)

    if (is.null(BPPARAM)) {
        BPPARAM <- matter::SnowfastParam(
            workers       = numberOfCores,
            tasks         = 0L,
            progressbar   = TRUE,
            force.GC      = TRUE,
            stop.on.error = FALSE)
    }

    started_here <- !BiocParallel::bpisup(BPPARAM)
    if (started_here) {
        BiocParallel::bpstart(BPPARAM)
        on.exit(BiocParallel::bpstop(BPPARAM), add = TRUE)
    }
    show_progress <- BiocParallel::bpprogressbar(BPPARAM)

    BiocParallel::bpprogressbar(BPPARAM) <- FALSE
    BiocParallel::bplapply(
        seq_len(BiocParallel::bpnworkers(BPPARAM)),
        .warmup_worker, BPPARAM = BPPARAM)
    BiocParallel::bpprogressbar(BPPARAM) <- show_progress

    next_slot <- 0L
    pack_next_bundle <- function() {
        if (next_slot >= num_proteins) return(NULL)
        bundle_end <- min(next_slot + bundle_size, num_proteins)
        bundle <- vector("list", bundle_end - next_slot)
        for (i in seq_along(bundle)) {
            slot_index <- next_slot + i
            bundle[[i]] <- .pack_protein_slot(
                input[protein_indices[[slot_index]], ], slot_index, all_runs)
        }
        next_slot <<- bundle_end
        bundle
    }

    results_by_bundle <- BiocParallel::bpiterate(
        pack_next_bundle, bundle_worker_fn, BPPARAM = BPPARAM)
    results <- unlist(results_by_bundle, recursive = FALSE)
    names(results) <- protein_ids

    worker_peaks <- NULL
    if (track_memory) {
        BiocParallel::bpprogressbar(BPPARAM) <- FALSE
        worker_peaks <- BiocParallel::bplapply(
            seq_len(BiocParallel::bpnworkers(BPPARAM)),
            .report_worker_peak, BPPARAM = BPPARAM)
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
