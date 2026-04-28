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
