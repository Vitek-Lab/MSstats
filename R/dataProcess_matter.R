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
