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

# Cross-platform peak-RSS reader. Reflects the true lifetime peak of the
# calling process, regardless of when you call it — no polling required.
.peakRSS_MB <- function() {
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
# checkpoints: named numeric vector of RSS snapshots (MB).
# worker_mems: numeric vector of per-worker peak RSS values (may be NA).
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





## ── V4: matrix-native TMP ─────────────────────────────────────────────────────
##
## .MSstatsSummarizeSingleTMPV2 takes the V3 packed double vector + meta list
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
.MSstatsSummarizeSingleTMPV2 <- function(packed, meta,
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
            MSstats:::.fitSurvival(fit_data, aft_iterations)
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
            predicted    = NA  # logical NA, matching MSstatsSummarizeSingleTMP's
            # bare `survival[, predicted := NA]` in its
            # equivalent non-imputed branch (dataProcess.R:586)
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
        
        tmp_fitted <- MSstats:::median_polish_summary(wide_mat)
        
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

#' Build the per-record worker closure for
#' \code{MSstatsSummarizeWithMultipleCoresV6}
#'
#' Defined at package top level — not nested via \code{local()} inside
#' \code{MSstatsSummarizeWithMultipleCoresV6} — so the returned closure's
#' enclosing environment chain is this factory's own (small) call frame plus
#' the package namespace. \code{BiocParallel}/\code{matter} serialize a
#' closure's entire enclosing environment chain to ship it to each socket
#' worker, not just the variables the closure body actually references. A
#' closure built via \code{local()} inside
#' \code{MSstatsSummarizeWithMultipleCoresV6} would have that function's own
#' evaluation frame in its chain — which holds \code{input},
#' \code{protein_records}, and other run-scale objects — so every worker
#' would receive and retain a serialized copy of them even though
#' \code{.worker} never touches them. Building the closure here instead keeps
#' its captured state to only the scalar run parameters.
#'
#' @keywords internal
.buildSummarizeWorkerV6 <- function(
        use_TMP, impute, censored_symbol, remove50missing,
        aft_iterations, equal_variance
) {
    reconstruct_     <- .reconstructProteinDTV3
    use_TMP_         <- use_TMP
    impute_          <- impute
    censored_symbol_ <- censored_symbol
    remove50missing_ <- remove50missing
    aft_iterations_  <- aft_iterations
    equal_variance_  <- equal_variance
    
    function(record) {
        packed <- record$packed
        meta   <- record$meta
        result <- if (use_TMP_) {
            .MSstatsSummarizeSingleTMPV2(
                packed, meta,
                impute_, censored_symbol_,
                remove50missing_, aft_iterations_)
        } else {
            dt <- reconstruct_(packed, meta)
            MSstatsSummarizeSingleLinear(
                dt, impute_, censored_symbol_,
                remove50missing_, aft_iterations_,
                equal_variances = equal_variance_)
        }
        # Normalize column types/levels at the V6 worker boundary rather than
        # inside .MSstatsSummarizeSingleTMPV2/MSstatsSummarizeSingleLinear:
        # those are shared with V3-V5 and dataProcess.R respectively, and the
        # TMP path in particular carries RUN as plain character (from
        # meta$runs) and cen as double (from the packed matrix). Re-leveling
        # against meta$runs (rather than a bare factor(RUN)) preserves the
        # numeric run order that levels(input$RUN) already established
        # upstream — a bare factor() call would instead re-sort alphabetically
        # ("1","10","11","2",... for >=10 runs). as.character() first strips
        # any existing level order before re-leveling, since .MSstatsSummarize
        # SingleTMPV2's survival table (result[[2L]]) already wraps RUN in a
        # bare factor() for the imputed branch — that call has the exact same
        # alphabetical-resort bug even though it looks pre-typed. droplevels()
        # then matches MSstatsSummarizeSingleTMP/SingleCore, whose factor(RUN)
        # only ever sees — and so only ever keeps — the runs present for that
        # protein, since meta$runs carries every run across the whole input.
        #
        # FEATURE gets the same treatment: .MSstatsSummarizeSingleTMPV2's
        # non-imputed survival branch leaves FEATURE as plain character (only
        # the imputed branch wraps it in a bare factor()), so result[[2L]]$
        # FEATURE isn't reliably a factor the way SingleCore's is (single_
        # protein[, FEATURE := factor(FEATURE)] before survival is sliced off
        # it).
        #
        # Deliberately NOT leveling against meta$feat_label_pep$FEATURE's
        # existing row order here: that order comes from .buildProteinSlotV3's
        # data.table::setorder(flp, FEATURE, LABEL), and data.table sorts
        # character columns in the C-locale (byte order: all uppercase before
        # any lowercase) for platform-independence — whereas SingleCore's
        # bare factor(FEATURE) goes through base R's factor()/sort(), which
        # use the session's collation locale (e.g. en_US.UTF-8, where case is
        # interleaved: "a" < "A" < "b" < "B"). Any FEATURE strings with mixed
        # case (e.g. modification tags) sort differently between the two, so
        # re-deriving the level order with a base R sort() reproduces
        # SingleCore's order instead of inheriting data.table's.
        feature_levels <- sort(unique(meta$feat_label_pep$FEATURE))
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

#' Per-worker peak-RSS query task for \code{MSstatsSummarizeWithMultipleCoresV6}
#'
#' Dispatched once per worker (via \code{seq_len(bpnworkers(BPPARAM))}) after
#' the main summarization \code{bplapply} call, while the persistent workers
#' are still alive, so each worker reports its own true lifetime-peak RSS
#' rather than a snapshot taken mid-run.
#'
#' @keywords internal
#' @noRd
.reportWorkerPeakV6 <- function(i) {
    list(worker = i, pid = Sys.getpid(), peak_mb = .peakRSS_MB())
}

# Per-worker warm-up task for MSstatsSummarizeWithMultipleCoresV6: loads
# MSstats once per persistent worker process and pins data.table to a single
# thread. Defined at top level rather than inline inside
# MSstatsSummarizeWithMultipleCoresV6: even though its body never references
# `input`/`protein_records`, a closure defined inline there would still carry
# that function's evaluation frame in its enclosing environment chain, and
# BiocParallel would serialize that whole frame — including those run-scale
# objects — to every worker just to ship this no-op task.
.warmupV6Worker <- function(i) {
    library(MSstats, quietly = TRUE, warn.conflicts = FALSE)
    data.table::setDTthreads(1)
    NULL
}


#' Feature-level data summarization via socket-dispatched protein records (V6)
#'
#' Fixes a hidden RAM cost in \code{MSstatsSummarizeWithMultipleCoresV5}: V5's
#' worker closure captured \code{meta_list} — the metadata for \emph{every}
#' protein — by reference, so each worker received and retained the full
#' metadata set for the whole run regardless of how many proteins were
#' actually assigned to it. Bounding the per-task \emph{packed-vector} payload
#' (e.g. via \code{SnowfastParam(tasks = ...)}) did nothing to reduce that,
#' since the metadata was baked into the closure once, not re-sliced per task.
#'
#' V6 pairs each protein's packed double vector with only its own metadata
#' into a single list element (\code{list(packed = ..., meta = ...)}) and
#' dispatches that combined list through \code{bplapply}. \code{BiocParallel}
#' then only serializes and sends the records actually assigned to a given
#' task, so a worker never holds metadata for proteins outside its own
#' task(s) — unlike V5, RAM scales with the batch actually being processed.
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
#' @return A named list with one element per protein slot, identical in
#'   structure to \code{MSstatsSummarizeWithMultipleCores}.
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
    
    t_start <- proc.time()[["elapsed"]]
    mem_log <- list()
    
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
                            paste0("V6: packing ", num_proteins, " proteins × ",
                                   length(all_runs), " runs into per-protein records"))
    
    # ── 2. Pack each protein into a (packed vector, own metadata) record ──────
    #
    # Unlike V5's parallel packed_list/meta_list arrays, each protein's
    # metadata travels bundled with its own packed vector. bplapply/BPPARAM
    # slice this combined list into tasks, so a worker's closure never needs
    # — and never receives — metadata for proteins outside its assigned
    # task(s).
    protein_records <- vector("list", num_proteins)
    for (k in seq_len(num_proteins)) {
        slot <- .buildProteinSlotV3(
            input[protein_indices[[k]], ], k, all_runs)
        protein_records[[k]] <- list(packed = slot$packed, meta = slot$meta)
    }
    
    payload_mb <- sum(vapply(protein_records,
                             function(r) length(r$packed), integer(1))) * 8 / 1024^2
    getOption("MSstatsLog")("INFO",
                            paste0("V6: dispatching via sockets (",
                                   format(round(payload_mb, 1)),
                                   " MB total packed payload; metadata sharded per task)"))
    
    # ── 3. Worker closure ─────────────────────────────────────────────────────
    #
    # No captured meta_list: each task's records already carry their own
    # metadata, so the closure only needs the scalar run parameters.
    #
    # Built via .buildSummarizeWorkerV6() (defined at package top level)
    # rather than local() here, so the closure's enclosing environment chain
    # never includes this function's own frame — which holds `input`,
    # `protein_records`, etc. — and BiocParallel doesn't serialize those
    # run-scale objects to every worker.
    use_TMP <- identical(method, "TMP")
    
    .worker <- .buildSummarizeWorkerV6(
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
                                paste0("V6: dispatching as ",
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
        .warmupV6Worker, BPPARAM = BPPARAM)
    BiocParallel::bpprogressbar(BPPARAM) <- show_progress
    
    results <- BiocParallel::bplapply(protein_records, .worker, BPPARAM = BPPARAM)
    names(results) <- protein_ids
    
    worker_peaks <- NULL
    if (track_memory) {
        BiocParallel::bpprogressbar(BPPARAM) <- FALSE
        worker_peaks <- BiocParallel::bplapply(
            seq_len(BiocParallel::bpnworkers(BPPARAM)),
            .reportWorkerPeakV6, BPPARAM = BPPARAM)   # must run BEFORE bpstop() while workers are alive
        BiocParallel::bpprogressbar(BPPARAM) <- show_progress
        mem_log[["parent peak (main)"]] <- .peakRSS_MB()
        worker_peak_mb <- vapply(worker_peaks, function(x) x$peak_mb, numeric(1L))
        .printMemReport(
            "MSstatsSummarizeWithMultipleCoresV6",
            mem_log, worker_peak_mb,
            elapsed = proc.time()[["elapsed"]] - t_start)
    }
    
    getOption("MSstatsLog")("INFO", "V6: summarization complete.")
    results
}
