# Tests that .fitSurvivalCG() - the conjugate-gradient alternative to
# .fitSurvival() - fits the same model, and agrees numerically with it.
#
# The scenarios below reuse the noiseless fixtures from
# test_utils_imputation.R purely to check that .fitSurvivalCG() selects the
# same predictors as .fitSurvival() (via the shared .buildAFTFormula()).
# For numeric agreement on the fitted values themselves, a Gaussian AFT
# model needs actual residual variation to estimate - a noiseless design
# has a degenerate (unbounded) scale MLE, so a second set of fixtures below
# adds realistic noise and left-censoring before comparing coefficients,
# scale, and predictions.

make_surv_labeled_single <- function() {
    runs <- paste0("R", 1:3)
    dt <- data.table::rbindlist(list(
        data.table::data.table(
            FEATURE      = factor(rep("F1", 9)),
            RUN          = factor(rep(runs, each = 3)),
            LABEL        = "H",
            newABUNDANCE = seq(10.1, by = 0.1, length.out = 9),
            cen          = 1L
        ),
        data.table::data.table(
            FEATURE      = factor(rep("F1", 9)),
            RUN          = factor(rep(runs, each = 3)),
            LABEL        = "L",
            newABUNDANCE = seq(14.1, by = 0.1, length.out = 9),
            cen          = 1L
        )
    ))
    ref_vals <- ifelse(dt$LABEL == "L", as.character(dt$RUN), "0")
    dt[["ref_covariate"]] <- factor(ref_vals, levels = c("0", runs))
    dt
}

make_surv_unlabeled_multi_welldetermined <- function() {
    dt <- data.table::CJ(
        FEATURE = paste0("F", 1:3),
        RUN     = paste0("R", 1:5)
    )
    dt[, FEATURE := factor(FEATURE)]
    dt[, RUN := factor(RUN)]
    dt[, LABEL := "L"]
    dt[, newABUNDANCE := seq(10, by = 0.5, length.out = .N)]
    dt[, cen := 1L]
    dt
}

# --- .fitSurvivalCG() selects the same predictors as .fitSurvival() -------

coef_names <- function(fit) names(coef(fit))

for (make_input in list(make_surv_labeled_single,
                        make_surv_unlabeled_multi_welldetermined)) {
    input <- make_input()
    chol_names <- sort(coef_names(MSstats:::.fitSurvival(input, 90)))
    cg_names <- sort(coef_names(MSstats:::.fitSurvivalCG(input, 90)))
    expect_equal(
        cg_names, chol_names,
        info = ".fitSurvivalCG must select the same predictors as .fitSurvival"
    )
}

# --- numeric agreement on realistic (noisy, censored) data ----------------

make_noisy_censored_input <- function(seed, is_labeled) {
    set.seed(seed)
    features <- paste0("F", 1:3)
    runs <- paste0("R", 1:4)
    labels <- if (is_labeled) c("H", "L") else "L"
    dt <- data.table::CJ(FEATURE = features, RUN = runs, LABEL = labels)
    dt[, FEATURE := factor(FEATURE)]
    dt[, RUN := factor(RUN)]
    dt[, newABUNDANCE :=
           10 + as.integer(FEATURE) + as.integer(RUN) * 0.5 +
           ifelse(LABEL == "L", 4, 0) + rnorm(.N, sd = 0.7)]
    dt[, cen := 1L]
    censoring_threshold <- stats::quantile(dt$newABUNDANCE, 0.2)
    dt[newABUNDANCE < censoring_threshold, cen := 0L]
    dt[cen == 0L, newABUNDANCE := censoring_threshold]
    if (is_labeled) {
        ref_vals <- ifelse(dt$LABEL == "L", as.character(dt$RUN), "0")
        dt[["ref_covariate"]] <- factor(ref_vals, levels = c("0", runs))
    }
    dt
}

check_solvers_agree <- function(input, tolerance, label) {
    fit_cholesky <- MSstats:::.fitSurvival(input, 90)
    fit_cg <- MSstats:::.fitSurvivalCG(input, 90)

    matched_names <- names(fit_cholesky$coefficients)
    expect_equal(
        fit_cg$coefficients[matched_names],
        fit_cholesky$coefficients,
        tolerance = tolerance, check.attributes = FALSE,
        info = paste(label, ": coefficients should match .fitSurvival")
    )
    expect_equal(
        fit_cg$scale, fit_cholesky$scale, tolerance = tolerance,
        check.attributes = FALSE,
        info = paste(label, ": scale should match .fitSurvival")
    )

    predictions_cholesky <- predict(fit_cholesky, newdata = input, se.fit = TRUE)
    predictions_cg <- predict(fit_cg, newdata = input, se.fit = TRUE)
    expect_equal(
        predictions_cg$fit, predictions_cholesky$fit,
        tolerance = tolerance, check.attributes = FALSE,
        info = paste(label, ": predicted values should match .fitSurvival")
    )
    expect_equal(
        predictions_cg$se.fit, predictions_cholesky$se.fit,
        tolerance = tolerance, check.attributes = FALSE,
        info = paste(label, ": prediction standard errors should match",
                     ".fitSurvival")
    )
}

check_solvers_agree(
    make_noisy_censored_input(seed = 1, is_labeled = TRUE),
    tolerance = 1e-4, label = "labeled, noisy, censored"
)
check_solvers_agree(
    make_noisy_censored_input(seed = 2, is_labeled = FALSE),
    tolerance = 1e-4, label = "unlabeled, noisy, censored"
)
