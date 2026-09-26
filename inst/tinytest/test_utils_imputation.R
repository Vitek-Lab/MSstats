add_ref_covariate <- function(dt) {
    dt[, ref_covariate := factor(ifelse(LABEL == "L", as.character(RUN), "0"),
                                 levels = c("0", levels(RUN)))]
}

make_surv_input <- function(n_features, n_runs, is_labeled, n_reps = 1L,
                            noise_sd = 0.1, censored_fraction = 0, seed = 1) {
    set.seed(seed)
    dt <- data.table::CJ(
        FEATURE   = paste0("F", seq_len(n_features)),
        RUN       = paste0("R", seq_len(n_runs)),
        LABEL     = if (is_labeled) c("H", "L") else "L",
        REPLICATE = seq_len(n_reps)
    )
    dt[, FEATURE := factor(FEATURE)]
    dt[, RUN := factor(RUN)]
    dt[, newABUNDANCE := 10 + as.integer(FEATURE) + as.integer(RUN) * 0.5 +
           ifelse(LABEL == "L", 4, 0) + (REPLICATE - 1) * 0.1 +
           rnorm(.N, sd = noise_sd)]
    dt[, REPLICATE := NULL]
    dt[, cen := 1L]
    if (censored_fraction > 0) {
        censoring_threshold <- stats::quantile(dt$newABUNDANCE, censored_fraction)
        dt[newABUNDANCE < censoring_threshold,
           `:=`(cen = 0L, newABUNDANCE = censoring_threshold)]
    }
    if (is_labeled) add_ref_covariate(dt)
    dt
}

make_surv_labeled_underdetermined <- function() {
    dt <- data.table::data.table(
        FEATURE      = factor(c(paste0("F", 1:8), "F1")),
        RUN          = factor(c(rep_len(paste0("R", 1:3), 8), "R1")),
        LABEL        = c(rep("H", 8), "L"),
        newABUNDANCE = c(seq(10, by = 0.5, length.out = 8), 14),
        cen          = 1L
    )
    add_ref_covariate(dt)
}

coef_names <- function(fit) names(coef(fit))

predictor_cases <- list(
    "labeled single-feature" = list(
        input = make_surv_input(1, 3, TRUE, n_reps = 3),
        has_ref_covariate = TRUE, has_feature = FALSE),
    "labeled multi well-determined" = list(
        input = make_surv_input(3, 4, TRUE),
        has_ref_covariate = TRUE, has_feature = TRUE),
    "labeled underdetermined" = list(
        input = make_surv_labeled_underdetermined(),
        has_ref_covariate = TRUE, has_feature = FALSE),
    "unlabeled single-feature" = list(
        input = make_surv_input(1, 5, FALSE, n_reps = 3),
        has_ref_covariate = FALSE, has_feature = FALSE),
    "unlabeled multi well-determined" = list(
        input = make_surv_input(3, 5, FALSE),
        has_ref_covariate = FALSE, has_feature = TRUE)
)

for (case_name in names(predictor_cases)) {
    case <- predictor_cases[[case_name]]
    chol_names <- coef_names(MSstats:::.fitSurvival(case$input, 90))
    expect_equal(
        any(grepl("ref_covariate", chol_names)), case$has_ref_covariate,
        info = paste(".fitSurvival", case_name, ": ref_covariate in coefficients")
    )
    expect_equal(
        any(grepl("^FEATURE", chol_names)), case$has_feature,
        info = paste(".fitSurvival", case_name, ": FEATURE in coefficients")
    )
    expect_equal(
        sort(coef_names(MSstats:::.fitSurvivalCG(case$input, 90))),
        sort(chol_names),
        info = paste(case_name, ": .fitSurvivalCG must select the same",
                     "predictors as .fitSurvival")
    )
}

check_solvers_agree <- function(input, label, use_jacobi_preconditioner) {
    fit_cholesky <- MSstats:::.fitSurvival(input, 90)
    fit_cg <- MSstats:::.fitSurvivalCG(
        input, 90, use_jacobi_preconditioner = use_jacobi_preconditioner)
    summarize_fit <- function(fit) {
        predictions <- predict(fit, newdata = input, se.fit = TRUE)
        list(
            coefficients = fit$coefficients[names(fit_cholesky$coefficients)],
            scale = fit$scale,
            `predicted values` = predictions$fit,
            `prediction standard errors` = predictions$se.fit
        )
    }
    expected <- summarize_fit(fit_cholesky)
    actual <- summarize_fit(fit_cg)
    for (quantity in names(expected)) {
        expect_equal(
            actual[[quantity]], expected[[quantity]],
            tolerance = 1e-4, check.attributes = FALSE,
            info = paste(label, ":", quantity, "should match .fitSurvival")
        )
    }
}

noisy_inputs <- list(
    labeled = make_surv_input(3, 4, TRUE, noise_sd = 0.7,
                              censored_fraction = 0.2, seed = 1),
    unlabeled = make_surv_input(3, 4, FALSE, noise_sd = 0.7,
                                censored_fraction = 0.2, seed = 2)
)
for (input_name in names(noisy_inputs)) {
    for (use_jacobi in c(FALSE, TRUE)) {
        check_solvers_agree(
            noisy_inputs[[input_name]],
            label = paste0(input_name, ", noisy, censored",
                           if (use_jacobi) ", jacobi-preconditioned"),
            use_jacobi_preconditioner = use_jacobi
        )
    }
}

noisy_input <- make_surv_input(3, 4, FALSE, noise_sd = 0.7,
                               censored_fraction = 0.2, seed = 3)

fit_cholesky <- MSstats:::.fitAFTModel(noisy_input, 90, "cholesky")
expect_inherits(
    fit_cholesky, "survreg",
    info = ".fitAFTModel(aft_solver = 'cholesky') should return a survreg fit"
)
expect_true(
    is.null(fit_cholesky$cg_diagnostics),
    info = "the cholesky path should not attach cg_diagnostics"
)
for (aft_solver in c("cg", "pcg")) {
    expect_false(
        is.null(MSstats:::.fitAFTModel(noisy_input, 90, aft_solver)$cg_diagnostics),
        info = paste0(".fitAFTModel(aft_solver = '", aft_solver,
                      "') should attach cg_diagnostics")
    )
}

expect_silent(
    MSstats:::.fitSurvivalCG(noisy_input, 90, verbose = FALSE)
)
expect_message(
    MSstats:::.fitSurvivalCG(noisy_input, 90, verbose = TRUE),
    pattern = "\\[AFT-CG\\] starting fit",
    info = "verbose = TRUE should report the problem size at the start of the fit"
)
expect_message(
    MSstats:::.fitSurvivalCG(noisy_input, 90, verbose = TRUE),
    pattern = "\\[AFT-CG\\] finished",
    info = "verbose = TRUE should report a summary once fitting finishes"
)

fit_with_diagnostics <- MSstats:::.fitSurvivalCG(noisy_input, 90)
expect_equal(
    nrow(fit_with_diagnostics$cg_diagnostics), fit_with_diagnostics$iter,
    info = "cg_diagnostics should have one row per Newton-Raphson iteration taken"
)
expect_true(
    all(fit_with_diagnostics$cg_diagnostics$cg_iterations >= 0),
    info = "cg_iterations should be a non-negative count for every Newton iteration"
)
