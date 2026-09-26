number_of_iterations <- 90

add_reference_covariate <- function(input) {
    input[, ref_covariate := factor(
        ifelse(LABEL == "L", as.character(RUN), "0"),
        levels = c("0", levels(RUN)))]
}

make_survival_input <- function(number_of_features, number_of_runs, is_labeled,
                                number_of_replicates = 1L,
                                noise_standard_deviation = 0.1,
                                censored_fraction = 0, seed = 1) {
    set.seed(seed)
    input <- data.table::CJ(
        FEATURE   = paste0("F", seq_len(number_of_features)),
        RUN       = paste0("R", seq_len(number_of_runs)),
        LABEL     = if (is_labeled) c("H", "L") else "L",
        REPLICATE = seq_len(number_of_replicates)
    )
    input[, FEATURE := factor(FEATURE)]
    input[, RUN := factor(RUN)]
    input[, newABUNDANCE := 10 + as.integer(FEATURE) + as.integer(RUN) * 0.5 +
              ifelse(LABEL == "L", 4, 0) + (REPLICATE - 1) * 0.1 +
              rnorm(.N, sd = noise_standard_deviation)]
    input[, REPLICATE := NULL]
    input[, cen := 1L]
    if (censored_fraction > 0) {
        censoring_threshold <- stats::quantile(input$newABUNDANCE,
                                               censored_fraction)
        input[newABUNDANCE < censoring_threshold,
              `:=`(cen = 0L, newABUNDANCE = censoring_threshold)]
    }
    if (is_labeled) add_reference_covariate(input)
    input
}

make_labeled_input_with_too_few_observations <- function() {
    input <- data.table::data.table(
        FEATURE      = factor(c(paste0("F", 1:8), "F1")),
        RUN          = factor(c(rep_len(paste0("R", 1:3), 8), "R1")),
        LABEL        = c(rep("H", 8), "L"),
        newABUNDANCE = c(seq(10, by = 0.5, length.out = 8), 14),
        cen          = 1L
    )
    add_reference_covariate(input)
}

coefficient_names <- function(fit) names(coef(fit))

predictor_cases <- list(
    "labeled single-feature" = list(
        input = make_survival_input(number_of_features = 1, number_of_runs = 3,
                                    is_labeled = TRUE, number_of_replicates = 3),
        has_reference_covariate = TRUE, has_feature = FALSE),
    "labeled multi-feature, enough observations" = list(
        input = make_survival_input(number_of_features = 3, number_of_runs = 4,
                                    is_labeled = TRUE),
        has_reference_covariate = TRUE, has_feature = TRUE),
    "labeled multi-feature, too few observations" = list(
        input = make_labeled_input_with_too_few_observations(),
        has_reference_covariate = TRUE, has_feature = FALSE),
    "unlabeled single-feature" = list(
        input = make_survival_input(number_of_features = 1, number_of_runs = 5,
                                    is_labeled = FALSE, number_of_replicates = 3),
        has_reference_covariate = FALSE, has_feature = FALSE),
    "unlabeled multi-feature, enough observations" = list(
        input = make_survival_input(number_of_features = 3, number_of_runs = 5,
                                    is_labeled = FALSE),
        has_reference_covariate = FALSE, has_feature = TRUE)
)

for (case_name in names(predictor_cases)) {
    case <- predictor_cases[[case_name]]
    cholesky_coefficient_names <- coefficient_names(
        MSstats:::.fitSurvival(case$input, number_of_iterations))
    expect_equal(
        any(grepl("ref_covariate", cholesky_coefficient_names)),
        case$has_reference_covariate,
        info = paste(".fitSurvival", case_name,
                     ": reference covariate in coefficients")
    )
    expect_equal(
        any(grepl("^FEATURE", cholesky_coefficient_names)), case$has_feature,
        info = paste(".fitSurvival", case_name, ": FEATURE in coefficients")
    )
    expect_equal(
        sort(coefficient_names(
            MSstats:::.fitSurvivalCG(case$input, number_of_iterations))),
        sort(cholesky_coefficient_names),
        info = paste(case_name, ": the conjugate gradient solver must select",
                     "the same predictors as the Cholesky solver")
    )
}

check_solvers_agree <- function(input, label, use_jacobi_preconditioner) {
    fit_cholesky <- MSstats:::.fitSurvival(input, number_of_iterations)
    fit_conjugate_gradient <- MSstats:::.fitSurvivalCG(
        input, number_of_iterations,
        use_jacobi_preconditioner = use_jacobi_preconditioner)
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
    actual <- summarize_fit(fit_conjugate_gradient)
    for (quantity in names(expected)) {
        expect_equal(
            actual[[quantity]], expected[[quantity]],
            tolerance = 1e-4, check.attributes = FALSE,
            info = paste(label, ":", quantity,
                         "should match the Cholesky solver")
        )
    }
}

make_noisy_censored_input <- function(is_labeled, seed) {
    make_survival_input(number_of_features = 3, number_of_runs = 4,
                        is_labeled = is_labeled,
                        noise_standard_deviation = 0.7,
                        censored_fraction = 0.2, seed = seed)
}

noisy_inputs <- list(
    labeled = make_noisy_censored_input(is_labeled = TRUE, seed = 1),
    unlabeled = make_noisy_censored_input(is_labeled = FALSE, seed = 2)
)
for (input_name in names(noisy_inputs)) {
    for (use_jacobi_preconditioner in c(FALSE, TRUE)) {
        check_solvers_agree(
            noisy_inputs[[input_name]],
            label = paste0(input_name, ", noisy, censored",
                           if (use_jacobi_preconditioner)
                               ", with Jacobi preconditioner"),
            use_jacobi_preconditioner = use_jacobi_preconditioner
        )
    }
}

noisy_input <- make_noisy_censored_input(is_labeled = FALSE, seed = 3)
fit_cholesky <- MSstats:::.fitAFTModel(noisy_input, number_of_iterations,
                                       "cholesky")
expect_inherits(
    fit_cholesky, "survreg",
    info = ".fitAFTModel(aft_solver = 'cholesky') should return a survreg fit"
)
expect_true(
    is.null(fit_cholesky$cg_diagnostics),
    info = "the Cholesky solver should not attach cg_diagnostics"
)
for (aft_solver in c("cg", "pcg")) {
    fit <- MSstats:::.fitAFTModel(noisy_input, number_of_iterations, aft_solver)
    expect_false(
        is.null(fit$cg_diagnostics),
        info = paste0(".fitAFTModel(aft_solver = '", aft_solver,
                      "') should attach cg_diagnostics")
    )
}
expect_error(
    MSstats:::.fitAFTModel(noisy_input, number_of_iterations, "cgp"),
    pattern = "aft_solver",
    info = ".fitAFTModel should reject an unsupported aft_solver instead of falling back to Cholesky"
)

expect_silent(
    MSstats:::.fitSurvivalCG(noisy_input, number_of_iterations, verbose = FALSE)
)
expect_message(
    MSstats:::.fitSurvivalCG(noisy_input, number_of_iterations, verbose = TRUE),
    pattern = "\\[AFT-CG\\] starting fit",
    info = "verbose = TRUE should report the problem size at the start of the fit"
)
expect_message(
    MSstats:::.fitSurvivalCG(noisy_input, number_of_iterations, verbose = TRUE),
    pattern = "\\[AFT-CG\\] finished",
    info = "verbose = TRUE should report a summary once fitting finishes"
)


fit_with_diagnostics <- MSstats:::.fitSurvivalCG(noisy_input,
                                                 number_of_iterations)
expect_equal(
    nrow(fit_with_diagnostics$cg_diagnostics), fit_with_diagnostics$iter,
    info = "cg_diagnostics should have one row per Newton-Raphson iteration taken"
)
expect_true(
    all(fit_with_diagnostics$cg_diagnostics$cg_iterations >= 0),
    info = paste("the conjugate gradient iteration count should be",
                 "non-negative for every Newton iteration")
)
