#' Decide which predictors go into a single protein's AFT imputation model
#'
#' MSstats fits an accelerated-failure-time (AFT) model per protein to
#' impute left-censored values, and predictors are chosen based on how much
#' information is actually available: whether this is a labeled (SRM)
#' experiment with a reference channel (\code{ref_covariate}), whether
#' there is more than one feature to estimate a \code{FEATURE} effect for,
#' and whether there are enough uncensored observations to estimate that
#' effect at all. Both \code{.fitSurvival} (Cholesky-based, via
#' \code{survival::survreg}) and \code{.fitSurvivalCG} (conjugate-gradient
#' based) share this selection logic, so the two solvers always fit the
#' same model and differ only in how the Newton step is solved.
#'
#' @param input data.table with columns \code{newABUNDANCE}, \code{cen},
#' \code{RUN}, \code{FEATURE}, \code{LABEL}, and (for labeled experiments)
#' \code{ref_covariate}.
#'
#' @return a formula whose left side is
#' \code{Surv(newABUNDANCE, cen, type = "left")}.
#'
#' @importFrom data.table uniqueN
#' @importFrom survival Surv
#' @keywords internal
.buildAFTFormula = function(input) {
    FEATURE = RUN = NULL

    missingness_filter = is.finite(input$newABUNDANCE)
    n_total = nrow(input[missingness_filter, ])
    n_features = data.table::uniqueN(input[missingness_filter, FEATURE])
    n_runs = data.table::uniqueN(input[missingness_filter, RUN])
    is_labeled = data.table::uniqueN(input$LABEL) > 1

    # With too few uncensored observations, there isn't enough information
    # left to also estimate a separate effect per feature.
    not_enough_data_for_feature_effect = n_total < n_features + n_runs - 1

    if (is_labeled) {
        if (length(unique(input$FEATURE)) == 1 ||
            not_enough_data_for_feature_effect) {
            # with a single feature (or too little data), a FEATURE term
            # either adds nothing or keeps the model from converging /
            # gives it the wrong intercept - need to check
            Surv(newABUNDANCE, cen, type = "left") ~ RUN + ref_covariate
        } else {
            Surv(newABUNDANCE, cen, type = "left") ~
                FEATURE + RUN + ref_covariate
        }
    } else {
        if (n_features == 1L || not_enough_data_for_feature_effect) {
            Surv(newABUNDANCE, cen, type = "left") ~ RUN
        } else {
            Surv(newABUNDANCE, cen, type = "left") ~ FEATURE + RUN
        }
    }
}

#' @importFrom survival survreg
#' @keywords internal
.fitSurvival = function(input, aft_iterations) {
    # TODO: set.seed here?
    set.seed(100)
    fit = survreg(.buildAFTFormula(input), data = input, dist = "gaussian",
                  control = list(maxiter = aft_iterations))
    fit$y = NULL
    fit$linear.predictors = NULL
    fit
}


#' Per-observation log-likelihood and derivatives for a Gaussian AFT model
#'
#' Computes what a Newton-Raphson step needs at the current parameter
#' guess: the log-likelihood, its first derivative with respect to the
#' linear predictor and to the log of the scale parameter, and the
#' corresponding second derivatives - all summed/assembled later into the
#' score vector and information matrix by \code{.fitSurvivalCG}. This only
#' covers the two cases MSstats' AFT imputation actually uses: an exact
#' (uncensored) observation, or one left-censored below a detection-limit
#' ceiling (\code{Surv(y, cen, type = "left")} with \code{cen == 0}).
#'
#' The formulas are transcribed term-for-term from \code{survival}'s own
#' C implementation (\code{survregc1.c}'s \code{gauss_d} function and its
#' "exact"/"left censored" cases) rather than re-derived by hand, since a
#' hand re-derivation is an easy place to introduce a sign error; this
#' function's correctness is instead checked against numerical
#' differentiation of the log-likelihood (see
#' \code{test_utils_imputation_cg.R}).
#'
#' @param linear_predictor current linear predictor
#' (\code{model_matrix \%*\% coefficients}).
#' @param log_scale current log of the scale parameter.
#' @param observed_value observed value (or, for censored rows, the
#' detection-limit ceiling substituted in by
#' \code{.setCensoredByThreshold}).
#' @param exact_indicator \code{1} for an exact/uncensored observation,
#' \code{0} for one left-censored below \code{observed_value}.
#'
#' @return a list with the total \code{log_likelihood}, and
#' per-observation vectors \code{gradient_wrt_linear_predictor},
#' \code{second_derivative_wrt_linear_predictor},
#' \code{gradient_wrt_log_scale}, \code{second_derivative_wrt_log_scale},
#' and \code{cross_derivative}
#' (d2 log_likelihood / d linear_predictor d log_scale).
#'
#' @importFrom stats dnorm pnorm
#' @keywords internal
.aftGaussianDerivatives = function(linear_predictor, log_scale,
                                    observed_value, exact_indicator) {
    scale = exp(log_scale)
    inverse_scale_squared = 1 / scale^2

    # How far the observation sits from its predicted value, in raw units
    # and in standard deviations.
    distance_from_prediction = observed_value - linear_predictor
    standardized_distance = distance_from_prediction / scale

    density_at_standardized_distance = dnorm(standardized_distance)
    cumulative_probability_at_standardized_distance =
        pnorm(standardized_distance)
    is_exact_observation = (exact_indicator == 1)

    # --- exact (uncensored) observations --------------------------------
    # log-likelihood contribution is log(density) - log(scale); what
    # follows is that expression's derivatives wrt linear_predictor and
    # log_scale.
    exact_log_likelihood =
        log(density_at_standardized_distance) - log_scale
    exact_gradient_wrt_linear_predictor = standardized_distance / scale
    exact_log_density_curvature =
        (standardized_distance^2 - 1) * inverse_scale_squared
    exact_second_derivative_wrt_linear_predictor =
        exact_log_density_curvature -
        exact_gradient_wrt_linear_predictor^2
    exact_gradient_wrt_log_scale_before_adjustment =
        exact_gradient_wrt_linear_predictor * distance_from_prediction
    exact_cross_derivative =
        distance_from_prediction * exact_log_density_curvature -
        exact_gradient_wrt_linear_predictor *
        (exact_gradient_wrt_log_scale_before_adjustment + 1)
    exact_second_derivative_wrt_log_scale =
        distance_from_prediction^2 * exact_log_density_curvature -
        exact_gradient_wrt_log_scale_before_adjustment *
        (1 + exact_gradient_wrt_log_scale_before_adjustment)
    exact_gradient_wrt_log_scale =
        exact_gradient_wrt_log_scale_before_adjustment - 1

    # Guard against the density underflowing to exactly zero (only
    # happens for astronomically large |standardized_distance|, e.g. from
    # a wild early Newton guess). Any reasonable derivative works here,
    # since the collapsed log-likelihood itself is what triggers
    # step-halving.
    exact_density_underflowed = density_at_standardized_distance <= 0
    exact_log_likelihood =
        ifelse(exact_density_underflowed, -200, exact_log_likelihood)
    exact_gradient_wrt_linear_predictor = ifelse(
        exact_density_underflowed, -standardized_distance / scale,
        exact_gradient_wrt_linear_predictor)
    exact_second_derivative_wrt_linear_predictor = ifelse(
        exact_density_underflowed, -1 / scale,
        exact_second_derivative_wrt_linear_predictor)
    exact_gradient_wrt_log_scale =
        ifelse(exact_density_underflowed, 0, exact_gradient_wrt_log_scale)
    exact_cross_derivative =
        ifelse(exact_density_underflowed, 0, exact_cross_derivative)
    exact_second_derivative_wrt_log_scale = ifelse(
        exact_density_underflowed, 0,
        exact_second_derivative_wrt_log_scale)

    # --- left-censored observations (true value <= the recorded ceiling) -
    # log-likelihood contribution is log(Phi(standardized_distance));
    # "censoring_hazard" plays the same role for these rows that the
    # density itself plays above.
    censored_log_likelihood =
        log(cumulative_probability_at_standardized_distance)
    censoring_hazard = density_at_standardized_distance /
        (cumulative_probability_at_standardized_distance * scale)
    censored_gradient_wrt_linear_predictor = -censoring_hazard
    censored_log_density_curvature =
        -standardized_distance * density_at_standardized_distance *
        inverse_scale_squared /
        cumulative_probability_at_standardized_distance
    censored_second_derivative_wrt_linear_predictor =
        censored_log_density_curvature -
        censored_gradient_wrt_linear_predictor^2
    censored_gradient_wrt_log_scale =
        censored_gradient_wrt_linear_predictor * distance_from_prediction
    censored_cross_derivative =
        distance_from_prediction * censored_log_density_curvature -
        censored_gradient_wrt_linear_predictor *
        (censored_gradient_wrt_log_scale + 1)
    censored_second_derivative_wrt_log_scale =
        distance_from_prediction^2 * censored_log_density_curvature -
        censored_gradient_wrt_log_scale * (1 + censored_gradient_wrt_log_scale)

    # Same underflow guard as above, triggered when the cumulative
    # probability collapses to zero (standardized_distance very
    # negative).
    censored_probability_underflowed =
        cumulative_probability_at_standardized_distance <= 0
    censored_log_likelihood = ifelse(
        censored_probability_underflowed, -200, censored_log_likelihood)
    censored_gradient_wrt_linear_predictor = ifelse(
        censored_probability_underflowed, -standardized_distance / scale,
        censored_gradient_wrt_linear_predictor)
    censored_second_derivative_wrt_linear_predictor = ifelse(
        censored_probability_underflowed, 0,
        censored_second_derivative_wrt_linear_predictor)
    censored_gradient_wrt_log_scale = ifelse(
        censored_probability_underflowed, 0, censored_gradient_wrt_log_scale)
    censored_cross_derivative = ifelse(
        censored_probability_underflowed, 0, censored_cross_derivative)
    censored_second_derivative_wrt_log_scale = ifelse(
        censored_probability_underflowed, 0,
        censored_second_derivative_wrt_log_scale)

    list(
        log_likelihood = sum(ifelse(
            is_exact_observation, exact_log_likelihood,
            censored_log_likelihood)),
        gradient_wrt_linear_predictor = ifelse(
            is_exact_observation, exact_gradient_wrt_linear_predictor,
            censored_gradient_wrt_linear_predictor),
        second_derivative_wrt_linear_predictor = ifelse(
            is_exact_observation,
            exact_second_derivative_wrt_linear_predictor,
            censored_second_derivative_wrt_linear_predictor),
        gradient_wrt_log_scale = ifelse(
            is_exact_observation, exact_gradient_wrt_log_scale,
            censored_gradient_wrt_log_scale),
        second_derivative_wrt_log_scale = ifelse(
            is_exact_observation, exact_second_derivative_wrt_log_scale,
            censored_second_derivative_wrt_log_scale),
        cross_derivative = ifelse(
            is_exact_observation, exact_cross_derivative,
            censored_cross_derivative)
    )
}

#' Fit a Gaussian, left-censored AFT model with a conjugate-gradient
#' Newton step
#'
#' An alternative to \code{.fitSurvival} for exactly the same imputation
#' model (Gaussian accelerated-failure-time regression, left-censoring
#' only, chosen by the same \code{.buildAFTFormula} both solvers share),
#' used when \code{aft_solver = "cg"}. It runs the same kind of
#' Newton-Raphson iteration \code{survival::survreg} does - repeatedly
#' solving \code{information_matrix \%*\% step = gradient} for the next
#' set of coefficients - but performs that linear solve with the
#' conjugate-gradient routine \code{.cgSolve} instead of the Cholesky
#' factorization \code{survreg} uses internally. The returned object is
#' classed \code{"survreg"} and carries the fields \code{predict.survreg}
#' needs, so it is a drop-in replacement anywhere \code{.fitSurvival}'s
#' result is used.
#'
#' @param input data.table, the same shape \code{.fitSurvival} expects.
#' @param aft_iterations maximum number of Newton-Raphson iterations.
#' @param convergence_tolerance stop once the change in log-likelihood
#' between iterations falls below this (matches the default
#' \code{rel.tolerance} in \code{survival::survreg.control}).
#'
#' @return a fitted model of class \code{"survreg"}.
#'
#' @importFrom stats model.frame model.matrix model.response lm.fit sd
#' @keywords internal
.fitSurvivalCG = function(input, aft_iterations,
                           convergence_tolerance = 1e-9) {
    model_frame = model.frame(.buildAFTFormula(input), data = input)
    model_terms = attr(model_frame, "terms")
    design_matrix = model.matrix(model_terms, model_frame)
    number_of_coefficients = ncol(design_matrix)

    response = model.response(model_frame)
    observed_value = response[, 1]
    exact_indicator = response[, 2]

    # Initial guess: an ordinary least-squares fit for the regression
    # coefficients (treating the detection-limit ceiling already
    # substituted into censored rows as if it were observed), and the
    # residual standard deviation for the scale parameter. A
    # rank-deficient design leaves some coefficients unidentified
    # (reported as NA by lm.fit); start those at zero.
    initial_fit = lm.fit(design_matrix, observed_value)
    coefficients = initial_fit$coefficients
    coefficients[!is.finite(coefficients)] = 0
    residual_standard_deviation = sd(initial_fit$residuals)
    log_scale = log(max(residual_standard_deviation, 1e-4))

    evaluate_log_likelihood_and_derivatives = function(coefficients,
                                                        log_scale) {
        .aftGaussianDerivatives(
            drop(design_matrix %*% coefficients), log_scale,
            observed_value, exact_indicator)
    }

    build_gradient = function(derivatives) {
        c(as.vector(crossprod(
              design_matrix, derivatives$gradient_wrt_linear_predictor)),
          sum(derivatives$gradient_wrt_log_scale))
    }

    build_information_matrix = function(derivatives) {
        # Regression block: -t(X) %*% diag(second_derivative) %*% X,
        # computed without forming the diagonal matrix explicitly.
        regression_block = -crossprod(
            design_matrix,
            design_matrix * derivatives$second_derivative_wrt_linear_predictor)
        cross_block = -as.vector(
            crossprod(design_matrix, derivatives$cross_derivative))
        scale_block = -sum(derivatives$second_derivative_wrt_log_scale)
        rbind(cbind(regression_block, cross_block),
              c(cross_block, scale_block))
    }

    is_finite_fit = function(derivatives) {
        is.finite(derivatives$log_likelihood) &&
            all(is.finite(derivatives$gradient_wrt_linear_predictor)) &&
            all(is.finite(derivatives$gradient_wrt_log_scale)) &&
            all(is.finite(derivatives$second_derivative_wrt_linear_predictor)) &&
            all(is.finite(derivatives$second_derivative_wrt_log_scale))
    }

    # A Newton step away from the optimum, the exact information matrix
    # is not guaranteed to be positive definite. survival::survreg falls
    # back, in that situation, to the sum of the outer products of each
    # observation's own contribution to the gradient - always positive
    # semi-definite by construction, and equal to the exact information
    # matrix in expectation (this is the classic Gauss-Newton / BHHH
    # approximation). Mirror that fallback here.
    build_gauss_newton_approximation = function(derivatives) {
        per_observation_gradient_contributions = cbind(
            design_matrix * derivatives$gradient_wrt_linear_predictor,
            derivatives$gradient_wrt_log_scale)
        crossprod(per_observation_gradient_contributions)
    }

    solve_newton_step = function(information_matrix, derivatives, gradient) {
        information_matrix_is_not_positive_definite = FALSE
        step = withCallingHandlers(
            .cgSolve(information_matrix, gradient),
            warning = function(w) {
                if (grepl("not positive definite", conditionMessage(w))) {
                    information_matrix_is_not_positive_definite <<- TRUE
                }
                invokeRestart("muffleWarning")
            })
        if (information_matrix_is_not_positive_definite) {
            step = .cgSolve(build_gauss_newton_approximation(derivatives),
                            gradient)
        }
        step
    }

    current_fit =
        evaluate_log_likelihood_and_derivatives(coefficients, log_scale)
    current_log_likelihood = current_fit$log_likelihood
    number_of_iterations_used = 0
    converged = FALSE

    for (iteration in seq_len(aft_iterations)) {
        number_of_iterations_used = iteration
        gradient = build_gradient(current_fit)
        information_matrix = build_information_matrix(current_fit)
        newton_step =
            solve_newton_step(information_matrix, current_fit, gradient)

        candidate_coefficients =
            coefficients + newton_step[seq_len(number_of_coefficients)]
        candidate_log_scale =
            log_scale + newton_step[number_of_coefficients + 1]

        # Step-halving: if the Newton step overshoots (a non-finite or
        # decreasing log-likelihood), back the trial point off toward the
        # last accepted one, mirroring survival::survreg's own recovery
        # strategy (survreg6.c) rather than simply rejecting the step
        # outright.
        number_of_halvings = 0
        halving_exhausted = FALSE
        repeat {
            candidate_fit = evaluate_log_likelihood_and_derivatives(
                candidate_coefficients, candidate_log_scale)
            candidate_improves = is_finite_fit(candidate_fit) &&
                candidate_fit$log_likelihood >= current_log_likelihood
            if (candidate_improves) {
                break
            }
            number_of_halvings = number_of_halvings + 1
            if (number_of_halvings > 30) {
                halving_exhausted = TRUE
                break
            }
            if (number_of_halvings == 1 &&
                (log_scale - candidate_log_scale) > 1.1) {
                # a single huge drop in scale is the most common cause of
                # a bad trial; keep the first back-off from cutting scale
                # by more than a factor of exp(1.1), same as survreg6.c
                candidate_log_scale = log_scale - 1.1
            }
            candidate_coefficients =
                (candidate_coefficients + 2 * coefficients) / 3
            candidate_log_scale = (candidate_log_scale + 2 * log_scale) / 3
        }

        if (halving_exhausted) {
            break
        }

        relative_change =
            abs(1 - current_log_likelihood / candidate_fit$log_likelihood)
        absolute_change =
            abs(candidate_fit$log_likelihood - current_log_likelihood)

        coefficients = candidate_coefficients
        log_scale = candidate_log_scale
        current_fit = candidate_fit
        current_log_likelihood = candidate_fit$log_likelihood

        if (relative_change <= convergence_tolerance ||
            absolute_change <= convergence_tolerance) {
            converged = TRUE
            break
        }
    }

    if (!converged) {
        warning("AFT model (CG solver) ran out of iterations and did not ",
                "converge")
    }

    final_information_matrix = build_information_matrix(current_fit)
    variance_covariance_matrix = tryCatch(
        solve(final_information_matrix),
        error = function(e) MASS::ginv(final_information_matrix))

    fitted_coefficients = coefficients
    names(fitted_coefficients) = colnames(design_matrix)

    is_factor_column = vapply(model_frame, is.factor, logical(1))

    fit = list(
        coefficients = fitted_coefficients,
        var = variance_covariance_matrix,
        scale = exp(log_scale),
        terms = model_terms,
        xlevels = lapply(model_frame[is_factor_column], levels),
        dist = "gaussian",
        iter = number_of_iterations_used,
        loglik = current_log_likelihood
    )
    class(fit) = "survreg"
    fit
}

#' Get predicted values from a survival model
#' @param input data.table
#' @return numeric vector of predictions
#' @importFrom stats predict
#' @keywords internal
.addSurvivalPredictions = function(input) {
    LABEL = NULL
    
    survival_fit = .fitSurvival(input[LABEL == "L", ])
    predict(survival_fit, newdata = input)
}
