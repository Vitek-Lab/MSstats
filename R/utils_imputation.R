#' Decide which predictors go into a single protein's AFT imputation model
#'
#' MSstats fits an accelerated-failure-time (AFT) model per protein to
#' impute left-censored values, and predictors are chosen based on how much
#' information is actually available: whether this is a labeled (SRM)
#' experiment (\code{ref_covariate}), whether there is more than one feature 
#' to estimate a \code{FEATURE} effect for, and whether there are enough 
#' uncensored observations to estimate that effect at all.
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
#' @noRd
.buildAFTFormula = function(input) {
    FEATURE = RUN = NULL

    missingness_filter = is.finite(input$newABUNDANCE)
    n_total = nrow(input[missingness_filter, ])
    n_features = data.table::uniqueN(input[missingness_filter, FEATURE])
    n_runs = data.table::uniqueN(input[missingness_filter, RUN])
    is_labeled = data.table::uniqueN(input$LABEL) > 1

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

#' Fit an AFT survival model with SurvReg dependency
#' 
#' @param input data.table with the columns \code{.buildAFTFormula} needs.
#' @param aft_iterations maximum number of iterations for AFT model fitting.
#' @param verbose if \code{TRUE}, \code{message()} the problem size
#' (observations and parameters) before fitting and the wall time the fit
#' took afterwards, mirroring what \code{.fitSurvivalCG}'s \code{verbose}
#' reports. Meant for comparing solvers, not for routine use (this fits
#' one protein at a time).
#'
#' @importFrom stats model.frame model.matrix
#' @importFrom survival survreg
#' @keywords internal
#' @noRd
.fitSurvival = function(input, aft_iterations, verbose = FALSE) {
    set.seed(100)
    aft_formula = .buildAFTFormula(input)
    if (verbose) {
        model_frame = model.frame(aft_formula, data = input)
        design_matrix = model.matrix(attr(model_frame, "terms"), model_frame)
        message(sprintf(
            "[AFT-Cholesky] starting fit: %d observations, %d parameters",
            nrow(design_matrix), ncol(design_matrix) + 1))
    }
    fit_start_time = Sys.time()
    fit = survreg(aft_formula, data = input, dist = "gaussian",
                  control = list(maxiter = aft_iterations))
    if (verbose) {
        message(sprintf(
            "[AFT-Cholesky] finished: %d iterations, %.4f sec",
            fit$iter[length(fit$iter)],
            as.numeric(Sys.time() - fit_start_time, units = "secs")))
    }
    fit$y = NULL
    fit$linear.predictors = NULL
    fit
}


#' Per-observation log-likelihood and derivatives for a Gaussian AFT model
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
#' @noRd
.aftGaussianDerivatives = function(linear_predictor, log_scale,
                                    observed_value, exact_indicator) {
    scale = exp(log_scale)
    inverse_scale_squared = 1 / scale^2
    distance_from_prediction = observed_value - linear_predictor
    standardized_distance = distance_from_prediction / scale

    density_at_standardized_distance = dnorm(standardized_distance)
    cumulative_probability_at_standardized_distance =
        pnorm(standardized_distance)
    is_exact_observation = (exact_indicator == 1)

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

#' Evaluate the Gaussian AFT log-likelihood and derivatives at a
#' parameter guess
#'
#' @param design_matrix model matrix of the AFT fit.
#' @param coefficients current regression coefficients.
#' @param log_scale current log of the scale parameter.
#' @param observed_value observed (or censoring-threshold) values.
#' @param exact_indicator \code{1} for exact rows, \code{0} for
#' left-censored rows.
#'
#' @return the list returned by \code{.aftGaussianDerivatives}.
#'
#' @keywords internal
#' @noRd
.evaluateAFTLogLikelihood = function(design_matrix, coefficients,
                                     log_scale, observed_value,
                                     exact_indicator) {
    .aftGaussianDerivatives(
        drop(design_matrix %*% coefficients), log_scale,
        observed_value, exact_indicator)
}

#' Assemble the AFT gradient vector
#'
#' @param design_matrix model matrix of the AFT fit.
#' @param derivatives output of \code{.aftGaussianDerivatives}.
#'
#' @return gradient of the log-likelihood with respect to the regression
#' coefficients followed by the log scale.
#'
#' @keywords internal
#' @noRd
.buildAFTGradient = function(design_matrix, derivatives) {
    c(as.vector(crossprod(
          design_matrix, derivatives$gradient_wrt_linear_predictor)),
      sum(derivatives$gradient_wrt_log_scale))
}

#' Assemble the negative Hessian of the AFT log-likelihood
#'
#' @param design_matrix model matrix of the AFT fit.
#' @param derivatives output of \code{.aftGaussianDerivatives}.
#'
#' @return negative Hessian of the log-likelihood over the regression
#' coefficients and the log scale.
#'
#' @keywords internal
#' @noRd
.buildAFTNegativeHessian = function(design_matrix, derivatives) {
    regression_block = -crossprod(
        design_matrix,
        design_matrix * derivatives$second_derivative_wrt_linear_predictor)
    cross_block = -as.vector(
        crossprod(design_matrix, derivatives$cross_derivative))
    scale_block = -sum(derivatives$second_derivative_wrt_log_scale)
    rbind(cbind(regression_block, cross_block),
          c(cross_block, scale_block))
}

#' Check that an AFT log-likelihood evaluation is finite
#'
#' @param derivatives output of \code{.aftGaussianDerivatives}.
#'
#' @return \code{TRUE} if the log-likelihood and all first/second
#' derivatives used by the Newton step are finite.
#'
#' @keywords internal
#' @noRd
.isFiniteAFTFit = function(derivatives) {
    is.finite(derivatives$log_likelihood) &&
        all(is.finite(derivatives$gradient_wrt_linear_predictor)) &&
        all(is.finite(derivatives$gradient_wrt_log_scale)) &&
        all(is.finite(derivatives$second_derivative_wrt_linear_predictor)) &&
        all(is.finite(derivatives$second_derivative_wrt_log_scale))
}

#' Gauss-Newton (outer-product-of-gradients) approximation to the AFT
#' negative Hessian
#'
#' A fallback when the
#' negative Hessian is not positive definite.
#'
#' @param design_matrix model matrix of the AFT fit.
#' @param derivatives output of \code{.aftGaussianDerivatives}.
#'
#' @return crossproduct of the per-observation gradient contributions.
#'
#' @keywords internal
#' @noRd
.buildGaussNewtonApproximation = function(design_matrix, derivatives) {
    per_observation_gradient_contributions = cbind(
        design_matrix * derivatives$gradient_wrt_linear_predictor,
        derivatives$gradient_wrt_log_scale)
    crossprod(per_observation_gradient_contributions)
}

#' Run \code{.cgSolve} with its "not positive definite" warning muffled
#'
#' @param ... passed to \code{.cgSolve}.
#'
#' @return the output of \code{.cgSolve}.
#'
#' @keywords internal
#' @noRd
.cgSolveMufflingPDWarning = function(...) {
    withCallingHandlers(
        .cgSolve(...),
        warning = function(w) {
            if (grepl("not positive definite", conditionMessage(w))) {
                invokeRestart("muffleWarning")
            }
        })
}

#' Solve for one AFT Newton-Raphson step with conjugate gradient
#'
#' @param design_matrix model matrix of the AFT fit.
#' @param negative_hessian output of \code{.buildAFTNegativeHessian}.
#' @param derivatives output of \code{.aftGaussianDerivatives}.
#' @param gradient output of \code{.buildAFTGradient}.
#' @param use_jacobi_preconditioner passed to \code{.cgSolve}.
#'
#' @return a list with the \code{step}, \code{primary_iterations},
#' \code{used_fallback}, and \code{fallback_iterations}.
#'
#' @keywords internal
#' @noRd
.solveAFTNewtonStep = function(design_matrix, negative_hessian,
                               derivatives, gradient,
                               use_jacobi_preconditioner) {
    primary_solve = .cgSolveMufflingPDWarning(
        negative_hessian, gradient,
        use_jacobi_preconditioner = use_jacobi_preconditioner)
    if (primary_solve$positive_definite) {
        list(step = primary_solve$solution,
            primary_iterations = primary_solve$iterations,
            used_fallback = FALSE, fallback_iterations = 0L)
    } else {
        fallback_solve = .cgSolve(
            .buildGaussNewtonApproximation(design_matrix, derivatives),
            gradient,
            use_jacobi_preconditioner = use_jacobi_preconditioner)
        list(step = fallback_solve$solution,
            primary_iterations = primary_solve$iterations,
            used_fallback = TRUE,
            fallback_iterations = fallback_solve$iterations)
    }
}

#' Fit a Gaussian, left-censored AFT model with a conjugate-gradient
#' Newton step
#'
#' An alternative to \code{.fitSurvival} for exactly the same imputation
#' model (Gaussian accelerated-failure-time regression, left-censoring
#' only, chosen by the same \code{.buildAFTFormula} both solvers share),
#' used when \code{aft_solver = "cg"}. It runs the same kind of
#' Newton-Raphson iteration \code{survival::survreg} does - repeatedly
#' solving \code{negative_hessian \%*\% step = gradient} for the next
#' set of coefficients - but performs that linear solve with the
#' conjugate-gradient routine \code{.cgSolve} instead of the Cholesky
#' factorization \code{survreg} uses internally. The returned object is
#' classed \code{"survreg"} and carries the fields \code{predict.survreg}
#' needs, so it is a drop-in replacement anywhere \code{.fitSurvival}'s
#' result is used.
#'
#' @section Maximum likelihood estimation:
#' \strong{The model.} Each row i has a log-intensity \code{y_i}, a row
#' \code{x_i} of \code{design_matrix}, and linear predictor
#' \code{mu_i = x_i' beta}. The Gaussian AFT model says
#' \preformatted{
#'     y_i = mu_i + sigma * eps_i,    eps_i ~ N(0, 1)
#' }
#' so each true log-intensity is normally distributed around its
#' prediction with a common standard deviation \code{sigma} (the fitted
#' \code{scale}). The unknowns are \code{theta = (beta, log sigma)};
#' \code{sigma} is estimated on the log scale so that Newton steps are
#' unconstrained and can never produce a negative standard deviation.
#' Write \code{z_i = (y_i - mu_i) / sigma} for the standardized distance
#' from the prediction, \code{phi} for the standard normal density
#' (\code{dnorm}), and \code{Phi} for its CDF (\code{pnorm}).
#'
#' \strong{The objective: density for observed rows, CDF for censored
#' rows.} Maximum likelihood picks the \code{theta} under which the data we
#' saw were most probable. What we "saw" differs by row type
#' (\code{exact_indicator}):
#' \itemize{
#'   \item An \emph{observed} (exact) row has a known value, so it
#'   contributes the normal density evaluated at that value:
#'   \code{L_i = (1 / sigma) phi(z_i)}.
#'   \item A \emph{censored} row is one whose intensity fell below the
#'   detection limit. Its true value is unknown; all we know is that it lies
#'   somewhere below the threshold \code{c_i} (which
#'   \code{.setCensoredByThreshold} has substituted in as \code{y_i}). The
#'   honest contribution is therefore the total probability of landing
#'   anywhere below that threshold - the normal CDF:
#'   \code{L_i = P(Y_i <= c_i) = Phi((c_i - mu_i) / sigma)}.
#' }
#' Taking logs and summing over rows gives the objective that is maximized:
#' \preformatted{
#'     l(theta) = sum_{observed} [ log phi(z_i) - log sigma ]
#'              + sum_{censored} log Phi(z_i)
#' }
#' (the first sum is, up to a constant, ordinary least squares; the second
#' is what pulls \code{mu_i} and \code{sigma} toward values that make the
#' censored rows plausibly low). If a censored row's \code{mu_i} is well
#' above its threshold, \code{Phi(z_i)} is tiny and \code{l} is heavily
#' penalized - so the fit, and the imputed values later predicted from it,
#' respect the information that those rows were below the limit, rather
#' than ignoring them or treating the threshold as an exact value.
#'
#' \strong{Setting the derivatives to zero.} The maximum is where the score
#' (gradient of \code{l}) is zero. By the chain rule through
#' \code{mu_i = x_i' beta}, each row only needs its derivatives with respect
#' to \code{mu_i} and \code{log sigma} (computed by
#' \code{.aftGaussianDerivatives}):
#' \preformatted{
#'     observed:  d l_i / d mu_i = z_i / sigma
#'     censored:  d l_i / d mu_i = -phi(z_i) / (sigma Phi(z_i))
#' }
#' The observed term is the usual least-squares residual pull; the censored
#' term (an inverse Mills ratio) always pushes \code{mu_i} down, strongly
#' when the prediction sits above the threshold and negligibly when it is
#' already well below. These assemble into the score
#' (\code{.buildAFTGradient}), with \code{d} the vector of
#' \code{d l_i / d mu_i}:
#' \preformatted{
#'     U(theta) = [ X' d                         ]   (gradient wrt beta)
#'                [ sum_i d l_i / d log sigma    ]   (gradient wrt log sigma)
#' }
#' There is no closed-form root because of the \code{Phi} terms, so the
#' root is found iteratively.
#'
#' \strong{Newton-Raphson.} Expanding the score to first order around the
#' current guess, \code{U(theta + step) ~ U(theta) + H(theta) step}, and
#' setting it to zero gives the Newton step
#' \preformatted{
#'     -H(theta) step = U(theta),    theta_new = theta + step
#' }
#' where \code{H = d^2 l / d theta d theta'} is the Hessian of the
#' log-likelihood, with blocks
#' \preformatted{
#'     H = [ X' W X        X' v   ]    W = diag(d^2 l_i / d mu_i^2)
#'         [ v' X      sum_i s_i  ]    v_i = d^2 l_i / d mu_i d log sigma
#'                                     s_i = d^2 l_i / d (log sigma)^2
#' }
#' The negative Hessian \code{-H} (\code{.buildAFTNegativeHessian}) is
#' also known as the observed information matrix.
#' This linear system is the part solved by \code{.cgSolve} (see its
#' documentation for the conjugate-gradient math). Starting values come
#' from ordinary least squares that ignores censoring. If a step fails to
#' increase \code{l} (or produces non-finite values), the candidate is
#' pulled toward the current guess, \code{(candidate + 2 * current) / 3},
#' cutting the step to a third each time until \code{l} improves; if \code{-H} is not positive definite (so the
#' Newton direction might not point uphill), the Gauss-Newton matrix
#' \code{sum_i g_i g_i'} of per-row score contributions is used instead,
#' which is positive semi-definite by construction. Iteration stops when \code{l} changes by less than
#' \code{convergence_tolerance}, and the inverse of the final negative
#' Hessian is returned as the coefficient variance-covariance matrix, as
#' \code{survreg} does.
#'
#' @param input data.table, the same shape \code{.fitSurvival} expects.
#' @param aft_iterations maximum number of log-likelihood evaluations the
#' fit may spend. Newton-Raphson iterations and the step-halvings used to
#' recover from an overshooting step share this one budget; once it is
#' exhausted fitting stops and the last accepted coefficients and scale
#' are returned (with a non-convergence warning), rather than failing.
#' @param convergence_tolerance stop once the change in log-likelihood
#' between iterations falls below this (matches the default
#' \code{rel.tolerance} in \code{survival::survreg.control}).
#' @param use_jacobi_preconditioner if \code{TRUE}, precondition every
#' conjugate-gradient solve with the inverse of the current negative
#' Hessian's own diagonal (see \code{.cgSolve}'s
#' \code{use_jacobi_preconditioner}). This is what \code{aft_solver =
#' "pcg"} enables, versus plain conjugate gradient for \code{"cg"}.
#' @param verbose if \code{TRUE}, \code{message()} a line per
#' Newton-Raphson iteration - conjugate-gradient iterations used, whether
#' the Gauss-Newton fallback (see below) was needed, elapsed time, and the
#' resulting log-likelihood - plus a one-line summary once fitting
#' finishes. Meant for evaluating how solver choice and problem size
#' trade off against iteration count and wall time, not for routine use
#' (this fits one protein at a time, so it is easy to generate a line per
#' protein across a whole \code{dataProcess()} run).
#'
#' @return a fitted model of class \code{"survreg"}, with one added field:
#' \code{cg_diagnostics}, a data.frame with one row per Newton-Raphson
#' iteration recording the conjugate-gradient iteration counts and timing
#' described above (populated regardless of \code{verbose}, so it can be
#' inspected/aggregated programmatically after the fact).
#'
#' @importFrom stats model.frame model.matrix model.response lm.fit sd
#' @keywords internal
#' @noRd
.fitSurvivalCG = function(input, aft_iterations,
                           convergence_tolerance = 1e-9,
                           use_jacobi_preconditioner = FALSE,
                           verbose = FALSE) {
    model_frame = model.frame(.buildAFTFormula(input), data = input)
    model_terms = attr(model_frame, "terms")
    design_matrix = model.matrix(model_terms, model_frame)
    number_of_coefficients = ncol(design_matrix)
    number_of_parameters = number_of_coefficients + 1
    number_of_observations = nrow(design_matrix)

    if (verbose) {
        message(sprintf(
            "[AFT-CG] starting fit: %d observations, %d parameters, preconditioner = %s",
            number_of_observations, number_of_parameters,
            if (use_jacobi_preconditioner) "jacobi" else "none"))
    }

    response = model.response(model_frame)
    observed_value = response[, 1]
    exact_indicator = response[, 2]

    initial_fit = lm.fit(design_matrix, observed_value)
    coefficients = initial_fit$coefficients
    coefficients[!is.finite(coefficients)] = 0
    residual_standard_deviation = sd(initial_fit$residuals)
    log_scale = log(max(residual_standard_deviation, 1e-4))

    current_fit =
        .evaluateAFTLogLikelihood(design_matrix, coefficients, log_scale,
                                  observed_value, exact_indicator)
    current_log_likelihood = current_fit$log_likelihood
    number_of_iterations_used = 0
    converged = FALSE
    iterations_remaining = aft_iterations
    cg_diagnostics = vector("list", aft_iterations)

    iteration = 0
    while (iterations_remaining > 0) {
        iteration = iteration + 1
        iterations_remaining = iterations_remaining - 1
        number_of_iterations_used = iteration
        iteration_start_time = Sys.time()

        gradient = .buildAFTGradient(design_matrix, current_fit)
        negative_hessian =
            .buildAFTNegativeHessian(design_matrix, current_fit)
        newton_step = .solveAFTNewtonStep(
            design_matrix, negative_hessian, current_fit, gradient,
            use_jacobi_preconditioner)

        elapsed_seconds =
            as.numeric(Sys.time() - iteration_start_time, units = "secs")
        cg_diagnostics[[iteration]] = data.frame(
            newton_iteration = iteration,
            cg_iterations = newton_step$primary_iterations +
                newton_step$fallback_iterations,
            used_gauss_newton_fallback = newton_step$used_fallback,
            elapsed_seconds = elapsed_seconds)
        if (verbose) {
            message(sprintf(
                "[AFT-CG] newton iter %d: cg iterations = %d%s, %.4f sec",
                iteration,
                newton_step$primary_iterations + newton_step$fallback_iterations,
                if (newton_step$used_fallback) " (Gauss-Newton fallback used)" else "",
                elapsed_seconds))
        }

        candidate_coefficients =
            coefficients + newton_step$step[seq_len(number_of_coefficients)]
        candidate_log_scale =
            log_scale + newton_step$step[number_of_coefficients + 1]

        number_of_halvings = 0
        repeat {
            candidate_fit = .evaluateAFTLogLikelihood(
                design_matrix, candidate_coefficients, candidate_log_scale,
                observed_value, exact_indicator)
            candidate_improves = .isFiniteAFTFit(candidate_fit) &&
                candidate_fit$log_likelihood >= current_log_likelihood
            if (candidate_improves || iterations_remaining <= 0) {
                break
            }
            iterations_remaining = iterations_remaining - 1
            number_of_halvings = number_of_halvings + 1
            if (number_of_halvings == 1 &&
                (log_scale - candidate_log_scale) > 1.1) {
                candidate_log_scale = log_scale - 1.1
            }
            candidate_coefficients =
                (candidate_coefficients + 2 * coefficients) / 3
            candidate_log_scale = (candidate_log_scale + 2 * log_scale) / 3
        }

        if (!candidate_improves) {
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
        warning("AFT model (CG solver) did not converge within its ",
                "iteration budget; returning the last accepted ",
                "coefficients")
    }

    cg_diagnostics = do.call(
        rbind, cg_diagnostics[seq_len(number_of_iterations_used)])

    if (verbose) {
        message(sprintf(
            paste0("[AFT-CG] finished: %d newton iterations, ",
                   "%d total cg iterations, %.4f sec total, converged = %s"),
            number_of_iterations_used, sum(cg_diagnostics$cg_iterations),
            sum(cg_diagnostics$elapsed_seconds), converged))
    }

    final_negative_hessian =
        .buildAFTNegativeHessian(design_matrix, current_fit)
    variance_covariance_matrix = tryCatch(
        solve(final_negative_hessian),
        error = function(e) MASS::ginv(final_negative_hessian))

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
        loglik = current_log_likelihood,
        cg_diagnostics = cg_diagnostics
    )
    class(fit) = "survreg"
    fit
}

#' Fit the AFT imputation model with the requested solver
#'
#' Shared dispatch used by both \code{MSstatsSummarizeSingleLinear} and
#' \code{MSstatsSummarizeSingleTMP} so the \code{aft_solver}/
#' \code{aft_verbose} logic lives in one place instead of being duplicated
#' at both call sites.
#'
#' @param input data.table, the same shape \code{.fitSurvival} expects.
#' @param aft_iterations maximum number of iterations for AFT model fitting.
#' @param aft_solver "cholesky" (default, via \code{survival::survreg}),
#' "cg" (conjugate gradient), or "pcg" (conjugate gradient with a
#' Jacobi/inverse-diagonal preconditioner).
#' @param aft_verbose passed through to the chosen solver's
#' \code{verbose}: \code{.fitSurvivalCG}'s for "cg"/"pcg",
#' \code{.fitSurvival}'s for "cholesky".
#'
#' @return a fitted model of class \code{"survreg"}.
#'
#' @keywords internal
#' @noRd
.fitAFTModel = function(input, aft_iterations, aft_solver = "cholesky",
                         aft_verbose = FALSE) {
    .checkAFTSolver(aft_solver)
    if (aft_solver == "pcg") {
        .fitSurvivalCG(input, aft_iterations,
                      use_jacobi_preconditioner = TRUE, verbose = aft_verbose)
    } else if (aft_solver == "cg") {
        .fitSurvivalCG(input, aft_iterations, verbose = aft_verbose)
    } else {
        .fitSurvival(input, aft_iterations, verbose = aft_verbose)
    }
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
