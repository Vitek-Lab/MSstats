#' Compute a numerically stable Givens rotation
#'
#' Given two numbers \code{a} and \code{b}, finds the rotation (a cosine and
#' a sine) that turns the pair \code{(a, b)} into \code{(length, 0)}. It
#' returns \code{cosine}, \code{sine} and \code{length} such that
#' \code{cosine * a + sine * b == length} and
#' \code{-sine * a + cosine * b == 0}. Used by \code{.solveLeastSquares} at
#' every iteration to zero out one entry of the bidiagonal matrix built up
#' during Golub-Kahan bidiagonalization. Dividing by whichever of
#' \code{a}/\code{b} is larger in magnitude avoids the overflow/precision
#' loss that the naive \code{a / sqrt(a^2 + b^2)} formula can suffer.
#' Ported from \code{_sym_ortho} in \code{scipy/sparse/linalg/_isolve/lsqr.py}.
#'
#' @param a,b numeric scalars
#' @return list with elements \code{cosine}, \code{sine}, \code{length}
#' @keywords internal
.computeGivensRotation = function(a, b) {
    if (b == 0) {
        list(cosine = sign(a), sine = 0, length = abs(a))
    } else if (a == 0) {
        list(cosine = 0, sine = sign(b), length = abs(b))
    } else if (abs(b) > abs(a)) {
        ratio = a / b
        sine = sign(b) / sqrt(1 + ratio ^ 2)
        cosine = sine * ratio
        list(cosine = cosine, sine = sine, length = b / sine)
    } else {
        ratio = b / a
        cosine = sign(a) / sqrt(1 + ratio ^ 2)
        sine = cosine * ratio
        list(cosine = cosine, sine = sine, length = a / cosine)
    }
}


#' Solve a linear least-squares problem without forming the matrix (LSQR)
#'
#' Finds the (least-squares) solution of \code{A x = target} using LSQR
#' (C. C. Paige and M. A. Saunders, "LSQR: An algorithm for sparse linear
#' equations and sparse least squares", ACM TOMS 8(1), 43-71, 1982), an
#' iterative method built on Golub-Kahan bidiagonalization. Ported from
#' \code{lsqr()} in \code{scipy/sparse/linalg/_isolve/lsqr.py}.
#'
#' The key property this function relies on: at no point does the
#' algorithm need \code{A} itself, or any factorization of it. Every step
#' only ever multiplies a vector by \code{A} or by its transpose. Those two
#' operations are passed in as plain R functions
#' (\code{multiply}/\code{multiply_transposed}) instead of \code{A} being a
#' matrix argument, mirroring scipy's \code{LinearOperator} abstraction.
#' This is what lets a caller solve a system whose matrix is sparse and
#' structured (like the one-hot design matrices built by
#' \code{.buildRunFeatureDesign}) without ever allocating the matrix in
#' memory.
#'
#' How it works: each iteration extends the bidiagonalization of \code{A}
#' by one step (calling \code{multiply} and then
#' \code{multiply_transposed}), then applies one Givens rotation to fold
#' the newly revealed matrix entry into a running upper-bidiagonal system.
#' That immediately gives the next update to the solution, along with
#' running estimates of the residual size, the size of \code{A}, and the
#' condition number of \code{A}, which are used to decide when to stop.
#'
#' @param multiply function(v) computing \code{A \%*\% v} for a
#' length-\code{num_unknowns} vector \code{v}, returning a vector the same
#' length as \code{target}
#' @param multiply_transposed function(u) computing \code{t(A) \%*\% u} for
#' a vector \code{u} the same length as \code{target}, returning a
#' length-\code{num_unknowns} vector
#' @param target numeric vector (the right-hand side of \code{A x = target})
#' @param num_unknowns number of unknowns (i.e. number of columns of the
#' implicit A)
#' @param damping damping coefficient (Tikhonov / ridge regularization),
#' default 0
#' @param absolute_tolerance,relative_tolerance stopping tolerances on the
#' backward error, default 1e-10 (\code{atol}/\code{btol} in scipy)
#' @param max_condition_number stop if the estimated condition number of A
#' exceeds this, default 1e8
#' @param max_iterations maximum number of iterations, default
#' \code{4 * num_unknowns}
#'
#' @return list with elements \code{solution} (solution vector),
#' \code{stop_reason} (a string saying why iteration stopped, see below)
#' and \code{iterations} (number of iterations actually performed). Possible
#' stop reasons, from most to least desirable:
#' \itemize{
#'   \item \code{"target is zero"}: \code{target} (or \code{t(A) \%*\% target})
#'     is zero, so the zero vector is the exact solution
#'   \item \code{"exact solution found"}: \code{A x = target} is solved to
#'     within tolerance
#'   \item \code{"least squares solution found"}: \code{x} minimizes
#'     \code{||A x - target||} to within tolerance
#'   \item \code{"reached machine precision"}: further progress is not
#'     possible in double precision
#'   \item \code{"matrix too ill-conditioned"}: the estimated condition
#'     number of \code{A} exceeded \code{max_condition_number} (or machine
#'     precision)
#'   \item \code{"reached max iterations"}
#' }
#' @keywords internal
.solveLeastSquares = function(multiply, multiply_transposed, target, num_unknowns,
                              damping = 0,
                              absolute_tolerance = 1e-10,
                              relative_tolerance = 1e-10,
                              max_condition_number = 1e8,
                              max_iterations = NULL) {
    if (is.null(max_iterations)) {
        max_iterations = 4L * num_unknowns
    }
    machine_epsilon = .Machine$double.eps
    min_inverse_condition = if (max_condition_number > 0) 1 / max_condition_number else 0
    damping_squared = damping ^ 2

    # Running estimates, updated incrementally every iteration rather than
    # recomputed from scratch (that's the whole point of the method: no
    # step ever looks back at A or at previous vectors).
    matrix_norm_estimate = 0         # estimate of norm(A) (with damping * I stacked below A)
    condition_number_estimate = 0    # estimate of cond(A)
    sum_sq_scaled_directions = 0     # running sum used to estimate norm(inverse(A)), hence cond(A)
    damping_residual_sq_total = 0    # part of the squared residual caused by the damping term
    solution_norm_estimate = 0       # estimate of norm(solution)
    solution_norm_sq_so_far = 0      # finalized part of norm(solution)^2
    solution_norm_component = 0      # latest finalized piece of norm(solution)
    norm_rotation_cosine = -1        # rotation used to track norm(solution)
    norm_rotation_sine = 0

    solution = rep(0, num_unknowns)

    # Set up the first pair of bidiagonalization vectors:
    # left_vector = target / norm(target) (since we start from solution = 0)
    # right_vector = t(A) %*% left_vector, normalized.
    # Each iteration below extends this by one more pair.
    left_vector = target
    target_norm = sqrt(sum(target ^ 2))
    left_vector_length = target_norm

    if (left_vector_length > 0) {
        left_vector = left_vector / left_vector_length
        right_vector = multiply_transposed(left_vector)
        right_vector_length = sqrt(sum(right_vector ^ 2))
    } else {
        right_vector = rep(0, num_unknowns)
        right_vector_length = 0
    }
    if (right_vector_length > 0) {
        right_vector = right_vector / right_vector_length
    }
    update_direction = right_vector

    # Diagonal entry of the bidiagonal system that has not yet been
    # finalized by a rotation.
    pending_diagonal = right_vector_length
    # Part of norm(target) that the solution does not explain yet.
    unexplained_residual = left_vector_length
    stop_reason = NA_character_
    iteration = 0

    # Estimate of norm(t(A) %*% residual). Zero here means target is zero
    # (or orthogonal to every column of A), so solution = 0 is exact.
    transposed_residual_norm = right_vector_length * left_vector_length
    if (transposed_residual_norm == 0) {
        return(list(solution = solution, stop_reason = "target is zero",
                    iterations = 0))
    }

    while (iteration < max_iterations) {
        iteration = iteration + 1

        # --- Extend the Golub-Kahan bidiagonalization by one step ---
        # Two-term recurrence producing the next orthonormal pair:
        #   left_length  * new_left  = A    %*% right - right_length * left
        #   right_length * new_right = t(A) %*% left  - left_length  * right
        # `right_vector` is kept separate from `update_direction` (the
        # vector actually used to update the solution below). They start
        # out equal but differ after the first iteration.
        left_vector = multiply(right_vector) - right_vector_length * left_vector
        left_vector_length = sqrt(sum(left_vector ^ 2))
        if (left_vector_length > 0) {
            left_vector = left_vector / left_vector_length
            matrix_norm_estimate = sqrt(matrix_norm_estimate ^ 2 +
                                            right_vector_length ^ 2 +
                                            left_vector_length ^ 2 +
                                            damping_squared)
            right_vector = multiply_transposed(left_vector) -
                left_vector_length * right_vector
            right_vector_length = sqrt(sum(right_vector ^ 2))
            if (right_vector_length > 0) {
                right_vector = right_vector / right_vector_length
            }
        }

        # --- Fold in the damping term (ridge regularization), if any ---
        # A rotation that absorbs damping * I into the current diagonal
        # entry before the main elimination step below.
        if (damping > 0) {
            damped_diagonal = sqrt(pending_diagonal ^ 2 + damping_squared)
            damping_cosine = pending_diagonal / damped_diagonal
            damping_sine = damping / damped_diagonal
            damping_residual = damping_sine * unexplained_residual
            unexplained_residual = damping_cosine * unexplained_residual
        } else {
            damped_diagonal = pending_diagonal
            damping_residual = 0
        }

        # --- Eliminate the newly revealed below-diagonal entry ---
        # This turns the lower-bidiagonal system built up so far into an
        # upper-bidiagonal one, one row at a time.
        rotation = .computeGivensRotation(damped_diagonal, left_vector_length)
        rotation_cosine = rotation$cosine
        rotation_sine = rotation$sine
        diagonal = rotation$length

        off_diagonal = rotation_sine * right_vector_length
        pending_diagonal = -rotation_cosine * right_vector_length
        explained_residual = rotation_cosine * unexplained_residual
        unexplained_residual = rotation_sine * unexplained_residual
        transposed_residual_factor = rotation_sine * explained_residual

        # --- Update the solution and the update direction ---
        step_size = explained_residual / diagonal
        direction_carryover = -off_diagonal / diagonal
        scaled_direction = update_direction / diagonal

        solution = solution + step_size * update_direction
        update_direction = right_vector + direction_carryover * update_direction
        sum_sq_scaled_directions = sum_sq_scaled_directions + sum(scaled_direction ^ 2)

        # --- Update the running estimate of norm(solution) ---
        # Uses one more rotation, so the norm never has to be recomputed
        # from the full solution vector.
        norm_off_diagonal = norm_rotation_sine * diagonal
        norm_pending_diagonal = -norm_rotation_cosine * diagonal
        norm_numerator = explained_residual - norm_off_diagonal * solution_norm_component
        tentative_norm_component = norm_numerator / norm_pending_diagonal
        solution_norm_estimate = sqrt(solution_norm_sq_so_far + tentative_norm_component ^ 2)
        norm_diagonal = sqrt(norm_pending_diagonal ^ 2 + off_diagonal ^ 2)
        norm_rotation_cosine = norm_pending_diagonal / norm_diagonal
        norm_rotation_sine = off_diagonal / norm_diagonal
        solution_norm_component = norm_numerator / norm_diagonal
        solution_norm_sq_so_far = solution_norm_sq_so_far + solution_norm_component ^ 2

        # --- Update size/condition estimates and decide whether to stop ---
        condition_number_estimate = matrix_norm_estimate * sqrt(sum_sq_scaled_directions)
        damping_residual_sq_total = damping_residual_sq_total + damping_residual ^ 2
        residual_norm = sqrt(unexplained_residual ^ 2 + damping_residual_sq_total)
        transposed_residual_norm = right_vector_length * abs(transposed_residual_factor)

        # How close is A x to target, relative to the size of target?
        relative_residual = residual_norm / target_norm
        # How close is x to minimizing ||A x - target|| (normal equations)?
        relative_least_squares_error = transposed_residual_norm /
            (matrix_norm_estimate * residual_norm + machine_epsilon)
        inverse_condition_estimate = 1 / (condition_number_estimate + machine_epsilon)
        solution_scale = matrix_norm_estimate * solution_norm_estimate / target_norm
        residual_vs_precision = relative_residual / (1 + solution_scale)
        residual_tolerance = relative_tolerance + absolute_tolerance * solution_scale

        # `1 + x <= 1` is true when x is too small to register next to 1 in
        # double precision.
        if (relative_residual <= residual_tolerance) {
            stop_reason = "exact solution found"
        } else if (relative_least_squares_error <= absolute_tolerance) {
            stop_reason = "least squares solution found"
        } else if (inverse_condition_estimate <= min_inverse_condition) {
            stop_reason = "matrix too ill-conditioned"
        } else if (1 + residual_vs_precision <= 1 ||
                   1 + relative_least_squares_error <= 1) {
            stop_reason = "reached machine precision"
        } else if (1 + inverse_condition_estimate <= 1) {
            stop_reason = "matrix too ill-conditioned"
        } else if (iteration >= max_iterations) {
            stop_reason = "reached max iterations"
        }

        if (!is.na(stop_reason)) break
    }

    list(solution = solution, stop_reason = stop_reason, iterations = iteration)
}


#' Build a matrix-free design for the model "1 + run + feature"
#'
#' The model \code{log2inty ~ run + feature} (treatment-contrast coded,
#' i.e. the first level of each factor is a "reference" level that gets no
#' dummy column of its own) has a design matrix where every row has at
#' most 3 nonzero entries: the intercept (always 1), at most one `run`
#' dummy, and at most one `feature` dummy. That means the matrix-vector
#' products \code{.solveLeastSquares} needs (\code{A \%*\% v} and
#' \code{t(A) \%*\% u}) can be computed directly from which run and
#' feature each row belongs to, without ever building the matrix \code{A}:
#' \itemize{
#'   \item \code{A \%*\% coefficients} (\code{multiply}) looks up, for every
#'     row, which run/feature coefficient (if any) applies and adds it to
#'     the intercept. This gives the fitted value for every row.
#'   \item \code{t(A) \%*\% values} (\code{multiply_transposed}) sums
#'     \code{values} within each run and within each feature (plus an
#'     overall sum for the intercept).
#' }
#' \code{run_levels}/\code{feature_levels} let the caller fix the full set
#' of levels (and which one is the reference) independently of which rows
#' are actually being fit. This matches how \code{stats::model.matrix}
#' derives factor levels from a column before any rows are dropped for
#' missing values, so a level with zero surviving observations still gets
#' a column (all zeros, with its coefficient kept at 0 because
#' \code{.solveLeastSquares} returns the minimum-norm solution) instead of
#' silently changing which coefficients the model has.
#'
#' @param run,feature vectors (same length), the two grouping variables
#' @param run_levels,feature_levels optional explicit level sets (else
#' taken from \code{factor()}'s default sorted ordering, matching
#' \code{model.matrix}'s default reference level)
#'
#' @return list with elements \code{multiply}, \code{multiply_transposed}
#' (functions, see \code{.solveLeastSquares}) and \code{num_coefficients}
#' @keywords internal
.buildRunFeatureDesign = function(run, feature, run_levels = NULL, feature_levels = NULL) {
    run_factor = if (is.null(run_levels)) factor(run) else factor(run, levels = run_levels)
    feature_factor = if (is.null(feature_levels)) factor(feature) else factor(feature, levels = feature_levels)

    num_rows = length(run)
    num_runs = nlevels(run_factor)
    num_features = nlevels(feature_factor)

    # 0 = reference level (no coefficient of its own, only the intercept
    # applies); 1..(num_runs - 1) numbers the non-reference run
    # coefficients, and likewise for feature.
    run_number = as.integer(run_factor) - 1L
    feature_number = as.integer(feature_factor) - 1L

    # Coefficient layout: [intercept, run coefficients..., feature coefficients...]
    run_start = 1L
    feature_start = 1L + (num_runs - 1L)
    num_coefficients = feature_start + (num_features - 1L)

    is_nonreference_run = run_number > 0L
    is_nonreference_feature = feature_number > 0L

    multiply = function(coefficients) {
        fitted = rep(coefficients[1], num_rows)
        fitted[is_nonreference_run] = fitted[is_nonreference_run] +
            coefficients[run_start + run_number[is_nonreference_run]]
        fitted[is_nonreference_feature] = fitted[is_nonreference_feature] +
            coefficients[feature_start + feature_number[is_nonreference_feature]]
        fitted
    }

    multiply_transposed = function(values) {
        totals = numeric(num_coefficients)
        totals[1] = sum(values)
        if (num_runs > 1L) {
            run_totals = rowsum(values, run_number)
            group_number = as.integer(rownames(run_totals))
            is_nonreference = group_number > 0L
            totals[run_start + group_number[is_nonreference]] = run_totals[is_nonreference, 1]
        }
        if (num_features > 1L) {
            feature_totals = rowsum(values, feature_number)
            group_number = as.integer(rownames(feature_totals))
            is_nonreference = group_number > 0L
            totals[feature_start + group_number[is_nonreference]] = feature_totals[is_nonreference, 1]
        }
        totals
    }

    list(multiply = multiply, multiply_transposed = multiply_transposed,
         num_coefficients = num_coefficients)
}


#' Apply observation weights to a matrix-free design
#'
#' Represents \code{diag(sqrt(weights)) \%*\% A} without ever forming that
#' product: \code{multiply} scales the design's output by
#' \code{sqrt(weights)}, and \code{multiply_transposed} scales its input by
#' \code{sqrt(weights)} first. Rebuilding just this thin wrapper each
#' reweighting iteration (rather than \code{.buildRunFeatureDesign} itself)
#' means the run/feature bookkeeping is done once per protein, not once
#' per iteration.
#'
#' @param design a matrix-free design (list with \code{multiply},
#' \code{multiply_transposed}, \code{num_coefficients}, see
#' \code{.buildRunFeatureDesign})
#' @param weights numeric vector of weights, one per observation
#' @return a matrix-free design of the same form for the weighted system
#' @keywords internal
.applyWeightsToDesign = function(design, weights) {
    sqrt_weights = sqrt(weights)
    list(
        multiply = function(coefficients) sqrt_weights * design$multiply(coefficients),
        multiply_transposed = function(values) design$multiply_transposed(sqrt_weights * values),
        num_coefficients = design$num_coefficients
    )
}


#' Fit a robust (Huber) regression of log2inty on run and feature without
#' building the design matrix
#'
#' Reimplements the "M" estimation / "Huber" scale-estimate branch of
#' \code{MASS::rlm.default} (see \code{MASS/R/rlm.R}) for the specific
#' model \code{log2inty ~ run + feature}. It uses iteratively reweighted
#' least squares: fit, downweight observations with large residuals,
#' refit, and repeat until the residuals stop changing. Every weighted
#' least-squares fit that \code{MASS::rlm} does with
#' \code{stats::lm.wfit(x, y, w, method = "qr")} is done here with
#' \code{.solveLeastSquares()} on \code{.buildRunFeatureDesign}'s
#' matrix-free design. Only the code path used by \code{.fitHuber} is
#' reproduced: \code{init = "ls"}, \code{psi = MASS::psi.huber},
#' \code{scale.est = "Huber"}, no case/prior weights.
#'
#' @param log2inty numeric response vector
#' @param run,feature grouping vectors, same length as \code{log2inty}
#' @param run_levels,feature_levels optional explicit level sets, see
#' \code{.buildRunFeatureDesign}
#' @param huber_threshold tuning constant for Huber's method: residuals
#' larger than this many scale units get downweighted (\code{k2} in
#' \code{MASS::rlm})
#' @param max_iterations maximum number of reweighting iterations
#' (\code{maxit} in \code{MASS::rlm})
#' @param convergence_tolerance stop once the relative change in residuals
#' between iterations is below this (\code{acc} in \code{MASS::rlm})
#'
#' @return list with elements \code{coefficients}, \code{residuals},
#' \code{fitted.values}, \code{scale}, \code{df.residual}, \code{rank} and
#' \code{converged} (named to match the corresponding \code{MASS::rlm}
#' fields)
#' @importFrom MASS psi.huber
#' @importFrom stats mad pnorm dnorm
#' @keywords internal
.fitRobustRunFeatureModel = function(log2inty, run, feature,
                                     run_levels = NULL, feature_levels = NULL,
                                     huber_threshold = 1.345,
                                     max_iterations = 20,
                                     convergence_tolerance = 1e-4) {
    num_observations = length(log2inty)
    design = .buildRunFeatureDesign(run, feature, run_levels, feature_levels)
    num_coefficients = design$num_coefficients
    residual_df = num_observations - num_coefficients

    # Start from the ordinary (unweighted) least-squares fit.
    coefficients = .solveLeastSquares(design$multiply, design$multiply_transposed,
                                      log2inty, num_coefficients)$solution
    fitted = design$multiply(coefficients)
    residuals = log2inty - fitted
    scale = stats::mad(residuals, center = 0)

    # Constant that makes Huber's scale estimate consistent for normally
    # distributed errors (Huber's "Proposal 2").
    prob_within_threshold = 2 * stats::pnorm(huber_threshold) - 1
    scale_correction = prob_within_threshold +
        huber_threshold ^ 2 * (1 - prob_within_threshold) -
        2 * huber_threshold * stats::dnorm(huber_threshold)

    converged = FALSE
    for (iteration in seq_len(max_iterations)) {
        previous_residuals = residuals

        # Re-estimate the scale, capping each residual at the threshold.
        capped_squared_residuals = pmin(residuals ^ 2, (huber_threshold * scale) ^ 2)
        scale = sqrt(sum(capped_squared_residuals) / (residual_df * scale_correction))
        if (scale == 0) {
            converged = TRUE
            break
        }

        # Downweight observations with large residuals, then refit.
        weights = MASS::psi.huber(residuals / scale)
        weighted_design = .applyWeightsToDesign(design, weights)
        coefficients = .solveLeastSquares(weighted_design$multiply,
                                          weighted_design$multiply_transposed,
                                          sqrt(weights) * log2inty,
                                          num_coefficients)$solution
        fitted = design$multiply(coefficients)
        residuals = log2inty - fitted

        relative_change = sqrt(sum((previous_residuals - residuals) ^ 2) /
                                   max(1e-20, sum(previous_residuals ^ 2)))
        converged = relative_change <= convergence_tolerance
        if (converged) break
    }

    list(
        coefficients = coefficients,
        residuals = residuals,
        fitted.values = fitted,
        scale = scale,
        df.residual = residual_df,
        rank = num_coefficients,
        converged = converged
    )
}
