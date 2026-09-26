#' Solve a symmetric positive (semi-)definite linear system via conjugate
#' gradient
#'
#' A minimal, single right-hand-side conjugate gradient solver, used as the
#' Newton-step linear solve in \code{.fitSurvivalCG}. Modeled on
#' \code{lfe::cgsolve}, stripped down to a single dense matrix and a single
#' right-hand-side vector (no multi-column batching, no \code{Matrix}-package
#' or operator/closure dispatch - neither is needed for the small, dense AFT
#' information matrices this is used on). Optionally applies a Jacobi
#' (inverse-diagonal) preconditioner, which \code{lfe::cgsolve} does not
#' support at all.
#'
#' Iteration always starts from the zero vector and stops once the true
#' (unpreconditioned) residual is below
#' \code{1e-8 * max(norm(right_hand_side), 1)} - relative for large
#' right-hand sides, absolute (\code{1e-8}) for small ones - so the
#' tolerance means the same thing whether or not
#' \code{use_jacobi_preconditioner} is set. The absolute floor is
#' intentional: the caller (\code{.fitSurvivalCG}) passes a gradient, so a
#' right-hand side this small means the Newton iteration has already
#' converged, and \code{.fitSurvivalCG} judges convergence by the change in
#' log-likelihood rather than by the step returned here. In exact arithmetic,
#' conjugate gradient converges within \code{nrow(coefficient_matrix)}
#' steps, but rounding error erodes that guarantee as the system grows, so
#' up to \code{10 * nrow(coefficient_matrix)} steps are allowed.
#'
#' @param coefficient_matrix symmetric positive (semi-)definite matrix,
#' e.g. the Hessian/information matrix from a Newton step.
#' @param right_hand_side vector the system is solved against, e.g. the
#' gradient/score vector from a Newton step.
#' @param use_jacobi_preconditioner if \code{TRUE}, precondition with the
#' inverse of \code{coefficient_matrix}'s own diagonal - cheap to apply,
#' and often enough to cut down the number of iterations needed when the
#' diagonal dominates (as it typically does for an AFT information matrix,
#' where each parameter's own curvature tends to be much larger than its
#' cross-terms with the other parameters). Defaults to \code{FALSE}, which
#' reduces exactly to plain (unpreconditioned) conjugate gradient.
#'
#' @return a list with: \code{solution}, the numeric vector solving
#' (approximately) \code{coefficient_matrix \%*\% solution =
#' right_hand_side}; \code{iterations}, how many conjugate-gradient steps
#' were actually taken; \code{converged}, whether the residual tolerance
#' was met; and \code{positive_definite}, whether \code{coefficient_matrix}
#' behaved as positive definite throughout (a caller can fall back to a
#' different matrix, e.g. a Gauss-Newton approximation, when this is
#' \code{FALSE}).
#'
#' @keywords internal
#' @noRd
.cgSolve = function(coefficient_matrix, right_hand_side,
                     use_jacobi_preconditioner = FALSE) {
    number_of_unknowns = nrow(coefficient_matrix)
    relative_tolerance = 1e-8
    max_iterations = 10 * number_of_unknowns
    solution = rep(0, number_of_unknowns)

    apply_preconditioner = if (use_jacobi_preconditioner) {
        diagonal = diag(coefficient_matrix)
        inverse_diagonal = ifelse(
            is.finite(diagonal) & diagonal > 0, 1 / diagonal, 1)
        function(vector) inverse_diagonal * vector
    } else {
        identity
    }

    residual = right_hand_side - drop(coefficient_matrix %*% solution)
    preconditioned_residual = apply_preconditioner(residual)
    search_direction = preconditioned_residual
    residual_size = sum(residual * residual)
    residual_dot_preconditioned_residual =
        sum(residual * preconditioned_residual)

    convergence_threshold =
        (relative_tolerance * max(sqrt(sum(right_hand_side^2)), 1))^2

    positive_definite = TRUE
    iterations_used = 0

    for (iteration in seq_len(max_iterations)) {
        if (residual_size <= convergence_threshold) {
            break
        }
        iterations_used = iteration
        matrix_times_search_direction =
            drop(coefficient_matrix %*% search_direction)
        curvature = sum(search_direction * matrix_times_search_direction)
        if (!is.finite(curvature) || curvature <= 0) {
            positive_definite = FALSE
            warning(".cgSolve: coefficient_matrix is not positive definite ",
                    "along the current search direction; returning the ",
                    "best iterate found so far")
            break
        }

        step_length = residual_dot_preconditioned_residual / curvature
        solution = solution + step_length * search_direction
        residual = residual - step_length * matrix_times_search_direction
        new_residual_size = sum(residual * residual)

        new_preconditioned_residual = apply_preconditioner(residual)
        new_residual_dot_preconditioned_residual =
            sum(residual * new_preconditioned_residual)
        search_direction = new_preconditioned_residual +
            (new_residual_dot_preconditioned_residual /
                 residual_dot_preconditioned_residual) * search_direction
        residual_size = new_residual_size
        residual_dot_preconditioned_residual =
            new_residual_dot_preconditioned_residual
    }

    converged = residual_size <= convergence_threshold
    if (!converged && positive_definite) {
        warning(".cgSolve: did not converge within max_iterations = ",
                max_iterations, " iterations; returning the best iterate ",
                "found so far")
    }

    list(solution = solution, iterations = iterations_used,
        converged = converged, positive_definite = positive_definite)
}
