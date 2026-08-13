#' Solve a symmetric positive (semi-)definite linear system via conjugate
#' gradient
#'
#' A minimal, single right-hand-side conjugate gradient solver, used as the
#' Newton-step linear solve in \code{.fitSurvivalCG}. Modeled on
#' \code{lfe::cgsolve}, stripped down to a single dense matrix and a single
#' right-hand-side vector (no multi-column batching, no \code{Matrix}-package
#' or operator/closure dispatch, no preconditioning - none of which are
#' needed for the small, dense AFT information matrices this is used on).
#'
#' @param coefficient_matrix symmetric positive (semi-)definite matrix,
#' e.g. the Hessian/information matrix from a Newton step.
#' @param right_hand_side vector the system is solved against, e.g. the
#' gradient/score vector from a Newton step.
#' @param initial_guess optional starting point for the iteration. Defaults
#' to the zero vector.
#' @param relative_tolerance how small the residual needs to shrink,
#' relative to the size of \code{right_hand_side}, before iteration stops.
#' @param max_iterations how many conjugate-gradient steps to try before
#' giving up. In exact arithmetic, conjugate gradient converges within
#' \code{nrow(coefficient_matrix)} steps, but rounding error erodes that
#' guarantee as the system grows, so the default allows for several times
#' that many steps.
#'
#' @return numeric vector solving (approximately)
#' \code{coefficient_matrix \%*\% solution = right_hand_side}.
#'
#' @keywords internal
.cgSolve = function(coefficient_matrix, right_hand_side, initial_guess = NULL,
                     relative_tolerance = 1e-8,
                     max_iterations = 10 * nrow(coefficient_matrix)) {
    number_of_unknowns = nrow(coefficient_matrix)
    solution = if (is.null(initial_guess)) {
        rep(0, number_of_unknowns)
    } else {
        initial_guess
    }

    # The residual measures how far the current guess is from solving the
    # system. Conjugate gradient starts out searching in that direction.
    residual = right_hand_side - drop(coefficient_matrix %*% solution)
    search_direction = residual
    residual_size = sum(residual * residual)
    smallest_residual_size_seen = residual_size

    # Stop once the residual has shrunk far enough, relative to the size of
    # the right-hand side (falling back to an absolute scale when that size
    # is tiny).
    convergence_threshold =
        (relative_tolerance * max(sqrt(sum(right_hand_side^2)), 1))^2

    for (iteration in seq_len(max_iterations)) {
        if (residual_size <= convergence_threshold) {
            break
        }

        # How far moving along the search direction changes things, as
        # measured through the matrix itself.
        matrix_times_search_direction =
            drop(coefficient_matrix %*% search_direction)
        curvature = sum(search_direction * matrix_times_search_direction)
        if (!is.finite(curvature) || curvature <= 0) {
            warning(".cgSolve: coefficient_matrix is not positive definite ",
                    "along the current search direction; returning the ",
                    "best iterate found so far")
            break
        }

        # Move as far as possible along the search direction without
        # overshooting the solution, then see how much residual remains.
        step_length = residual_size / curvature
        solution = solution + step_length * search_direction
        residual = residual - step_length * matrix_times_search_direction
        new_residual_size = sum(residual * residual)
        smallest_residual_size_seen =
            min(smallest_residual_size_seen, new_residual_size)

        # If the residual has grown far past its best value so far, the
        # iteration is diverging (e.g. because coefficient_matrix is
        # ill-conditioned) - give up and return what we have rather than
        # loop until max_iterations.
        if (iteration > 10 &&
            new_residual_size > 1e4 * smallest_residual_size_seen) {
            warning(".cgSolve: residual is diverging; returning the best ",
                    "iterate found so far")
            break
        }

        # Choose the next search direction so it doesn't undo the progress
        # made by earlier directions.
        search_direction = residual +
            (new_residual_size / residual_size) * search_direction
        residual_size = new_residual_size
    }

    if (residual_size > convergence_threshold) {
        warning(".cgSolve: did not converge within max_iterations = ",
                max_iterations, " iterations; returning the best iterate ",
                "found so far")
    }
    solution
}
