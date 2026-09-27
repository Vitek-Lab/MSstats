#' Solve a system of linear equations Ax = b with the conjugate gradient method.
#' 
#' @section Conjugate gradient is similar to gradient descent, except with
#' how it picks its search direction.
#' 
#' Conjugate gradient is an iterative method that solves for "x" in Ax=b. 
#' The method reformulates the linear system as an optimization problem where 
#' one attempts to minimize Ax^2 - bx. Similar to gradient descent, the solution 
#' is initialized at some arbitrary point, then the method iteratively descends 
#' toward the optimal solution.  But as opposed to gradient descent, which 
#' moves in the direction of steepest descent, conjugate gradient picks a 
#' different search direction.
#' 
#' @section Conjugate gradient search direction is determined based on a linear
#' transformation of the previous iteration's search direction.
#' 
#' \strong{Pseudocode for one iteration, written as matrix operations.} 
#' Starting from \code{x_0 = 0}, \code{r_0 = b}, \code{z_0 = M^-1 r_0},
#' \code{p_0 = z_0}, each step computes
#' \preformatted{
#'     q_k     = A p_k                        (the only matrix-vector
#'                                             product per iteration)
#'     alpha_k = (r_k' z_k) / (p_k' q_k)      step length; p_k' A p_k is
#'                                             the "curvature"
#'     x_k+1   = x_k + alpha_k p_k
#'     r_k+1   = r_k - alpha_k q_k            (= b - A x_k+1, updated
#'                                             without a new product)
#'     z_k+1   = M^-1 r_k+1                   (elementwise scaling for
#'                                             Jacobi)
#'     beta_k  = (r_k+1' z_k+1) / (r_k' z_k)  (guarantees p_k A p_k+1 = 0,
#'                                             i.e. A-conjugacy)
#'     p_k+1   = z_k+1 + beta_k p_k
#' }
#' 
#' @section If A is a (p x p) matrix, conjugate gradient is guaranteed to 
#' converge in p iterations
#' 
#' \strong{Residuals are orthogonal to the step taken.} Given a direction 
#' \code{p_k}, CG picks the step size that minimizes the objective along the
#' search direction by setting the derivative to zero with respect to alpha:
#' \preformatted{
#'     d/d alpha  phi(x_k + alpha p_k)
#'         = p_k' (A x_k + alpha A p_k - b)
#'         = -p_k' r_k + alpha p_k' A p_k  = 0
#'     =>  alpha_k = (p_k' r_k) / (p_k' A p_k)
#'     
#'     x_k+1 = x_k + alpha_k p_k
#'     r_k+1 = b - A x_k+1 = r_k - alpha_k A p_k
#'     p_k' r_k+1 = p_k' r_k - alpha_k p_k' A p_k
#'                = p_k' r_k - (p_k' r_k / p_k' A p_k) p_k' A p_k
#'                = p_k' r_k - p_k' r_k
#'                = 0
#' }
#'
#' \strong{Residuals are orthogonal to all other residuals.} 
#' 
#' Dotting \code{r_k+1 = r_k - alpha_k A p_k} with \code{p_k-1} instead:
#' 
#' \preformatted{
#'     p_k-1' r_k = 0
#'     p_k-1' A p_k = 0 (by design of beta to ensure A-conjugacy)
#'     p_k-1' r_k+1 = p_k-1' r_k - alpha_k p_k-1' A p_k = 0
#'     
#'     r_k = p_k - beta_k-1 p_k-1
#'     r_k' r_k+1 = (p_k - beta_k-1 p_k-1)' r_k+1
#'                = p_k' r_k+1 - beta_k-1 (p_k-1' r_k+1)
#'                = 0 - beta_k-1 (0)
#'                = 0
#' }
#' 
#' So each new residual is orthogonal to the previous one as well.  Since in
#' a p-dimensional system, there can only be p orthogonal residuals, CG
#' is guaranteed to converge in p-iterations.
#' 
#' @section Time complexity is empirically linear with respect to the number
#' of entries in the coefficient matrix
#' 
#' The time complexity of CG is O(kn), where n is the number of entries in 
#' matrix A, and k is the number of iterations.  As shown earlier, k is 
#' guaranteed to be at most the number of rows in matrix A, i.e. n^0.5.
#' 
#' But we can ensure k is small and near a constant given certain data
#' constraints.  In the case of MSstats, the following intuition leads
#' k to be empirically small:
#' 
#' 1. MSstats Hessian matrix contains zeros on 
#'    (run x run) or (feature x feature) entries
#' 2. Sparsity makes it more likely diagonals of a matrix dominate
#' 3. Because diagonals more likely dominate, scaling the matrix by the inverse 
#' of the diagonals (Jacobi preconditioning) makes the eigenvalues of the 
#' resultant matrix to cluster well
#' 4. Well-clustered eigenvalues reduces the number of iterations needed
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
