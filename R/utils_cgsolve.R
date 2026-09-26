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
#' @section The math behind conjugate gradient:
#' Notation (with the variable each symbol lives in): \code{A} is
#' \code{coefficient_matrix} (n x n, symmetric positive definite), \code{b}
#' is \code{right_hand_side}, \code{x_k} is \code{solution} after k steps,
#' \code{r_k = b - A x_k} is \code{residual}, \code{M^-1} is the
#' preconditioner (\code{diag(A)^-1} for Jacobi, the identity otherwise),
#' \code{z_k = M^-1 r_k} is \code{preconditioned_residual}, and \code{p_k}
#' is \code{search_direction}. A prime (') denotes transpose, so
#' \code{u'v} is a dot product, computed here as \code{sum(u * v)}.
#'
#' \strong{1. Solving Ax = b as a minimization.} Because \code{A} is
#' symmetric positive definite, the quadratic
#' \preformatted{
#'     phi(x) = (1/2) x'A x - b'x
#' }
#' is a bowl with a unique minimum. Its gradient is
#' \code{grad phi(x) = A x - b = -r}, so setting the derivative to zero
#' gives exactly \code{A x = b}: the minimizer of \code{phi} is the
#' solution, and the residual is the negative gradient (the steepest
#' downhill direction).
#'
#' \strong{2. One iteration, written as matrix operations.} Starting from
#' \code{x_0 = 0}, \code{r_0 = b}, \code{z_0 = M^-1 r_0},
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
#'     beta_k  = (r_k+1' z_k+1) / (r_k' z_k)
#'     p_k+1   = z_k+1 + beta_k p_k
#' }
#' In the code, \code{q_k} is \code{matrix_times_search_direction},
#' \code{p_k' q_k} is \code{curvature}, \code{alpha_k} is
#' \code{step_length}, and \code{r_k' z_k} is
#' \code{residual_dot_preconditioned_residual}. The per-step cost is one
#' n x n matrix-vector product (O(n^2)) plus a handful of O(n) dot products
#' and vector updates, and only the current \code{x}, \code{r}, \code{z},
#' and \code{p} are kept - no earlier directions need to be stored.
#'
#' \strong{3. Where alpha comes from: residuals orthogonal to the step
#' taken.} Given a direction \code{p_k}, pick the step that minimizes
#' \code{phi} along that line by setting the derivative to zero:
#' \preformatted{
#'     d/d alpha  phi(x_k + alpha p_k)
#'         = p_k' (A x_k + alpha A p_k - b)
#'         = -p_k' r_k + alpha p_k' A p_k  = 0
#'     =>  alpha_k = (p_k' r_k) / (p_k' A p_k)
#' }
#' The zero-derivative condition is the same statement as
#' \code{p_k' r_k+1 = 0}: the new residual (the new downhill direction) is
#' orthogonal to the direction just searched, so there is nothing left to
#' gain along \code{p_k}. Since \code{p_k = z_k + beta_k-1 p_k-1} and
#' \code{p_k-1' r_k = 0} by the same argument one step earlier,
#' \code{p_k' r_k = z_k' r_k}, which is the numerator used in the code.
#'
#' \strong{4. Where beta comes from: A-conjugate directions.} Steepest
#' descent (always stepping along \code{r_k}) also does an exact line
#' search, but a later step can undo progress made along an earlier
#' direction, so it zig-zags across narrow valleys. Conjugate gradient
#' avoids this by choosing each new direction to be \emph{A-conjugate} to
#' the previous one:
#' \preformatted{
#'     p_k+1' A p_k = 0
#'     =>  (z_k+1 + beta_k p_k)' A p_k = 0
#'     =>  beta_k = -(z_k+1' A p_k) / (p_k' A p_k)
#' }
#' Substituting \code{A p_k = (r_k - r_k+1) / alpha_k} (from the residual
#' update) together with the orthogonality \code{z_k' r_k+1 = 0} collapses
#' this to \code{beta_k = (r_k+1' z_k+1) / (r_k' z_k)}, a ratio of two
#' dot products that are already on hand. \code{p_i' A p_j = 0} says the
#' directions are orthogonal under the inner product
#' \code{<u, v>_A = u' A v} - i.e. they would be ordinary perpendicular
#' vectors if space were stretched by \code{A} so the elliptical contours
#' of \code{phi} became circles. By induction the short recurrence gives
#' this for all pairs, not just neighbours:
#' \preformatted{
#'     p_i' A p_j   = 0   for all i != j     (directions A-conjugate)
#'     r_i' M^-1 r_j = 0   for all i != j     (residuals orthogonal; plain
#'                                             r_i' r_j = 0 when M = I)
#' }
#'
#' \strong{5. Verifying optimality with derivatives, and why this is
#' fast.} After k steps \code{x_k = alpha_0 p_0 + ... + alpha_k-1 p_k-1}.
#' Take the derivative of \code{phi} with respect to the coefficient on
#' any earlier direction \code{p_j} (j < k):
#' \preformatted{
#'     d phi / d c_j = p_j' (A x_k - b) = -p_j' r_k = 0
#' }
#' Every one of these partial derivatives is zero, so \code{x_k} is not
#' just the best point on the latest line but the exact minimizer of
#' \code{phi} over the whole subspace \code{span(p_0, ..., p_k-1)} - which
#' equals the Krylov subspace \code{span(r_0, A r_0, ..., A^(k-1) r_0)}
#' (with \code{M^-1} folded in when preconditioned). Conjugacy is what
#' makes this possible: writing the exact answer as
#' \code{x* = sum_j alpha_j p_j} and multiplying \code{A x* = b} on the
#' left by \code{p_j'} kills every cross term, leaving
#' \code{alpha_j = (p_j' b) / (p_j' A p_j)} - each coefficient is
#' determined independently, so each direction is solved once and never
#' revisited. With n independent directions available, exact arithmetic
#' reaches the solution in at most n steps.
#'
#' \strong{6. Why clustered eigenvalues make it faster still.} Because
#' \code{x_k} is optimal over the Krylov subspace, its error
#' \code{e_k = x* - x_k} can be written \code{e_k = P_k(A) e_0} for the
#' degree-k polynomial \code{P_k} with \code{P_k(0) = 1} that minimizes the
#' A-norm of the error. Expanding \code{e_0 = sum_i c_i v_i} in the
#' eigenvectors \code{v_i} of \code{A} (eigenvalues \code{lambda_i}):
#' \preformatted{
#'     ||e_k||_A^2 = sum_i lambda_i c_i^2 P_k(lambda_i)^2
#'     ||e_k||_A / ||e_0||_A <= min over P_k of  max_i |P_k(lambda_i)|
#' }
#' so convergence depends on how many \emph{distinct groups} of
#' eigenvalues there are, not on n. If the eigenvalues fall into m tight
#' clusters, a degree-m polynomial with one root near the centre of each
#' cluster is small at every eigenvalue, and CG is essentially done after
#' about m steps. Intuitively, within a cluster \code{A} acts almost like a
#' single scalar, so the curvature \code{p' A p / p' p} - and hence the
#' step length \code{alpha} - is nearly the same for every eigen-direction
#' in that cluster; one step of that size removes the error along all of
#' those directions at once instead of one at a time (in the extreme case
#' \code{A = c I}, a single step solves the system exactly). When the
#' spectrum is spread out instead, the standard bound
#' \preformatted{
#'     ||e_k||_A <= 2 ((sqrt(kappa) - 1) / (sqrt(kappa) + 1))^k ||e_0||_A,
#'     kappa = lambda_max / lambda_min
#' }
#' applies. This is the motivation for the Jacobi preconditioner: CG is
#' effectively run on \code{M^-1 A}, and for a diagonally dominant AFT
#' information matrix dividing by the diagonal pulls the eigenvalues
#' toward 1 - clustering them and shrinking \code{kappa}.
#'
#' Finally, the curvature \code{p_k' A p_k} must be positive for
#' \code{phi} to have a minimum along \code{p_k}; if it is zero, negative,
#' or non-finite, \code{A} is not positive definite along that direction,
#' the line search has no minimizer, and the loop stops and reports
#' \code{positive_definite = FALSE}.
#'
#' @section Computational cost:
#' In this section, n is the number of entries in
#' \code{coefficient_matrix} - for m unknowns, \code{n = m^2} - not the
#' number of unknowns used in the sections above.
#'
#' \strong{A solve costs O(k * n): linear in the size of the matrix.} Each
#' iteration does exactly one product \code{A p}. Computing
#' \code{coefficient_matrix \%*\% search_direction} visits every entry of
#' \code{A} once (one multiply and one add each), so it costs O(n). The
#' rest of the iteration is a fixed number of length-m dot products and
#' vector updates, O(m) = O(sqrt(n)), which the product dominates. Over k
#' iterations the total is
#' \preformatted{
#'     O(k * (n + m)) = O(k * n)
#' }
#' \code{A} is never factorized or modified; it is only ever read, once per
#' iteration. By comparison, a Cholesky factorization costs O(m^3) =
#' O(n^1.5) regardless of how quickly the problem could converge, so CG
#' wins whenever k is small relative to m.
#'
#' \strong{k stays small on MSstats data.} Section 6 above shows k is
#' governed by the number of eigenvalue clusters, not by the number of
#' unknowns. The AFT information matrix is dominated by its diagonal (each
#' parameter's own curvature is much larger than its cross-terms), so
#' Jacobi preconditioning pulls most of the spectrum of \code{M^-1 A} to
#' 1, with only a few outlying eigenvalues. For example, in a simulated
#' 12-feature x 10-run protein (m = 22 unknowns, 20\% censored), 14 of the
#' 22 preconditioned eigenvalues lie in [0.9, 1.1], the condition number
#' drops from about 234 to about 103, and CG reaches tolerance in 12
#' iterations instead of 16. Because the bulk cluster is handled in a few
#' steps, and the remaining iterations are spent on the few outliers, k
#' grows much more slowly than m.
#'
#' \strong{Symmetry and sparsity: the real work is in the upper-right
#' block.} For the \code{~ FEATURE + RUN} model, every row of the design
#' matrix \code{X} has exactly one feature indicator and one run indicator.
#' In \code{A = -X' W X} (plus the log-scale row and column), this means:
#' \itemize{
#'   \item the feature-feature block is \emph{diagonal} - two different
#'   features never appear in the same row, so their cross-term is zero;
#'   \item the run-run block is likewise \emph{diagonal};
#'   \item the only dense off-diagonal content is the feature x run block
#'   \code{C} in the upper-right corner (one entry per feature/run cell
#'   that has data), together with the intercept and log-scale border
#'   rows;
#'   \item by symmetry the lower-left block is \code{C'}, so it carries no
#'   new information.
#' }
#' Ignoring the border, the product therefore reduces to
#' \preformatted{
#'     A = [ D_F   C  ]        A p = [ D_F p_F + C  p_R ]
#'         [ C'   D_R ]              [ C' p_F + D_R p_R ]
#' }
#' with \code{D_F} and \code{D_R} diagonal. The only real computation in
#' each product is two elementwise scalings plus one pass over \code{C},
#' used once as-is and once transposed. Everything else in the n entries is
#' either zero or a mirror image of \code{C}, so the number of distinct
#' nonzero values is about \code{F + R + (number of feature/run cells)},
#' well below n. The dense product still visits all n entries, which is
#' what the O(k * n) bound counts. A product written against this
#' structure could skip the zeros and reuse \code{C} for both halves, but
#' at per-protein sizes the dense product is already cheap.
#'
#' The same structure is why the Jacobi preconditioner works well here.
#' The diagonal blocks are exactly diagonal, so scaling by \code{diag(A)}
#' turns them into identity blocks, and the preconditioned matrix is the
#' identity plus only the scaled \code{C} coupling (and the border). Its
#' eigenvalues therefore sit near 1, spread only as far as that coupling
#' pushes them - the clustering that keeps k small.
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
