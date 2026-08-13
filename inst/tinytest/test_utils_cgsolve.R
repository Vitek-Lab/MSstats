# Tests for .cgSolve(), the vendored conjugate-gradient linear solver used
# as the Newton-step solve in .fitSurvivalCG(). Returns a list:
# solution/iterations/converged/positive_definite.

make_random_spd_matrix <- function(size, seed, ridge = 0.01) {
    set.seed(seed)
    random_factor <- matrix(rnorm(size * size), size, size)
    random_factor %*% t(random_factor) + diag(size) * ridge
}

for (size in c(2, 5, 10, 30, 80)) {
    coefficient_matrix <- make_random_spd_matrix(size, seed = size)
    set.seed(size + 1000)
    right_hand_side <- rnorm(size)

    cg_result <- MSstats:::.cgSolve(coefficient_matrix, right_hand_side)
    exact_solution <- solve(coefficient_matrix, right_hand_side)

    expect_equal(
        cg_result$solution, exact_solution, tolerance = 1e-6,
        info = paste0(".cgSolve should match solve() on a random SPD ",
                      "system of size ", size)
    )
    expect_true(
        cg_result$converged && cg_result$positive_definite,
        info = paste0("A well-conditioned SPD system of size ", size,
                      " should report converged/positive_definite = TRUE")
    )
    expect_true(
        cg_result$iterations >= 1 && cg_result$iterations <= size * 10,
        info = "iterations should be a small positive count, not the default cap"
    )
}

# --- near-singular system: still returns a finite result, with a warning ---

near_singular_matrix <- make_random_spd_matrix(10, seed = 42)
near_singular_matrix[1, ] <- 0
near_singular_matrix[, 1] <- 0
set.seed(43)
right_hand_side <- rnorm(10)

expect_warning(
    singular_result <- MSstats:::.cgSolve(near_singular_matrix, right_hand_side),
    info = paste("A singular coefficient_matrix should warn rather than",
                "error or hang")
)
expect_true(
    all(is.finite(singular_result$solution)),
    info = "A singular system should still return a finite (partial) solution"
)
expect_false(
    singular_result$converged && singular_result$positive_definite,
    info = paste("A singular coefficient_matrix should signal trouble via",
                "converged = FALSE and/or positive_definite = FALSE")
)

# --- an initial guess that is already the solution converges immediately ---

exact_matrix <- make_random_spd_matrix(6, seed = 7)
set.seed(8)
exact_rhs <- rnorm(6)
exact_answer <- solve(exact_matrix, exact_rhs)

result_from_exact_start <- MSstats:::.cgSolve(
    exact_matrix, exact_rhs, initial_guess = exact_answer)
expect_equal(
    result_from_exact_start$solution, exact_answer, tolerance = 1e-8,
    info = "Starting from the exact solution should return it unchanged"
)
expect_equal(
    result_from_exact_start$iterations, 0,
    info = "Starting from the exact solution should take zero iterations"
)

# --- relative_tolerance controls how tightly the system is solved ---

loose_matrix <- make_random_spd_matrix(20, seed = 99)
set.seed(100)
loose_rhs <- rnorm(20)
exact_loose_answer <- solve(loose_matrix, loose_rhs)

loose_result <- MSstats:::.cgSolve(
    loose_matrix, loose_rhs, relative_tolerance = 1e-2)
tight_result <- MSstats:::.cgSolve(
    loose_matrix, loose_rhs, relative_tolerance = 1e-10)

loose_error <- max(abs(loose_result$solution - exact_loose_answer))
tight_error <- max(abs(tight_result$solution - exact_loose_answer))
expect_true(
    tight_error < loose_error,
    info = "A tighter relative_tolerance should produce a more accurate solution"
)

# --- Jacobi preconditioner: same answer, fewer or equal iterations --------
# on a diagonally-dominant system (where a diagonal preconditioner is most
# effective), preconditioned CG should converge in no more iterations than
# plain CG, and to the same solution.

make_diagonally_dominant_matrix <- function(size, seed) {
    set.seed(seed)
    matrix_off_diagonal <- matrix(runif(size * size, -0.1, 0.1), size, size)
    matrix_off_diagonal <- (matrix_off_diagonal + t(matrix_off_diagonal)) / 2
    diag(matrix_off_diagonal) <- 0
    diag(size) * runif(size, 5, 10) + matrix_off_diagonal
}

dominant_matrix <- make_diagonally_dominant_matrix(40, seed = 11)
set.seed(12)
dominant_rhs <- rnorm(40)
exact_dominant_answer <- solve(dominant_matrix, dominant_rhs)

plain_cg_result <- MSstats:::.cgSolve(dominant_matrix, dominant_rhs)
preconditioned_result <- MSstats:::.cgSolve(
    dominant_matrix, dominant_rhs, use_jacobi_preconditioner = TRUE)

expect_equal(
    preconditioned_result$solution, exact_dominant_answer, tolerance = 1e-6,
    info = "Preconditioned CG should still match solve() on a diagonally dominant system"
)
expect_true(
    preconditioned_result$iterations <= plain_cg_result$iterations,
    info = paste("Jacobi preconditioning should not need more iterations",
                "than plain CG on a diagonally dominant system (plain =",
                plain_cg_result$iterations, ", preconditioned =",
                preconditioned_result$iterations, ")")
)

# A degenerate (all-zero) diagonal entry should not blow up the
# preconditioner (falls back to an identity-like scale of 1 for that entry).
degenerate_diagonal_matrix <- make_random_spd_matrix(8, seed = 55)
degenerate_diagonal_matrix[3, 3] <- 0
set.seed(56)
degenerate_rhs <- rnorm(8)
expect_true(
    all(is.finite(suppressWarnings(MSstats:::.cgSolve(
        degenerate_diagonal_matrix, degenerate_rhs,
        use_jacobi_preconditioner = TRUE))$solution)),
    info = "A zero diagonal entry should not produce a non-finite preconditioned solution"
)
