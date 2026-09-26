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
