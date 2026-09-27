make_random_solvable_matrix <- function(size, seed, diagonal_boost = 0.01) {
    # Builds a random symmetric positive-definite matrix
    set.seed(seed)
    random_matrix <- matrix(rnorm(size * size), size, size)
    random_matrix %*% t(random_matrix) + diag(size) * diagonal_boost
}

for (size in c(2, 5, 10, 30, 80)) {
    coefficient_matrix <- make_random_solvable_matrix(size, seed = size)
    set.seed(size + 1000)
    right_hand_side <- rnorm(size)

    iterative_result <- MSstats:::.cgSolve(coefficient_matrix, right_hand_side)
    exact_solution <- solve(coefficient_matrix, right_hand_side)

    expect_equal(
        iterative_result$solution, exact_solution, tolerance = 1e-6,
        info = paste0(".cgSolve should give the same answer as solve() on a ",
                      "random solvable symmetric system of size ", size)
    )
    expect_true(
        iterative_result$converged && iterative_result$positive_definite,
        info = paste0("A well-behaved symmetric system of size ", size,
                      " should report converged = TRUE and ",
                      "positive_definite = TRUE")
    )
    expect_true(
        iterative_result$iterations >= 1 &&
            iterative_result$iterations <= size * 10,
        info = paste("The number of steps taken should be at least one and",
                     "should not go over the maximum allowed")
    )
}

make_random_unsolvable_matrix <- function(size, seed) {
    unsolvable_matrix <- make_random_solvable_matrix(size, seed)
    unsolvable_matrix[1, ] <- 0
    unsolvable_matrix[, 1] <- 0
    unsolvable_matrix
}

unsolvable_matrix <- make_random_unsolvable_matrix(10, seed = 42)
set.seed(43)
right_hand_side <- rnorm(10)

expect_warning(
    unsolvable_result <- MSstats:::.cgSolve(unsolvable_matrix, right_hand_side),
    info = paste("An unsolvable matrix should produce a warning rather than",
                 "an error or an endless loop")
)
expect_true(
    all(is.finite(unsolvable_result$solution)),
    info = paste("An unsolvable system should still return a partial answer",
                 "made of ordinary finite numbers")
)
expect_false(
    unsolvable_result$converged && unsolvable_result$positive_definite,
    info = paste("An unsolvable matrix should be flagged by reporting",
                 "converged = FALSE, positive_definite = FALSE, or both")
)

make_large_diagonal_matrix <- function(size, seed) {
    set.seed(seed)
    small_off_diagonal_entries <-
        matrix(runif(size * size, -0.1, 0.1), size, size)
    small_off_diagonal_entries <-
        (small_off_diagonal_entries + t(small_off_diagonal_entries)) / 2
    diag(small_off_diagonal_entries) <- 0
    diag(size) * runif(size, 5, 10) + small_off_diagonal_entries
}

large_diagonal_matrix <- make_large_diagonal_matrix(40, seed = 11)
set.seed(12)
large_diagonal_right_hand_side <- rnorm(40)
large_diagonal_exact_solution <-
    solve(large_diagonal_matrix, large_diagonal_right_hand_side)

result_without_scaling <- MSstats:::.cgSolve(
    large_diagonal_matrix, large_diagonal_right_hand_side)
result_with_scaling <- MSstats:::.cgSolve(
    large_diagonal_matrix, large_diagonal_right_hand_side,
    use_jacobi_preconditioner = TRUE)

expect_equal(
    result_with_scaling$solution, large_diagonal_exact_solution,
    tolerance = 1e-6,
    info = paste("With diagonal scaling turned on, .cgSolve should still",
                 "give the same answer as solve()")
)
expect_true(
    result_with_scaling$iterations <= result_without_scaling$iterations,
    info = paste("Diagonal scaling should not need more steps than",
                 "no scaling when the diagonal entries are large (without",
                 "scaling =", result_without_scaling$iterations,
                 ", with scaling =", result_with_scaling$iterations, ")")
)

zero_diagonal_matrix <- make_random_solvable_matrix(8, seed = 55)
zero_diagonal_matrix[3, 3] <- 0
set.seed(56)
zero_diagonal_right_hand_side <- rnorm(8)
expect_true(
    all(is.finite(suppressWarnings(MSstats:::.cgSolve(
        zero_diagonal_matrix, zero_diagonal_right_hand_side,
        use_jacobi_preconditioner = TRUE))$solution)),
    info = paste("A zero on the diagonal should not cause diagonal scaling",
                 "to return infinite or missing values")
)
