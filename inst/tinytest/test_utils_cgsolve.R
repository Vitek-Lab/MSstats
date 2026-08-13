# Tests for .cgSolve(), the vendored conjugate-gradient linear solver used
# as the Newton-step solve in .fitSurvivalCG().

make_random_spd_matrix <- function(size, seed, ridge = 0.01) {
    set.seed(seed)
    random_factor <- matrix(rnorm(size * size), size, size)
    random_factor %*% t(random_factor) + diag(size) * ridge
}

for (size in c(2, 5, 10, 30, 80)) {
    coefficient_matrix <- make_random_spd_matrix(size, seed = size)
    set.seed(size + 1000)
    right_hand_side <- rnorm(size)

    cg_solution <- MSstats:::.cgSolve(coefficient_matrix, right_hand_side)
    exact_solution <- solve(coefficient_matrix, right_hand_side)

    expect_equal(
        cg_solution, exact_solution, tolerance = 1e-6,
        info = paste0(".cgSolve should match solve() on a random SPD ",
                      "system of size ", size)
    )
}

# --- near-singular system: still returns a finite result, with a warning ---

near_singular_matrix <- make_random_spd_matrix(10, seed = 42)
near_singular_matrix[1, ] <- 0
near_singular_matrix[, 1] <- 0
set.seed(43)
right_hand_side <- rnorm(10)

expect_warning(
    solution <- MSstats:::.cgSolve(near_singular_matrix, right_hand_side),
    info = paste("A singular coefficient_matrix should warn rather than",
                "error or hang")
)
expect_true(
    all(is.finite(solution)),
    info = "A singular system should still return a finite (partial) solution"
)

# --- an initial guess that is already the solution converges immediately ---

exact_matrix <- make_random_spd_matrix(6, seed = 7)
set.seed(8)
exact_rhs <- rnorm(6)
exact_answer <- solve(exact_matrix, exact_rhs)

solution_from_exact_start <- MSstats:::.cgSolve(
    exact_matrix, exact_rhs, initial_guess = exact_answer)
expect_equal(
    solution_from_exact_start, exact_answer, tolerance = 1e-8,
    info = "Starting from the exact solution should return it unchanged"
)

# --- relative_tolerance controls how tightly the system is solved ---

loose_matrix <- make_random_spd_matrix(20, seed = 99)
set.seed(100)
loose_rhs <- rnorm(20)
exact_loose_answer <- solve(loose_matrix, loose_rhs)

loose_solution <- MSstats:::.cgSolve(loose_matrix, loose_rhs, relative_tolerance = 1e-2)
tight_solution <- MSstats:::.cgSolve(loose_matrix, loose_rhs, relative_tolerance = 1e-10)

loose_error <- max(abs(loose_solution - exact_loose_answer))
tight_error <- max(abs(tight_solution - exact_loose_answer))
expect_true(
    tight_error < loose_error,
    info = "A tighter relative_tolerance should produce a more accurate solution"
)
