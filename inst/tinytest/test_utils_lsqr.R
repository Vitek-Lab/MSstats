# Setup ----------------------------------------------------------------------
# .buildRunFeatureDesign/.fitRobustRunFeatureModel model a "1 + run + feature"
# design, the same formula .fitHuber fits via MASS::rlm(log2inty ~ run +
# feature). Ground truth for the design itself is an explicit
# stats::model.matrix() build of that same design.
set.seed(20260815)
n = 40
run = sample(paste0("run", 1:5), n, replace = TRUE)
feature = sample(paste0("feat", 1:4), n, replace = TRUE)
design_df = data.frame(run = run, feature = feature)
X = stats::model.matrix(~ run + feature, design_df)

# Test 1: .buildRunFeatureDesign matches an explicit model.matrix build ------
design = MSstats:::.buildRunFeatureDesign(run, feature)
expect_equal(design$num_coefficients, ncol(X))

set.seed(1)
v = rnorm(design$num_coefficients)
u = rnorm(n)
expect_equal(design$multiply(v), as.vector(X %*% v))
expect_equal(design$multiply_transposed(u), as.vector(crossprod(X, u)))

# Test 2: single-level run/feature collapses to intercept-only design -------
design_single = MSstats:::.buildRunFeatureDesign(rep("onlyrun", 10), rep("onlyfeat", 10))
expect_equal(design_single$num_coefficients, 1L)
expect_equal(design_single$multiply(3), rep(3, 10))

# Test 3: .solveLeastSquares matches direct solvers on a well-conditioned system
set.seed(2)
A = matrix(rnorm(30 * 5), 30, 5)
b = rnorm(30)
multiply = function(v) as.vector(A %*% v)
multiply_transposed = function(u) as.vector(crossprod(A, u))
lsqr_fit = MSstats:::.solveLeastSquares(multiply, multiply_transposed, b, ncol(A))
direct_fit = qr.solve(A, b)
expect_equal(lsqr_fit$solution, direct_fit, tolerance = 1e-6)
expect_true(lsqr_fit$stop_reason %in% c("exact solution found",
                                        "least squares solution found"))

# Test 4: .solveLeastSquares degrades gracefully (minimum-norm solution, no
# error) on a rank-deficient system, unlike MASS::rlm's explicit singularity
# check -----------------------------------------------------------------------
A_deficient = cbind(A, A[, 1] + A[, 2])
multiply_d = function(v) as.vector(A_deficient %*% v)
multiply_transposed_d = function(u) as.vector(crossprod(A_deficient, u))
lsqr_deficient = MSstats:::.solveLeastSquares(multiply_d, multiply_transposed_d,
                                              b, ncol(A_deficient))
fitted_lsqr = as.vector(A_deficient %*% lsqr_deficient$solution)
fitted_lm = fitted(lm.fit(A_deficient, b))
expect_equal(fitted_lsqr, fitted_lm, tolerance = 1e-6)

# Test 5: .fitRobustRunFeatureModel agrees with MASS::rlm(scale.est = "Huber")
set.seed(3)
run2 = sample(paste0("run", 1:6), 60, replace = TRUE)
feature2 = sample(paste0("feat", 1:5), 60, replace = TRUE)
run_effect = c(0, 0.5, -0.3, 0.8, -0.6, 0.2)[as.integer(factor(run2))]
feature_effect = c(0, 0.4, -0.2, 0.3, -0.5)[as.integer(factor(feature2))]
y2 = 10 + run_effect + feature_effect + rnorm(60, sd = 0.2)
y2[1:3] = y2[1:3] + 5 # a few outliers, exercising the Huber downweighting

fit_mass = MASS::rlm(y2 ~ run2 + feature2, scale.est = "Huber")
fit_lsqr = MSstats:::.fitRobustRunFeatureModel(y2, run2, feature2)

expect_true(fit_lsqr$converged)
expect_equal(unname(coef(fit_mass)), unname(fit_lsqr$coefficients), tolerance = 1e-4)
expect_equal(fit_mass$s, fit_lsqr$scale, tolerance = 1e-4)
expect_equal(unname(residuals(fit_mass)), unname(fit_lsqr$residuals), tolerance = 1e-4)
expect_equal(fit_lsqr$df.residual, nrow(model.matrix(fit_mass)) - ncol(model.matrix(fit_mass)))

# Test 6: a run/feature level with zero surviving observations (e.g. all its
# rows have NA log2inty for this protein) is handled gracefully, matching how
# .fitHuberLSQR calls .fitRobustRunFeatureModel with levels fixed from the
# full column ------------------------------------------------------------------
na_rows = which(run2 == "run6")
y3 = y2
y3[na_rows] = NA
keep = !is.na(y3)
fit_edge = MSstats:::.fitRobustRunFeatureModel(y3[keep], run2[keep], feature2[keep],
                                               run_levels = sort(unique(run2)),
                                               feature_levels = sort(unique(feature2)))
expect_true(fit_edge$converged)
expect_false(anyNA(fit_edge$coefficients))
full_design_df = data.frame(run2 = run2, feature2 = feature2)
expect_equal(fit_edge$rank, ncol(model.matrix(~ run2 + feature2, full_design_df)))
