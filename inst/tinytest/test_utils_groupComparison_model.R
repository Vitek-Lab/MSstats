# Test for .fitModelForGroupComparison function
test_input_with_variance = data.frame(
    ABUNDANCE = c(1.2, 1.5, 1.8, 2.1, 1.9, 2.3),
    GROUP = factor(c("A", "A", "A", "B", "B", "B")),
    SUBJECT = factor(c("S1", "S2", "S3", "S1", "S2", "S3")),
    Variance = c(0.1, 0.15, 0.08, 0.12, 0.09, 0.11)
)
result = MSstats:::.fitModelForGroupComparison(
    input = test_input_with_variance,
    repeated = FALSE,
    is_single_subject = FALSE,
    has_tech_replicates = FALSE
)
expect_false(is.null(result$full_fit$weights),
             "Fitted model should contain weights when Variance is provided")
expected_weights = 1/test_input_with_variance$Variance
expect_equal(result$full_fit$weights, expected_weights,
             info = "Model weights should match 1/Variance")
