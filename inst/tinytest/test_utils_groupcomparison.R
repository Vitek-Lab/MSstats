library(data.table)

# Test 1: Basic functionality with no missing values
# Setup
contrast_matrix <- matrix(c(1, -1, 0, 
                            0, 1, -1), nrow = 2, ncol = 3)
colnames(contrast_matrix) <- c("Group1", "Group2", "Group3")

summarized <- data.table(
    GROUP = c("Group1", "Group1", "Group2", "Group2", "Group3", "Group3"),
    TotalGroupMeasurements = c(100, 100, 100, 100, 100, 100),
    NumMeasuredFeature = c(50, 50, 50, 50, 50, 50),
    NumImputedFeature = c(0, 0, 0, 0, 0, 0)
)

result <- list()
samples_info <- data.table(GROUP = c("Group1", "Group2", "Group3"), NumRuns = c(10, 10, 10))

# Execute
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, FALSE)

# Assert
expect_equal(length(output$MissingPercentage), 2, info = "Basic functionality: MissingPercentage length")
expect_equal(output$MissingPercentage, c(0, 0), info = "Basic functionality: No missing values")
expect_true(is.null(output$ImputationPercentage), info = "Basic functionality: No imputation when has_imputed = FALSE")

# Test 2: With imputed values
# Setup
contrast_matrix <- matrix(c(1, -1), nrow = 1, ncol = 2)
colnames(contrast_matrix) <- c("Group1", "Group2")

summarized <- data.table(
    GROUP = c("Group1", "Group2"),
    TotalGroupMeasurements = c(100, 100),
    NumMeasuredFeature = c(80, 70),
    NumImputedFeature = c(10, 20)
)

result <- list()
samples_info <- data.table(GROUP = c("Group1", "Group2"), NumRuns = c(10, 10))

# Execute
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)

# Assert
expected_missing <- 1 - (80 + 70) / (100 + 100) # 0.25
expected_imputed <- (10 + 20) / (100 + 100) # 0.15

expect_equal(output$MissingPercentage[1], expected_missing, info = "Imputed values: Missing percentage calculation")
expect_equal(output$ImputationPercentage[1], expected_imputed, info = "Imputed values: Imputation percentage calculation")

# Test 3: With empty conditions (groups not in summarized data)
# Setup
contrast_matrix <- matrix(c(1, -1, 0), nrow = 1, ncol = 3)
colnames(contrast_matrix) <- c("Group1", "Group2", "Group3")

summarized <- data.table(
    GROUP = c("Group1"),
    TotalGroupMeasurements = c(100),
    NumMeasuredFeature = c(80),
    NumImputedFeature = c(0)
)

result <- list()
samples_info <- data.table(GROUP = c("Group1", "Group2", "Group3"), NumRuns = c(10, 10, 10))

# Execute
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, FALSE)

# Assert
expect_equal(length(output$MissingPercentage), 1, info = "Empty conditions: MissingPercentage length")
expect_true(is.numeric(output$MissingPercentage), info = "Empty conditions: Numeric output")

# Test 4: Multiple contrasts with different missing patterns
# Setup
contrast_matrix <- matrix(c(1, -1, 0, 0, 1, -1, 1, 0, -1), nrow = 3, ncol = 3)
colnames(contrast_matrix) <- c("Group1", "Group2", "Group3")

summarized <- data.table(
    GROUP = c("Group1", "Group2", "Group3"),
    TotalGroupMeasurements = c(100, 100, 100),
    NumMeasuredFeature = c(90, 80, 70),
    NumImputedFeature = c(5, 10, 15)
)

result <- list()
samples_info <- data.table(GROUP = c("Group1", "Group2", "Group3"), NumRuns = c(10, 10, 10))

# Execute
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)

# Assert
expect_equal(length(output$MissingPercentage), 3, info = "Multiple contrasts: MissingPercentage length")
expect_equal(length(output$ImputationPercentage), 3, info = "Multiple contrasts: ImputationPercentage length")
expect_true(all(output$MissingPercentage >= 0), info = "Multiple contrasts: MissingPercentage non-negative")
expect_true(all(output$MissingPercentage <= 1), info = "Multiple contrasts: MissingPercentage <= 1")

# Test 5: Edge case with all values missing in one group
# Setup
contrast_matrix <- matrix(c(1, -1), nrow = 1, ncol = 2)
colnames(contrast_matrix) <- c("Group1", "Group2")

summarized <- data.table(
    GROUP = c("Group1", "Group2"),
    TotalGroupMeasurements = c(100, 100),
    NumMeasuredFeature = c(0, 100),
    NumImputedFeature = c(0, 0)
)

result <- list()
samples_info <- data.table(GROUP = c("Group1", "Group2"), NumRuns = c(10, 10))

# Execute
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, FALSE)

# Assert
expected_missing <- 1 - (0 + 100) / (100 + 100) # 0.5
expect_equal(output$MissingPercentage[1], expected_missing, info = "Complete missing group: Missing percentage calculation")

# Test 6: Single group contrast
# Setup
contrast_matrix <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
colnames(contrast_matrix) <- c("Group1", "Group2", "Group3")

summarized <- data.table(
    GROUP = c("Group1", "Group2", "Group3"),
    TotalGroupMeasurements = c(100, 100, 100),
    NumMeasuredFeature = c(75, 80, 85),
    NumImputedFeature = c(10, 5, 0)
)

result <- list()
samples_info <- data.table(GROUP = c("Group1", "Group2", "Group3"), NumRuns = c(10, 10, 10))

# Execute
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)

# Assert
expected_missing <- 1 - 75 / 100 # Only Group1 is involved
expected_imputed <- 10 / 100

expect_equal(output$MissingPercentage[1], expected_missing, info = "Single group contrast: Missing percentage")
expect_equal(output$ImputationPercentage[1], expected_imputed, info = "Single group contrast: Imputation percentage")

# Test 7: Test with NA values in summarized data
# Setup
contrast_matrix <- matrix(c(1, -1), nrow = 1, ncol = 2)
colnames(contrast_matrix) <- c("Group1", "Group2")

summarized <- data.table(
    GROUP = c("Group1", "Group2"),
    TotalGroupMeasurements = c(100, 100),
    NumMeasuredFeature = c(NA, 80),
    NumImputedFeature = c(0, NA)
)

result <- list()
samples_info <- data.table(GROUP = c("Group1", "Group2"), NumRuns = c(10, 10))

# Execute
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)

# Assert - should handle NA values with na.rm = TRUE
expect_true(is.numeric(output$MissingPercentage), info = "NA values: Numeric MissingPercentage")
expect_true(is.numeric(output$ImputationPercentage), info = "NA values: Numeric ImputationPercentage")
expect_equal(length(output$MissingPercentage), 1, info = "NA values: Correct length")

# Test 8: Test preserving existing result data
# Setup
contrast_matrix <- matrix(c(1, -1), nrow = 1, ncol = 2)
colnames(contrast_matrix) <- c("Group1", "Group2")

summarized <- data.table(
    GROUP = c("Group1", "Group2"),
    TotalGroupMeasurements = c(100, 100),
    NumMeasuredFeature = c(90, 80),
    NumImputedFeature = c(5, 10)
)

result <- list(existing_data = "should_be_preserved", other_field = 123)
samples_info <- data.table(GROUP = c("Group1", "Group2"), NumRuns = c(10, 10))

# Execute
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, FALSE)

# Assert
expect_equal(output$existing_data, "should_be_preserved", info = "Preserve existing: existing_data field")
expect_equal(output$other_field, 123, info = "Preserve existing: other_field")
expect_true("MissingPercentage" %in% names(output), info = "Preserve existing: MissingPercentage added")

# Test 9: Test with complex contrast matrix (multiple comparisons)
# Complex contrast matrix with multiple non-zero values

# Setup - contrast involving multiple groups with different weights
contrast_matrix <- matrix(c(1, 1, -2, 0.5, -0.5, 0), nrow = 2, ncol = 3)
colnames(contrast_matrix) <- c("Group1", "Group2", "Group3")

summarized <- data.table(
    GROUP = c("Group1", "Group2", "Group3"),
    TotalGroupMeasurements = c(200, 150, 100),
    NumMeasuredFeature = c(180, 120, 80),
    NumImputedFeature = c(10, 15, 5)
)

result <- list()
samples_info <- data.table(GROUP = c("Group1", "Group2", "Group3"), NumRuns = c(20, 15, 10))

# Execute
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)

# Assert
expect_equal(length(output$MissingPercentage), 2)
expect_equal(length(output$ImputationPercentage), 2)
expect_true(all(output$MissingPercentage >= 0))
expect_true(all(output$ImputationPercentage >= 0))

# Test 10: Edge case with zero total measurements
# Handles zero total measurements edge case

# Setup
contrast_matrix <- matrix(c(1, -1), nrow = 1, ncol = 2)
colnames(contrast_matrix) <- c("Group1", "Group2")

summarized <- data.table(
    GROUP = c("Group1", "Group2"),
    TotalGroupMeasurements = c(0, 100),
    NumMeasuredFeature = c(0, 80),
    NumImputedFeature = c(0, 10)
)

result <- list()
samples_info <- data.table(GROUP = c("Group1", "Group2"), NumRuns = c(10, 10))

# Execute
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)

# Assert - should handle division by potentially small numbers
expect_true(is.numeric(output$MissingPercentage))
expect_true(is.numeric(output$ImputationPercentage))
expect_false(is.na(output$MissingPercentage[1]))

