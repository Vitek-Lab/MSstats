library(data.table)

# Test suite for .countMissingPercentage function
## Test 1: Basic functionality with no missing values
contrast_matrix <- matrix(c(1, -1, 0, 
                            0, 1, -1), nrow = 2, ncol = 3, byrow = TRUE)
colnames(contrast_matrix) <- c("Group1", "Group2", "Group3")

summarized <- data.table(
    GROUP = c("Group1", "Group1", "Group2", "Group2", "Group3", "Group3"),
    TotalGroupMeasurements = c(100, 100, 100, 100, 100, 100),
    NumMeasuredFeature = c(50, 50, 50, 50, 50, 50),
    NumImputedFeature = c(0, 0, 0, 0, 0, 0)
)

result <- data.table(
    logFC = c(6.154384, 6.154384),
    SE = c(0.2917031, 0.2917031),
    Tvalue = c(21.09811, 21.09811),
    DF = c(4, 4),
    pvalue = c(0.0000381, 0.0000381),
    Protein = c("IDHC", "IDHC"),
    Label = c("Group1 - Group2", "Group2 - Group3"),
    issue = c(NA, NA)
)
samples_info <- data.table(GROUP = c("Group1", "Group2", "Group3"), NumRuns = c(2, 2, 2))
output <- MSstats:::.countMissingPercentage(
    contrast_matrix, summarized, result, samples_info, FALSE
)
expect_equal(length(output$MissingPercentage), 2, info = "Basic functionality: MissingPercentage length")
expect_equal(output$MissingPercentage, c(0, 0), info = "Basic functionality: No missing values")
expect_true(is.null(output$ImputationPercentage), info = "Basic functionality: No imputation when has_imputed = FALSE")
expect_true(all(names(output) %in% c(names(result), "MissingPercentage")), info = "Basic functionality: Preserve existing result data")

## Test 2: With imputed values
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
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)
expected_missing <- 1 - (80 + 70) / (100 + 100) # 0.25
expected_imputed <- (10 + 20) / (100 + 100) # 0.15
expect_equal(output$MissingPercentage[1], expected_missing, info = "Imputed values: Missing percentage calculation")
expect_equal(output$ImputationPercentage[1], expected_imputed, info = "Imputed values: Imputation percentage calculation")

## Test 3: With empty conditions (groups not in summarized data)

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


output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, FALSE)


expect_equal(length(output$MissingPercentage), 1, info = "Empty conditions: MissingPercentage length")
expect_true(is.numeric(output$MissingPercentage), info = "Empty conditions: Numeric output")

## Test 4: Multiple contrasts with different missing patterns
## THIS TEST IS FAILING, NEED TO MAKE SURE IT DOES NOT FAIL

contrast_matrix <- matrix(c(1, -1, 0, 
                            0, 1, -1, 
                            1, 0, -1), nrow = 3, ncol = 3, byrow = TRUE)
colnames(contrast_matrix) <- c("Group3", "Group2", "Group1")

summarized <- data.table(
    GROUP = c("Group1", "Group2", "Group3"),
    TotalGroupMeasurements = c(100, 100, 100),
    NumMeasuredFeature = c(90, 80, 70),
    NumImputedFeature = c(5, 10, 15)
)

result <- list()
samples_info <- data.table(
    GROUP = c("Group3", "Group2", "Group1"), 
    NumRuns = c(1, 1, 1)
)

# Execute
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)

# Assert - Verify specific calculations for each contrast
# Contrast 1: Group3 vs Group2 (1*Group3 + -1*Group2 + 0*Group1)
# Groups involved: Group3 (70 measured, 100 total), Group2 (80 measured, 100 total)
expected_missing_1 <- 1 - (70 + 80) / (100 + 100)  # 1 - 150/200 = 0.25
expected_imputed_1 <- (15 + 10) / (100 + 100)      # 25/200 = 0.125

# Contrast 2: Group2 vs Group1 (0*Group3 + 1*Group2 + -1*Group1)  
# Groups involved: Group2 (80 measured, 100 total), Group1 (90 measured, 100 total)
expected_missing_2 <- 1 - (80 + 90) / (100 + 100)  # 1 - 170/200 = 0.15
expected_imputed_2 <- (10 + 5) / (100 + 100)       # 15/200 = 0.075

# Contrast 3: Group3 vs Group1 (1*Group3 + 0*Group2 + -1*Group1)
# Groups involved: Group3 (70 measured, 100 total), Group1 (90 measured, 100 total)  
expected_missing_3 <- 1 - (70 + 90) / (100 + 100)  # 1 - 160/200 = 0.20
expected_imputed_3 <- (15 + 5) / (100 + 100)       # 20/200 = 0.10

expect_equal(length(output$MissingPercentage), 3, info = "Column ordering: MissingPercentage length")
expect_equal(length(output$ImputationPercentage), 3, info = "Column ordering: ImputationPercentage length")

# Verify exact calculations regardless of column order
expect_equal(output$MissingPercentage[1], expected_missing_1, info = "Column ordering: Contrast 1 missing percentage (Group3 vs Group2)")
expect_equal(output$ImputationPercentage[1], expected_imputed_1, info = "Column ordering: Contrast 1 imputation percentage")

expect_equal(output$MissingPercentage[2], expected_missing_2, info = "Column ordering: Contrast 2 missing percentage (Group2 vs Group1)")
expect_equal(output$ImputationPercentage[2], expected_imputed_2, info = "Column ordering: Contrast 2 imputation percentage")

expect_equal(output$MissingPercentage[3], expected_missing_3, info = "Column ordering: Contrast 3 missing percentage (Group3 vs Group1)")
expect_equal(output$ImputationPercentage[3], expected_imputed_3, info = "Column ordering: Contrast 3 imputation percentage")

## Test 5: Edge case with all values missing in one group

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


output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, FALSE)


expected_missing <- 1 - (0 + 100) / (100 + 100) # 0.5
expect_equal(output$MissingPercentage[1], expected_missing, info = "Complete missing group: Missing percentage calculation")

## Test 6: Single group contrast

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


output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)


expected_missing <- 1 - 75 / 100 # Only Group1 is involved
expected_imputed <- 10 / 100

expect_equal(output$MissingPercentage[1], expected_missing, info = "Single group contrast: Missing percentage")
expect_equal(output$ImputationPercentage[1], expected_imputed, info = "Single group contrast: Imputation percentage")

## Test 7: Test with NA values in summarized data

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


output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)


expect_true(is.numeric(output$MissingPercentage), info = "NA values: Numeric MissingPercentage")
expect_true(is.numeric(output$ImputationPercentage), info = "NA values: Numeric ImputationPercentage")
expect_equal(length(output$MissingPercentage), 1, info = "NA values: Correct length")

## Test 8: Test preserving existing result data

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


output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, FALSE)


expect_equal(output$existing_data, "should_be_preserved", info = "Preserve existing: existing_data field")
expect_equal(output$other_field, 123, info = "Preserve existing: other_field")
expect_true("MissingPercentage" %in% names(output), info = "Preserve existing: MissingPercentage added")

## Test 9: Test with complex contrast matrix (multiple comparisons)
# Complex contrast matrix with multiple non-zero values


contrast_matrix <- matrix(c(0.5, 0.5, -1, 0.5, -0.5, 0), nrow = 2, ncol = 3, byrow = TRUE)
colnames(contrast_matrix) <- c("Group1", "Group2", "Group3")

summarized <- data.table(
    GROUP = c("Group1", "Group2", "Group3"),
    TotalGroupMeasurements = c(200, 150, 100),
    NumMeasuredFeature = c(180, 120, 80),
    NumImputedFeature = c(10, 15, 5)
)

result <- list()
samples_info <- data.table(GROUP = c("Group1", "Group2", "Group3"), NumRuns = c(20, 15, 10))


output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)


expect_equal(length(output$MissingPercentage), 2)
expect_equal(length(output$ImputationPercentage), 2)
expect_true(all(output$MissingPercentage >= 0))
expect_true(all(output$ImputationPercentage >= 0))

## Test 10: Edge case with zero total measurements
# Handles zero total measurements edge case


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


output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)

expect_true(is.numeric(output$MissingPercentage))
expect_true(is.numeric(output$ImputationPercentage))
expect_false(is.na(output$MissingPercentage[1]))

