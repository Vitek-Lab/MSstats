# Test suite for .countMissingPercentage function
## Test 1: Basic functionality with no missing values
contrast_matrix <- matrix(c(1, -1, 0, 
                            0, 1, -1), nrow = 2, ncol = 3, byrow = TRUE)
colnames(contrast_matrix) <- c("Group1", "Group2", "Group3")
summarized <- data.table::data.table(
    GROUP = c("Group1", "Group1", "Group2", "Group2", "Group3", "Group3"),
    TotalGroupMeasurements = c(100, 100, 100, 100, 100, 100),
    NumMeasuredFeature = c(50, 50, 50, 50, 50, 50),
    NumImputedFeature = c(0, 0, 0, 0, 0, 0)
)
result <- data.table::data.table(
    logFC = c(6.154384, 6.154384),
    SE = c(0.2917031, 0.2917031),
    Tvalue = c(21.09811, 21.09811),
    DF = c(4, 4),
    pvalue = c(0.0000381, 0.0000381),
    Protein = c("IDHC", "IDHC"),
    Label = c("Group1 - Group2", "Group2 - Group3"),
    issue = c(NA, NA)
)
samples_info <- data.table::data.table(GROUP = c("Group1", "Group2", "Group3"), NumRuns = c(2, 2, 2))
output <- MSstats:::.countMissingPercentage(
    contrast_matrix, summarized, result, samples_info, FALSE
)
expect_equal(length(output$MissingPercentage), 2, info = "Basic functionality: MissingPercentage length")
expect_equal(output$MissingPercentage, c(0, 0), info = "Basic functionality: No missing values")
expect_true(is.null(output$ImputationPercentage), info = "Basic functionality: No imputation when has_imputed = FALSE")
expect_true(all(names(result) %in% names(output)), 
            info = "Basic functionality: Preserve existing result columns")
expect_true(setequal(names(output), c(names(result), "MissingPercentage")), 
            info = "Basic functionality: No extraneous columns added")

## Test 2: With imputed values
contrast_matrix <- matrix(c(1, -1), nrow = 1, ncol = 2)
colnames(contrast_matrix) <- c("Group1", "Group2")
summarized <- data.table::data.table(
    GROUP = c("Group1", "Group2"),
    TotalGroupMeasurements = c(100, 100),
    NumMeasuredFeature = c(80, 70),
    NumImputedFeature = c(10, 20)
)
result <- list()
samples_info <- data.table::data.table(GROUP = c("Group1", "Group2"), NumRuns = c(10, 10))
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)
expected_missing <- 1 - (80 + 70) / (100 + 100) # 0.25
expected_imputed <- (10 + 20) / (100 + 100) # 0.15
expect_equal(output$MissingPercentage[1], expected_missing, info = "Imputed values: Missing percentage calculation")
expect_equal(output$ImputationPercentage[1], expected_imputed, info = "Imputed values: Imputation percentage calculation")

## Test 3: With empty conditions (groups not in summarized data)
contrast_matrix <- matrix(c(1, -1, 0), nrow = 1, ncol = 3)
colnames(contrast_matrix) <- c("Group1", "Group2", "Group3")
summarized <- data.table::data.table(
    GROUP = c("Group1"),
    TotalGroupMeasurements = c(100),
    NumMeasuredFeature = c(80),
    NumImputedFeature = c(0)
)

result <- list()
samples_info <- data.table::data.table(GROUP = c("Group1", "Group2", "Group3"), NumRuns = c(10, 10, 10))
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, FALSE)
expect_equal(length(output$MissingPercentage), 1, info = "Empty conditions: MissingPercentage length")
expect_true(is.numeric(output$MissingPercentage), info = "Empty conditions: Numeric output")

## Test 4: Multiple contrasts with different missing patterns
contrast_matrix <- matrix(c(1, -1, 0, 
                            0, 1, -1, 
                            1, 0, -1), nrow = 3, ncol = 3, byrow = TRUE)
colnames(contrast_matrix) <- c("Group3", "Group2", "Group1")
summarized <- data.table::data.table(
    GROUP = c("Group1", "Group2", "Group3"),
    TotalGroupMeasurements = c(100, 100, 100),
    NumMeasuredFeature = c(90, 80, 70),
    NumImputedFeature = c(5, 10, 15)
)
result <- list()
samples_info <- data.table::data.table(
    GROUP = c("Group3", "Group2", "Group1"), 
    NumRuns = c(1, 1, 1)
)
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)

expected_missing_1 <- 1 - (70 + 80) / (100 + 100)  # 1 - 150/200 = 0.25
expected_imputed_1 <- (15 + 10) / (100 + 100)      # 25/200 = 0.125
expected_missing_2 <- 1 - (80 + 90) / (100 + 100)  # 1 - 170/200 = 0.15
expected_imputed_2 <- (10 + 5) / (100 + 100)       # 15/200 = 0.075
expected_missing_3 <- 1 - (70 + 90) / (100 + 100)  # 1 - 160/200 = 0.20
expected_imputed_3 <- (15 + 5) / (100 + 100)       # 20/200 = 0.10

expect_equal(length(output$MissingPercentage), 3, info = "Column ordering: MissingPercentage length")
expect_equal(length(output$ImputationPercentage), 3, info = "Column ordering: ImputationPercentage length")
expect_equal(output$MissingPercentage[1], expected_missing_1, info = "Column ordering: Contrast 1 missing percentage (Group3 vs Group2)")
expect_equal(output$ImputationPercentage[1], expected_imputed_1, info = "Column ordering: Contrast 1 imputation percentage")
expect_equal(output$MissingPercentage[2], expected_missing_2, info = "Column ordering: Contrast 2 missing percentage (Group2 vs Group1)")
expect_equal(output$ImputationPercentage[2], expected_imputed_2, info = "Column ordering: Contrast 2 imputation percentage")
expect_equal(output$MissingPercentage[3], expected_missing_3, info = "Column ordering: Contrast 3 missing percentage (Group3 vs Group1)")
expect_equal(output$ImputationPercentage[3], expected_imputed_3, info = "Column ordering: Contrast 3 imputation percentage")

## Test 5: Edge case with all values missing in one group
contrast_matrix <- matrix(c(1, -1), nrow = 1, ncol = 2)
colnames(contrast_matrix) <- c("Group1", "Group2")
summarized <- data.table::data.table(
    GROUP = c("Group1", "Group2"),
    TotalGroupMeasurements = c(0, 100),
    NumMeasuredFeature = c(0, 80),
    NumImputedFeature = c(0, 20)
)
result <- list()
samples_info <- data.table::data.table(GROUP = c("Group1", "Group2"), NumRuns = c(10, 10))
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, FALSE)
expected_missing <- 1 - (0 + 80) / (0 + 100) # 0.2
expect_equal(output$MissingPercentage[1], expected_missing, info = "Complete missing group: Missing percentage calculation")

## Test 6: Test with complex contrast matrix (multiple comparisons)
contrast_matrix <- matrix(c(0.5, 0.5, -1, 1, -1, 0), nrow = 2, ncol = 3, byrow = TRUE)
colnames(contrast_matrix) <- c("Group1", "Group2", "Group3")
summarized <- data.table::data.table(
    GROUP = c("Group1", "Group2", "Group3"),
    TotalGroupMeasurements = c(200, 150, 100),
    NumMeasuredFeature = c(180, 120, 80),
    NumImputedFeature = c(10, 15, 5)
)
result <- list()
samples_info <- data.table::data.table(GROUP = c("Group1", "Group2", "Group3"), NumRuns = c(20, 15, 10))
output <- MSstats:::.countMissingPercentage(contrast_matrix, summarized, result, samples_info, TRUE)
expected_missing_1 <- 1 - 380 / 450  
expected_imputed_1 <- 30 / 450       
expected_missing_2 <- 1 - 300 / 350  # 1 - 0.8571 = 0.1429 (approximately)
expected_imputed_2 <- 25 / 350       # 0.0714 (approximately)
expect_equal(length(output$MissingPercentage), 2, info = "Complex contrast: MissingPercentage length")
expect_equal(length(output$ImputationPercentage), 2, info = "Complex contrast: ImputationPercentage length")
expect_equal(output$MissingPercentage[1], expected_missing_1, tolerance = 1e-10, 
             info = "Complex contrast: Contrast 1 missing percentage (0.5*Group1 + 0.5*Group2 - Group3)")
expect_equal(output$ImputationPercentage[1], expected_imputed_1, tolerance = 1e-10,
             info = "Complex contrast: Contrast 1 imputation percentage")
expect_equal(output$MissingPercentage[2], expected_missing_2, tolerance = 1e-10,
             info = "Complex contrast: Contrast 2 missing percentage (Group1 - Group2)")  
expect_equal(output$ImputationPercentage[2], expected_imputed_2, tolerance = 1e-10,
             info = "Complex contrast: Contrast 2 imputation percentage")