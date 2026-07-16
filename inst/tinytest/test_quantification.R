# Setup ------------------------------------------------------------------
QuantData = dataProcess(SRMRawData, use_log_file = FALSE)

# Test 1: Sample quantification, matrix format (default) -------------------
sample_matrix = quantification(QuantData, use_log_file = FALSE)
expect_true(is.data.frame(sample_matrix))
expect_true("Protein" %in% colnames(sample_matrix))
expect_equal(nrow(sample_matrix), length(unique(QuantData$ProteinLevelData$Protein)))

# Test 2: Sample quantification, long format --------------------------------
sample_long = quantification(QuantData, type = "Sample", format = "long",
                              use_log_file = FALSE)
expect_true(is.data.frame(sample_long))
expect_true(all(c("Protein", "Group_Subject", "LogIntensity") %in% colnames(sample_long)))

# Test 3: Group quantification, matrix format -------------------------------
group_matrix = quantification(QuantData, type = "Group", use_log_file = FALSE)
expect_true(is.data.frame(group_matrix))
expect_true("Protein" %in% colnames(group_matrix))
expect_equal(nrow(group_matrix), length(unique(QuantData$ProteinLevelData$Protein)))

# Test 4: Group quantification, long format ---------------------------------
group_long = quantification(QuantData, type = "Group", format = "long",
                             use_log_file = FALSE)
expect_true(is.data.frame(group_long))
expect_true(all(c("Protein", "Group", "LogIntensity") %in% colnames(group_long)))

# Test 5: type is case-insensitive ------------------------------------------
sample_lower = quantification(QuantData, type = "sample", use_log_file = FALSE)
expect_equal(dim(sample_lower), dim(sample_matrix))

# Test 6: invalid type errors -----------------------------------------------
expect_error(quantification(QuantData, type = "Invalid", use_log_file = FALSE))

# Test 7: invalid format errors ----------------------------------------------
expect_error(quantification(QuantData, format = "Invalid", use_log_file = FALSE))
