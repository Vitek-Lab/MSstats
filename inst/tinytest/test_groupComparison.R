# Setup ------------------------------------------------------------------
QuantData = dataProcess(SRMRawData, use_log_file = FALSE)
comparison = matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
row.names(comparison) = "T7-T1"
groups = levels(QuantData$ProteinLevelData$GROUP)
colnames(comparison) = groups[order(as.numeric(groups))]

# Test groupComparison with default parameters ---------------------------
testResultDefaultComparison = groupComparison(contrast.matrix=comparison,
                                           data=QuantData,
                                           use_log_file = FALSE)

# Test groupComparison with numberOfCores parameter ----------------------
testResultParallelComparison = groupComparison(contrast.matrix=comparison,
                                           data=QuantData,
                                           use_log_file = FALSE,
                                           numberOfCores = 2)

expect_equal(nrow(testResultDefaultComparison$ComparisonResult),
             nrow(testResultParallelComparison$ComparisonResult))