# Test dataProcess with default parameters ------------------------------------
QuantDataDefault = dataProcess(SRMRawData, use_log_file = FALSE)
QuantDataDefaultLinear = dataProcess(DDARawData, use_log_file = FALSE, 
                                     summaryMethod = "linear")

# Test dataProcess with numberOfCores parameter ----------------------
QuantDataParallel = dataProcess(SRMRawData, use_log_file = FALSE,
                                               numberOfCores = 2)
QuantDataParallelLinear = dataProcess(DDARawData, use_log_file = FALSE, 
                                     summaryMethod = "linear", numberOfCores = 2)

expect_equal(nrow(QuantDataDefault$FeatureLevelData),
             nrow(QuantDataParallel$FeatureLevelData))

expect_equal(nrow(QuantDataDefaultLinear$FeatureLevelData),
             nrow(QuantDataParallelLinear$FeatureLevelData))