# Test dataProcess with default parameters ------------------------------------
QuantDataDefault <- dataProcess(SRMRawData, use_log_file = FALSE)
QuantDataDefaultLinear <- dataProcess(DDARawData,
  use_log_file = FALSE,
  summaryMethod = "linear"
)

# Test dataProcess with numberOfCores parameter ----------------------
QuantDataParallel <- dataProcess(SRMRawData,
  use_log_file = FALSE,
  numberOfCores = 2
)
QuantDataParallelLinear <- dataProcess(DDARawData,
  use_log_file = FALSE,
  summaryMethod = "linear", numberOfCores = 2
)

expect_equal(
  nrow(QuantDataDefault$FeatureLevelData),
  nrow(QuantDataParallel$FeatureLevelData)
)

expect_equal(
  nrow(QuantDataDefaultLinear$FeatureLevelData),
  nrow(QuantDataParallelLinear$FeatureLevelData)
)


# Test dataProcess with technical replicates & fractions ------------------
msstats_input_fractions_techreps <- data.table::fread(
  system.file("tinytest/processed_data/input_techreps_fractions.csv",
    package = "MSstats"
  )
)
QuantDataTechRepsFractions <- dataProcess(msstats_input_fractions_techreps,
  use_log_file = FALSE
)
expect_true(!is.null(QuantDataTechRepsFractions))
expect_true(nrow(QuantDataTechRepsFractions$FeatureLevelData) > 0)
