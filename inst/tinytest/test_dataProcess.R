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


# Test dataProcess with technical replicates & fractions ------------------
msstats_input_fractions_techreps = data.table::fread(
    system.file("tinytest/processed_data/input_techreps_fractions.csv",
                package = "MSstats")
)
QuantDataTechRepsFractions = dataProcess(msstats_input_fractions_techreps, 
                                         use_log_file = FALSE)
expect_true(!is.null(QuantDataTechRepsFractions))
expect_true(nrow(QuantDataTechRepsFractions$FeatureLevelData) > 0)

# Test dataProcess with linear summary method and anomaly scores --------------
msstats_input = data.table::data.table(
    ProteinName = c(rep("Q9UFW8", 10), rep("Q96S19", 15)),
    PeptideSequence = c(rep("AEFEEQNVR", 5), rep("TALYVTPLDR", 5),
                        rep("AFPLAEWQPSDVDQR", 5), rep("ASGLLLER", 5),
                        rep("LowAbundancePeptide", 5)),
    PrecursorCharge = rep(2, 25),
    FragmentIon = rep("y3", 25),
    ProductCharge = rep(1, 25),
    IsotopeLabelType = rep("L", 25),
    Condition = rep(c("A", "A", "A", "B", "B"), 5),
    BioReplicate = rep(seq(1:5), 5),
    Run = rep(paste0("Run", seq(1:5)), 5),
    Intensity = c(1000, 1500, 2000, 2500, 3000,
                  1100, 1600, 2100, 2600, 3100,
                  1200, 1700, 2200, 2000, 1700,
                  1300, 1800, 1200, 2800, 1800,
                  100, 200, 300, 400, 500),
    AnomalyScores = c(rep(0.03, 10), c(0.01,0.01,0.01,0.9,0.8), c(0.01,0.01,0.01,0.7,0.5), rep(0.10, 5))
)

# Must run with linear summarization
expect_error(
    dataProcess(msstats_input, summaryMethod="TMP"),
    pattern = "AnomalyScores column detected in your input columns.*set summaryMethod to 'linear'",
    info = "Should error when AnomalyScores is present but summaryMethod is not 'linear'"
)

# Process data with linear summarization
QuantDataLinearAnomaly = dataProcess(msstats_input,
                                     normalization=FALSE, # no normalization
                                     MBimpute=TRUE, # Use imputation
                                     summaryMethod="linear" # Key MSstats+ parameter
)

expect_true("Variance" %in% colnames(QuantDataLinearAnomaly$ProteinLevelData),
            info = "Variance column should be present in ProteinLevelData")
protein_data = QuantDataLinearAnomaly$ProteinLevelData
expect_true(all(protein_data$Variance > 0, na.rm = TRUE),
            info = "All variance values should be positive")
expect_true(all(is.finite(protein_data$Variance), na.rm = TRUE),
            info = "All variance values should be finite")
variance_values = protein_data[protein_data$Protein == "Q9UFW8"]$Variance
expect_true(all(variance_values == variance_values[1]),
            info = paste("Not all variance values are equal for Q9UFW8"))

# Create comparison dataset with low anomaly scores for contrast
msstats_input_low_anomaly = data.table::data.table(
    ProteinName = c(rep("Q9UFW8", 10), rep("Q96S19", 15)),
    PeptideSequence = c(rep("AEFEEQNVR", 5), rep("TALYVTPLDR", 5),
                        rep("AFPLAEWQPSDVDQR", 5), rep("ASGLLLER", 5),
                        rep("LowAbundancePeptide", 5)),
    PrecursorCharge = rep(2, 25),
    FragmentIon = rep("y3", 25),
    ProductCharge = rep(1, 25),
    IsotopeLabelType = rep("L", 25),
    Condition = rep(c("A", "A", "A", "B", "B"), 5),
    BioReplicate = rep(seq(1:5), 5),
    Run = rep(paste0("Run", seq(1:5)), 5),
    Intensity = c(1000, 1500, 2000, 2500, 3000,
                  1100, 1600, 2100, 2600, 3100,
                  1200, 1700, 2200, 2000, 1700,
                  1300, 1800, 1200, 2800, 1800,
                  100, 200, 300, 400, 500),
    AnomalyScores = rep(0.01, 25) # All low anomaly scores
)

QuantDataLowAnomaly = dataProcess(msstats_input_low_anomaly,
                                  normalization=FALSE,
                                  MBimpute=TRUE,
                                  summaryMethod="linear")
high_anomaly_variance = protein_data$Variance
low_anomaly_variance = QuantDataLowAnomaly$ProteinLevelData$Variance
expect_true(mean(high_anomaly_variance, na.rm = TRUE) > mean(low_anomaly_variance, na.rm = TRUE),
            info = "Mean variance should be higher for dataset with high anomaly scores")


# Verify that without AnomalyScores, different behavior occurs
msstats_input_no_anomaly = msstats_input[, !c("AnomalyScores")]
QuantDataNoAnomaly = dataProcess(msstats_input_no_anomaly,
                                 normalization=FALSE,
                                 MBimpute=TRUE,
                                 summaryMethod="linear")
no_anomaly_variance = QuantDataNoAnomaly$ProteinLevelData$Variance
expect_false(identical(high_anomaly_variance, no_anomaly_variance),
             info = "Variance should be different when AnomalyScores are present vs absent")