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

# SRMRawData is a label-based experiment: heavy ("H") rows must be preserved
# in FeatureLevelData after dataProcess
expect_true(
    "H" %in% QuantDataDefault$FeatureLevelData$LABEL,
    info = "FeatureLevelData should contain heavy-label (H) rows for label-based SRM data"
)
expect_true(
    "L" %in% QuantDataDefault$FeatureLevelData$LABEL,
    info = "SRMRawData FeatureLevelData must contain L rows"
)
expect_true(
    nrow(QuantDataDefault$ProteinLevelData) > 0,
    info = "SRMRawData must produce non-empty ProteinLevelData"
)
expect_true(
    "LogIntensities" %in% colnames(QuantDataDefault$ProteinLevelData),
    info = "regression: ProteinLevelData must have LogIntensities column"
)
expect_true(
    "GROUP" %in% colnames(QuantDataDefault$ProteinLevelData),
    info = "regression: ProteinLevelData must have GROUP column"
)


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



# Test MSstatsSummarizeSingleLinear function ------------------------------

single_protein = data.table::data.table(
    PROTEIN = c( "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19" ),
    PEPTIDE = c( "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2" ),
    TRANSITION = c( "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1" ),
    FEATURE = c( "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1" ),
    LABEL = c( "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L" ),
    GROUP_ORIGINAL = c( "A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B" ),
    SUBJECT_ORIGINAL = c( 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5 ),
    RUN = c( "1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4", "5", "5", "5" ),
    GROUP = c( "A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B" ),
    SUBJECT = c( "1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4", "5", "5", "5" ),
    FRACTION = c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ),
    INTENSITY = c( 1200, 1300, 100, 1700, 1800, 200, 2200, 1200, 300, 2000, 2800, 400, 1700, 1800, 500 ),
    ANOMALYSCORES = c( 0.01, 0.01, 0.1, 0.01, 0.01, 0.1, 0.01, 0.01, 0.1, 0.9, 0.3, 0.1, 0.8, 0.2, 0.1 ),
    ABUNDANCE = c( 10.229, 10.344, 6.644, 10.731, 10.814, 7.644, 11.103, 10.229, 8.229, 10.966, 11.451, 8.644, 10.731, 10.814, 8.966 ),
    originalRUN = c( "Run1", "Run1", "Run1", "Run2", "Run2", "Run2", "Run3", "Run3", "Run3", "Run4", "Run4", "Run4", "Run5", "Run5", "Run5" ),
    censored = c( FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE ),
    newABUNDANCE = c( 10.229, 10.344, NA, 10.731, 10.814, NA, 11.103, 10.229, NA, 10.966, 11.451, NA, 10.731, 10.814, NA ),
    cen = c( 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0 ),
    nonmissing = c( TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE ),
    n_obs = c( 5, 5, 0, 5, 5, 0, 5, 5, 0, 5, 5, 0, 5, 5, 0 ),
    n_obs_run = c( 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ),
    total_features = c( 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ),
    prop_features = c( 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667 ),
    nonmissing_all = c( TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE ),
    ABUNDANCE_cut = c( 10.127, 10.127, NA, 10.127, 10.127, NA, 10.127, 10.127, NA, 10.127, 10.127, NA, 10.127, 10.127, NA ),
    any_censored = c( FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE )
)

summarized = MSstats::MSstatsSummarizeSingleLinear(
    single_protein,
    impute = FALSE,
    censored_symbol = "NA", 
    remove50missing = FALSE, 
    aft_iterations = 90,
    equal_variances = TRUE
)

protein_level_summary = summarized[[1]]
feature_level_summary = summarized[[2]]

## Basic tests
expect_equal(nrow(protein_level_summary), 5)
expect_equal(ncol(protein_level_summary), 4)
expect_true(all(protein_level_summary$Protein == "Q96S19"))
expect_equal(as.numeric(protein_level_summary$RUN), 1:5)
expect_equal(nrow(feature_level_summary), 10)
expect_equal(ncol(feature_level_summary), 5)
expect_equal(length(unique(feature_level_summary$FEATURE)), 2)
expect_equal(length(unique(feature_level_summary$RUN)), 5)

## Verify that variance is higher in runs 4-5 with high anomaly scores
variance_runs_1_3 = protein_level_summary[RUN %in% 1:3, Variance]
variance_runs_4_5 = protein_level_summary[RUN %in% 4:5, Variance]

for (var_45 in variance_runs_4_5) {
    for (var_13 in variance_runs_1_3) {
        expect_true(var_45 > var_13, 
                  info = paste("Variance in runs 4-5 should be > variance in runs 1-3:",
                               var_45, ">", var_13))
    }
}

mean_variance_1_3 = mean(variance_runs_1_3)
mean_variance_4_5 = mean(variance_runs_4_5)

expect_true(mean_variance_4_5 > 10 * mean_variance_1_3,
          info = paste("Mean variance of runs 4-5 should be at least 10x higher:",
                       mean_variance_4_5, "vs", 10 * mean_variance_1_3))

## Verify that protein levels in runs 4-5 are closer to ASGLLLER than AFPLAEWQPSDVDQR
protein_runs_4_5 = protein_level_summary[RUN %in% 4:5]
afpla_runs_4_5 = feature_level_summary[FEATURE == "AFPLAEWQPSDVDQR_2_y3_1" & RUN %in% 4:5]
asgll_runs_4_5 = feature_level_summary[FEATURE == "ASGLLLER_2_y3_1" & RUN %in% 4:5]
for (i in 1:nrow(protein_runs_4_5)) {
    run_num = protein_runs_4_5$RUN[i]
    protein_logint = protein_runs_4_5$LogIntensities[i]
    
    afpla_abundance = afpla_runs_4_5[RUN == run_num, newABUNDANCE]
    asgll_abundance = asgll_runs_4_5[RUN == run_num, newABUNDANCE]
    
    dist_to_afpla = abs(protein_logint - afpla_abundance)
    dist_to_asgll = abs(protein_logint - asgll_abundance)
    
    expect_true(dist_to_asgll < dist_to_afpla,
              info = paste("Run", run_num, ": Protein LogIntensity", protein_logint,
                           "should be closer to ASGLLLER", asgll_abundance,
                           "(dist =", round(dist_to_asgll, 4), ")",
                           "than to AFPLAEWQPSDVDQR", afpla_abundance,
                           "(dist =", round(dist_to_afpla, 4), ")"))
}
