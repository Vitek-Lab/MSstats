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

# Test groupComparison with weights
protein_level_data_with_anomalies = data.table::data.table(
    RUN = c( "1", "2", "3", "4", "5", "1", "2", "3", "4", "5" ),
    Protein = c( "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q9UFW8", "Q9UFW8", "Q9UFW8", "Q9UFW8", "Q9UFW8" ),
    LogIntensities = c(10.28, 10.39, 10.63, 11.26, 11.25, 10.03, 10.59, 11.00, 11.31, 11.57),
    Variance = c( 0.0020, 0.0020, 0.15, 0.13, 0.0020, rep(0.00032, 5)),
    originalRUN = c( "Run1", "Run2", "Run3", "Run4", "Run5", "Run1", "Run2", "Run3", "Run4", "Run5" ),
    GROUP = c( "A", "A", "A", "B", "B", "A", "A", "A", "B", "B" ),
    SUBJECT = c( 1, 2, 3, 4, 5, 1, 2, 3, 4, 5 ),
    TotalGroupMeasurements = c( 9, 9, 9, 6, 6, 6, 6, 6, 4, 4 ),
    NumMeasuredFeature = c( 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ),
    MissingPercentage = c( 0.333333333333333, 0.333333333333333, 0.333333333333333, 0.333333333333333, 0.333333333333333, 0, 0, 0, 0, 0 ),
    more50missing = c( FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE ),
    NumImputedFeature = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
)
protein_level_data_no_anomalies = protein_level_data_with_anomalies
protein_level_data_no_anomalies$Variance = NULL
QuantDataNoAnomaly = list(ProteinLevelData = protein_level_data_no_anomalies)
QuantDataWithAnomaly = list(ProteinLevelData = protein_level_data_with_anomalies)
model_no_anomaly = groupComparison("pairwise", QuantDataNoAnomaly,
                        use_log_file = FALSE)
model_with_anomaly = groupComparison("pairwise", QuantDataWithAnomaly,
                         use_log_file = FALSE)

no_anomaly_Q9UFW8 = model_no_anomaly$ComparisonResult[model_no_anomaly$ComparisonResult$Protein == "Q9UFW8", ]
with_anomaly_Q9UFW8 = model_with_anomaly$ComparisonResult[model_with_anomaly$ComparisonResult$Protein == "Q9UFW8", ]

expect_true(nrow(no_anomaly_Q9UFW8) > 0, info = "Q9UFW8 not found in model_no_anomaly")
expect_true(nrow(with_anomaly_Q9UFW8) > 0, info = "Q9UFW8 not found in model_with_anomaly")
expect_equal(no_anomaly_Q9UFW8$SE, with_anomaly_Q9UFW8$SE, 
             info = "SE values for Q9UFW8 should be the same")
expect_equal(no_anomaly_Q9UFW8$pvalue, with_anomaly_Q9UFW8$pvalue,
             info = "pvalue values for Q9UFW8 should be the same")
no_anomaly_Q96S19 = model_no_anomaly$ComparisonResult[model_no_anomaly$ComparisonResult$Protein == "Q96S19", ]
with_anomaly_Q96S19 = model_with_anomaly$ComparisonResult[model_with_anomaly$ComparisonResult$Protein == "Q96S19", ]
expect_true(nrow(no_anomaly_Q96S19) > 0, info = "Q96S19 not found in model_no_anomaly")
expect_true(nrow(with_anomaly_Q96S19) > 0, info = "Q96S19 not found in model_with_anomaly")
expect_true(with_anomaly_Q96S19$SE < no_anomaly_Q96S19$SE,
            info = "SE for Q96S19 should be lower in model_with_anomaly")
expect_true(with_anomaly_Q96S19$pvalue < no_anomaly_Q96S19$pvalue,
            info = "pvalue for Q96S19 should be lower in model_with_anomaly")

