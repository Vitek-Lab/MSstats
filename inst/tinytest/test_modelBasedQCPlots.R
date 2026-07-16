# These tests exercise the real modelBasedQCPlots()/groupComparisonQCPlots()
# plotting logic. The existing test_groupComparisonQCPlots.R mocks out
# modelBasedQCPlots entirely, so none of this logic was previously covered.

QuantData = dataProcess(SRMRawData, use_log_file = FALSE)
comparison = matrix(c(-1,0,0,0,0,0,1,0,0,0), nrow=1)
row.names(comparison) = "T7-T1"
colnames(comparison) = unique(QuantData$ProteinLevelData$GROUP)
testResultOneComparison = groupComparison(contrast.matrix = comparison, data = QuantData,
                                           use_log_file = FALSE)
protein_name = as.character(unique(testResultOneComparison$ComparisonResult$Protein))[1]

# Test 1: invalid type errors -------------------------------------------------
expect_error(
    modelBasedQCPlots(data = testResultOneComparison, type = "INVALID",
                       which.Protein = protein_name, address = FALSE)
)

tmp_dir = tempfile("msstats_modelbasedqc_")
dir.create(tmp_dir)
address_prefix = paste0(tmp_dir, "/")

# Test 2: QQPlots creates a pdf file ------------------------------------------
invisible(capture.output(
    modelBasedQCPlots(data = testResultOneComparison, type = "QQPlots",
                       which.Protein = protein_name, address = address_prefix)
))
expect_true(file.exists(paste0(address_prefix, "QQPlot.pdf")))

# Test 3: ResidualPlots creates a pdf file ------------------------------------
invisible(capture.output(
    modelBasedQCPlots(data = testResultOneComparison, type = "ResidualPlots",
                       which.Protein = protein_name, address = address_prefix)
))
expect_true(file.exists(paste0(address_prefix, "ResidualPlot.pdf")))

unlink(tmp_dir, recursive = TRUE)

# Test 4: groupComparisonQCPlots() wrapper actually calls through to
# modelBasedQCPlots and produces the same output ------------------------------
tmp_dir2 = tempfile("msstats_groupcomparisonqcplots_")
dir.create(tmp_dir2)
address_prefix2 = paste0(tmp_dir2, "/")
invisible(capture.output(
    groupComparisonQCPlots(data = testResultOneComparison, type = "QQPlots",
                            which.Protein = protein_name, address = address_prefix2)
))
expect_true(file.exists(paste0(address_prefix2, "QQPlot.pdf")))
unlink(tmp_dir2, recursive = TRUE)
