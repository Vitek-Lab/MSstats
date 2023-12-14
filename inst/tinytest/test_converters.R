# Test MetamorpheustoMSstatsFormat ---------------------------

input_file_path = system.file("tinytest/raw_data/Metamorpheus/AllQuantifiedPeaks.tsv", package="MSstatsConvert")
annotation_file_path = system.file("tinytest/raw_data/Metamorpheus/Annotation.tsv", package = "MSstatsConvert")
input = data.table::fread(input_file_path)
annotation = data.table::fread(annotation_file_path)
output = MSstats:::MetamorpheustoMSstatsFormat(input, annotation = annotation)
expect_equal(ncol(output), 11)
expect_true(nrow(output) > 0)
expected_column_names = c(
    "Run",
    "ProteinName",
    "PeptideSequence",
    "PrecursorCharge",
    "Intensity",
    "FragmentIon",
    "ProductCharge",
    "IsotopeLabelType",
    "Condition",
    "BioReplicate",
    "Fraction"
)
missing_columns = setdiff(expected_column_names, colnames(output))
expect_equal(length(missing_columns), 0)