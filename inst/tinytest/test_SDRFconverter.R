# Test SDRFtoAnnotation ------------------------------------------------------

# Test 1: default column names on the built-in example_SDRF dataset
annot = SDRFtoAnnotation(example_SDRF)
expect_true(is.data.frame(annot))
expect_equal(colnames(annot), c("Run", "Condition", "BioReplicate"))
expect_equal(nrow(annot), nrow(example_SDRF))

# Test 2: fraction column can be added when provided
annot_fraction = SDRFtoAnnotation(
    example_SDRF,
    fraction = "comment[fraction identifier]"
)
expect_true("Fraction" %in% colnames(annot_fraction))
expect_equal(nrow(annot_fraction), nrow(example_SDRF))

# Test 3: missing column name errors
expect_error(
    SDRFtoAnnotation(example_SDRF, run_name = "not a real column")
)

# Test SDRF/extractSDRF -------------------------------------------------------

msstats_format_data = data.table::data.table(
    Condition = c("A", "A", "B", "B"),
    BioReplicate = c("1", "2", "3", "4"),
    Run = c("run1", "run2", "run3", "run4"),
    Fraction = c(1, 1, 1, 1),
    Protein = c("P1", "P1", "P1", "P1"),
    Intensity = c(10, 20, 30, 40)
)

# Test 4: default call keeps Fraction column, and values land in the
# correctly-named columns (Run -> run_name, Condition -> condition_name,
# BioReplicate -> biological_replicate), not just column names matching
sdrf_with_fraction = extractSDRF(msstats_format_data, fraction = "Fraction")
expect_true(is.data.frame(sdrf_with_fraction))
expect_equal(nrow(sdrf_with_fraction), 4)
expect_true("Fraction" %in% colnames(sdrf_with_fraction))
expect_true(all(sdrf_with_fraction[["comment[data file]"]] %in%
                msstats_format_data$Run))
expect_true(all(sdrf_with_fraction[["characteristics[disease]"]] %in%
                msstats_format_data$Condition))
expect_true(all(sdrf_with_fraction[["characteristics[biological replicate]"]] %in%
                msstats_format_data$BioReplicate))

# Test 5: fraction = NULL drops the Fraction column and renames the others
sdrf_no_fraction = extractSDRF(msstats_format_data)
expect_false("Fraction" %in% colnames(sdrf_no_fraction))
expect_true(all(c("comment[data file]", "characteristics[disease]",
                  "characteristics[biological replicate]") %in%
                colnames(sdrf_no_fraction)))
expect_true(all(sdrf_no_fraction[["comment[data file]"]] %in%
                msstats_format_data$Run))

# Test 6: meta_data is merged in by run name
meta = data.frame(
    `comment[data file]` = c("run1", "run2", "run3", "run4"),
    instrument = c("QE", "QE", "QE", "QE"),
    check.names = FALSE
)
sdrf_with_meta = extractSDRF(msstats_format_data, meta_data = meta)
expect_true("instrument" %in% colnames(sdrf_with_meta))
expect_equal(nrow(sdrf_with_meta), 4)
