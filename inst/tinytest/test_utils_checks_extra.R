# Additional coverage for utils_checks.R: validateAnnotation,
# makePeptidesDictionary, .checkExperimentDesign, .preProcessIntensities.

# validateAnnotation ---------------------------------------------------------

# Test 1: valid group-comparison annotation (each BioReplicate in one Condition)
gc_annot = data.table::data.table(
    BioReplicate = c("s1", "s2", "s3", "s4"),
    Condition = c("H", "H", "D", "D")
)
expect_true(validateAnnotation(gc_annot, design_type = "group comparison"))

# Test 2: only one condition -> error (relative quantification needs >1)
one_cond = data.table::data.table(
    BioReplicate = c("s1", "s2"),
    Condition = c("H", "H")
)
expect_error(validateAnnotation(one_cond))

# Test 3: group comparison with a BioReplicate spanning conditions -> error
gc_bad = data.table::data.table(
    BioReplicate = c("s1", "s1", "s2"),
    Condition = c("H", "D", "D")
)
expect_error(validateAnnotation(gc_bad, design_type = "group comparison"))

# Test 4: valid repeated-measures annotation (BioReplicate across conditions)
rm_annot = data.table::data.table(
    BioReplicate = c("s1", "s1", "s2", "s2"),
    Condition = c("H", "D", "H", "D")
)
expect_true(validateAnnotation(rm_annot, design_type = "repeated measures"))

# Test 5: repeated measures but each BioReplicate has a single condition -> error
rm_bad = data.table::data.table(
    BioReplicate = c("s1", "s2", "s3", "s4"),
    Condition = c("H", "H", "D", "D")
)
expect_error(validateAnnotation(rm_bad, design_type = "repeated measures"))

# Test 6: unrecognized design type -> error
expect_error(validateAnnotation(gc_annot, design_type = "nonsense"))

# makePeptidesDictionary -----------------------------------------------------

dda = data.table::as.data.table(DDARawData)

# Test 7: GLOBALSTANDARDS builds a dictionary with a PEPTIDE key column
peptides_dict = makePeptidesDictionary(dda, "GLOBALSTANDARDS")
expect_true(data.table::is.data.table(peptides_dict))
expect_true("PEPTIDE" %in% colnames(peptides_dict))
expect_true(all(grepl("_", peptides_dict$PEPTIDE)))

# Test 8: any other normalization returns NULL
expect_true(is.null(makePeptidesDictionary(dda, "EQUALIZEMEDIANS")))

# .checkExperimentDesign -----------------------------------------------------

design_ok = data.table::data.table(Condition = c("A", "B"))
# Test 9: complete column passes silently (returns NULL invisibly, no error)
expect_silent(MSstats:::.checkExperimentDesign(design_ok, "Condition"))

# Test 10: NA in the checked column errors
design_bad = data.table::data.table(Condition = c("A", NA))
expect_error(MSstats:::.checkExperimentDesign(design_bad, "Condition"))

# .preProcessIntensities -----------------------------------------------------
MSstatsConvert::MSstatsLogsSettings(FALSE)

# Test 11: intensities below 1 are floored to 1 before log transform, and an
#          ABUNDANCE column is created
intens = data.table::data.table(INTENSITY = c(0.5, 100, 1000, NA))
MSstats:::.preProcessIntensities(intens, log_base = 2)
expect_true("ABUNDANCE" %in% colnames(intens))
expect_equal(intens$INTENSITY[1], 1)        # 0.5 floored to 1
expect_equal(intens$ABUNDANCE[1], 0)        # log2(1) == 0
expect_equal(intens$ABUNDANCE[2], log(100, 2))
expect_true(is.na(intens$ABUNDANCE[4]))     # NA stays NA

# Test 12: with no sub-1 intensities the flooring branch is skipped
intens2 = data.table::data.table(INTENSITY = c(10, 100, 1000))
MSstats:::.preProcessIntensities(intens2, log_base = 10)
expect_equal(intens2$ABUNDANCE, log(c(10, 100, 1000), 10))
