# Group comparison edge cases: proteins measured in only some conditions.
# These exercise .getEmptyComparison (protein present in a single group) and
# .handleEmptyConditions (a contrast references a condition with no
# measurements for that protein) in utils_groupcomparison.R.

# Build synthetic summarized (ProteinLevelData) with three conditions A/B/C.
# SUBJECT ids are scoped to their group so the design is case-control
# (checkRepeatedDesign == FALSE).
mk = function(prot, groups_vec) {
    n = length(groups_vec)
    subj = paste0(groups_vec, ave(seq_along(groups_vec), groups_vec, FUN = seq_along))
    data.table::data.table(
        RUN = paste0(prot, "_", seq_len(n)),
        Protein = prot,
        LogIntensities = 10 + seq_len(n) / 10,
        originalRUN = paste0(prot, "_Run", seq_len(n)),
        GROUP = groups_vec,
        SUBJECT = subj,
        TotalGroupMeasurements = 1, NumMeasuredFeature = 2,
        MissingPercentage = 0, more50missing = FALSE, NumImputedFeature = 0
    )
}
pld = data.table::rbindlist(list(
    mk("Pall", rep(c("A", "B", "C"), each = 3)),       # all 3 conditions
    mk("PsingleGroup", rep("A", 6)),                   # only condition A
    mk("PmissingC", rep(c("A", "B"), each = 3))        # A and B, C missing
))
pld[, GROUP := factor(GROUP, levels = c("A", "B", "C"))]
data_in = list(ProteinLevelData = pld,
               FeatureLevelData = data.table::data.table(Label = "L"))

comparison = matrix(c(1, -1, 0, 1, 0, -1), nrow = 2, byrow = TRUE)
rownames(comparison) = c("A-B", "A-C")
colnames(comparison) = c("A", "B", "C")

res = suppressWarnings(groupComparison(comparison, data_in, use_log_file = FALSE))
cr = res$ComparisonResult

# Test 1: all three proteins produce results for both contrasts
expect_equal(nrow(cr), 6)
expect_true(all(c("Pall", "PsingleGroup", "PmissingC") %in% cr$Protein))

# Test 2: fully-measured protein has finite, non-issue log-fold changes
pall = cr[cr$Protein == "Pall", ]
expect_true(all(is.finite(pall$log2FC)))
expect_true(all(is.na(pall$issue)))

# Test 3: single-group protein (.getEmptyComparison) flags oneConditionMissing
#         with an infinite fold change for both contrasts
psingle = cr[cr$Protein == "PsingleGroup", ]
expect_true(all(psingle$issue == "oneConditionMissing"))
expect_true(all(is.infinite(psingle$log2FC)))

# Test 4: protein missing condition C (.handleEmptyConditions) - the A-B
#         contrast is estimable, the A-C contrast is flagged
pmissing = cr[cr$Protein == "PmissingC", ]
ab = pmissing[pmissing$Label == "A-B", ]
ac = pmissing[pmissing$Label == "A-C", ]
expect_true(is.finite(ab$log2FC))
expect_true(is.na(ab$issue))
expect_equal(ac$issue, "oneConditionMissing")
expect_true(is.infinite(ac$log2FC))

# Test 5: .getEmptyComparison completeMissing branch - a contrast where BOTH
#         referenced conditions are absent for the single-group protein
comparison_bc = matrix(c(0, 1, -1), nrow = 1)
rownames(comparison_bc) = "B-C"
colnames(comparison_bc) = c("A", "B", "C")
res_bc = suppressWarnings(groupComparison(comparison_bc, data_in, use_log_file = FALSE))
psingle_bc = res_bc$ComparisonResult[res_bc$ComparisonResult$Protein == "PsingleGroup", ]
expect_equal(psingle_bc$issue, "completeMissing")
expect_true(is.na(psingle_bc$log2FC))
