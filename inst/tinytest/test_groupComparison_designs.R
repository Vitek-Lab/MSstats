# Coverage for the model-fitting branches of .fitModelForGroupComparison in
# utils_groupcomparison.R. The existing test_groupComparison.R only exercises
# the case-control / no-tech-replicate path (a plain lm). Here we drive the
# mixed-model branches: repeated measures, technical replicates, and both.

set.seed(101)

build = function(proteins, group_vec, subject_vec) {
    data.table::rbindlist(lapply(proteins, function(p) {
        n = length(group_vec)
        data.table::data.table(
            RUN = paste0(p, "_", seq_len(n)),
            Protein = p,
            LogIntensities = 10 + rnorm(n, 0, 0.3) + (group_vec == "B") * 0.6,
            originalRUN = paste0(p, "_Run", seq_len(n)),
            GROUP = group_vec,
            SUBJECT = subject_vec,
            TotalGroupMeasurements = n / 2, NumMeasuredFeature = 2,
            MissingPercentage = 0, more50missing = FALSE, NumImputedFeature = 0
        )
    }))
}

run_gc = function(pld) {
    pld[, GROUP := factor(GROUP, levels = c("A", "B"))]
    data_in = list(ProteinLevelData = pld,
                   FeatureLevelData = data.table::data.table(Label = "L"))
    cm = matrix(c(1, -1), nrow = 1)
    rownames(cm) = "A-B"; colnames(cm) = c("A", "B")
    suppressWarnings(groupComparison(cm, data_in, use_log_file = FALSE))
}

# Test 1: repeated measures, no technical replicates (subjects across groups)
#         -> lmer(ABUNDANCE ~ GROUP + (1|SUBJECT))
pld_rep = build(c("P1", "P2"),
                group_vec = rep(c("A", "B"), each = 3),
                subject_vec = rep(c("s1", "s2", "s3"), 2))
expect_true(checkRepeatedDesign(list(ProteinLevelData = pld_rep)))
res_rep = run_gc(pld_rep)
expect_equal(nrow(res_rep$ComparisonResult), 2)
expect_true(all(is.finite(res_rep$ComparisonResult$log2FC)))
# a mixed model was fitted (lmerMod), not a plain lm
expect_true(inherits(res_rep$FittedModel[[1]], "lmerMod"))

# Test 2: case-control WITH technical replicates (subject repeated within group)
#         -> lmer(ABUNDANCE ~ GROUP + (1|SUBJECT))
pld_tech = build(c("P1", "P2"),
                 group_vec = rep(c("A", "B"), each = 4),
                 subject_vec = c("s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"))
expect_false(checkRepeatedDesign(list(ProteinLevelData = pld_tech)))
expect_true(MSstats:::.checkTechReplicate(pld_tech))
res_tech = run_gc(pld_tech)
expect_equal(nrow(res_tech$ComparisonResult), 2)
expect_true(inherits(res_tech$FittedModel[[1]], "lmerMod"))

# Test 3: repeated measures WITH technical replicates
#         -> lmer(ABUNDANCE ~ GROUP + (1|SUBJECT) + (1|GROUP:SUBJECT))
pld_both = build(c("P1", "P2"),
                 group_vec = rep(c("A", "B"), each = 6),
                 subject_vec = rep(c("s1", "s1", "s2", "s2", "s3", "s3"), 2))
expect_true(checkRepeatedDesign(list(ProteinLevelData = pld_both)))
expect_true(MSstats:::.checkTechReplicate(pld_both))
res_both = run_gc(pld_both)
expect_equal(nrow(res_both$ComparisonResult), 2)
expect_true(all(is.finite(res_both$ComparisonResult$log2FC)))
expect_true(inherits(res_both$FittedModel[[1]], "lmerMod"))
