# Characterization tests for MSstatsMergeFractions --------------------------
#
# MSstatsMergeFractions collapses the multiple fraction-runs of a sample into a
# single merged run and drops feature/fraction combinations that were never
# observed. The function is hard to read, so these tests pin its observable
# behavior so the in-place-update rewrite (and any future simplification) stays
# equivalent. They focus on the non-TECHREPLICATE branch, which is the part this
# change rewrote; the TECHREPLICATE branch is unchanged from the previous code.

MSstatsConvert::MSstatsLogsSettings(FALSE)

# Two subjects, two fractions each. f2 in fraction 2 is all-zero abundance, so
# that (FEATURE, FRACTION) combination has no positive observation.
make_fraction_input <- function(f2_frac2 = c(0.0, 0.0)) {
    input <- data.table::data.table(
        PROTEIN = "P1",
        FEATURE = c("f1", "f1", "f2", "f2", "f1", "f1", "f2", "f2"),
        GROUP_ORIGINAL = "G",
        SUBJECT_ORIGINAL = c("S1", "S1", "S1", "S1", "S2", "S2", "S2", "S2"),
        RUN = factor(c(1, 2, 1, 2, 3, 4, 3, 4)),
        originalRUN = c("a1", "a2", "a1", "a2", "b1", "b2", "b1", "b2"),
        FRACTION = c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L),
        INTENSITY = c(10, 20, 30, 40, 50, 60, 70, 80),
        ABUNDANCE = c(3.0, 4.0, 2.0, 0.0, 3.5, 4.5, 2.5, 0.0)
    )
    # The two f2 / fraction-2 rows (positions 4 and 8) take the parametrized
    # abundance. Assigned after construction to avoid referencing a local in a
    # data.table() column expression.
    input$ABUNDANCE[c(4L, 8L)] <- f2_frac2
    input
}

# Single fraction -> nothing to merge, returned untouched (besides clamping) ---
single_fraction = data.table::data.table(
    PROTEIN = "P1", FEATURE = c("f1", "f2"),
    GROUP_ORIGINAL = "G", SUBJECT_ORIGINAL = "S1",
    RUN = factor(c(1, 1)), originalRUN = c("a1", "a1"), FRACTION = c(1L, 1L),
    INTENSITY = c(10, 20), ABUNDANCE = c(3.0, 4.0))
res_single = MSstatsMergeFractions(data.table::copy(single_fraction))
expect_equal(nrow(res_single), 2L,
             info = "single-fraction input keeps all rows")
expect_equal(res_single$originalRUN, c("a1", "a1"),
             info = "single-fraction input leaves originalRUN untouched")

# Abundance clamping happens before the fraction check ------------------------
# ABUNDANCE < 0 becomes 0; an INTENSITY of exactly 1 forces ABUNDANCE to 0.
clamp_input = data.table::data.table(
    PROTEIN = "P1", FEATURE = c("f1", "f2"),
    GROUP_ORIGINAL = "G", SUBJECT_ORIGINAL = "S1",
    RUN = factor(c(1, 1)), originalRUN = c("a1", "a1"), FRACTION = c(1L, 1L),
    INTENSITY = c(1, 50), ABUNDANCE = c(2.5, -3.0))
res_clamp = MSstatsMergeFractions(data.table::copy(clamp_input))
expect_equal(res_clamp$ABUNDANCE, c(0, 0),
             info = "INTENSITY == 1 and negative ABUNDANCE are both clamped to 0")

# Non-TECHREPLICATE: fractions merge and unobserved combos drop ---------------
res = MSstatsMergeFractions(make_fraction_input())

# f2/fraction-2 was all-zero, so those 2 rows are dropped; 6 remain.
expect_equal(nrow(res), 6L,
             info = "feature/fraction with no positive abundance is removed")
expect_equal(nrow(res[FEATURE == "f2" & FRACTION == 2L]), 0L,
             info = "the unobserved f2/fraction-2 rows are gone")

# Each sample's fraction-runs collapse to one '<group>_<subject>_<run>_merged'.
expect_equal(sort(unique(res$originalRUN)),
             c("G_S1_a1_merged", "G_S2_b1_merged"),
             info = "originalRUN is the merged run name, one per sample")

# RUN is refactored to one level per merged run (2 samples -> 2 levels).
expect_true(is.factor(res$RUN), info = "RUN stays a factor")
expect_equal(nlevels(res$RUN), 2L,
             info = "two merged samples give two RUN levels")

# The temporary newRun column must not leak into the output.
expect_false("newRun" %in% colnames(res),
             info = "newRun helper column is removed")

# Non-TECHREPLICATE: fully observed features are all kept ---------------------
res_keep = MSstatsMergeFractions(make_fraction_input(f2_frac2 = c(1.0, 1.5)))
expect_equal(nrow(res_keep), 8L,
             info = "nothing is dropped when every feature/fraction is observed")
expect_equal(sort(unique(paste(res_keep$FEATURE, res_keep$FRACTION))),
             c("f1 1", "f1 2", "f2 1", "f2 2"),
             info = "all feature/fraction combinations survive")
