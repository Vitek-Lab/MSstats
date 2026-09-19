# Condition label layout helpers.
#
# These decide what is drawn in place of a condition name that does not fit the
# horizontal room it is given. They are pure functions of the names and the
# canvas geometry, so they are tested directly rather than through a rendered
# plot.

strip = MSstats:::.stripCommonAffix
slot_chars = MSstats:::.conditionSlotChars
wrap = MSstats:::.wrapConditionLabels
layout_labels = MSstats:::.layoutConditionLabels

# Test .stripCommonAffix ----------------------------------------------------

# Test 1: the shared stem is removed and reported
result = strip(c("Study_Tissue_Timepoint_0hr", "Study_Tissue_Timepoint_12hrs"))
expect_equal(result$labels, c("0hr", "12hrs"))
expect_equal(result$prefix, "Study_Tissue_Timepoint_")

# Test 2: names sharing nothing are left alone
result = strip(c("Alpha", "Beta"))
expect_equal(result$labels, c("Alpha", "Beta"))
expect_equal(result$prefix, "")

# Test 3: a shared stem that is not on a separator boundary is not split.
# "Control" and "Contrast" share "Cont", but chopping mid-token would leave
# labels that do not correspond to anything in the data.
result = strip(c("Control_1", "Contrast_1"))
expect_equal(result$prefix, "")

# Test 4: identical names are left alone rather than reduced to nothing
result = strip(c("same", "same"))
expect_equal(result$labels, c("same", "same"))
expect_equal(result$prefix, "")

# Test 5: a name is never consumed entirely. Every name here starts with the
# whole of the first, so stripping greedily would leave an empty label.
result = strip(c("A_B", "A_B_C"))
expect_true(all(nzchar(result$labels)))

# Test 6: a single condition has no shared stem to speak of
result = strip("OnlyOne")
expect_equal(result$labels, "OnlyOne")
expect_equal(result$prefix, "")

# Test 7: separators other than underscore are honoured
expect_equal(strip(c("run.a", "run.b"))$prefix, "run.")
expect_equal(strip(c("run a", "run b"))$prefix, "run ")
expect_equal(strip(c("run-a", "run-b"))$prefix, "run-")

# Test .conditionSlotChars --------------------------------------------------

# Test 8: more conditions in the same canvas means fewer characters each
expect_true(slot_chars(20, 1, 1400, 4) < slot_chars(5, 1, 1400, 4))

# Test 9: a wider canvas means more characters
expect_true(slot_chars(10, 1, 1400, 4) > slot_chars(10, 1, 800, 4))

# Test 10: splitting the canvas across facets means fewer characters
expect_true(slot_chars(10, 2, 1400, 4) < slot_chars(10, 1, 1400, 4))

# Test 11: a larger font means fewer characters
expect_true(slot_chars(10, 1, 1400, 8) < slot_chars(10, 1, 1400, 4))

# Test 12: never returns less than one character, however cramped
expect_true(slot_chars(500, 4, 200, 12) >= 1L)

# Test 13: a canvas of unknown width imposes no limit, so nothing is shortened
expect_equal(slot_chars(10, 1, NA, 4), .Machine$integer.max)
expect_equal(slot_chars(10, 1, 0, 4), .Machine$integer.max)

# Test .wrapConditionLabels -------------------------------------------------

# Test 14: names that already fit are returned untouched
expect_equal(wrap(c("0hr", "12hrs"), 10), c("0hr", "12hrs"))

# Test 15: wrapping happens at separators, not mid-token
expect_equal(wrap("aaaa_bbbb_cccc", 6), "aaaa_\nbbbb_\ncccc")

# Test 16: a single token wider than the slot cannot be broken, so it is
# shortened to exactly the slot width, keeping both ends. Head-only truncation
# would drop the tail, which is the part that tells two conditions apart.
result = wrap("ABCDEFGHIJKLMNOP", 6)
expect_equal(result, "A...OP")
expect_equal(nchar(result), 6L)

# Test 17: every wrapped line respects the limit
lines = unlist(strsplit(wrap("alpha_beta_gamma_delta", 8), "\n", fixed = TRUE))
expect_true(all(nchar(lines) <= 8L))

# Test .layoutConditionLabels -----------------------------------------------

short = c("1", "2", "3")
long = paste0("Study_Tissue_Timepoint_", c("0hr", "12hrs", "168hrs"))

# Test 18: labels that already fit are returned unchanged, at the caller's font
# size, on one line. This is the path every dataset that renders correctly today
# takes.
result = layout_labels(short, 1, 1400, 4)
expect_equal(result$labels, short)
expect_equal(result$size, 4)
expect_equal(result$n_lines, 1L)

# Test 19: labels that do not fit are shortened. The x-axis title stays the
# standard "MS runs" whatever the layout does, so it is not asserted on here.
result = layout_labels(long, 2, 800, 4)
expect_equal(result$labels, c("0hr", "12hrs", "168hrs"))

# Test 20: the layout takes no text.angle. It is only ever computed for the
# Plotly output, which ggplotly() draws horizontally whatever the caller asked
# for, so rotation cannot be a mitigation and cannot suppress one either.
expect_false("text.angle" %in% names(formals(layout_labels)))

# Test 21: a single condition cannot collide with anything
result = layout_labels("OnlyOneVeryLongConditionName", 1, 400, 4)
expect_equal(result$labels, "OnlyOneVeryLongConditionName")

# Test 22: when stripping cannot help, the font shrinks rather than giving up,
# but not below the legibility floor
no_stem = c("AlphaHepatocyteBaseline", "BetaRenalCortexStimulated",
            "GammaCardiacTissue")
result = layout_labels(no_stem, 2, 800, 4)
expect_true(result$size < 4)
expect_true(result$size >= 2.5)

# Test 23: the drawn labels stay distinguishable from one another even in that
# worst case, which is the whole point of the exercise
expect_equal(length(unique(result$labels)), length(no_stem))

# Test 24: n_lines reports the tallest label, so the caller knows how much
# headroom to add above the data
result = layout_labels(c("alpha_beta_gamma", "delta_epsilon_zeta"), 1, 300, 4)
expect_equal(result$n_lines,
             max(lengths(strsplit(result$labels, "\n", fixed = TRUE))))

# Test 25: wrapping stops at three lines however cramped the canvas gets.
# Uncapped, this fixture reached 4 lines at 900px and 8 at 400px.
long_condition_names = c(
    "0hr_0hr_20240101_XX_Sample_ctrl_f1_merged",
    "12hrs_12hrs_20240101_Sample_Tissue_12h_f1_merged",
    "168hrs_168hrs_202401011_XX_Sample_168h_f1_merged",
    "1hr_1hr_20240101_XX_Sample_1h_f1_merged",
    "24hrs_24hrs_20240101_XX_Sample_24h_f1_merged",
    "48hrs_48hrs_20240101_XX_Sample_48h_f1_merged",
    "4hr_4hr_20240101_XX_Sample_4h_f1_merged",
    "96hrs_96hrs_20240101_XX_Sample_96h_f1_merged")
for (canvas in c(1400, 900, 600, 400)) {
    result = layout_labels(long_condition_names, 1, canvas, 4)
    expect_true(result$n_lines <= 3L)
}

# Test 26: and the conditions stay tellable apart at every one of those widths
for (canvas in c(1400, 900, 600, 400)) {
    result = layout_labels(long_condition_names, 1, canvas, 4)
    expect_equal(length(unique(result$labels)), length(long_condition_names))
}

# Test 27: names that differ only in their tail survive the fold onto the last
# line. Keeping the first three lines and dropping the rest would render these
# two conditions as the same string.
shared_head = c("Cohort_Baseline_Liver_Replicate_Alpha_Treated",
                "Cohort_Baseline_Liver_Replicate_Alpha_Control")
expect_equal(length(unique(wrap(shared_head, 10))), 2L)

# Test 28: when no amount of shortening keeps the conditions distinct, the full
# names are drawn instead. A crowded axis is recoverable; two conditions sharing
# one label is not. These share a stem with no separator to break on, so the
# wrapper alone collapses them below eight characters.
covariates = c("DiseaseGroupMale", "DiseaseGroupFemale")
expect_equal(length(unique(wrap(covariates, 8))), 1L)
for (canvas in c(800, 400, 200, 120)) {
    result = layout_labels(covariates, 1, canvas, 4)
    expect_equal(length(unique(result$labels)), 2L)
}

# Covariate designs ---------------------------------------------------------
# "Condition_Gender" is a very common way to encode a covariate, and the
# condition half of the name must survive.

covariate_design = c("Disease_Male", "Disease_Female",
                     "Control_Male", "Control_Female")

# Test 29: no single stem is shared by every name here -- Disease_ and Control_
# each cover only half -- so nothing is dropped.
expect_equal(strip(covariate_design)$prefix, "")
expect_equal(strip(covariate_design)$labels, covariate_design)

# Test 30: and that holds through the whole layout at any canvas width. The
# conditions stay distinct and every label still names its condition, even at
# widths cramped enough to force wrapping and truncation.
for (canvas in c(1400, 800, 500, 300, 200)) {
    result = layout_labels(covariate_design, 1, canvas, 4)
    expect_equal(length(unique(result$labels)), 4L)
    expect_true(all(grepl("^(Dis|Con)", result$labels)))
}

# Test 31: a third factor does not make the strip loop over-consume. "Week1" is
# shared by every name but is not a leading token, so it stays put.
expect_equal(strip(paste0(covariate_design, "_Week1"))$prefix, "")

# Test 32: when every condition genuinely does share a leading stem it is
# dropped, covariate or not. Only reachable once the labels no longer fit, and
# only in the Plotly output, where the hover still carries the full name.
two_level = c("Disease_Male", "Disease_Female")
expect_equal(layout_labels(two_level, 1, 1400, 4)$labels, two_level)
expect_equal(layout_labels(two_level, 1, 300, 4)$labels, c("Male", "Female"))
