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

# Test 3: a shared stem that is not on a separator boundary is not split
result = strip(c("Control_1", "Contrast_1"))
expect_equal(result$prefix, "")

# Test 4: identical names are left alone rather than reduced to nothing
result = strip(c("same", "same"))
expect_equal(result$labels, c("same", "same"))
expect_equal(result$prefix, "")

# Test 5: a name is never consumed entirely
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

# Test 16: a token wider than the slot is shortened keeping both ends
result = wrap("ABCDEFGHIJKLMNOP", 6)
expect_equal(result, "A...OP")
expect_equal(nchar(result), 6L)

# Test 17: every wrapped line respects the limit
lines = unlist(strsplit(wrap("alpha_beta_gamma_delta", 8), "\n", fixed = TRUE))
expect_true(all(nchar(lines) <= 8L))

# Test .layoutConditionLabels -----------------------------------------------

short = c("1", "2", "3")
long = paste0("Study_Tissue_Timepoint_", c("0hr", "12hrs", "168hrs"))

# Test 18: labels that already fit are returned unchanged, on one line
result = layout_labels(short, 1, 1400, 4)
expect_equal(result$labels, short)
expect_equal(result$size, 4)
expect_equal(result$n_lines, 1L)

# Test 19: labels that do not fit are shortened
result = layout_labels(long, 2, 800, 4)
expect_equal(result$labels, c("0hr", "12hrs", "168hrs"))

# Test 20: the layout takes no text.angle
expect_false("text.angle" %in% names(formals(layout_labels)))

# Test 21: a single condition cannot collide with anything
result = layout_labels("OnlyOneVeryLongConditionName", 1, 400, 4)
expect_equal(result$labels, "OnlyOneVeryLongConditionName")

# Test 22: when stripping cannot help, the font shrinks but not below the floor
no_stem = c("AlphaHepatocyteBaseline", "BetaRenalCortexStimulated",
            "GammaCardiacTissue")
result = layout_labels(no_stem, 2, 800, 4)
expect_true(result$size < 4)
expect_true(result$size >= 2.5)

# Test 23: the drawn labels stay distinct even in that worst case
expect_equal(length(unique(result$labels)), length(no_stem))

# Test 24: n_lines reports the tallest label
result = layout_labels(c("alpha_beta_gamma", "delta_epsilon_zeta"), 1, 300, 4)
expect_equal(result$n_lines,
             max(lengths(strsplit(result$labels, "\n", fixed = TRUE))))

# Test 25: wrapping stops at three lines however narrow the canvas
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

# Test 26: the conditions stay distinct at every one of those widths
for (canvas in c(1400, 900, 600, 400)) {
    result = layout_labels(long_condition_names, 1, canvas, 4)
    expect_equal(length(unique(result$labels)), length(long_condition_names))
}

# Test 27: names that differ only in their tail stay distinct when wrapped
shared_head = c("Cohort_Baseline_Liver_Replicate_Alpha_Treated",
                "Cohort_Baseline_Liver_Replicate_Alpha_Control")
expect_equal(length(unique(wrap(shared_head, 10))), 2L)

# Test 28: full names are drawn when shortening cannot keep them distinct
covariates = c("DiseaseGroupMale", "DiseaseGroupFemale")
expect_equal(length(unique(wrap(covariates, 8))), 1L)
for (canvas in c(800, 400, 200, 120)) {
    result = layout_labels(covariates, 1, canvas, 4)
    expect_equal(length(unique(result$labels)), 2L)
}

# Covariate designs ---------------------------------------------------------

covariate_design = c("Disease_Male", "Disease_Female",
                     "Control_Male", "Control_Female")

# Test 29: nothing is stripped when no stem is shared by every name
expect_equal(strip(covariate_design)$prefix, "")
expect_equal(strip(covariate_design)$labels, covariate_design)

# Test 30: covariate labels stay distinct and keep their condition at any width
for (canvas in c(1400, 800, 500, 300, 200)) {
    result = layout_labels(covariate_design, 1, canvas, 4)
    expect_equal(length(unique(result$labels)), 4L)
    expect_true(all(grepl("^(Dis|Con)", result$labels)))
}

# Test 31: a shared token that is not leading is not stripped
expect_equal(strip(paste0(covariate_design, "_Week1"))$prefix, "")

# Test 32: a leading stem shared by every name is dropped once labels do not fit
two_level = c("Disease_Male", "Disease_Female")
expect_equal(layout_labels(two_level, 1, 1400, 4)$labels, two_level)
expect_equal(layout_labels(two_level, 1, 300, 4)$labels, c("Male", "Female"))
