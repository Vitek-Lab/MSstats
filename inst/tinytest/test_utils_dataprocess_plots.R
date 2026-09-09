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
result = strip(c("Cyno_Colon_Timepoint_0hr", "Cyno_Colon_Timepoint_12hrs"))
expect_equal(result$labels, c("0hr", "12hrs"))
expect_equal(result$prefix, "Cyno_Colon_Timepoint_")

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
# truncated to exactly the slot width
result = wrap("ABCDEFGHIJKLMNOP", 6)
expect_equal(result, "ABC...")
expect_equal(nchar(result), 6L)

# Test 17: every wrapped line respects the limit
lines = unlist(strsplit(wrap("alpha_beta_gamma_delta", 8), "\n", fixed = TRUE))
expect_true(all(nchar(lines) <= 8L))

# Test .layoutConditionLabels -----------------------------------------------

short = c("1", "2", "3")
long = paste0("Cyno_Colon_Timepoint_", c("0hr", "12hrs", "168hrs"))

# Test 18: labels that already fit are returned unchanged, at the caller's font
# size, on one line, under the plain axis title. This is the path every dataset
# that renders correctly today takes.
result = layout_labels(short, 1, 1400, 4)
expect_equal(result$labels, short)
expect_equal(result$size, 4)
expect_equal(result$n_lines, 1L)
expect_equal(result$xaxis, "MS runs")

# Test 19: labels that do not fit are shortened, and the stem moves to the axis
# title so it is still reported
result = layout_labels(long, 2, 800, 4)
expect_equal(result$labels, c("0hr", "12hrs", "168hrs"))
expect_true(grepl("Cyno_Colon_Timepoint_", result$xaxis, fixed = TRUE))

# Test 20: a non-zero text.angle is a deliberate choice by the caller, so the
# layout is left alone. Note that ggplotly() does not carry rotation through,
# which is why rotation is not used as a mitigation here.
result = layout_labels(long, 2, 800, 4, text.angle = 90)
expect_equal(result$labels, long)
expect_equal(result$xaxis, "MS runs")

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
