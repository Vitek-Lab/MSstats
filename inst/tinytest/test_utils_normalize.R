# Test for .normalizeGlobalStandards function ----------------------------------

# Setup test data ---------------------------------------------------------------

create_peptide_dictionary <- function() {
    data.table::data.table(
        PeptideSequence = c("AAAAAAAAAAAAAAGAGAGAK", "AAAAAAAAAAAAAAGAGAGAK", "AAAAAAAAAAAAVSRD"),
        PrecursorCharge = c(3, 2, 3),
        PEPTIDE = c("AAAAAAAAAAAAAAGAGAGAK_3", "AAAAAAAAAAAAAAGAGAGAK_2", "AAAAAAAAAAAAVSRD_3")
    )
}

create_test_input <- function(standard_intensities, peptide2_intensities, peptide3_intensities) {
    # Constants
    n_proteins <- 3
    n_runs <- 48
    n_subjects <- 4
    n_fractions <- 6
    
    # Create base structure
    input <- data.table::data.table(
        PROTEIN = rep(c("P1", "P1", "P3"), each = n_runs),
        PEPTIDE = rep(c("AAAAAAAAAAAAAAGAGAGAK_3", "AAAAAAAAAAAAAAGAGAGAK_2", "AAAAAAAAAAAAVSRD_3"), 
                      each = n_runs),
        TRANSITION = rep("NA_NA", n_proteins * n_runs),
        FEATURE = rep(c("AAAAAAAAAAAAAAGAGAGAK_3_NA_NA", 
                        "AAAAAAAAAAAAAAGAGAGAK_2_NA_NA", 
                        "AAAAAAAAAAAAVSRD_3_NA_NA"), 
                      each = n_runs),
        LABEL = rep("L", n_proteins * n_runs),
        GROUP_ORIGINAL = rep(rep(c("Control", "Treatment"), each = n_runs/2), n_proteins),
        SUBJECT_ORIGINAL = rep(paste0("Subject", rep(1:n_subjects, each = n_fractions)), 
                               n_proteins * 2),
        RUN = rep(1:n_runs, n_proteins),
        GROUP = rep(rep(c("Control", "Treatment"), each = n_runs/2), n_proteins),
        SUBJECT = rep(paste0("Subject", rep(1:n_subjects, each = n_fractions)), 
                      n_proteins * 2),
        FRACTION = rep(rep(1:n_fractions, n_subjects), n_proteins * 2),
        INTENSITY = c(standard_intensities, peptide2_intensities, peptide3_intensities),
        ANOMALYSCORES = rep(NA, n_proteins * n_runs),
        originalRUN = rep(paste0("Run", 1:n_runs), n_proteins)
    )
    
    input[, ABUNDANCE := log2(INTENSITY)]
    return(input)
}

# Test 1: Standards with different intensities between groups -------------------
test_different_group_intensities <- function() {
    peptide_dict <- create_peptide_dictionary()
    standards <- c("AAAAAAAAAAAAAAGAGAGAK")
    
    # Control group: 262144, Treatment group: 524288
    standard_intensities <- c(
        rep(262144, 24),  # Control (runs 1-24)
        rep(524288, 24)   # Treatment (runs 25-48)
    )
    
    # Non-standard peptides: all 262144
    peptide3_intensities <- rep(262144, 48)
    
    input <- create_test_input(standard_intensities, standard_intensities, peptide3_intensities)
    output <- MSstats:::.normalizeGlobalStandards(input, peptide_dict, standards)
    
    # Verify normalization: Control runs should be shifted up, Treatment runs shifted down
    control_runs <- 1:24
    treatment_runs <- 25:48
    
    # Check Control group (shifted up to match treatment standard)
    control_abundance <- output[RUN %in% control_runs & 
                                    !is.na(ABUNDANCE) & 
                                    !grepl(standards, PEPTIDE)]$ABUNDANCE
    expect_true(all(abs(control_abundance - 18.5) < 1e-10),
                info = "Control group non-standard peptides should be normalized to 18.5")
    
    # Check Treatment group (shifted down to match control standard)
    treatment_abundance <- output[RUN %in% treatment_runs & 
                                      !is.na(ABUNDANCE) & 
                                      !grepl(standards, PEPTIDE)]$ABUNDANCE
    expect_true(all(abs(treatment_abundance - 17.5) < 1e-10),
                info = "Treatment group non-standard peptides should be normalized to 17.5")
}

# Test 2: Standards with alternating intensities within fractions ---------------
test_alternating_intensities_within_fractions <- function() {
    peptide_dict <- create_peptide_dictionary()
    standards <- c("AAAAAAAAAAAAAAGAGAGAK")
    
    # Standard alternates between 262144 and 524288 within each fraction
    standard_intensities <- rep(c(262144, 524288), 24)
    
    # Non-standard peptides: all 262144
    peptide3_intensities <- rep(262144, 48)
    
    input <- create_test_input(standard_intensities, standard_intensities, peptide3_intensities)
    output <- MSstats:::.normalizeGlobalStandards(input, peptide_dict, standards)
    
    # When standards vary within fractions but average to same level,
    # no net normalization should occur
    all_runs <- 1:48
    normalized_abundance <- output[RUN %in% all_runs & 
                                       !is.na(ABUNDANCE) & 
                                       !grepl(standards, PEPTIDE)]$ABUNDANCE
    
    expect_true(all(abs(normalized_abundance - 18) < 1e-10),
                info = "No normalization should occur when standard averages are equal across fractions")
}

# Run tests ---------------------------------------------------------------------
test_different_group_intensities()
test_alternating_intensities_within_fractions()


# Tests for MSstatsNormalize ---------------------------------------------------

create_normalize_input <- function(n_runs = 4) {
    # Two features across n_runs runs, single fraction
    # Feature A: intensities 100, 200, 400, 800
    # Feature B: intensities 50, 100, 200, 400
    intensities <- c(100, 200, 400, 800, 50, 100, 200, 400)
    n_features <- 2
    data.table::data.table(
        PROTEIN = rep("P1", n_features * n_runs),
        PEPTIDE = rep(c("PEP_A", "PEP_B"), each = n_runs),
        TRANSITION = rep("NA_NA", n_features * n_runs),
        FEATURE = rep(c("PEP_A_NA_NA", "PEP_B_NA_NA"), each = n_runs),
        LABEL = rep("L", n_features * n_runs),
        GROUP_ORIGINAL = rep(c("G1", "G1", "G2", "G2"), n_features),
        SUBJECT_ORIGINAL = rep(c("S1", "S2", "S3", "S4"), n_features),
        RUN = rep(1:n_runs, n_features),
        GROUP = rep(c("G1", "G1", "G2", "G2"), n_features),
        SUBJECT = rep(c("S1", "S2", "S3", "S4"), n_features),
        FRACTION = rep(1L, n_features * n_runs),
        INTENSITY = intensities,
        ANOMALYSCORES = rep(NA_real_, n_features * n_runs),
        originalRUN = rep(paste0("Run", 1:n_runs), n_features),
        ABUNDANCE = log2(intensities)
    )
}

# Test: NONE / FALSE returns input unchanged
input_raw <- create_normalize_input()
out_none <- MSstatsNormalize(input_raw, "NONE")
expect_identical(out_none, input_raw,
                 info = "NONE normalization should return input unchanged")

out_false <- MSstatsNormalize(input_raw, "FALSE")
expect_identical(out_false, input_raw,
                 info = "FALSE normalization should return input unchanged")

# Test: EQUALIZEMEDIANS shifts all runs to the same median
out_median <- MSstatsNormalize(input_raw, "EQUALIZEMEDIANS")

# After median normalization, median ABUNDANCE per run should be equal
run_medians <- out_median[, list(med = median(ABUNDANCE, na.rm = TRUE)), by = "RUN"]
expect_true(
    max(run_medians$med) - min(run_medians$med) < 1e-10,
    info = "EQUALIZEMEDIANS: all run medians should be equal after normalization"
)

# Overall number of rows should be unchanged
expect_equal(nrow(out_median), nrow(input_raw),
             info = "EQUALIZEMEDIANS: row count should be unchanged")

# Relative differences within each run should be preserved
# Feature A - Feature B difference should stay the same within each run before vs after
for (run_id in 1:4) {
    diff_before <- input_raw[RUN == run_id & FEATURE == "PEP_A_NA_NA", ABUNDANCE] -
                   input_raw[RUN == run_id & FEATURE == "PEP_B_NA_NA", ABUNDANCE]
    diff_after  <- out_median[RUN == run_id & FEATURE == "PEP_A_NA_NA", ABUNDANCE] -
                   out_median[RUN == run_id & FEATURE == "PEP_B_NA_NA", ABUNDANCE]
    expect_true(abs(diff_before - diff_after) < 1e-10,
                info = paste("EQUALIZEMEDIANS: within-run relative differences preserved for run", run_id))
}


# Tests for .normalizeMedian ---------------------------------------------------

# Test: single-fraction, single-label data
input_median <- create_normalize_input()

# Run medians before normalization differ (runs have different scales)
before_medians <- input_median[, list(med = median(ABUNDANCE, na.rm = TRUE)), by = "RUN"]
expect_true(
    max(before_medians$med) - min(before_medians$med) > 0.1,
    info = ".normalizeMedian: run medians should differ before normalization"
)

out <- MSstats:::.normalizeMedian(input_median)

# All run medians should converge to the same value
after_medians <- out[, list(med = median(ABUNDANCE, na.rm = TRUE)), by = "RUN"]
expect_true(
    max(after_medians$med) - min(after_medians$med) < 1e-10,
    info = ".normalizeMedian: run medians should be equal after normalization"
)

# Temporary columns should be cleaned up
expect_false("ABUNDANCE_RUN" %in% colnames(out),
             info = ".normalizeMedian: ABUNDANCE_RUN column should be removed")
expect_false("ABUNDANCE_FRACTION" %in% colnames(out),
             info = ".normalizeMedian: ABUNDANCE_FRACTION column should be removed")

# Test: all-missing run is not affected and does not corrupt other runs
input_with_na <- data.table::copy(input_median)
input_with_na[RUN == 1, ABUNDANCE := NA_real_]
out_na <- MSstats:::.normalizeMedian(input_with_na)
non_na_medians <- out_na[RUN != 1, list(med = median(ABUNDANCE, na.rm = TRUE)), by = "RUN"]
expect_true(
    max(non_na_medians$med) - min(non_na_medians$med) < 1e-10,
    info = ".normalizeMedian: all-NA run should not corrupt normalization of other runs"
)


# Tests for .normalizeMedian with heavy labels ---------------------------------

create_labeled_input <- function(n_runs = 4) {
    # Both L and H rows per run. H abundances vary across runs (the reference
    # used for normalization). L abundances are constant so any shift applied
    # to them comes entirely from the H-derived normalization factor.
    #
    # H abundances per run: 16, 17, 18, 19  (log2 scale)
    # L abundances per run: all 16
    n_features <- 2
    runs <- rep(1:n_runs, n_features * 2)
    data.table::data.table(
        PROTEIN   = rep("P1", n_features * n_runs * 2),
        PEPTIDE   = rep(rep(c("PEP_A", "PEP_B"), each = n_runs), 2),
        TRANSITION = rep("NA_NA", n_features * n_runs * 2),
        FEATURE   = rep(rep(c("PEP_A_NA_NA", "PEP_B_NA_NA"), each = n_runs), 2),
        LABEL     = rep(c("L", "H"), each = n_features * n_runs),
        GROUP_ORIGINAL = rep(c("G1", "G1", "G2", "G2"), n_features * 2),
        SUBJECT_ORIGINAL = rep(c("S1", "S2", "S3", "S4"), n_features * 2),
        RUN       = runs,
        GROUP     = rep(c("G1", "G1", "G2", "G2"), n_features * 2),
        SUBJECT   = rep(c("S1", "S2", "S3", "S4"), n_features * 2),
        FRACTION  = rep(1L, n_features * n_runs * 2),
        INTENSITY = rep(1L, n_features * n_runs * 2),  # placeholder
        ANOMALYSCORES = rep(NA_real_, n_features * n_runs * 2),
        originalRUN = rep(paste0("Run", 1:n_runs), n_features * 2),
        ABUNDANCE = c(
            rep(16, n_features * n_runs),          # L: constant across all runs
            rep(c(16, 17, 18, 19), n_features)     # H: increasing across runs
        )
    )
}

input_labeled <- create_labeled_input()
input_labeled_copy <- create_labeled_input()
out_labeled <- MSstats:::.normalizeMedian(input_labeled_copy)

# Normalization factor is derived from H abundances. H run medians are
# 16, 17, 18, 19 and their overall median is 17.5, so the per-run shifts
# are +1.5, +0.5, -0.5, -1.5.

# After normalization, H run medians should be equalized
h_medians_after <- out_labeled[LABEL == "H",
                                list(med = median(ABUNDANCE, na.rm = TRUE)),
                                by = "RUN"]
expect_true(
    max(h_medians_after$med) - min(h_medians_after$med) < 1e-10,
    info = ".normalizeMedian (labeled): H run medians should be equal after normalization"
)

# The same shift must also be applied to L rows (not just H rows)
l_medians_after <- out_labeled[LABEL == "L",
                                list(med = median(ABUNDANCE, na.rm = TRUE)),
                                by = "RUN"]
expect_true(
    max(l_medians_after$med) - min(l_medians_after$med) > 0.1,
    info = ".normalizeMedian (labeled): L run medians should not be equal anymore"
)

# The shift applied to L rows should equal the shift applied to H rows
# for the same run (both rows move by the same delta)
for (run_id in 1:4) {
    h_shift <- out_labeled[LABEL == "H" & RUN == run_id, ABUNDANCE[1]] -
               input_labeled[LABEL == "H" & RUN == run_id, ABUNDANCE[1]]
    l_shift <- out_labeled[LABEL == "L" & RUN == run_id, ABUNDANCE[1]] -
               input_labeled[LABEL == "L" & RUN == run_id, ABUNDANCE[1]]
    expect_true(
        abs(h_shift - l_shift) < 1e-10,
        info = paste(".normalizeMedian (labeled): L and H rows in run", run_id,
                     "should receive the same abundance shift")
    )
}

# L run medians before normalization are all equal (constant input),
# confirming that any equalization of L medians is driven by H, not L
l_medians_before <- input_labeled[LABEL == "L",
                                   list(med = median(ABUNDANCE, na.rm = TRUE)),
                                   by = "RUN"]
expect_true(
    max(l_medians_before$med) - min(l_medians_before$med) < 1e-10,
    info = ".normalizeMedian (labeled): L medians are already equal before normalization (sanity check)"
)

# Temporary columns should be cleaned up in the labeled case too
expect_false("ABUNDANCE_RUN" %in% colnames(out_labeled),
             info = ".normalizeMedian (labeled): ABUNDANCE_RUN column should be removed")
expect_false("ABUNDANCE_FRACTION" %in% colnames(out_labeled),
             info = ".normalizeMedian (labeled): ABUNDANCE_FRACTION column should be removed")
