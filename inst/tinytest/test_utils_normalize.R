# Test for .normalizeGlobalStandards function ----------------------------------

# Setup test data ---------------------------------------------------------------

create_peptide_dictionary <- function() {
    data.table::data.table(
        PeptideSequence = c("AAAAAAAAAAAAAAGAGAGAK", "AAAAAAAAAAAAVSR", "AAAAAAAAAAAAVSRD"),
        PrecursorCharge = c(3, 2, 3),
        PEPTIDE = c("AAAAAAAAAAAAAAGAGAGAK_3", "AAAAAAAAAAAAVSR_2", "AAAAAAAAAAAAVSRD_3")
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
        PROTEIN = rep(c("P1", "P2", "P3"), each = n_runs),
        PEPTIDE = rep(c("AAAAAAAAAAAAAAGAGAGAK_3", "AAAAAAAAAAAAVSR_2", "AAAAAAAAAAAAVSRD_3"), 
                      each = n_runs),
        TRANSITION = rep("NA_NA", n_proteins * n_runs),
        FEATURE = rep(c("AAAAAAAAAAAAAAGAGAGAK_3_NA_NA", 
                        "AAAAAAAAAAAAVSR_2_NA_NA", 
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
    peptide2_intensities <- rep(262144, 48)
    peptide3_intensities <- rep(262144, 48)
    
    input <- create_test_input(standard_intensities, peptide2_intensities, peptide3_intensities)
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
    peptide2_intensities <- rep(262144, 48)
    peptide3_intensities <- rep(262144, 48)
    
    input <- create_test_input(standard_intensities, peptide2_intensities, peptide3_intensities)
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