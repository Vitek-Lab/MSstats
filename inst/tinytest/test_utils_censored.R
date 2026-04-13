# Test .setCensoredByThreshold
dt_na <- data.table::data.table(
    PROTEIN = c("P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2"),
    FEATURE = c("F1", "F1", "F2", "F2", "F1", "F2", "F1", "F1", "F2", "F2", "F1", "F2", "F3", "F3", "F4", "F4", "F3", "F4", "F3", "F3", "F4", "F4", "F3", "F4"),
    LABEL = c("L", "L", "L", "L", "H", "H", "L", "L", "L", "L", "H", "H", "L", "L", "L", "L", "H", "H", "L", "L", "L", "L", "H", "H"),
    RUN = c("R1", "R2", "R1", "R2", "R1", "R2", "R3", "R4", "R3", "R4", "R3", "R4", "R1", "R2", "R1", "R2", "R1", "R2", "R3", "R4", "R3", "R4", "R3", "R4"),
    newABUNDANCE = c(1.5, NA, 2, 2.2, 1.4, 1.6, 1.5, 1.9, 2, 2.2, 1.4, 1.6, 0, 4, 2.5, NA, 3, 3.2, NA, 4, NA, NA, 3, 3.2),
    censored = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE),
    nonmissing = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
    n_obs = c(3, 3, 4, 4, 3, 4, 3, 3, 4, 4, 3, 4, 3, 3, 1, 1, 3, 1, 3, 3, 1, 1, 3, 1),
    n_obs_run = c(2, 1, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1),
    total_features = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
)

# === Run NA-based test ===
MSstats:::.setCensoredByThreshold(dt_na, censored_symbol = "NA", remove50missing = FALSE)

# Check imputation for P1-F1-L (should be 0.99 * 1.5)
expected_val_p1 <- 0.99 * 1.5
imputed_val_p1 <- dt_na[
  PROTEIN == "P1" & FEATURE == "F1" & LABEL == "L" & RUN == "R2",
  newABUNDANCE
]
expect_equal(imputed_val_p1, expected_val_p1)
non_imputed_val_p2 <- dt_na[
    PROTEIN == "P2" & FEATURE == "F3" & LABEL == "L" & RUN == "R3",
    newABUNDANCE
]
expect_true(is.na(non_imputed_val_p2))
non_imputed_val_p2_f4 <- dt_na[
    PROTEIN == "P2" & FEATURE == "F4" & LABEL == "L" & RUN == "R2",
    newABUNDANCE
]
expect_true(is.na(non_imputed_val_p2_f4))

dt_zero <- data.table::data.table(
  PROTEIN = c("P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", 
              "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2"),
  FEATURE = c("F1", "F1", "F2", "F2", "F1", "F2", "F1", "F1", "F2", "F2", "F1", "F2", 
              "F1", "F1", "F2", "F2", "F1", "F2", "F1", "F1", "F2", "F2", "F1", "F2"),
  LABEL = c("L", "L", "L", "L", "H", "H", "L", "L", "L", "L", "H", "H", 
            "L", "L", "L", "L", "H", "H", "L", "L", "L", "L", "H", "H"),
  RUN = c("R1", "R2", "R1", "R2", "R1", "R2", "R3", "R4", "R3", "R4", "R3", "R4", 
          "R1", "R2", "R1", "R2", "R1", "R2", "R3", "R4", "R3", "R4", "R3", "R4"),
  newABUNDANCE = c(1.5, NA, 2.0, 2.2, 1.4, 1.6, 1.5, 1.9, 2.0, 2.2, 1.4, 1.6, 
                   0.0, 4.0, 2.5, 2.7, 3.0, 3.2, 1.7, 4.0, 2.5, 2.7, 3.0, 3.2),
  censored = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
               TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  nonmissing = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, 
                 FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  n_obs = c(3, 3, 4, 4, 3, 4, 3, 3, 4, 4, 3, 4, 
            3, 3, 4, 4, 3, 4, 3, 3, 4, 4, 3, 4),
  n_obs_run = c(2, 1, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 
                1, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2),
  total_features = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
                     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
)
MSstats:::.setCensoredByThreshold(dt_zero, censored_symbol = "0", remove50missing = FALSE)

# Check imputation for P2-F1-L
expected_val_p2 <- 0.99 * 1.7
imputed_val_p2 <- dt_zero[
  PROTEIN == "P2" & FEATURE == "F1" & LABEL == "L" & RUN == "R1",
  newABUNDANCE
]
expect_equal(imputed_val_p2, expected_val_p2)


# Test MSstatsHandleMissing
make_cens_input <- function() {
  data.table::data.table(
    PROTEIN   = rep("P1", 8),
    FEATURE   = rep("F1", 8),
    LABEL     = c("L","L","L","L","H","H","H","H"),
    RUN       = rep(c("R1","R2","R3","R4"), 2),
    GROUP     = rep(c("G1","G1","G2","G2"), 2),
    FRACTION  = rep(1L, 8),
    INTENSITY = c(1L, 1L, 100L, 200L,   # L: first two will be censored
                  1L, 1L, 100L, 200L),  # H: same intensities but must NOT be censored
    ABUNDANCE = c(0, 0, log2(100), log2(200),
                  0, 0, log2(100), log2(200)),
    ref       = c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE)
  )
}

cens_in <- make_cens_input()
out_cens <- MSstatsHandleMissing(
  cens_in,
  summary_method  = "TMP",
  impute          = TRUE,
  missing_symbol  = "0",
  censored_cutoff = NULL
)

# H rows must never be flagged as censored regardless of intensity
expect_false(
  any(out_cens[LABEL == "H", censored]),
  info = "H rows (ref=TRUE) must never be flagged censored (use_for_analysis=FALSE)"
)

# L rows with intensity==1 must be flagged censored
expect_true(
  any(out_cens[LABEL == "L" & INTENSITY == 1L, censored]),
  info = "L rows with intensity=1 (use_for_analysis=TRUE) must be flagged censored"
)

# L rows with intensity>1 must NOT be flagged censored
expect_false(
  any(out_cens[LABEL == "L" & INTENSITY > 1L, censored]),
  info = "L rows with high intensity must not be flagged censored"
)


# Tests for .getNonMissingFilter -----------------------------------------------

dt_nmf_no_ref <- data.table::data.table(
    newABUNDANCE = c(1.5, 0.0, NA_real_, 2.0)
)
result_no_ref <- MSstats:::.getNonMissingFilter(dt_nmf_no_ref, impute = FALSE, censored_symbol = "0")
expect_equal(
    result_no_ref, c(TRUE, FALSE, FALSE, TRUE),
    info = ".getNonMissingFilter: without 'ref' column all rows have use_for_analysis=TRUE; zero and NA excluded"
)

dt_nmf_with_ref <- data.table::data.table(
    newABUNDANCE = c(1.5, 2.0, 3.0, 1.0),
    ref          = c(TRUE, FALSE, FALSE, TRUE)
)
result_with_ref <- MSstats:::.getNonMissingFilter(dt_nmf_with_ref, impute = FALSE, censored_symbol = "0")
expect_equal(
    result_with_ref, c(FALSE, TRUE, TRUE, FALSE),
    info = ".getNonMissingFilter: ref=TRUE rows must be excluded (use_for_analysis=FALSE) regardless of abundance"
)

dt_nmf_ref_valid <- data.table::data.table(
    newABUNDANCE = c(5.0, 3.0),
    ref          = c(TRUE, FALSE)
)
result_ref_valid <- MSstats:::.getNonMissingFilter(dt_nmf_ref_valid, impute = FALSE, censored_symbol = "0")
expect_false(
    result_ref_valid[1],
    info = ".getNonMissingFilter: ref=TRUE row with valid abundance must still yield FALSE"
)
expect_true(
    result_ref_valid[2],
    info = ".getNonMissingFilter: ref=FALSE row with valid abundance must yield TRUE"
)

dt_nmf_impute_na <- data.table::data.table(
    newABUNDANCE = c(1.5, 0.0, NA_real_, 2.0)
)
result_impute_na <- MSstats:::.getNonMissingFilter(dt_nmf_impute_na, impute = TRUE, censored_symbol = "NA")
expect_equal(
    result_impute_na, c(TRUE, TRUE, FALSE, TRUE),
    info = ".getNonMissingFilter: impute=TRUE, censored_symbol='NA' treats zero as non-missing"
)

dt_nmf_impute_zero <- data.table::data.table(
    newABUNDANCE = c(1.5, 0.0, NA_real_, 2.0)
)
result_impute_zero <- MSstats:::.getNonMissingFilter(dt_nmf_impute_zero, impute = TRUE, censored_symbol = "0")
expect_equal(
    result_impute_zero, c(TRUE, FALSE, FALSE, TRUE),
    info = ".getNonMissingFilter: impute=TRUE, censored_symbol='0' treats zero as missing"
)

dt_nmf_ref_impute_na <- data.table::data.table(
    newABUNDANCE = c(1.5, 2.0, 0.0, 3.0),
    ref          = c(TRUE, FALSE, FALSE, TRUE)
)
result_ref_impute_na <- MSstats:::.getNonMissingFilter(dt_nmf_ref_impute_na, impute = TRUE, censored_symbol = "NA")
expect_equal(
    result_ref_impute_na, c(FALSE, TRUE, TRUE, FALSE),
    info = ".getNonMissingFilter: ref=TRUE rows excluded even when impute=TRUE, censored_symbol='NA'"
)