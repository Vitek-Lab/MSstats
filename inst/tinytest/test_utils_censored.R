# Test .setCensoredByThreshold

# Construct dataset ensuring n_obs == 2 for all PROTEIN-FEATURE combinations
dt <- data.table::data.table(
  PROTEIN = rep(c("P1", "P2"), each = 12),
  FEATURE = rep(c("F1", "F1", "F2", "F2", "F1", "F2"), 4),
  LABEL   = rep(c("L", "L", "L", "L", "H", "H"), 4),
  RUN     = rep(c("R1", "R2", "R1", "R2", "R1", "R2",
                  "R3", "R4", "R3", "R4", "R3", "R4"), 2),
  newABUNDANCE = c(
    # P1
    1.5, NA, 2.0, 2.2, 1.4, 1.6, 1.5, 1.9, 2.0, 2.2, 1.4, 1.6,
    # P2
    0, 4.0, 2.5, 2.7, 3.0, 3.2, 1.7, 4.0, 2.5, 2.7, 3.0, 3.2
  ),
  censored = c(
    # P1
    FALSE, TRUE, FALSE, FALSE, FALSE, FALSE,FALSE, TRUE, FALSE, FALSE, FALSE, FALSE,
    # P2
    TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,TRUE, FALSE, FALSE, FALSE, FALSE, FALSE
  )
)

# === Compute n_obs and friends: censored_symbol = "NA" case ===
dt[, nonmissing := LABEL == "L" & !is.na(newABUNDANCE)]
dt[, n_obs := sum(nonmissing), by = .(PROTEIN, FEATURE)]
dt[, nonmissing := ifelse(n_obs <= 1, FALSE, nonmissing)]
dt[, n_obs_run := sum(nonmissing), by = .(PROTEIN, RUN)]
dt[, total_features := uniqueN(FEATURE), by = "PROTEIN"]

# === Run NA-based test ===
dt_na <- copy(dt)
MSstats:::.setCensoredByThreshold(dt_na, censored_symbol = "NA", remove50missing = FALSE)

# Check imputation for P1-F1-L (should be 0.99 * 1.5)
expected_val_p1 <- 0.99 * 1.5
imputed_val_p1 <- dt_na[
  PROTEIN == "P1" & FEATURE == "F1" & LABEL == "L" & RUN == "R2",
  newABUNDANCE
]
expect_equal(imputed_val_p1, expected_val_p1)

# === Compute n_obs and friends: censored_symbol = "0" case ===
dt[, nonmissing := LABEL == "L" & !is.na(newABUNDANCE) & newABUNDANCE != 0]
dt[, n_obs := sum(nonmissing), by = .(PROTEIN, FEATURE)]
dt[, nonmissing := ifelse(n_obs <= 1, FALSE, nonmissing)]
dt[, n_obs_run := sum(nonmissing), by = .(PROTEIN, RUN)]
dt[, total_features := uniqueN(FEATURE), by = "PROTEIN"]

# === Run 0-based test ===
dt_zero <- copy(dt)
MSstats:::.setCensoredByThreshold(dt_zero, censored_symbol = "0", remove50missing = FALSE)

# Check imputation for P2-F1-L
expected_val_p2 <- 0.99 * 1.7
imputed_val_p2 <- dt_zero[
  PROTEIN == "P2" & FEATURE == "F1" & LABEL == "L" & RUN == "R1",
  newABUNDANCE
]
expect_equal(imputed_val_p2, expected_val_p2)