# Test .setCensoredByThreshold
dt <- data.table(
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
                 TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  n_obs = c(3, 3, 4, 4, 3, 4, 3, 3, 4, 4, 3, 4, 
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4),
  n_obs_run = c(2, 1, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
  total_features = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
                     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
)

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

dt <- data.table(
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
dt_zero <- copy(dt)
MSstats:::.setCensoredByThreshold(dt_zero, censored_symbol = "0", remove50missing = FALSE)

# Check imputation for P2-F1-L
expected_val_p2 <- 0.99 * 1.7
imputed_val_p2 <- dt_zero[
  PROTEIN == "P2" & FEATURE == "F1" & LABEL == "L" & RUN == "R1",
  newABUNDANCE
]
expect_equal(imputed_val_p2, expected_val_p2)
