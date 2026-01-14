library(data.table)

# Peptide dictionary
peptide_dict <- data.table(
  PeptideSequence = c("AAAAAAAAAAAAAAGAGAGAK", "AAAAAAAAAAAAVSR", "AAAAAAAAAAAAVSRD"),
  PrecursorCharge = c(3, 2, 3),
  PEPTIDE = c("AAAAAAAAAAAAAAGAGAGAK_3", "AAAAAAAAAAAAVSR_2", "AAAAAAAAAAAAVSRD_3")
)

# Standards (matching peptide sequences)
standards = c("AAAAAAAAAAAAAAGAGAGAK", "AAAAAAAAAAAAVSR", "AAAAAAAAAAAAVSRD")

# Create input data with 4 runs, varying fractions (1-3), 2 groups
set.seed(123)
input <- data.table(
  PROTEIN = rep(c("P1", "P2", "P3"), each = 48),  # 3 proteins
  PEPTIDE = rep(c("AAAAAAAAAAAAAAGAGAGAK_3", "AAAAAAAAAAAAVSR_2", "AAAAAAAAAAAAVSRD_3"), each = 48),
  TRANSITION = rep("NA_NA", 144),
  FEATURE = rep(c("AAAAAAAAAAAAAAGAGAGAK_3_NA_NA", "AAAAAAAAAAAAVSR_2_NA_NA", "AAAAAAAAAAAAVSRD_3_NA_NA"), each = 48),
  LABEL = rep("L", 144),
  GROUP_ORIGINAL = rep(rep(c("Control", "Treatment"), each = 24), 3),
  SUBJECT_ORIGINAL = rep(rep(c("Control", "Treatment"), each = 24), 3),
  RUN = rep(rep(1:4, each = 6), 6),  # 4 runs
  GROUP = rep(rep(c("Control", "Treatment"), each = 24), 3),
  SUBJECT = rep(rep(c("Control", "Treatment"), each = 24), 3),
  FRACTION = rep(rep(c(1, 1, 2, 2, 3, 3), 4), 6),  # Fractions 1-3
  INTENSITY = c(
    # Standard 1 (AAAAAAAAAAAAAAGAGAGAK_3)
    150000, 148000, 152000, 149000, 151000, 147000,  # Run 1
    145000, 143000, 147000, 144000, 146000, 142000,  # Run 2
    200000, 198000, 202000, 199000, 201000, 197000,  # Run 3
    195000, 193000, 197000, 194000, 196000, 192000,  # Run 4
    160000, 158000, 162000, 159000, 161000, 157000,  # Run 1 treatment
    155000, 153000, 157000, 154000, 156000, 152000,  # Run 2 treatment
    210000, 208000, 212000, 209000, 211000, 207000,  # Run 3 treatment
    205000, 203000, 207000, 204000, 206000, 202000,  # Run 4 treatment
    # Standard 2 (AAAAAAAAAAAAVSR_2)
    500000, 498000, 502000, 499000, 501000, 497000,
    495000, 493000, 497000, 494000, 496000, 492000,
    550000, 548000, 552000, 549000, 551000, 547000,
    545000, 543000, 547000, 544000, 546000, 542000,
    510000, 508000, 512000, 509000, 511000, 507000,
    505000, 503000, 507000, 504000, 506000, 502000,
    560000, 558000, 562000, 559000, 561000, 557000,
    555000, 553000, 557000, 554000, 556000, 552000,
    # Standard 3 (AAAAAAAAAAAAVSRD_3)
    250000, 248000, 252000, 249000, 251000, 247000,
    245000, 243000, 247000, 244000, 246000, 242000,
    300000, 298000, 302000, 299000, 301000, 297000,
    295000, 293000, 297000, 294000, 296000, 292000,
    260000, 258000, 262000, 259000, 261000, 257000,
    255000, 253000, 257000, 254000, 256000, 252000,
    310000, 308000, 312000, 309000, 311000, 307000,
    305000, 303000, 307000, 304000, 306000, 302000
  ),
  ANOMALYSCORES = rep(NA, 144),
  originalRUN = rep(paste0("Run", rep(1:4, each = 6)), 6)
)

# Calculate ABUNDANCE as log2(INTENSITY)
input[, ABUNDANCE := log2(INTENSITY)]

# Add some missing values for realism
input[c(5, 12, 23, 34, 45, 67, 89, 111), ABUNDANCE := NA]
input[c(5, 12, 23, 34, 45, 67, 89, 111), INTENSITY := NA]

# Test the function
output <- MSstats:::.normalizeGlobalStandards(input, peptide_dict, standards)
