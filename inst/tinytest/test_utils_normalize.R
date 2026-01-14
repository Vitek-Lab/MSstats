peptide_dict <- data.table::data.table(
  PeptideSequence = c( "AAAAAAAAAAAAAAGAGAGAK", "AAAAAAAAAAAAVSR", "AAAAAAAAAAAAVSRD"),
  PrecursorCharge = c( 3, 2, 3),
  PEPTIDE = c( "AAAAAAAAAAAAAAGAGAGAK_3", "AAAAAAAAAAAAVSR_2", "AAAAAAAAAAAAVSRD_3")
)

input <- data.table(
  PROTEIN = c( "P1", "P1", "P1", "P1", "P1", "P1" ),
  PEPTIDE = c( "PEP1", "PEP1", "PEP1", "PEP1", "PEP1", "PEP1" ),
  TRANSITION = c( "NA_NA", "NA_NA", "NA_NA", "NA_NA", "NA_NA", "NA_NA" ),
  FEATURE = c( "PEP1_2_NA_NA", "PEP1_3_NA_NA", "PEP1_2_NA_NA", "PEP1_2_NA_NA", "PEP1_3_NA_NA", "PEP1_2_NA_NA" ),
  LABEL = c( "L", "L", "L", "L", "L", "L" ),
  GROUP_ORIGINAL = c( "0hr", "0hr", "0hr", "0hr", "0hr", "0hr" ),
  SUBJECT_ORIGINAL = c( "0hr", "0hr", "0hr", "0hr", "0hr", "0hr" ),
  RUN = c( 1, 1, 1, 1, 1, 1 ),
  GROUP = c( "0hr", "0hr", "0hr", "0hr", "0hr", "0hr" ),
  SUBJECT = c( "0hr", "0hr", "0hr", "0hr", "0hr", "0hr" ),
  FRACTION = c( 1, 1, 1, 1, 1, 1 ),
  INTENSITY = c( 180697888, 674297.25, NA, 267224.25, NA, NA ),
  ANOMALYSCORES = c( NA, NA, NA, NA, NA, NA ),
  ABUNDANCE = c( 27.429, 19.363, NA, 18.028, NA, NA ),
  originalRUN = c( "Run1", "Run1", "Run1", "Run1", "Run1", "Run1" )
)

standards = c("AAAAAAAAAAAAAAGAGAGAK", "AAAAAAAAAAAAVSR", "AAAAAAAAAAAAVSRD")
output = MSstats:::.normalizeGlobalStandards(input, peptide_dict, standards)
