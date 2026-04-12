# .prepareForDataProcess tests ---------------------------------------------
input <- data.frame(
    PeptideModifiedSequence = c("Apeptide", "BPEPTIDE"),
    PrecursorCharge         = c(2L, 3L),
    FragmentIon             = c("b2", "y3"),
    ProductCharge           = c(1L, 1L),
    IsotopeLabelType        = c("H", "L")
)
result_peptideModifiedSequence <- MSstats:::.prepareForDataProcess(input)
expect_true(
    "PEPTIDESEQUENCE" %in% colnames(result_peptideModifiedSequence),
    info = "PEPTIDESEQUENCE should exist after renaming PEPTIDEMODIFIEDSEQUENCE"
)
expect_false(
    "PEPTIDEMODIFIEDSEQUENCE" %in% colnames(result_peptideModifiedSequence),
    info = "PEPTIDEMODIFIEDSEQUENCE should no longer exist after renaming"
)


make_input <- function(labels) {
    dt <- data.table::data.table(
        PeptideSequence = rep("PEPT", length(labels)),
        PrecursorCharge = rep(2L, length(labels)),
        FragmentIon     = rep("y3", length(labels)),
        ProductCharge   = rep(1L, length(labels)),
        IsotopeLabelType = labels
    )
    return(dt)
}

result <- MSstats:::.prepareForDataProcess(make_input(c("heavy", "light")))
expect_true(
    "H" %in% as.character(result$ISOTOPELABELTYPE),
    info = "label_map: 'heavy' must map to factor level 'H'"
)
expect_true(
    "L" %in% as.character(result$ISOTOPELABELTYPE),
    info = "label_map: 'light' must map to factor level 'L'"
)
expect_equal(
    levels(result$ISOTOPELABELTYPE), c("H", "L"),
    info = "label_map: factor levels must be exactly c('H', 'L')"
)


result_light <- MSstats:::.prepareForDataProcess(make_input(rep("light", 4)))
expect_equal(
    unique(as.character(result_light$ISOTOPELABELTYPE)), "L",
    info = "label_map: all 'light' input must produce all 'L' output"
)

result_hl <- MSstats:::.prepareForDataProcess(make_input(c("H", "L", "H", "L")))
expect_equal(
    as.character(result_hl$ISOTOPELABELTYPE), c("H", "L", "H", "L"),
    info = "label_map: 'H' / 'L' strings must still pass through unchanged"
)
expect_equal(
    levels(result_hl$ISOTOPELABELTYPE), c("H", "L"),
    info = "label_map: factor levels for H/L input must still be c('H', 'L')"
)

result_ll <- MSstats:::.prepareForDataProcess(make_input(rep("L", 5)))
expect_equal(
    unique(as.character(result_ll$ISOTOPELABELTYPE)), "L",
    info = "label_map: all 'L' input must produce all 'L' output"
)

result_other <- MSstats:::.prepareForDataProcess(make_input(rep("test", 5)))
expect_true(
    all(is.na(result_other$ISOTOPELABELTYPE)),
    info = "Other IsotopeLabelType maps to NA"
)
expect_equal(
    levels(result_other$ISOTOPELABELTYPE), c("H", "L"),
    info = "label_map: factor levels must be exactly c('H', 'L')"
)

result_na <- MSstats:::.prepareForDataProcess(make_input(rep(NA, 5)))
expect_true(
    all(is.na(result_na$ISOTOPELABELTYPE)),
    info = "NA IsotopeLabelType maps to NA"
)
expect_equal(
    levels(result_na$ISOTOPELABELTYPE), c("H", "L"),
    info = "label_map: factor levels must be exactly c('H', 'L')"
)
