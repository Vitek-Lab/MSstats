
if (requireNamespace("tinytest", quietly = TRUE)) {
    MSstatsConvert::MSstatsLogsSettings(FALSE, FALSE, TRUE)
    tinytest::test_package("MSstatsConvert")
}

