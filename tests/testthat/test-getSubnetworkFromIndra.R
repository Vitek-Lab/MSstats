test_that("getSubnetworkFromIndra works correctly", {
    input = data.table::fread(
        system.file("tinytest/processed_data/groupComparisonModel.csv", package = "MSstats")
    )
    local_mocked_bindings(.callIndraCogexApi = function(x) {
        return(readRDS(system.file("tinytest/processed_data/indraResponse.rds", package = "MSstats")))
    })
    subnetwork = getSubnetworkFromIndra(input)
    expect_equal(nrow(subnetwork$nodes), 4)
    expect_equal(nrow(subnetwork$edges), 16)
})

test_that("getSubnetworkFromIndra with pvalue filter works correctly", {
    input = data.table::fread(
        system.file("tinytest/processed_data/groupComparisonModel.csv", package = "MSstats")
    )
    local_mocked_bindings(.callIndraCogexApi = function(x) {
        return(readRDS(system.file("tinytest/processed_data/indraResponse.rds", package = "MSstats")))
    })
    subnetwork = getSubnetworkFromIndra(input, pvalue_cutoff = 0.05)
    expect_equal(nrow(subnetwork$nodes), 3)
    expect_equal(nrow(subnetwork$edges), 16)
})