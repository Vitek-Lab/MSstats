library(mockery)
library(testthat)
library(MSstats)

QuantData <- dataProcess(SRMRawData, use_log_file = FALSE)

comparison <- matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
row.names(comparison) <- "T7-T1"
colnames(comparison) <- unique(QuantData$ProteinLevelData$GROUP)

# Tests for differentially abundant proteins with models:
# label-based SRM experiment with expanded scope of biological replication.

testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData,use_log_file = FALSE)



test_that('Validate groupComparisonQCPlots call', {
  
  mock_modelBasedQCPlots <- mock()
  
  stub(groupComparisonQCPlots, 'modelBasedQCPlots', mock_modelBasedQCPlots)
  
  groupComparisonQCPlots(
    data = testResultOneComparison, 
    type = "QQPlots", 
    axis.size = 12, 
    dot.size = 3, 
    width = 7, 
    height = 5, 
    which.Protein = "P0A894", 
    address = FALSE
  )
  
  expect_called(mock_modelBasedQCPlots, 1)
  
  args <- mock_args(mock_modelBasedQCPlots)
  
  expect_equal(args[[1]][[1]], testResultOneComparison)   # data
  expect_equal(args[[1]][[2]], "QQPlots")                 # type
  expect_equal(args[[1]][[3]], 12)                        # axis.size
  expect_equal(args[[1]][[4]], 3)                         # dot.size
  expect_equal(args[[1]][[5]], 7)                         # width
  expect_equal(args[[1]][[6]], 5)                         # height
  expect_equal(args[[1]][[7]], "P0A894")                  # which.Protein
  expect_equal(args[[1]][[8]], FALSE)                     # address
})
