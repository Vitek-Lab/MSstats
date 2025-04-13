
library(MSstats)
library(data.table)

msstats_input = fread("C:\\Users\\Devon Kohler\\Downloads\\test_input.csv")

msstats_input$Intensity = ifelse(msstats_input$Intensity < 2, NA, msstats_input$Intensity)

summarized = dataProcess(msstats_input[msstats_input$ProteinName =="Q9BW30", ], 
                         normalization=FALSE,
                         featureSubset = "topN",
                         n_top_feature = 100,
                         MBimpute=TRUE,
                         summaryMethod="linear",
                         numberOfCores = 1)