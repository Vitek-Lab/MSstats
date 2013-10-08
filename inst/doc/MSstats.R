### R code from vignette source 'MSstats.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library(MSstats)


###################################################
### code chunk number 2: preliminaries (eval = FALSE)
###################################################
## ?MSstats


###################################################
### code chunk number 3: preliminaries (eval = FALSE)
###################################################
## ?RawData


###################################################
### code chunk number 4: MSstats.Rnw:384-385
###################################################
head(RawData)


###################################################
### code chunk number 5: MSstats.Rnw:410-411 (eval = FALSE)
###################################################
## ?dataProcess


###################################################
### code chunk number 6: MSstats.Rnw:416-417
###################################################
QuantData<-dataProcess(RawData)


###################################################
### code chunk number 7: MSstats.Rnw:422-423
###################################################
head(QuantData)


###################################################
### code chunk number 8: MSstats.Rnw:446-447 (eval = FALSE)
###################################################
## ?dataProcessPlots


###################################################
### code chunk number 9: QCplot (eval = FALSE)
###################################################
## dataProcessPlots(data=QuantData,type="QCPlot")


###################################################
### code chunk number 10: Profileplot (eval = FALSE)
###################################################
## dataProcessPlots(data=QuantData,type="ProfilePlot")


###################################################
### code chunk number 11: Conditionplot (eval = FALSE)
###################################################
## dataProcessPlots(data=QuantData,type="ConditionPlot")


###################################################
### code chunk number 12: MSstats.Rnw:536-537 (eval = FALSE)
###################################################
## ?groupComparison


###################################################
### code chunk number 13: MSstats.Rnw:547-548
###################################################
levels(QuantData$GROUP_ORIGINAL)


###################################################
### code chunk number 14: MSstats.Rnw:551-553
###################################################
comparison<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
row.names(comparison)<-"T7-T1"


###################################################
### code chunk number 15: MSstats.Rnw:556-558
###################################################
testResultOneComparison<-groupComparison(contrast.matrix=comparison, data=QuantData)
testResultOneComparison$ComparisonResult


###################################################
### code chunk number 16: MSstats.Rnw:568-575
###################################################
comparison1<-matrix(c(-1,0,1,0,0,0,0,0,0,0),nrow=1)
comparison2<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
comparison3<-matrix(c(-1,0,0,0,0,0,0,0,1,0),nrow=1)
comparison<-rbind(comparison1,comparison2, comparison3)
row.names(comparison)<-c("T3-T1","T7-T1","T9-T1")
testResultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=QuantData)
testResultMultiComparisons$ComparisonResult 


###################################################
### code chunk number 17: MSstats.Rnw:593-594 (eval = FALSE)
###################################################
## ?modelBasedQCPlots


###################################################
### code chunk number 18: MSstats.Rnw:597-599 (eval = FALSE)
###################################################
## modelBasedQCPlots(data=testResultMultiComparisons$ModelQC, type="ResidualPlots", 
##                   which.Protein="PMG2", address=FALSE)


###################################################
### code chunk number 19: MSstats.Rnw:602-604 (eval = FALSE)
###################################################
## modelBasedQCPlots(data=testResultMultiComparisons$ModelQC,type="QQPlots",
##                   which.Protein="PMG2", address=FALSE)


###################################################
### code chunk number 20: MSstats.Rnw:631-632 (eval = FALSE)
###################################################
## ?groupComparisonPlots


###################################################
### code chunk number 21: MSstats.Rnw:650-651
###################################################
head(testResultMultiComparisons$ComparisonResult)


###################################################
### code chunk number 22: groupComparisoneplot1 (eval = FALSE)
###################################################
## groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult,
##                      type="VolcanoPlot",which.Comparison=c("T7-T1"), address=FALSE)


###################################################
### code chunk number 23: groupComparisoneplot2 (eval = FALSE)
###################################################
## groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult,
##                      type="VolcanoPlot",FCcutoff=70, ylimUp=100, 
##                      ProteinName=FALSE, which.Comparison=c("T7-T1"), address=FALSE)


###################################################
### code chunk number 24: groupComparisoneplot3 (eval = FALSE)
###################################################
## groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult,
##                      type="Heatmap",address=FALSE)


###################################################
### code chunk number 25: groupComparisoneplot4 (eval = FALSE)
###################################################
## groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult,
##                      type="Heatmap",FCcutoff=70,address=FALSE)


###################################################
### code chunk number 26: groupComparisoneplot5 (eval = FALSE)
###################################################
## groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult,
##                      type="ComparisonPlot")


###################################################
### code chunk number 27: MSstats.Rnw:784-785 (eval = FALSE)
###################################################
## ?designSampleSize


###################################################
### code chunk number 28: MSstats.Rnw:794-796
###################################################
designSampleSize(data=QuantData,numSample=TRUE,numPep=3,numTran=4,power=0.8,
                 desiredFC=c(1.25,1.75),FDR=0.05)


###################################################
### code chunk number 29: MSstats.Rnw:802-804
###################################################
designSampleSize(data=QuantData,numSample=2,numPep=3,numTran=4,power=TRUE,
                 desiredFC=c(1.25,1.75),FDR=0.05)


###################################################
### code chunk number 30: MSstats.Rnw:817-818 (eval = FALSE)
###################################################
## ?designSampleSizePlots


###################################################
### code chunk number 31: MSstats.Rnw:822-826 (eval = FALSE)
###################################################
## # Minimal number of biological replicates per condition
## result.sample<-designSampleSize(data=QuantData,numSample=TRUE,numPep=3,numTran=4,power=0.8,
##                                 desiredFC=c(1.25,1.75),FDR=0.05)
## designSampleSizePlots(data=result.sample)


###################################################
### code chunk number 32: MSstats.Rnw:836-840 (eval = FALSE)
###################################################
## # Power
## result.power<-designSampleSize(data=QuantData,numSample=2,numPep=3,numTran=4,power=TRUE,
##                                desiredFC=c(1.25,1.75),FDR=0.05)
## designSampleSizePlots(data=result.power)


###################################################
### code chunk number 33: MSstats.Rnw:864-865 (eval = FALSE)
###################################################
## ?quantification


###################################################
### code chunk number 34: MSstats.Rnw:872-875
###################################################
#  (1): Sample quantification
subQuant<-quantification(QuantData)
head(subQuant)


###################################################
### code chunk number 35: MSstats.Rnw:878-881
###################################################
# (2): Group quantification
groupQuant<-quantification(QuantData, type="Group")
head(groupQuant)


