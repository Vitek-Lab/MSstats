#############################################
#############################################
## Part 0 : Transform data
#############################################
#############################################

#######################################
#### MSnSet -> input for MSstats
#######################################

transformMSnSetToMSstats <- function(ProteinName,PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Bioreplicate,Run, Condition, data) {

    if (!inherits(data, "MSnSet"))
        stop("only MSnSet class can be converted to input format for MSstats.")

    ## if the data is an expression set, extract the different data components and reshape to "long" format 

    ## extract the three data components
    sampleData <- pData(data)
    featureData <- fData(data)
    expressionMatrix <- as.data.frame(exprs(data))
    
    ## extract the variable names 
    sample.variables <- dimnames(sampleData)[[2]]
    feature.variables <- dimnames(featureData)[[2]] ## ProteinAccession, PeptideSequence, charge
    
    
    ## create a patient variable by which to merge with intensities 
    sampleData$pData_rownames <- rownames(sampleData)
    
    ## creat a feature variable by which to merge with intensities
    featureData$fData_rownames <- rownames(featureData)
    
    ## transform the matrix of intensities so that each row is a sample and feature combination
    long.abundances <- reshape(expressionMatrix, idvar = "fData_rownames", ids = rownames(expressionMatrix), times = colnames(expressionMatrix), timevar = "pData_rownames", varying = list(dimnames(sampleData)[[1]]), v.names = "ABUNDANCE", direction = "long")
    rownames(long.abundances) <- NULL

    ## merge with featureData and patientData to get the feature -> protein and the bio.rep -> group mappings
    long.abundances$ABUNDANCE <- as.numeric(as.character(long.abundances$ABUNDANCE))
    long.abundances.2 <- merge(long.abundances, featureData, by = "fData_rownames")
    final.data <- merge(long.abundances.2, sampleData, by = "pData_rownames")

#### extract required information
    ## default : Protein="ProteinAccession",PeptideSequence="PeptideSequence", PrecursorCharge="charge", FragmentIon, ProductCharge ,IsotopeLabelType="mz", Bioreplicate=NA,Run=NA,

    if(missing(ProteinName)) ProteinName<-"ProteinAccession"
    if(missing(PeptideSequence)) PeptideSequence<-"PeptideSequence"
    if(missing(PrecursorCharge)) PrecursorCharge<-"charge"
    if(missing(FragmentIon)){
	final.data$FragmentIon<-NA
	FragmentIon<-"FragmentIon"
    }
    if(missing(ProductCharge)){
	final.data$ProductCharge<-NA
	ProductCharge<-"ProductCharge"
    } 
    if(missing(IsotopeLabelType)) IsotopeLabelType<-"mz"
    if(missing(Bioreplicate)) Bioreplicate<-"pData_rownames"
    if(missing(Run)) Run<-"file"
    if(missing(Condition)) stop("The condition arguments must be specified")


    ## if there are any missing variable name, warn it and stop
    check.name<-c(ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,Bioreplicate,Run,Condition)

    diff.name<-setdiff(check.name, colnames(final.data))
    if(length(diff.name)>0)
	stop(paste("error : Please check the variable name. The provided variable name",paste(diff.name,collapse=","), "is not present in the data set.",sep=" "))


    ## define the group column 

    if(length(Condition) > 1){
	group.variable <- final.data[, which(colnames(final.data) == Condition[1])]
	
	for(i in 2:length(Condition)){
            var <- final.data[, which(colnames(final.data) == Condition[i])]
            group.variable <- paste(group.variable, var, sep = ".")
	}
	final.data$group.column <- group.variable
    } else final.data$group.column <- final.data[, which(colnames(final.data) == Condition)]

    ## combine into a data frame
    final.data.2 <- final.data[,c(ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,"group.column",Bioreplicate,Run,"ABUNDANCE")]

    colnames(final.data.2)<-c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Run","Intensity")
    rownames(final.data.2) <- NULL

    return(final.data.2)
}


#######################################
#### MSnSet -> input for MSstats
#######################################
transformMSstatsToMSnSet<-function(data){
	data<-droplevels(data)
	
	colnames(data)<-toupper(colnames(data))
	
	## the expression matrix
	xx<-with(data, by(INTENSITY, list(ISOTOPELABELTYPE, BIOREPLICATE, CONDITION,RUN), c))
	xx<-do.call(cbind, xx)
	
	## sample annotation
	pd<-data[!duplicated(data[, c("ISOTOPELABELTYPE", "CONDITION", "BIOREPLICATE", "RUN")]), c("ISOTOPELABELTYPE", "CONDITION", "BIOREPLICATE", "RUN")]

	## feature annotation
	fd<-data[!duplicated(data[, c("PROTEINNAME", "PEPTIDESEQUENCE", "PRECURSORCHARGE", "FRAGMENTION", "PRODUCTCHARGE")]), c("PROTEINNAME", "PEPTIDESEQUENCE", "PRECURSORCHARGE", "FRAGMENTION", "PRODUCTCHARGE")]


### need to make as MSnSet class
e <- MSnSet(xx, fd, pd)
	
sampleNames(e)<-paste(e$ISOTOPELABELTYPE, e$CONDITION, e$BIOREPLICATE,e$RUN, sep = ".")
featureNames(e)<-paste(fData(e)$PEPTIDESEQUENCE,fData(e)$PRECURSORCHARGE, fData(e)$FRAGMENTION,fData(e)$PRODUCTCHARGE, sep = ".")

        if (validObject(e))
            return(e)
	
}




#############################################
#############################################
#############################################
#############################################
# Part 1 dataProcess
#############################################
#############################################
#############################################
#############################################

dataProcess<-function(raw,logTrans=2, normalization="constant",nameStandards=NULL, betweenRunInterferenceScore=FALSE, address=""){
	
	## save process output in each step
	processout<-NULL
	
	# check right input
	requiredinput<-c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Run","Intensity")
	
	if(length(setdiff(toupper(requiredinput),toupper(colnames(raw))))!=0){
		processout<-rbind(processout,c(paste("The required input : ",setdiff(toupper(requiredinput),toupper(colnames(raw)))," are not provided in input - stop")))
		colnames(processout)<-c("Process-dataProcess function")
		write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)
		
		stop(message("Please check the required input. The required input needs (ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity)"))
	}else{
		processout<-rbind(processout,c("The required input : provided - okay"))
		colnames(processout)<-c("Process-dataProcess function")
		write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)
	}
	
	#  check whether class of intensity is factor or chaterer, if yes, neec to chage as numeric
    if(is.factor(raw$Intensity) | is.character(raw$Intensity)){	
		suppressWarnings(raw$Intensity<-as.numeric(as.character(raw$Intensity)))
	}
	
	# check whether the intensity has 0 value or negative value
	if(length(which(raw$Intensity<=0))>0){
		processout<-rbind(processout,c("There are some intensities which are zero or negative values. need to change them. - stop"))
		write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)
		
		stop(message("Intensity has 0 or negative values.\n"))
	}	
	
#######################################	
#### here, need to get standard protein name
## column name : standardtype..
## what value it has, normzalition, unique(proteinname)
## if normalition== "standard" & no normalizaion selection, error message
#######################################	

	
### For Skyline
### required cols : ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, and Condition, BioReplicate, Run, Intensity
	
	## make letters case-insensitive
	colnames(raw)<-toupper(colnames(raw))
	raw.temp<-raw[,c("PROTEINNAME", "PEPTIDESEQUENCE", "PRECURSORCHARGE", "FRAGMENTION", "PRODUCTCHARGE", "ISOTOPELABELTYPE", "CONDITION", "BIOREPLICATE", "RUN", "INTENSITY")]

	rm(raw)

	
	## assign peptide, transition
	raw.temp<-data.frame(raw.temp,PEPTIDE=paste(raw.temp$PEPTIDESEQUENCE,raw.temp$PRECURSORCHARGE,sep="_"), TRANSITION=paste(raw.temp$FRAGMENTION, raw.temp$PRODUCTCHARGE,sep="_"))
	
	if(length(unique(raw.temp$ISOTOPELABELTYPE))>2){
		processout<-rbind(processout,c("There are more than two levels of labeling. So far, only label-free or reference-labeled experiment are supported. - stop"))
		write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)
		
		stop(message("Statistical tools in MSstats are only proper for label-free or with reference peptide experiments."))
	}

	## change light, heavy -> L,H
	raw.temp$ISOTOPELABELTYPE<-factor(raw.temp$ISOTOPELABELTYPE)
	if(nlevels(raw.temp$ISOTOPELABELTYPE)==2){
	levels(raw.temp$ISOTOPELABELTYPE)<-c("H","L")
	}
	if(nlevels(raw.temp$ISOTOPELABELTYPE)==1){
	levels(raw.temp$ISOTOPELABELTYPE)<-c("L")
	}

	raw.temp<-raw.temp[,c("PROTEINNAME","PEPTIDE","TRANSITION","ISOTOPELABELTYPE","CONDITION","BIOREPLICATE","RUN","INTENSITY")]
	
	colnames(raw.temp)<-c("Protein","Peptide","Transition","Label","Condition","Sample","Run","Intensity")




##############################
## create work data for quant analysis

work<-data.frame(PROTEIN=raw.temp[,"Protein"],PEPTIDE=raw.temp[,"Peptide"],TRANSITION=raw.temp[,"Transition"],FEATURE=paste(raw.temp[,"Peptide"],raw.temp[,"Transition"],sep="_"),LABEL=raw.temp[,"Label"],GROUP_ORIGINAL=raw.temp[,"Condition"],SUBJECT_ORIGINAL=raw.temp[,"Sample"],RUN=raw.temp[,"Run"],GROUP=0,SUBJECT=0)

work$GROUP_ORIGINAL<-factor(work$GROUP_ORIGINAL)
work$SUBJECT_ORIGINAL<-factor(work$SUBJECT_ORIGINAL,levels=unique(work$SUBJECT_ORIGINAL))

work[work$LABEL=="L","GROUP"]<-work[work$LABEL=="L","GROUP_ORIGINAL"]
work[work$LABEL=="L","SUBJECT"]<-work[work$LABEL=="L","SUBJECT_ORIGINAL"]

work<-data.frame(work,SUBJECT_NESTED=paste(work[,"GROUP"],work[,"SUBJECT"],sep="."))


processout<-rbind(processout,c("New input format : made new columns for analysis - okay"))
write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)

############## log transformation
## check logTrans is 2,10 or not
if(logTrans!=2 & logTrans!=10){
	processout<-rbind(processout,c("Logarithm transformation : log2 or log10 only - stop"))
write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)
	
	stop(message("only log2 or log10 are posssible.\n"))
}

## based on logTrans option, assign log transformation
# remove log2 or log10 intensity
if(logTrans==2){
work<-data.frame(work,INTENSITY=raw.temp$Intensity,ABUNDANCE=log2(raw.temp$Intensity))
	} else if(logTrans==10){
work<-data.frame(work,INTENSITY=raw.temp$Intensity,ABUNDANCE=log10(raw.temp$Intensity))
		} 	

processout<-rbind(processout,c(paste("Logarithm transformation: log",logTrans," transformation is done - okay", sep="")))
write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)


##############################
## check no value for some feature : at lease NA is needed.

work.l<-work[work$LABEL=="L",]

structure = tapply ( work.l$ABUNDANCE, list ( work.l$FEATURE, work.l$SUBJECT_ORIGINAL ) , function ( x ) length ( x ) ) 

flag = any ( apply ( structure, 1, function ( x ) sum ( is.na ( x ) ) ) > 0 ) 

if ( flag ){
	processout<-rbind(processout,c("Missing feature intensities should be given by 'NA' in the abundance column - stop"))
write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)

	 stop ( message("Missing feature intensities should be given by 'NA' in the abundance column." ))
}else{
	processout<-rbind(processout,c("Balanced data format with NA for missing feature intensities - okay"))
write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)

}



##############################
## factorize GROUP, SUBJECT, GROUP_ORIGINAL, SUBJECT_ORIGINAL, SUBJECT_ORIGINAL_NESTED, FEATURE, RUN


work$PROTEIN<-factor(work$PROTEIN)
work$PEPTIDE<-factor(work$PEPTIDE)
work$TRANSITION<-factor(work$TRANSITION)

work<-work[with(work,order(LABEL,GROUP_ORIGINAL,SUBJECT_ORIGINAL,RUN,PROTEIN,PEPTIDE,TRANSITION)),]

work$GROUP<-factor(work$GROUP)
work$SUBJECT<-factor(work$SUBJECT)
# SUBJECT_ORIGINAL_NESTED will sorted as GROUP_ORIGINAL, SUBJECT_ORIGINAL

work$SUBJECT_NESTED<-factor(work$SUBJECT_NESTED,levels=unique(work$SUBJECT_NESTED))

# FEATURE will sorted as PROTEIN, PEPTIDE, TRANSITION
work$FEATURE<-factor(work$FEATURE,levels=unique(work$FEATURE))
# RUN will sorted as GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN
work$RUN<-factor(work$RUN,levels=unique(work$RUN),labels=seq(1,length(unique(work$RUN))))


processout<-rbind(processout,c("Factorize in columns(GROUP, SUBJECT, GROUP_ORIGINAL, SUBJECT_ORIGINAL, SUBJECT_ORIGINAL_NESTED, FEATURE, RUN) - okay"))
write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)


############################################################
## Normalization : option 0. none
if(normalization=="None" | normalization==FALSE){
	processout<-rbind(processout,c("Normalization : no normalization - okay"))
	write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)
}

############################################################
## Normalization : option 1. constant normalization
if(normalization=="constant"){
	if(nlevels(work$LABEL)==1){
		### Constant normalization by heavy standard
		median.run<-tapply(work$ABUNDANCE,work[,"RUN"], function(x) median(x,na.rm=TRUE))
		median.all<-median(work$ABUNDANCE, na.rm=TRUE)
		
		for (i in 1:length(unique(work$RUN))){
			## ABUNDANCE is normalized
			work[work$RUN==i,"ABUNDANCE"]<-work[work$RUN==i,"ABUNDANCE"]-median.run[i]+median.all
		}
	}
		
	if(nlevels(work$LABEL)==2){
		### Constant normalization by heavy standard
		h<-work[work$LABEL=="H",]
		median.run<-tapply(h$ABUNDANCE,h[,"RUN"], function(x) median(x,na.rm=TRUE))
		median.all<-median(h$ABUNDANCE, na.rm=TRUE)
		
		for (i in 1:length(unique(work$RUN))){
			## ABUNDANCE is normalized
			work[work$RUN==i,"ABUNDANCE"]<-work[work$RUN==i,"ABUNDANCE"]-median.run[i]+median.all
		}		
	}
	
	processout<-rbind(processout,c("Normalization : Constant normalization - okay"))
	write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)
}
	

############################################################
## Normalization : option 2. quantile normalization
if(normalization=="quantile"){
	
	if(nlevels(work$LABEL)==1){
		
		## for label-free, just use endogenous
		worktemp<-work[!is.na(work$INTENSITY),]
		quantiletemp<-as.matrix(xtabs(ABUNDANCE~FEATURE+RUN, data=worktemp))
		
		## need to put NA for missing value in endogenous
		quantiletemp[quantiletemp==0]<-NA
		
		## using preprocessCore library
		quantiledone<-normalize.quantiles(quantiletemp)
		rownames(quantiledone)<-rownames(quantiletemp)
		colnames(quantiledone)<-colnames(quantiletemp)
		
		### get quantiled to long format for apply difference endogenous
		quantilelong<-melt(quantiledone, id=rownames(quantiledone))
		colnames(quantilelong)<-c("FEATURE","RUN","ABUNDANCE_quantile")

		## merge with original data
		work<-merge(work,quantilelong, by=c("FEATURE","RUN"))
		
		## reorder
		work<-data.frame("PROTEIN"=work$PROTEIN, "PEPTIDE"=work$PEPTIDE, "TRANSITION"=work$TRANSITION, "FEATURE"=work$FEATURE, "LABEL"=work$LABEL, "GROUP_ORIGINAL"=work$GROUP_ORIGINAL, "SUBJECT_ORIGINAL"=work$SUBJECT_ORIGINAL, "RUN"=work$RUN, "GROUP"=work$GROUP, "SUBJECT"=work$SUBJECT, "SUBJECT_NESTED"=work$SUBJECT_NESTED, "INTENSITY"=work$INTENSITY, "ABUNDANCE"=work$ABUNDANCE_quantile)

		work<-work[with(work,order(LABEL,GROUP_ORIGINAL,SUBJECT_ORIGINAL,RUN,PROTEIN,PEPTIDE,TRANSITION)),]

	}
		
	if(nlevels(work$LABEL)==2){
		
		##### for label-based, make quantile normalization for reference
		worktemp<-work[work$LABEL=="H" & !is.na(work$INTENSITY),]
		quantiletemp<-as.matrix(xtabs(ABUNDANCE~FEATURE+RUN, data=worktemp))
		
		## need to put NA for missing value in endogenous
		quantiletemp[quantiletemp==0]<-NA
		
		## using preprocessCore library
		quantiledone<-normalize.quantiles(quantiletemp)
		rownames(quantiledone)<-rownames(quantiletemp)
		colnames(quantiledone)<-colnames(quantiletemp)
		
		## get quantiled to long format for apply difference endogenous
		quantilelong.h<-melt(quantiledone, id=rownames(quantiledone))
		colnames(quantilelong.h)<-c("FEATURE","RUN","ABUNDANCE_quantile")
		quantilelong.h<-data.frame(quantilelong.h, LABEL="H")
		
		
		##### endogenous, in order to applying
		worktemp.l<-work[work$LABEL=="L" & !is.na(work$INTENSITY),]
		quantiletemp.l<-as.matrix(xtabs(ABUNDANCE~FEATURE+RUN, data=worktemp.l))
		
		## need to put NA for missing value in endogenous
		quantiletemp.l[quantiletemp.l==0]<-NA
		
		## apply the difference from reference
		quantiledone.l<-quantiletemp.l-(quantiletemp-quantiledone)
		
		## get quantiled to long format for apply difference endogenous
		quantilelong.l<-melt(quantiledone.l, id=rownames(quantiledone.l))
		colnames(quantilelong.l)<-c("FEATURE","RUN","ABUNDANCE_quantile")
		quantilelong.l<-data.frame(quantilelong.l, LABEL="L")
		
		quantilelong<-rbind(quantilelong.h,quantilelong.l)
		
		##### merge with original data
		work<-merge(work,quantilelong, by=c("FEATURE","RUN","LABEL"))
		
		## reorder
		work<-data.frame("PROTEIN"=work$PROTEIN, "PEPTIDE"=work$PEPTIDE, "TRANSITION"=work$TRANSITION, "FEATURE"=work$FEATURE, "LABEL"=work$LABEL, "GROUP_ORIGINAL"=work$GROUP_ORIGINAL, "SUBJECT_ORIGINAL"=work$SUBJECT_ORIGINAL, "RUN"=work$RUN, "GROUP"=work$GROUP, "SUBJECT"=work$SUBJECT, "SUBJECT_NESTED"=work$SUBJECT_NESTED, "INTENSITY"=work$INTENSITY, "ABUNDANCE"=work$ABUNDANCE_quantile)

		work<-work[with(work,order(LABEL,GROUP_ORIGINAL,SUBJECT_ORIGINAL,RUN,PROTEIN,PEPTIDE,TRANSITION)),]

	}
	
	processout<-rbind(processout,c("Normalization : Quantile normalization - okay"))
	write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)
}
	
	
	
############################################################
## Normalization : option 3. global standards - for endogenous
if(normalization=="globalStandards"){
	
	## need the names of global standards
	if(is.null(nameStandards)) stop ( message("For normalization with global standards, the names of global standards are needed. Please add 'nameStandards' input." ))
	
	combine<-NULL
	
	for(i in 1:length(nameStandards)){
		
		tempStandard<-work[work$PROTEIN==nameStandards[i],]

		### here, by RUN, but need to check !!!
		tempStandard<-tempStandard[tempStandard$GROUP!="0",]
		tempStandard$RUN<-factor(tempStandard$RUN)

		tempStandard<-tempStandard[!is.na(tempStandard$ABUNDANCE),]
		meanStandard<-tapply(tempStandard$ABUNDANCE, tempStandard$RUN, mean)

		combine<-cbind(combine, meanStandard)
	}
	
	## get mean among global standards
	if(ncol(combine)!=1){
		allmean<-apply(combine,1,mean)
	}else{
		allmean<-meanStandard
	}

	## median among all 
	median.est<-median(allmean[!is.na(allmean)])

	#### adjust
	for (i in 1:length(unique(work$RUN))){

		if(!is.na(allmean[i])) work[work$RUN==i & work$LABEL=="L","ABUNDANCE"]<-work[work$RUN==i & work$LABEL=="L","ABUNDANCE"]-allmean[i]+median.est

	}	

	processout<-rbind(processout,c("Normalization : normalization with global standards protein - okay"))
	write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)

}
	
	
	
##############################
## BetweenRunInterferenceScore
## need to make new function
##############################

if(betweenRunInterferenceScore){

	# only output light

	l<-subset(work,LABEL=="L")

	# add ProtFeature and ProtPeptide, because the shared peptides appear in multiple proteinss
	l$ProtFeature<-paste(l$PROTEIN,l$FEATURE,sep="/")	
	l$ProtPeptide<-paste(l$PROTEIN,l$PEPTIDE,sep="/")	


	temp<-tapply(l$ABUNDANCE,l[,c("RUN","ProtPeptide")],function(x) mean(x,na.rm=TRUE))
	temp1<-data.frame(ProtPeptide=rep(colnames(temp),each=dim(temp)[1]),RUN=rep(rownames(temp),dim(temp)[2]),meanPEPTIDE=as.numeric(unlist(temp)))

	temp2<-merge(l[,c("PROTEIN","PEPTIDE","FEATURE","ProtPeptide","ProtFeature","RUN","ABUNDANCE")],temp1,by=c("ProtPeptide","RUN"))

	temp3<-temp2[!is.na(temp2$ABUNDANCE),]

	temp4<-tapply(rownames(temp3),temp3[,c("ProtFeature")], function(x) cor(temp3[x,"ABUNDANCE"],temp3[x,"meanPEPTIDE"]))

	names<-unique(temp2[,c("PROTEIN","PEPTIDE","FEATURE","ProtFeature")])
	names<-names[with(names,order(ProtFeature)),]
	BetweenRunInterferenceFile<-data.frame(names[,c("PROTEIN","PEPTIDE","FEATURE")],BetweenRunInterferenceScore=temp4)

	BetweenRunInterferenceFile<-BetweenRunInterferenceFile[with(BetweenRunInterferenceFile,order(PROTEIN,PEPTIDE,FEATURE)),]

	write.table(BetweenRunInterferenceFile,file=paste(address,"BetweenRunInterferenceFile.txt",sep=""))
	
	processout<-rbind(processout,c("Between Run Interference Score is calculated and saved in .csv file - okay"))
	write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)

}else{
	processout<-rbind(processout,c("Between Run Interference Score is not calculated."))
	write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)
}

####### check missingness 
####### transitions are completely missing in one condition : missingness #######
if (nlevels(work$LABEL)==1){
	all.work<-work	
	test<-tapply(is.na(work[,"ABUNDANCE"]),work[,c("GROUP_ORIGINAL","FEATURE")],function(x) sum(x,na.rm=TRUE))
	numObs<-tapply(work[,"ABUNDANCE"],work[,c("GROUP_ORIGINAL","FEATURE")],function(x) length(x))
	test1<-test==numObs
	test2<-apply(test1,2,function(x) sum(x,na.rm=TRUE))
	filterList<-names(test2)[test2>0]
	final.decision<-ifelse(test2>0,1,0)
}	

if (nlevels(work$LABEL)==2){
	## first, remove NA
	all.work<-work   # with all NA observations
	work.miss<-na.omit(work)

	## draw table
	light<-subset(work.miss,LABEL=="L")
	heavy<-subset(work.miss,LABEL=="H")


	## use FEATURE because the name of transition can be used in other peptide
	count.light<-xtabs(~FEATURE+GROUP_ORIGINAL, light)
	count.heavy<-xtabs(~FEATURE+GROUP_ORIGINAL, heavy)

	count.light<-count.light==0
	count.heavy<-count.heavy==0

	count.light<-as.data.frame(count.light)
	count.heavy<-as.data.frame(count.heavy)

	## summary of missingness
	decision<-count.light
	decision[]<-0

	for(i in 1:ncol(decision)){
		for(j in 1:nrow(decision)){
		
		# either light or heavy has no obs -> subject to filter
			if(count.light[j,i]==TRUE || count.heavy[j,i]==TRUE) { decision[j,i]<-1 }
		}
	}
	
	final.decision<-apply(decision,1,sum)

## assign "subject to filter" column
work<-data.frame(work, "SuggestToFilter"=0)

	for(i in 1:length(final.decision)){
		# # assign subject_to_filter=1 for entire transition
		if(final.decision[i]!=0) work[work$FEATURE==names(final.decision[i]), "SuggestToFilter"]<-1
	}	
}


#### output : summary #####

	## output for label
	processout<-rbind(processout,c(paste(length(unique(work$LABEL)), " level of Isotope type labeling in this experiment",sep="")))
	write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)

temp<-data.frame("Summary of Features :")
colnames(temp)<-" "
rownames(temp)<-" "
print(temp)

summary.f<-matrix(NA,nrow=3)
summary.f[1]<-nlevels(work$PROTEIN)

temp<-unique(work[,c("PROTEIN","PEPTIDE")])
temp1<-xtabs(~PROTEIN,data=temp)
temp2<-summary(as.numeric(temp1))
summary.f[2]<-paste(temp2["Min."],temp2["Max."],sep="-")

temp<-unique(work[,c("PEPTIDE","FEATURE")])
temp1<-xtabs(~PEPTIDE,data=temp)
temp2<-summary(as.numeric(temp1))
summary.f[3]<-paste(temp2["Min."],temp2["Max."],sep="-")

colnames(summary.f)<-"count"
rownames(summary.f)<-c("# of Protein","# of Peptides/Protein", "# of Transitions/Peptide")

print(as.data.frame(summary.f))

	## output for process
	processout<-rbind(processout,c("Summary of Features :"))
	processout<-rbind(processout,c(paste(rownames(summary.f)[1]," : ",summary.f[1],sep="")))
	processout<-rbind(processout,c(paste(rownames(summary.f)[2]," : ",summary.f[2],sep="")))
	processout<-rbind(processout,c(paste(rownames(summary.f)[3]," : ",summary.f[3],sep="")))

	write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)

#### protein list with 1 feature
temp<-unique(work[,c("PROTEIN","FEATURE")])
temp1<-xtabs(~PROTEIN,data=temp)
temp2<-as.data.frame(temp1[temp1==1])
if(nrow(temp2)>0) message("\n","** Protein (",rownames(temp2),") has only single transition : Consider excluding this protein from the dataset.", "\n")

temp<-data.frame("Summary of Samples :")
colnames(temp)<-" "
rownames(temp)<-" "
print(temp)
	
summary.s<-matrix(NA,ncol=nlevels(work$GROUP_ORIGINAL),nrow=3)

## # of MS runs
temp<-unique(work[,c("GROUP_ORIGINAL","RUN")])
temp1<-xtabs(~GROUP_ORIGINAL,data=temp)
summary.s[1,]<-temp1

## # of biological replicates
temp<-unique(work[,c("GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
temp1<-xtabs(~GROUP_ORIGINAL,data=temp)
summary.s[2,]<-temp1

## # of tehcnical replicates
c.tech<-round(summary.s[1,]/summary.s[2,])
summary.s[3,]<-ifelse(c.tech==1,0,c.tech)
	
colnames(summary.s)<-unique(work$GROUP_ORIGINAL)
rownames(summary.s)<-c("# of MS runs","# of Biological Replicates", "# of Technical Replicates")

print(summary.s)

	
	
message("\n Summary of Missingness :\n" )
message("  # transitions are completely missing in one condition: ", sum(final.decision!=0), "\n")
if(sum(final.decision!=0)!=0) message("    -> ", paste(names(final.decision[final.decision!=0]),collapse = ", "))

without<-xtabs(~RUN,work)
withall<-xtabs(~RUN,all.work)
run.missing<-without/withall
message("  # run with 75% missing observations: ", sum(run.missing<0.25), "\n")
if(sum(run.missing<0.25)!=0) message("    -> ", paste("RUN",names(without[run.missing<0.25]),sep=" "))
	
	## output process
	processout<-rbind(processout,c("Summary of Missingness :"))
	processout<-rbind(processout,c(paste("  # transitions are completely missing in one condition: ", sum(final.decision!=0), sep="")))
	if(sum(final.decision!=0)!=0) processout<-rbind(processout,"    -> ", paste(names(final.decision[final.decision!=0]),collapse = ", "))

	processout<-rbind(processout,c(paste("  # run with 75% missing observations: ", sum(run.missing<0.25), sep="")))
	if(sum(run.missing<0.25)!=0) processout<-rbind(processout,"    -> ", paste("RUN",names(without[run.missing<0.25]),sep=" "))

	write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)


#### check any protein has only light for labeled-experiment
if (nlevels(work$LABEL)==2){
	temp<-unique(work[,c("PROTEIN","LABEL")])
	temp1<-xtabs(~PROTEIN,data=temp)
	
	if(any(temp1!=2)){
		# check that is L or H
		namepro<-names(temp1[temp1!=2])
		for(j in 1:length(namepro)){
			if(unique(work[work$PROTEIN==namepro[j],"LABEL"])=="L")
				message("\n *** ",namepro[j]," has only endogeneous intensities in label-based experiment. Please check this protein or remove it.")
			
			if(unique(work[work$PROTEIN==namepro[j],"LABEL"])=="H")
				message("\n *** ",namepro[j]," has only reference intensities in label-based experiment. Please check this protein or remove it.")
		}
	}	
}

		###
		processout<-rbind(processout,c("Processing data for analysis is done. - okay"))
		write.table(processout, file="Process_dataProcess.txt",row.names=FALSE)

## return work data.frame	
return(work)

}




#############################################
#############################################
# Part 2 dataProcessPlots
#############################################
#############################################


dataProcessPlots<-function(data=data,type=type,featureName="Transition",ylimUp=FALSE,ylimDown=FALSE,scale=FALSE,interval="SE",axis.size=10,text.size=4,text.angle=0,legend.size=7,width=10,height=10, which.Protein="all", address=""){
	
	data$PROTEIN<-factor(data$PROTEIN)	

	#################
	## Profile plot
	#################
	if (type=="ProfilePlot"){
		
		#### choose Proteins or not
		if(which.Protein!="all"){
			## check which.Protein is name of Protein
			if(is.character(which.Protein)){
				temp.name<-which.Protein
			
				## message if name of Protein is wrong.
				if(length(setdiff(temp.name,unique(data$PROTEIN)))>0)
					stop(message(paste("error : Please check protein name. Data set does not have this protein. -", paste(temp.name, collapse=", "),sep=" ")))
			}
		
			## check which.Protein is order number of Protein
			if(is.numeric(which.Protein)){
				temp.name<-levels(data$PROTEIN)[which.Protein]
				
				## message if name of Protein is wrong.
				if(length(levels(data$PROTEIN))<max(which.Protein))
					stop(message(paste("error : Please check your selection of proteins. There are ", length(levels(data$PROTEIN))," proteins in this dataset.",sep=" ")))
			}
			
			## use only assigned proteins
			data<-data[which(data$PROTEIN %in% temp.name),]
			data$PROTEIN<-factor(data$PROTEIN)
		}
		
		#### save the plots as pdf or not
		if(address!=FALSE) pdf(paste(address,"ProfilePlot.pdf",sep=""), width=width, height=height)
	

		# assign upper or lower limit
		## ylimUp
		y.limup<-30
		if(is.numeric(ylimUp)) y.limup<-ylimUp 

		## ylimDown
		y.limdown=-1
		if(is.numeric(ylimDown)) y.limdown<-ylimDown 

		data<-data[with(data,order(GROUP_ORIGINAL,SUBJECT_ORIGINAL,LABEL)),]
		data$RUN<-factor(data$RUN,levels=unique(data$RUN),labels=seq(1,length(unique(data$RUN))))
		data$RUN<-as.numeric(data$RUN)
		tempGroupName<-unique(data[,c("GROUP_ORIGINAL","RUN")])

		groupAxis<-as.numeric(xtabs(~GROUP_ORIGINAL,tempGroupName))
		cumGroupAxis<-cumsum(groupAxis)
		lineNameAxis<-cumGroupAxis[-nlevels(data$GROUP_ORIGINAL)]

		groupName<-data.frame(RUN=c(0,lineNameAxis)+groupAxis/2+0.5,y=rep(y.limup-1,length(groupAxis)),Name=levels(data$GROUP_ORIGINAL))

		for (i in 1:nlevels(data$PROTEIN)){	
			sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
			sub$FEATURE<-factor(as.character(sub$FEATURE))	
			sub$SUBJECT<-factor(sub$SUBJECT)	
			sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
			sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)
			sub$PEPTIDE<-factor(as.character(sub$PEPTIDE))

			if(length(unique(data$LABEL))==2){
				sub$LABEL<-factor(sub$LABEL,labels=c("Reference","Endogenous"))	
			}else{
				if(unique(data$LABEL)=="L"){
					sub$LABEL<-factor(sub$LABEL,labels=c("Endogenous"))	
				}
				if(unique(data$LABEL)=="H"){
					sub$LABEL<-factor(sub$LABEL,labels=c("Reference"))
				}
			}

			#sub<- sub[with(sub, order(LABEL,RUN,FEATURE)), ]

			# seq for peptide and transition
			b<-unique(sub[,c("PEPTIDE","FEATURE")])
			b<-b[with(b,order(PEPTIDE,FEATURE)),] ## add because if there are missing value, orders are different.
			
			temp1<-xtabs(~b[,1])
			ss<-NULL
			s<-NULL
			
			for(j in 1:length(temp1)){
				temp3<-rep(j,temp1[j])
				s<-c(s,temp3)
				temp2<-seq(1,temp1[j])
				ss<-c(ss,temp2)	
			}

			#options(show.error.messages = FALSE)

			if(featureName=="Transition"){
	
            	ptemp<-ggplot(aes_string(x='RUN', y='ABUNDANCE', color='FEATURE',linetype='FEATURE'), data=sub)+facet_grid(~LABEL)+geom_point(size=1.3)+geom_line(size=0.5)+scale_colour_manual(values=s)+scale_linetype_manual(values=ss,guide="none")+scale_x_continuous('MS runs',breaks=cumGroupAxis)+annotate("text",x=groupName$RUN,y=groupName$y,label=groupName$Name,size=text.size,angle=text.angle)+geom_vline(xintercept=lineNameAxis+0.5,colour="grey",linetype="longdash")+labs(title=unique(sub$PROTEIN))+theme(
				panel.background=element_rect(fill='white', colour="black"),
				legend.key=element_rect(fill='white',colour='white'),
				panel.grid.minor = element_blank(),
				strip.background=element_rect(fill='gray95'),
				strip.text.x=element_text(colour=c("#00B0F6"),size=14),
				axis.text.x=element_text(size=axis.size,colour="black"),
				axis.text.y=element_text(size=axis.size,colour="black"),
				axis.ticks=element_line(colour="black"),
				axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
				axis.title.y=element_text(size=axis.size+5,vjust=0.3),
				title=element_text(size=axis.size+8,vjust=1.5),
				legend.position="top",
				legend.text=element_text(size=legend.size))+guides(color=guide_legend(title=NULL,ncol=3,override.aes=list(linetype=ss,size=0.5)))

				## y-axis labeling
				temp<-sub[!is.na(sub[,"ABUNDANCE"]),]
				temptest<-abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])

				if(temptest){
					ptemp<-ptemp+scale_y_continuous('Log2-intensities',limit=c(y.limdown, y.limup))
				}else{
					ptemp<-ptemp+scale_y_continuous('Log10-intensities',limit=c(y.limdown, y.limup))
				}	

				print(ptemp)
				
				message(paste("Drew the Profile plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))

			}

			if(featureName=="Peptide"){

				ptemp<-ggplot(aes_string(x='RUN', y='ABUNDANCE', color='PEPTIDE',linetype='FEATURE'), data=sub)+facet_grid(~LABEL)+geom_point(size=1.3)+geom_line(size=0.5)+scale_colour_manual(values=unique(s))+scale_linetype_manual(values=ss,guide="none")+scale_x_continuous('MS runs',breaks=cumGroupAxis)+annotate("text",x=groupName$RUN,y=groupName$y,label=groupName$Name,size=text.size,angle=text.angle)+geom_vline(xintercept=lineNameAxis+0.5,colour="grey",linetype="longdash")+labs(title=unique(sub$PROTEIN))+theme(
				panel.background=element_rect(fill='white', colour="black"),
				legend.key=element_rect(fill='white',colour='white'),
				panel.grid.minor = element_blank(),
				strip.background=element_rect(fill='gray95'),	
				strip.text.x=element_text(colour=c("#00B0F6"),size=14),
				axis.text.x=element_text(size=axis.size,colour="black"),
				axis.text.y=element_text(size=axis.size,colour="black"),
				axis.ticks=element_line(colour="black"),
				axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
				axis.title.y=element_text(size=axis.size+5,vjust=0.3),
				title=element_text(size=axis.size+8,vjust=1.5),
				legend.position="top",
				legend.text=element_text(size=legend.size))+guides(color=guide_legend(title=paste("# peptide:",nlevels(sub$PEPTIDE)),ncol=3))

				## y-axis labeling
				temp<-sub[!is.na(sub[,"ABUNDANCE"]),]
				temptest<-abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])

				if(temptest){
					ptemp<-ptemp+scale_y_continuous('Log2-intensities',limit=c(y.limdown, y.limup))	
				}else{
					ptemp<-ptemp+scale_y_continuous('Log10-intensities',limit=c(y.limdown, y.limup))
				}	

				print(ptemp)
				
				message(paste("Drew the Profile plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))

			}

			if(featureName=="NA"){
	
				ptemp<-ggplot(aes_string(x='RUN', y='ABUNDANCE', color='PEPTIDE',linetype='FEATURE'), data=sub)+facet_grid(~LABEL)+geom_point(size=1.3)+geom_line(size=0.5)+scale_colour_manual(values=unique(s))+scale_linetype_manual(values=ss,guide="none")+scale_x_continuous('MS runs',breaks=cumGroupAxis)+annotate("text",x=groupName$RUN,y=groupName$y,label=groupName$Name,size=text.size,angle=text.angle)+geom_vline(xintercept=lineNameAxis+0.5,colour="grey",linetype="longdash")+labs(title=unique(sub$PROTEIN))+theme(
				panel.background=element_rect(fill='white', colour="black"),
				legend.key=element_rect(fill='white',colour='white'),
				panel.grid.minor = element_blank(),
				strip.background=element_rect(fill='gray95'),	
				strip.text.x=element_text(colour=c("#00B0F6"),size=14),
				axis.text.x=element_text(size=axis.size,colour="black"),
				axis.text.y=element_text(size=axis.size,colour="black"),
				axis.ticks=element_line(colour="black"),
				axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
				axis.title.y=element_text(size=axis.size+5,vjust=0.3),
				title=element_text(size=axis.size+8,vjust=1.5),
				legend.position="none")

				## y-axis labeling
				temp<-sub[!is.na(sub[,"ABUNDANCE"]),]
				temptest<-abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])

				if(temptest){
					ptemp<-ptemp+scale_y_continuous('Log2-intensities',limit=c(y.limdown, y.limup))
				}else{
					ptemp<-ptemp+scale_y_continuous('Log10-intensities',limit=c(y.limdown, y.limup))
				}	

				print(ptemp)
				
				message(paste("Drew the Profile plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))

			}
		} ## end-loop for each protein
		if(address!=FALSE) dev.off() 
	} ## end Profile plot	


#################### 
## QC plot (Quality control plot)
#################### 
	if (type=="QCPlot"){
		
		#### save the plots as pdf or not
		if(address!=FALSE) pdf(paste(address,"QCPlot.pdf",sep=""), width=width, height=height)

		## options(warn=-1)

		# assign upper or lower limit
 		## ylimUp
		y.limup<-30
		if(is.numeric(ylimUp)) y.limup<-ylimUp 

		## ylimDown
		y.limdown=-1
		if(is.numeric(ylimDown)) y.limdown<-ylimDown 

		# relabel the Run (make it sorted by group first)
		data<-data[with(data,order(GROUP_ORIGINAL,SUBJECT_ORIGINAL)),]
		data$RUN<-factor(data$RUN,levels=unique(data$RUN),labels=seq(1,length(unique(data$RUN))))

		if(length(unique(data$LABEL))==2){
			data$LABEL<-factor(data$LABEL,labels=c("Reference","Endogenous"))	
			label.color<-c("darkseagreen1","lightblue")
		}else{
			if(unique(data$LABEL)=="L"){
				data$LABEL<-factor(data$LABEL,labels=c("Endogenous"))
				label.color<-c("lightblue")	
			}
			if(unique(data$LABEL)=="H"){
				data$LABEL<-factor(data$LABEL,labels=c("Reference"))
				label.color<-c("darkseagreen1")
			}
		}

		tempGroupName<-unique(data[,c("GROUP_ORIGINAL","RUN")])
		data<-data[with(data,order(LABEL,GROUP_ORIGINAL,SUBJECT_ORIGINAL)),]

		groupAxis<-as.numeric(xtabs(~GROUP_ORIGINAL,tempGroupName))
		cumGroupAxis<-cumsum(groupAxis)
		lineNameAxis<-cumGroupAxis[-nlevels(data$GROUP_ORIGINAL)]

		groupName<-data.frame(RUN=c(0,lineNameAxis)+groupAxis/2+0.5,y=rep(y.limup-1,length(groupAxis)),Name=levels(data$GROUP_ORIGINAL))

		#### all protein
		ptemp<-ggplot(aes_string(x='RUN', y='ABUNDANCE'), data=data)+facet_grid(~LABEL)+geom_boxplot(aes_string(fill='LABEL'),outlier.shape=1,outlier.size=1.5)+scale_fill_manual(values=label.color, guide="none")+scale_x_discrete('MS runs',breaks=cumGroupAxis)+annotate("text",x=groupName$RUN,y=groupName$y,label=groupName$Name,size=text.size,angle=text.angle)+geom_vline(xintercept=lineNameAxis+0.5,colour="grey",linetype="longdash")+labs(title="All proteins")+theme(
		panel.background=element_rect(fill='white', colour="black"),
		legend.key=element_rect(fill='white',colour='white'),
		panel.grid.minor = element_blank(),
		strip.background=element_rect(fill='gray95'),	
		strip.text.x=element_text(colour=c("#00B0F6"),size=14),
		axis.text.x=element_text(size=axis.size,colour="black"),
		axis.text.y=element_text(size=axis.size,colour="black"),
		axis.ticks=element_line(colour="black"),
		axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
		axis.title.y=element_text(size=axis.size+5,vjust=0.3),
		title=element_text(size=axis.size+8,vjust=1.5))

		## y-axis labeling
		temp<-data[!is.na(data[,"ABUNDANCE"]),]
		temptest<-abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])

		if(temptest){
			ptemp<-ptemp+scale_y_continuous('Log2-intensities',limit=c(y.limdown, y.limup))
		}else{
			ptemp<-ptemp+scale_y_continuous('Log10-intensities',limit=c(y.limdown, y.limup))
		}	

		print(ptemp)
		
		message("Drew the Quality Contol plot(boxplot) for all proteins.")

		#### each protein
		
		#### choose Proteins or not
		if(which.Protein!="all"){
			## check which.Protein is name of Protein
			if(is.character(which.Protein)){
				
				temp.name<-which.Protein
			
				## message if name of Protein is wrong.
				if(length(setdiff(temp.name,unique(data$PROTEIN)))>0){
					message(paste("error : Please check protein name. Data set does not have this protein. -", paste(temp.name, collapse=", "),sep=" "))
					dev.off()
					stop()
				}
			}
				
			## check which.Protein is order number of Protein
			if(is.numeric(which.Protein)){
				temp.name<-levels(data$PROTEIN)[which.Protein]
			
				## message if name of Protein is wrong.
				if(length(levels(data$PROTEIN))<max(which.Protein)){
					message(paste("error : Please check your selection of proteins. There are ", length(levels(data$PROTEIN))," proteins in this dataset.",sep=" "))
					dev.off()
					stop()
				}
			}
			
			## use only assigned proteins
			data<-data[which(data$PROTEIN %in% temp.name),]
			data$PROTEIN<-factor(data$PROTEIN)
		}
				
		for (i in 1:nlevels(data$PROTEIN)){	
			sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
			subTemp<-sub[!is.na(sub$ABUNDANCE),]
			subTemp$LABEL<-factor(subTemp$LABEL)	
			sub<-sub[with(sub, order(LABEL,RUN)),]

			if(length(unique(subTemp$LABEL))==2){
				label.color<-c("darkseagreen1","lightblue")
			}else{
				if(unique(subTemp$LABEL)=="Endogenous"){
					label.color<-c("lightblue")	
				}
				if(unique(subTemp$LABEL)=="Reference"){
					label.color<-c("darkseagreen1")
				}
			}

			##options(show.error.messages = FALSE)

			ptemp<-ggplot(aes_string(x='RUN', y='ABUNDANCE'), data=sub)+facet_grid(~LABEL)+geom_boxplot(aes_string(fill='LABEL'),outlier.shape=1,outlier.size=1.5)+scale_fill_manual(values=label.color, guide="none")+scale_x_discrete('MS runs',breaks=cumGroupAxis)+annotate("text",x=groupName$RUN,y=groupName$y,label=groupName$Name,size=text.size,angle=text.angle)+geom_vline(xintercept=lineNameAxis+0.5,colour="grey",linetype="longdash")+labs(title=unique(sub$PROTEIN))+theme(
			panel.background=element_rect(fill='white', colour="black"),
			legend.key=element_rect(fill='white',colour='white'),
			panel.grid.minor = element_blank(),
			strip.background=element_rect(fill='gray95'),	
			strip.text.x=element_text(colour=c("#00B0F6"),size=14),
			axis.text.x=element_text(size=axis.size,colour="black"),
			axis.text.y=element_text(size=axis.size,colour="black"),
			axis.ticks=element_line(colour="black"),
			axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
			axis.title.y=element_text(size=axis.size+5,vjust=0.3),
			title=element_text(size=axis.size+8,vjust=1.5))

			## y-axis labeling
			if(temptest){
				ptemp<-ptemp+scale_y_continuous('Log2-intensities',limit=c(y.limdown, y.limup))
			}else{
				ptemp<-ptemp+scale_y_continuous('Log10-intensities',limit=c(y.limdown, y.limup))
			}	

			print(ptemp)
			
			message(paste("Drew the Quality Contol plot(boxplot) for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))

		} ## end-loop
		if(address!=FALSE) dev.off()
	} # end QC plot	


#################### 
### Condition plot
#################### 
	if(type=="ConditionPlot"){
		
		#### choose Proteins or not
		if(which.Protein!="all"){
			## check which.Protein is name of Protein
			if(is.character(which.Protein)){
				
				temp.name<-which.Protein
			
				## message if name of Protein is wrong.
				if(length(setdiff(temp.name,unique(data$PROTEIN)))>0)
					stop(message(paste("error : Please check protein name. Dataset does not have this protein. -", paste(temp.name, collapse=", "),sep=" ")))
			}
		
			## check which.Protein is order number of Protein
			if(is.numeric(which.Protein)){
				
				temp.name<-levels(data$PROTEIN)[which.Protein]
				
				## message if name of Protein is wrong.
				if(length(levels(data$PROTEIN))<max(which.Protein))
					stop(message(paste("error : Please check your selection of proteins. There are ", length(levels(data$PROTEIN))," proteins in this dataset.",sep=" ")))
			}
			
			## use only assigned proteins
			data<-data[which(data$PROTEIN %in% temp.name),]
			data$PROTEIN<-factor(data$PROTEIN)
		}
		
		#### save the plots as pdf or not
		if(address!=FALSE) pdf(paste(address,"ConditionPlot.pdf",sep=""), width=width, height=height)

		if(nlevels(data$LABEL)==1){
			for (i in 1:nlevels(data$PROTEIN)){	
				sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
				sub<-na.omit(sub)	
				sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
				sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)	
				sub$FEATURE<-factor(sub$FEATURE)	
	
				######## statistics
				sub.mean<-by(sub$ABUNDANCE, sub$GROUP_ORIGINAL, function(x) mean(x,na.rm=TRUE))
				sub.sd<-by(sub$ABUNDANCE, sub$GROUP_ORIGINAL, sd)
				sub.len<-by(sub$ABUNDANCE, sub$GROUP_ORIGINAL, length)
				if(interval=="SE") ciw<-qt(0.975,sub.len)*sub.sd/sqrt(sub.len)
				if(interval=="SD") ciw<-sub.sd
				
				if(sum(is.na(ciw))>=1) ciw[is.na(ciw)]<-0

				# assign upper or lower limit
 				## ylimUp
				y.limup<-ceiling(max(sub.mean+ciw))
				if(is.numeric(ylimUp)) y.limup<-ylimUp 

				## ylimDown
				y.limdown<-floor(min(sub.mean-ciw))
				if(is.numeric(ylimDown)) y.limdown<-ylimDown 

				if(!scale){  ## scale: false
		
					## reformat as data.frame
					tempsummary<-data.frame(Label=unique(sub$GROUP_ORIGINAL),mean=as.vector(sub.mean),ciw=as.vector(ciw))

 					ptemp<-ggplot(aes_string(x='Label', y='mean'), data=tempsummary)+geom_errorbar(aes(ymax = mean + ciw, ymin=mean - ciw),data=tempsummary, width=0.1,colour="red")+geom_point(size=3,colour="darkred")+scale_x_discrete('Condition')+geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+labs(title=unique(sub$PROTEIN))+theme(
 					panel.background=element_rect(fill='white', colour="black"),
 					panel.grid.major.y = element_line(colour="grey95"),
 					panel.grid.minor.y = element_blank(),
 					axis.text.x=element_text(size=axis.size,colour="black",angle=text.angle),
 					axis.text.y=element_text(size=axis.size,colour="black"),
 					axis.ticks=element_line(colour="black"),
 					axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
 					axis.title.y=element_text(size=axis.size+5,vjust=0.3),
 					title=element_text(size=axis.size+8,vjust=1.5))

				}else{
					#extract numeric value, don't use levels (because T1,T10,T3,...)
					# ?? sort??
	
					## reformat as data.frame
					tempsummary<-data.frame(Label=as.numeric(gsub("\\D","",unique(sub$GROUP_ORIGINAL))),mean=as.vector(sub.mean),ciw=as.vector(ciw))

 					ptemp<-ggplot(aes_string(x='Label', y='mean'), data=tempsummary)+geom_errorbar(aes(ymax = mean + ciw, ymin=mean - ciw),data=tempsummary, width=0.1,colour="red")+geom_point(size=3,colour="darkred")+scale_x_continuous('Condition')+geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+labs(title=unique(sub$PROTEIN))+theme(
 					panel.background=element_rect(fill='white', colour="black"),
 					panel.grid.major.y = element_line(colour="grey95"),
 					panel.grid.minor.y = element_blank(),
 					axis.text.x=element_text(size=axis.size,colour="black",angle=text.angle),
 					axis.text.y=element_text(size=axis.size,colour="black"),
 					axis.ticks=element_line(colour="black"),
 					axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
 					axis.title.y=element_text(size=axis.size+5,vjust=0.3),
 					title=element_text(size=axis.size+8,vjust=1.5))

				} ## scale : true
				
				## y-axis labeling, find log 2 or log 10
				temp<-sub[!is.na(sub[,"ABUNDANCE"]),]
				temptest<-abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])

				if(temptest){
					ptemp<-ptemp+scale_y_continuous("Log2-intensities",limit=c(y.limdown, y.limup))
				}else{
					ptemp<-ptemp+scale_y_continuous("Log10-intensities",limit=c(y.limdown, y.limup))
				}
					
				print(ptemp)
				
				message(paste("Drew the condition plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))

			} ## end-loop
		} ## label-free

		if(nlevels(data$LABEL)==2){
			for (i in 1:nlevels(data$PROTEIN)){	

				sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
				sub<-na.omit(sub)	
				sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
				sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)	
				sub$FEATURE<-factor(sub$FEATURE)	
	
				sub.h<-subset(sub, LABEL=="H")
				sub.l<-subset(sub, LABEL=="L")
	
				sub.l<-data.frame(sub.l,"HEAVY"=0)
		
				## matching heavy and light

				for(j in 1:nrow(sub.l)){
					if(length(sub.h[sub.h$FEATURE==sub.l[j,"FEATURE"] & sub.h$GROUP_ORIGINAL==sub.l[j,"GROUP_ORIGINAL"] & sub.h$SUBJECT_ORIGINAL==sub.l[j,"SUBJECT_ORIGINAL"] & sub.h$RUN==sub.l[j,"RUN"],"ABUNDANCE"])!=0)
						sub.l[j,"HEAVY"]<-sub.h[sub.h$FEATURE==sub.l[j,"FEATURE"] & sub.h$GROUP_ORIGINAL==sub.l[j,"GROUP_ORIGINAL"] & sub.h$SUBJECT_ORIGINAL==sub.l[j,"SUBJECT_ORIGINAL"] & sub.h$RUN==sub.l[j,"RUN"],"ABUNDANCE"] 
				}

  				sub.l[sub.l$HEAVY==0,"HEAVY"]<-mean(sub.h[sub.h$GROUP_ORIGINAL==sub.l[i,"GROUP_ORIGINAL"] & sub.h$SUBJECT_ORIGINAL==sub.l[i,"SUBJECT_ORIGINAL"] & sub.h$RUN==sub.l[i,"RUN"],"ABUNDANCE"])
  	
				sub.l$ratio<-sub.l$ABUNDANCE - sub.l$HEAVY  ## log(L/H)

				sub.mean<-by(sub.l$ratio, sub.l$GROUP_ORIGINAL, function(x) mean(x,na.rm=TRUE))
				sub.sd<-by(sub.l$ratio, sub.l$GROUP_ORIGINAL, sd)
				sub.len<-by(sub.l$ratio, sub.l$GROUP_ORIGINAL, length)
				if(interval=="SE") ciw<-qt(0.975,sub.len)*sub.sd/sqrt(sub.len)
				if(interval=="SD") ciw<-sub.sd
				
				if(sum(is.na(ciw))>=1) ciw[is.na(ciw)]<-0
				
				# assign upper or lower limit
 				## ylimUp
				y.limup<-ceiling(max(sub.mean+ciw))
				if(is.numeric(ylimUp)) y.limup<-ylimUp 

				## ylimDown
				y.limdown<-floor(min(sub.mean-ciw))
				if(is.numeric(ylimDown)) y.limdown<-ylimDown 
  
				if(!scale){
	
	  				## reformat as data.frame
					tempsummary<-data.frame(Label=unique(sub.l$GROUP_ORIGINAL),mean=as.vector(sub.mean),ciw=as.vector(ciw))
	
  					ptemp<-ggplot(aes_string(x='Label', y='mean'), data=tempsummary)+geom_errorbar(aes(ymax = mean + ciw, ymin=mean - ciw),data=tempsummary, width=0.1,colour="red")+geom_point(size=3,colour="darkred")+scale_x_discrete('Condition')+geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+labs(title=unique(sub$PROTEIN))+theme(
  					panel.background=element_rect(fill='white', colour="black"),
  					panel.grid.major.y = element_line(colour="grey95"),
  					panel.grid.minor.y = element_blank(),
  					axis.text.x=element_text(size=axis.size,colour="black",angle=text.angle),
  					axis.text.y=element_text(size=axis.size,colour="black"),
  					axis.ticks=element_line(colour="black"),
  					axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
  					axis.title.y=element_text(size=axis.size+5,vjust=0.3),
  					title=element_text(size=axis.size+8,vjust=1.5))
  				}else{
	  				## reformat as data.frame
					tempsummary<-data.frame(Label=as.numeric(gsub("\\D","",unique(sub.l$GROUP_ORIGINAL))),mean=as.vector(sub.mean),ciw=as.vector(ciw))

   					ptemp<-ggplot(aes_string(x='Label', y='mean'), data=tempsummary)+geom_errorbar(aes(ymax = mean + ciw, ymin=mean - ciw),data=tempsummary, width=0.1,colour="red")+geom_point(size=3,colour="darkred")+scale_x_continuous('Condition')+geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+labs(title=unique(sub$PROTEIN))+theme(
   					panel.background=element_rect(fill='white', colour="black"),
   					panel.grid.major.y = element_line(colour="grey95"),
   					panel.grid.minor.y = element_blank(),
   					axis.text.x=element_text(size=axis.size,colour="black",angle=text.angle),
   					axis.text.y=element_text(size=axis.size,colour="black"),
   					axis.ticks=element_line(colour="black"),
   					axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
   					axis.title.y=element_text(size=axis.size+5,vjust=0.3),
   					title=element_text(size=axis.size+8,vjust=1.5))
  				} ## scale : true
  				
  				## y-axis labeling, find log 2 or log 10
				temp<-sub[!is.na(sub[,"ABUNDANCE"]),]
				temptest<-abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])

				if(temptest){
					ptemp<-ptemp+scale_y_continuous("Log2-Ratio(L/H)",limit=c(y.limdown, y.limup))
				}else{
					ptemp<-ptemp+scale_y_continuous("Log10-Ratio(L/H)",limit=c(y.limdown, y.limup))
				}
			
				print(ptemp)
				message(paste("Drew the condition plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
			} # end-loop
		} ## label-based
		
		
		if(address!=FALSE) dev.off()
	} ## Condition plot
}




#############################################
#############################################
#############################################
#############################################
# Part 1 detect differential expressed proteins
#############################################
#############################################
#############################################
#############################################

groupComparison<-function(contrast.matrix=contrast.matrix,data=data,labeled=TRUE, scopeOfBioReplication="restricted", scopeOfTechReplication="expanded", interference=TRUE,featureVar=FALSE,missing.action = "nointeraction"){

	## save process output in each step
	processout<-NULL
	
	
	## check input is correct
	## data format
	rawinput<-c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Run","Intensity")
	
	if(length(setdiff(toupper(rawinput),toupper(colnames(data))))==0){
		processout<-rbind(processout,c(paste("The required input - data : did not process from dataProcess function. - stop")))
		write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)

		stop(message("Please use 'dataProcess' first. Then use output of dataProcess function as input in groupComparison."))
	}
	
	
	## contrast. matrix
	if(ncol(contrast.matrix)!=length(unique(data$GROUP_ORIGINAL))){
		processout<-rbind(processout,c(paste("The required input - contrast.matrix: the number of column and the number of group are not the same. - stop")))
		colnames(processout)<-c("Process-groupComparison function")
		write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)
		
		stop(message("Please check contrast matrix. The number of group in data set is different with columns of contrast.matrix."))
	}
		
	# check whether row.names of contrast.matrix.sub exists or not
	if(sum(is.null(row.names(contrast.matrix)))>0){
		processout<-rbind(processout,c(paste("The required input - contrast.matrix: need row names of contrast.matrix . - stop")))
		write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)
		
		stop(message("No row.names of comparison exist.\n"))
	}
	



	## labeled value
	if(!(labeled==TRUE | labeled==FALSE) | !is.logical(labeled)){
		processout<-rbind(processout,c(paste("The required input - labeled : 'labeled' value is wrong. - stop")))
		write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)

		 stop(message("'labeled' must be one of TRUE or FALSE as a logical value."))
	}
	
	if(sum(unique(data$LABEL) %in% "H")==0 & labeled==TRUE){
		processout<-rbind(processout,c(paste("The required input - labeled : 'labeled' value is wrong. labeled should be 'FALSE' - stop")))
		write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)

		 stop(message("Dataset is based on label-free experiment. Please set \"labeled\"=FALSE."))
	}


	## other option value
	if(!(scopeOfTechReplication=="restricted" | scopeOfTechReplication=="expanded")){
		processout<-rbind(processout,c(paste("The required input - scopeOfTechReplication : 'scopeOfTechReplication' value is wrong. - stop")))
		write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)
		
		stop(message("'scopeOfTechReplication' must be one of \"restricted\" or \"expanded\".
	"))
	}
	
		
	if(!(scopeOfBioReplication=="restricted" | scopeOfBioReplication=="expanded")){
		processout<-rbind(processout,c(paste("The required input - scopeOfBioReplication : 'scopeOfBioReplication' value is wrong. - stop")))
		write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)
		
		stop(message("'scopeOfBioReplication' must be one of \"restricted\" or \"expanded\".
	"))
	} 
	
	if(!(interference==TRUE | interference==FALSE) | !is.logical(interference)){
		processout<-rbind(processout,c(paste("The required input - interference : 'interference' value is wrong. - stop")))
		write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)

		 stop(message("'interference' must be one of TRUE or FALSE as a logical value."))
	}
	
	if(!(featureVar==TRUE | featureVar==FALSE) | !is.logical(featureVar)){
		processout<-rbind(processout,c(paste("The required input - featureVar : 'featureVar' value is wrong. - stop")))
		write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)

		stop(message("'featureVar' must be one of TRUE or FALSE as a logical value."))
	}
	
	if(!(missing.action %in% c("nointeraction", "impute", "remove"))){
		processout<-rbind(processout,c(paste("The required input - missing.action : 'missing.action' value is wrong. - stop")))
		write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)

		 stop(message("'missing.action' must be one of \"nointeraction\", \"impute\", or \"remove\"."))
	}

	## all input
		processout<-rbind(processout,c(paste("labeled = ",labeled,sep="")))
		processout<-rbind(processout,c(paste("scopeOfBioReplication = ",scopeOfBioReplication,sep="")))
		processout<-rbind(processout,c(paste("scopeOfTechReplication = ", scopeOfTechReplication,sep="")))
		processout<-rbind(processout,c(paste("interference = ",interference,sep="")))
		processout<-rbind(processout,c(paste("featureVar = ",featureVar,sep="")))
		colnames(processout)<-c("Process-groupComparison function")

		write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)


## check whether case-control(FALSE) or time-course(TRUE)
repeated<-.checkRepeated(data)

	if(repeated){ 
		processout<-rbind(processout,c(paste("Time course design of experiment - okay")))
	}else{
		processout<-rbind(processout,c(paste("Case control design of experiment - okay")))

	}
	write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)


## since case-control(FALSE) with fixed subject and random run will fit the crossed model
# we need to set subject_original to non-unique

if(repeated==FALSE&scopeOfBioReplication=="restricted"&scopeOfTechReplication=="expanded"){
	test<-unique(data[,c("GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
	test<- test[with(test, order(GROUP_ORIGINAL,SUBJECT_ORIGINAL)),]
	test$GROUP_ORIGINAL<-factor(test$GROUP_ORIGINAL)
	test1<-as.matrix(xtabs(~test[,1]))
	test$SUBJECT_ORIGINAL_NOUNIQUE<-as.numeric(unlist(apply(test1,1,function(x) seq(x))))

	data$SUBJECT_ORIGINAL<-as.character(data$SUBJECT_ORIGINAL)

	for (i in 1:length(unique(test$SUBJECT_ORIGINAL_NOUNIQUE))){
		list<-test$SUBJECT_ORIGINAL[test$SUBJECT_ORIGINAL_NOUNIQUE==i]
		data$SUBJECT_ORIGINAL[data$SUBJECT_ORIGINAL%in%list]<-i	
	}

	data$SUBJECT_ORIGINAL<-factor(data$SUBJECT_ORIGINAL)
}


					
data$PROTEIN<-factor(data$PROTEIN)	
out<-NULL
outsummary<-NULL
dataafterfit<-NULL

#################################
### how to handle missingness for endogenous

data.l<-data[data$LABEL=="L",]
data.h<-data[data$LABEL=="H",]

missingPeptides<-.checkMissFeature(data.l)

protein.list = tapply ( data.l$FEATURE, data.l$PROTEIN, function ( x ) unique ( as.character ( x ) ) )

missing.results = sapply ( protein.list, function ( x ) any ( x %in% missingPeptides ) )



## Impute for missing endogenous intensity
if ( missing.action == "impute" ) { 
	if(length(missingPeptides) > 0){
	
		dataBySample <- tapply(data.l$ABUNDANCE, list(data.l$SUBJECT_ORIGINAL, data.l$FEATURE), function(x) min(x, na.rm = TRUE))
		dataBySample <- dataBySample[!(apply(dataBySample, 1, function(x) sum(is.na(x))) == ncol(dataBySample)), ]

		if(!is.numeric(dataBySample)){ ## only one subject
			minValuePerSample <- apply(dataBySample, 1, function(x) min(x, na.rm = TRUE))
		}else{
			minValuePerSample<-min(dataBySample, na.rm=TRUE)
		}
	
		imputeValue <- mean(minValuePerSample[!is.infinite(minValuePerSample)], na.rm = TRUE)

		for(i in 1:length(missingPeptides)){
			sub <- data.l[data.l$FEATURE == missingPeptides[i], ]
			t <- tapply(sub$ABUNDANCE, sub$GROUP, function(x) sum(x > 0, na.rm = TRUE))
			missingConds <- names(t)[which(t == 0 | is.na(t))]		
			data.l[data.l$FEATURE %in% missingPeptides[i] & data.l$GROUP %in% missingConds, ]$ABUNDANCE <- imputeValue
		}

		if(length(missingPeptides) > 0){
			number.missing <- length(missingPeptides)
			message(paste("There are ", number.missing, " features (", paste(missingPeptides, collapse = ", "), ") that are missing intensities for an entire condition.  Intensities in these conditions have been imputed with the minimum intensity across all samples.", sep = ""))		
		}
	}
	
	data<-rbind(data.l,data.h)
}


## even though it is not the case, user can do with no interaction ( no else command)


## 
if ( missing.action == "remove" ){
	data<-data[-which(data$FEATURE %in% missingPeptides),]

	message("The features that are missing intensities for an entire condition in Protein will be removed for fitting the model.")

}


## for assigning interference
missing.results = sapply ( protein.list, function ( x ) any ( x %in% missingPeptides ) )

##
	processout<-rbind(processout,c(paste("missing.action : ",missing.action," - okay",sep="")))
	write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)


###==================================================
### start to analyze by protein ID

	### need original group information
	origGroup<-unique(data$GROUP_ORIGINAL)
	
for (i in 1:nlevels(data$PROTEIN)){

	sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]

	# it is important to remove NA first, before we have the correct structure for each factor
	sub<-sub[!is.na(sub$ABUNDANCE),]

	sub$GROUP<-factor(sub$GROUP)
	sub$SUBJECT<-factor(sub$SUBJECT)
	sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
	sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)
	sub$SUBJECT_NESTED<-factor(sub$SUBJECT_NESTED)
	sub$FEATURE<-factor(sub$FEATURE)	
	sub$RUN<-factor(sub$RUN)


	singleFeature<-.checkSingleFeature(sub)
	singleSubject<-.checkSingleSubject(sub)
	TechReplicate<-.checkTechReplicate(sub) ## use for label-free model

	MissGroupByFeature<-.checkMissGroupByFeature(sub)
	MissRunByFeature<-.checkMissRunByFeature(sub)
	MissSubjectByGroup<-.checkRunbyFeature(sub)
	UnequalSubject<-.checkUnequalSubject(sub)
	
	
	## need to assigning whether interaction term is included or not
	remove.interaction<-interference
	
if ( missing.action == "nointeraction" & missing.results [ i ]){
	remove.interaction = FALSE
	
	message(paste(levels(data$PROTEIN)[i]," has some features that are missing intensities for an entire condition in Protein, The additive model (without interaction) will be fitted.", sep = ""))
}

	message(paste("Testing a comparison for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
	
	if(singleFeature){
		message("Protein ",levels(data$PROTEIN)[i]," has only single transition, the simplied model is fitted")
		temp<-try(.fit.model.single(contrast.matrix,sub,labeled,scopeOfTechReplication,scopeOfBioReplication,TechReplicate,repeated,origGroup),silent=TRUE)

	}	

	
	if(!singleFeature){

# not sure it is correct and what it means
#	noRunFeature<-.checkRunbyFeature(sub)
#	if(noRunFeature) unbalanced=TRUE
	
	temp<-try(.fit.model(contrast.matrix,sub,labeled, scopeOfBioReplication,scopeOfTechReplication,interference=remove.interaction,repeated,MissGroupByFeature,MissRunByFeature,UnequalSubject,singleSubject,TechReplicate,featureVar,origGroup),silent=TRUE)

	}


## fix(apr 16)
	if(class(temp)=="try-error") {
		message("error : can't analyze ", levels(data$PROTEIN)[i], " for comparison.")
		tempresult<-list(result=NULL,valueresid=NULL, valuefitted=NULL)
		for(k in 1:nrow(contrast.matrix)){	
			tempresult$result<-rbind(tempresult$result, data.frame(Protein=levels(data$PROTEIN)[i],Label=row.names(contrast.matrix)[k], logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA))
		}
	}else{
		tempresult<-temp
	}

	out<-rbind(out,tempresult$result)

## add residual and fitted after fitting the model
	if(class(temp)=="try-error") {
		sub$residuals<-NA
		sub$fitted<-NA
	}else{
		sub$residuals<-temp$valueresid
		sub$fitted<-temp$valuefitted
	}

	
	## order concerned
#	residuals<-data.frame(temp$valueresid)
#	fitted<-data.frame(temp$valuefitted)
#	sub<-merge(sub,residuals,by="row.names",all=T)
#	rownames(sub)<-sub$Row.names
#	sub<-merge(sub, fitted, by="row.names",all=T)
#	rownames(sub)<-data$Row.names

	dataafterfit<-rbind(dataafterfit,sub)
	
##
	processout<-rbind(processout,c(paste("Finished a comparison for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")")))
	write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)


} ### end protein loop

##
	processout<-rbind(processout,c("Comparisons for all proteins are done.- okay"))
	write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)


##### finalize result
## need to FDR per comparison
out.all<-NULL

out$Label<-as.character(out$Label)
for(i in 1:nrow(contrast.matrix)){
	outsub<-out[out$Label==row.names(contrast.matrix)[i],]
	outsub<-data.frame(outsub,adj.pvalue=p.adjust(outsub$pvalue,method="BH"))
	out.all<-rbind(out.all, outsub)
}
out.all$Label<-factor(out.all$Label)


##
	processout<-rbind(processout,c("Adjust p-values per comparison are calculated - okay."))
	write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)


temp<-data[!is.na(data[,"ABUNDANCE"]),]

if(abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<
abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])){
	colnames(out.all)[3]<-"log2FC"	
}

if(abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])>
abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])){
	colnames(out.all)[3]<-"log10FC"	
}

##
	processout<-rbind(processout,c("Group comparison is done. - okay"))
	write.table(processout, file="Process_groupComparison.txt",row.names=FALSE)


finalout<-list(ComparisonResult=out.all,ModelQC=dataafterfit)
return(finalout)	

}


#############################################
#############################################
### Model-based quality control plot 
#############################################
#############################################

modelBasedQCPlots<-function(data,type,axis.size=10,point.size=3,text.size=7,width=10, height=10,featureName=TRUE,feature.QQPlot="all",which.Protein="all",address=""){
	
	data<-data[!is.na(data$fitted),]
	data<-data[!is.na(data$residuals),]
	
	data$PROTEIN<-factor(data$PROTEIN)	
	
	#### choose Proteins or not
	if(which.Protein!="all"){
		## check which.Protein is name of Protein
		if(is.character(which.Protein)){
			
			temp.name<-which.Protein
		
			## message if name of Protein is wrong.
			if(length(setdiff(temp.name,unique(data$PROTEIN)))>0)
				stop(message(paste("error : Please check protein name. Dataset does not have this protein. -", paste(temp.name, collapse=", "),sep=" ")))
		}

		## check which.Protein is order number of Protein		
		if(is.numeric(which.Protein)){
			
			temp.name<-levels(data$PROTEIN)[which.Protein]
		
			## message if name of Protein is wrong.
			if(length(levels(data$PROTEIN))<max(which.Protein))
				stop(message(paste("error : Please check your selection of proteins. There are ", length(levels(data$PROTEIN))," proteins in this dataset.",sep=" ")))
		}
			
		## use only assigned proteins
		data<-data[which(data$PROTEIN %in% temp.name),]
		data$PROTEIN<-factor(data$PROTEIN)
	}
	
#############################################
### normality, QQ plot
#############################################
	if(type=="QQPlots"){
	
	### test feature.QQPlot is something wrong.
	
	
	#### one QQplot per protein, all features together
		if(feature.QQPlot=="all"){
		
			if(address!=FALSE) pdf(paste(address,"QQPlot_allFeatures.pdf",sep=""))

			for (i in 1:nlevels(data$PROTEIN)){	

				sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]

				qqnorm(sub$residual,main=paste("Normal Q-Q Plot (",unique(sub$PROTEIN),")"))
				qqline(sub$residual, col=2)
			} ## end loop
		
			if(address!=FALSE) dev.off()
		} ## end QQplot by all feature
	
	
	#### panel for each feature,
		if(feature.QQPlot=="byFeature"){
		
			if(address!=FALSE) pdf(paste(address,"QQPlot_byFeatures.pdf",sep=""))

			for (i in 1:nlevels(data$PROTEIN)){	

				sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
			
				## label-free
				if(length(unique(sub$LABEL))==1){
					ptemp<-ggplot(sub, aes_string(sample='residuals',color='FEATURE'))+geom_point(stat="qq",alpha=0.8)+facet_wrap(~FEATURE)+scale_y_continuous('Sample Quantiles')+scale_x_continuous('Theoretical Quantiles')+labs(title=paste("Normal Q-Q Plot (",unique(sub$PROTEIN),")"))+theme(
			panel.background=element_rect(fill='white', colour="black"),
			panel.grid.major = element_line(colour="grey95"),
			panel.grid.minor =element_blank(),
			axis.text.x=element_text(size=axis.size,colour="black"),
			axis.text.y=element_text(size=axis.size,colour="black"),
			axis.ticks=element_line(colour="black"),
			axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
			axis.title.y=element_text(size=axis.size+5,vjust=0.3),
			title=element_text(size=axis.size+8,vjust=1.5),
			strip.text.x=element_text(size=text.size),
			legend.position="none")
					
					print(ptemp)	
				}
			
				## label-based : seperate endogenous and reference
				if(length(unique(sub$LABEL))==2){
				
					### endogenous intensities
					sub.l<-sub[sub$LABEL=="L",]
					ptemp<-ggplot(sub.l, aes_string(sample='residuals',color='FEATURE'))+geom_point(stat="qq",alpha=0.8)+facet_wrap(~FEATURE)+scale_y_continuous('Sample Quantiles')+scale_x_continuous('Theoretical Quantiles')+labs(title=paste("Normal Q-Q Plot (",unique(sub$PROTEIN),") - Endogenous Intensities"))+theme(
			panel.background=element_rect(fill='white', colour="black"),
			panel.grid.major = element_line(colour="grey95"),
			panel.grid.minor =element_blank(),
			axis.text.x=element_text(size=axis.size,colour="black"),
			axis.text.y=element_text(size=axis.size,colour="black"),
			axis.ticks=element_line(colour="black"),
			axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
			axis.title.y=element_text(size=axis.size+5,vjust=0.3),
			title=element_text(size=axis.size+8,vjust=1.5),
			strip.text.x=element_text(size=text.size),
			legend.position="none")
					
					print(ptemp)
				
					### reference intensities
					sub.h<-sub[sub$LABEL=="H",]
					ptemp<-ggplot(sub.h, aes_string(sample='residuals',color='FEATURE'))+geom_point(stat="qq",alpha=0.8)+facet_wrap(~FEATURE)+scale_y_continuous('Sample Quantiles')+scale_x_continuous('Theoretical Quantiles')+labs(title=paste("Normal Q-Q Plot (",unique(sub$PROTEIN),") - Reference Intensities"))+theme(
			panel.background=element_rect(fill='white', colour="black"),
			panel.grid.major = element_line(colour="grey95"),
			panel.grid.minor =element_blank(),
			axis.text.x=element_text(size=axis.size,colour="black"),
			axis.text.y=element_text(size=axis.size,colour="black"),
			axis.ticks=element_line(colour="black"),
			axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
			axis.title.y=element_text(size=axis.size+5,vjust=0.3),
			title=element_text(size=axis.size+8,vjust=1.5),
			strip.text.x=element_text(size=text.size),
			legend.position="none")
					
					print(ptemp)	
				}
			
			} ## end loop
		
			if(address!=FALSE) dev.off()	
		} ### end for QQ by FEATURE
		
	} ### end QQplots


#############################################
#### Residual plot
#############################################
	if(type=="ResidualPlots"){

		y.limdown<-min(data$residuals,na.rm=TRUE)
		y.limup<-max(data$residuals,na.rm=TRUE)

		x.limdown<-min(data$fitted,na.rm=TRUE)
		x.limup<-max(data$fitted,na.rm=TRUE)

		# output separate pdf plot
		if(address!=FALSE) pdf(paste(address,"ResidualPlot.pdf",sep=""))

		for (i in 1:nlevels(data$PROTEIN)){	

			sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
			sub$PEPTIDE<-factor(sub$PEPTIDE)
			sub$FEATURE<-factor(sub$FEATURE)

			ptemp<-ggplot(aes_string(x='fitted', y='residuals', color='FEATURE', shape='LABEL'), data=sub)+geom_point(size=point.size,alpha=0.5)+geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+scale_y_continuous('Residuals',limit=c(y.limdown,y.limup))+scale_x_continuous('Predicted Abundance',limit=c(x.limdown,x.limup))+labs(title=levels(data$PROTEIN)[i])
			
			if(length(unique(sub$LABEL))==2){
				ptemp<-ptemp+scale_shape_manual(values=c(2,19),name="",labels=c("Reference","Endogenous"))
			}else{
				ptemp<-ptemp+scale_shape_manual(values=c(19),name="",labels=c("Endogenous"))
			}
			
			
			if(featureName){
				ptemp<-ptemp+theme(
		panel.background=element_rect(fill='white', colour="black"),
		legend.key=element_rect(fill='white',colour='white'),
		panel.grid.major = element_line(colour="grey95"),
		panel.grid.minor = element_blank(),
		axis.text.x=element_text(size=axis.size,colour="black"),
		axis.text.y=element_text(size=axis.size,colour="black"),
		axis.ticks=element_line(colour="black"),
		axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
		axis.title.y=element_text(size=axis.size+5,vjust=0.3),
		title=element_text(size=axis.size+8,vjust=1.5),
		legend.position=c("top"))+guides(color=guide_legend(ncol=2))
			}else{
				ptemp<-ptemp+theme(
		panel.background=element_rect(fill='white', colour="black"),
		legend.key=element_rect(fill='white',colour='white'),
		panel.grid.major = element_line(colour="grey95"),
		panel.grid.minor = element_blank(),
		axis.text.x=element_text(size=axis.size,colour="black"),
		axis.text.y=element_text(size=axis.size,colour="black"),
		axis.ticks=element_line(colour="black"),
		axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
		axis.title.y=element_text(size=axis.size+5,vjust=0.3),
		title=element_text(size=axis.size+8,vjust=1.5),
		legend.position="none")
			}
			
		print(ptemp)
		}
		
		if(address!=FALSE) dev.off()	
	} ## end residualplots
	
}







#############################################
#############################################
#############################################
#############################################
# Part 2 fit.model
#############################################
#############################################
#############################################
#############################################

#############################################
# check if measurements are missing for entire group
# if yes, length of group and length of contrast won't agree
#############################################

.chechGroupComparisonAgreement<-function(sub1,contrast.matrix){
	tempSub<-as.numeric(as.character(unique(sub1[,c("GROUP")])))
	positionMiss<-setdiff(seq(1,length(contrast.matrix)),tempSub)
	contrast.matrix.sub1<-contrast.matrix[tempSub]
	# either one of the groups for the comparison of interest is not present 

	return(list(sign=length(setdiff(contrast.matrix[tempSub],0))<2,positionMiss=positionMiss))
}


#############################################
# check repeated (case-control? or time-course?)
#############################################

.checkRepeated<-function(data){
	
	data.light<-data[data$LABEL=="L",]	
	subjectByGroup<-table(data.light$SUBJECT_ORIGINAL, data.light$GROUP_ORIGINAL)
	subjectAppearances<-apply(subjectByGroup, 1, function(x) sum(x>0))
	crossedIndicator<-any(subjectAppearances > 1)
	
	return(crossedIndicator)
}

#############################################
# check single subject for both case-control and time-course?
#############################################

.checkSingleSubject<- function(data){
		
	temp<-unique(data[,c("GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
	temp$GROUP_ORIGINAL<-factor(temp$GROUP_ORIGINAL)
	temp1<-xtabs(~GROUP_ORIGINAL,data=temp)
	singleSubject<-all(temp1=="1")
	
	return(singleSubject)	
}

#############################################
# check .checkSingleFeature
#############################################

.checkSingleFeature<- function(data){
	
	sigleFeature<-length(unique(data$FEATURE))<2
	
	return(sigleFeature)	
}

#############################################
# check .checkTechReplicate
#############################################

.checkTechReplicate<- function(data){

	data.light<-data[data$LABEL=="L",]
	temp<-unique(data.light[,c("SUBJECT_NESTED","RUN")])
	temp$SUBJECT_NESTED<-factor(temp$SUBJECT_NESTED)
	temp1<-xtabs(~SUBJECT_NESTED,data=temp)
	TechReplicate<-all(temp1!="1")
	
	return(TechReplicate)	
}


#############################################
# check .checkRunByFeature
#############################################
# it might not be right
.checkRunbyFeature<-function(data){
	
	data.light<-data[data$LABEL=="L",]	
	RunByFeature<-table(data.light$RUN, data.light$FEATURE)
	emptyRow <- apply(RunByFeature,1, sum)
	noRunFeature <- any(emptyRow == 0)
	
	return(noRunFeature)	
}


#############################################
# check .checkMissGroupByFeature
#############################################
.checkMissGroupByFeature<-function(data){
	
	temp<-unique(data[,c("GROUP","FEATURE")])
	temp1<-xtabs(~GROUP,data=temp)
	
	return(any(temp1!=temp1[1]))
}


#############################################
# check .checkMissRunByFeature
#############################################

.checkMissRunByFeature<-function(data){
	
	##temp<-unique(data[,c("RUN","FEATURE")])
	temp<-unique(data[data$LABEL=="L",c("RUN","FEATURE")])
	temp1<-xtabs(~RUN,data=temp)
	#return(any(temp1!=temp1[1]))
	
	return(any(temp1!=length(unique(data$FEATURE))))
}


#############################################
# check .checkMissFeature for label-free missingness
#############################################
.checkMissFeature<-function(data){
	
	dataByPeptide <- tapply(as.numeric(data$ABUNDANCE), list(data$FEATURE, data$GROUP_ORIGINAL), function(x) sum(x > 0, na.rm = TRUE))
	
	missPeptideInd <- apply(dataByPeptide, 1, function(x) any(x == 0 | is.na(x)))
	missingPeptides <- names(missPeptideInd)[missPeptideInd == TRUE]

	return(missingPeptides)
}


#############################################
# check .checkUnequalSubject
#############################################

.checkUnequalSubject<-function(data){

	#temp<-unique(data[,c("GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
	temp<-unique(data[data$LABEL=="L",c("GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
	temp1<-xtabs(~GROUP_ORIGINAL,data=temp)
	
	return(any(temp1!=temp1[1]))
}

########################################################
.fit.model.single<-function(contrast.matrix,data,labeled,scopeOfTechReplication,scopeOfBioReplication,TechReplicate,repeated, origGroup){
	
	if(labeled){ ## labeled==TRUE
	
# when subject is fixed, it is ok using lm function. When we use lm, model with nested or crossed samples have the same results. In the future, if we want to have sample quantification. We need to separate the cases of nested or crossed samples.
		
		### (1) fixed Subject, fixed Run
		if(scopeOfTechReplication=="restricted" & scopeOfBioReplication=="restricted"){
			
			# case-control
			if (!repeated){	
				if(TechReplicate){
					fit.full<-lm(ABUNDANCE ~ SUBJECT_NESTED + GROUP + RUN , data = data)

				}else{
					fit.full<-lm(ABUNDANCE ~ GROUP + RUN , data = data)
				}
			}else{ ### time-course
				if(TechReplicate){
					fit.full<-lm(ABUNDANCE ~ SUBJECT_NESTED + GROUP + RUN , data = data)

				}else{
					fit.full<-lm(ABUNDANCE ~ SUBJECT + GROUP + RUN , data = data)
				}
			}
			

			## get parameter from model
			fixedPara<-.getParameterFixed(fit.full)
	
			#### each comparison
			allout<-NULL
	
			for(k in 1:nrow(contrast.matrix)){
	
				# choose each comparison
				contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
				row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
				GroupComparisonAgreement<-.chechGroupComparisonAgreement(data,contrast.matrix.sub)

				if(GroupComparisonAgreement$sign){
					message("error : results of Protein ", levels(data$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

					out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
				}else{
					contrast<-.make.contrast.based(fit.full,contrast.matrix.sub,data)
					out<-.estimableFixedRandom(fixedPara,contrast)
	
					## any error for out, just NA
					if(is.null(out)){
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
					}else{
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
					}
				}
	
				allout<-rbind(allout, out)
			} # end loop for comparion

			#finalout<-list(result=allout,summary=fit.full)		
			#return(finalout)
		}	### tech replication=restricted


		### (2) random Subject, fixed Run
		if(scopeOfTechReplication=="restricted" & scopeOfBioReplication=="expanded"){

			# case-control
			if (!repeated){	
				if(TechReplicate){
					fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT_NESTED) + GROUP + RUN , data = data)
					# calculate DF of SUBJECT_NESTED
					temp<-unique(data[,c("GROUP","SUBJECT")])
					temp1<-xtabs(~GROUP,data=temp)
					df.full<-sum(temp1[-1]-1)	
				}else{
					message(paste(unique(data$PROTEIN),"This protein has single feature and no technical replicate. Therefore, we can't analysis with expanded scope of biolofical replication."))
				}
			}else{ ### time-course
				if(TechReplicate){
					fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + RUN , data = data)
					# calculate DF of SUBJECT:GROUP
					temp<-unique(data[,c("GROUP","SUBJECT")])
					temp1<-xtabs(~GROUP,data=temp)
					df.full<-sum(temp1[-c(1,2)]-1)
				}else{
					fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT) + GROUP + RUN , data = data)
					df.full<-lm(ABUNDANCE ~ SUBJECT + GROUP + RUN , data = data)$df.residual
				}
			}

			
			## get parameter from model
			randomPara<-.getParameterRandom(fit.full,df.full)
	
			#### each comparison
			allout<-NULL
			for(k in 1:nrow(contrast.matrix)){
	
				# choose each comparison
				contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
				row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
				GroupComparisonAgreement<-.chechGroupComparisonAgreement(data,contrast.matrix.sub)

				if(GroupComparisonAgreement$sign==TRUE){
					message("error : results of Protein ", levels(data$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

					out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
				}else{
					contrast<-.make.contrast.based(fit.full,contrast.matrix.sub,data)
					out<-.estimableFixedRandom(randomPara,contrast)
	
					## any error for out, just NA
					if(is.null(out)){
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
					}else{
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
					}
				}
	
				allout<-rbind(allout, out)
			} # end loop for comparion
	
		} # random subject, fixed run
		
		### (3) fixed subject, random Run
		if(scopeOfTechReplication=="expanded" & scopeOfBioReplication=="restricted"){

			# case-control
			if (!repeated){	
				if(TechReplicate){
					fit.full<-lmer(ABUNDANCE ~ SUBJECT_ORIGINAL + GROUP + (1|RUN) , data = data)
					df.full<-lm(ABUNDANCE ~ SUBJECT_ORIGINAL  + GROUP + RUN , data = data)$df.residual

				}else{
					fit.full<-lmer(ABUNDANCE ~ GROUP + (1|RUN) , data = data)
					df.full<-lm(ABUNDANCE ~ GROUP + RUN , data = data)$df.residual
				}
			}else{ ### time-course
				if(TechReplicate){
					fit.full<-lmer(ABUNDANCE ~ SUBJECT_ORIGINAL + SUBJECT_ORIGINAL:GROUP + GROUP + (1|RUN) , data = data)
					df.full<-lm(ABUNDANCE ~ SUBJECT_ORIGINAL + SUBJECT_ORIGINAL:GROUP + GROUP + RUN , data = data)$df.residual
	
				}else{
					fit.full<-lmer(ABUNDANCE ~ SUBJECT_ORIGINAL + GROUP + (1|RUN) , data = data)
					df.full<-lm(ABUNDANCE ~ SUBJECT_ORIGINAL + GROUP + RUN , data = data)$df.residual
				}
			}
			
			## get parameter from model
			randomPara<-.getParameterRandom(fit.full,df.full)
	
			#### each comparison
			allout<-NULL
			for(k in 1:nrow(contrast.matrix)){
	
				# choose each comparison
				contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
				row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
				GroupComparisonAgreement<-.chechGroupComparisonAgreement(data,contrast.matrix.sub)

				if(GroupComparisonAgreement$sign==TRUE){
					message("error : results of Protein ", levels(data$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

					out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
				}else{
					contrast<-.make.contrast.based(fit.full,contrast.matrix.sub,data)
					out<-.estimableFixedRandom(randomPara,contrast)
	
					## any error for out, just NA
					if(is.null(out)){
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
					}else{
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
					}
				}
	
				allout<-rbind(allout, out)
			} # end loop for comparion
	
		} # random run, fixed subject	
		
		### (4) random subject, random Run
		if(scopeOfTechReplication=="expanded" & scopeOfBioReplication=="expanded"){

			# case-control
			if (!repeated){	
				if(TechReplicate){
					fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT_NESTED) + GROUP + (1|RUN) , data = data)
					# calculate DF of SUBJECT_NESTED
					temp<-unique(data[,c("GROUP","SUBJECT")])
					temp1<-xtabs(~GROUP,data=temp)
					df.full<-sum(temp1[-1]-1)
				}else{
					fit.full<-lmer(ABUNDANCE ~ GROUP + (1|RUN) , data = data)
					df.full<-lm(ABUNDANCE ~ GROUP + RUN , data = data)$df.residual
				}
			}else{ ### time-course
				if(TechReplicate){
					fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + (1|RUN) , data = data)
					# calculate DF of SUBJECT:GROUP
					temp<-unique(data[,c("GROUP","SUBJECT")])
					temp1<-xtabs(~GROUP,data=temp)
					df.full<-sum(temp1[-c(1,2)]-1)	
				}else{
					fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT) + GROUP + (1|RUN) , data = data)
					# calculate DF of SUBJECT:GROUP
					temp<-unique(data[,c("GROUP","SUBJECT")])
					temp1<-xtabs(~GROUP,data=temp)
					df.full<-sum(temp1[-c(1,2)]-1)
				}
			}
			
			## get parameter from model
			randomPara<-.getParameterRandom(fit.full,df.full)
	
			#### each comparison
			allout<-NULL
			for(k in 1:nrow(contrast.matrix)){
	
				# choose each comparison
				contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
				row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
				GroupComparisonAgreement<-.chechGroupComparisonAgreement(data,contrast.matrix.sub)

				if(GroupComparisonAgreement$sign==TRUE){
					message("error : results of Protein ", levels(data$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

					out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
				}else{
					contrast<-.make.contrast.based(fit.full,contrast.matrix.sub,data)
					out<-.estimableFixedRandom(randomPara,contrast)
	
					## any error for out, just NA
					if(is.null(out)){
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
					}else{
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
					}
				}
	
				allout<-rbind(allout, out)
			} # end loop for comparion
	
		} # random run, fixed subject	
		
	} ## label-based

	
## label-free : use only light intensity, need to change contrast.matrix without Group0	
	if(!labeled){
		## use only light
		data2<-subset(data,LABEL=="L")
		data2$GROUP<-factor(data2$GROUP)
		data2$SUBJECT<-factor(data2$SUBJECT)

# when subject is fixed, it is ok using lm function.

# when single feature, consider technical replicates for time-course.

		if (scopeOfBioReplication=="restricted"){
	
			# case-control: impossible for GROUP+SUBJECT for case-control (parameter for the last SUBJECT is NA even though with technical replicates)
			if (!repeated){ ## case-control
				if(!TechReplicate){
					fit.full<-lm(ABUNDANCE ~ GROUP , data = data2)
				}else{
					fit.full<-lm(ABUNDANCE ~ GROUP + SUBJECT , data = data2)

				}
			}else{ ### repeated==TRUE, time-course
				if(!TechReplicate){
					fit.full<-lm(ABUNDANCE ~ GROUP+SUBJECT , data = data2)
				}else{
					fit.full<-lm(ABUNDANCE ~ GROUP+SUBJECT+GROUP:SUBJECT, data = data2) ## SUBJECT==SUBJECT_NESTED here
				}
			} ## time-course


			## get parameter from model
			fixedPara<-.getParameterFixed(fit.full)
	
			#### each comparison
			allout<-NULL
			
			for(k in 1:nrow(contrast.matrix)){
	
				# choose each comparison
				contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
				row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
				GroupComparisonAgreement<-.chechGroupComparisonAgreement(data2,contrast.matrix.sub)

				if(GroupComparisonAgreement$sign==TRUE){
					message("error : results of Protein ", levels(data2$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.\n")
					out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
				}else{
					contrast<-.make.contrast.free.single(fit.full,contrast.matrix.sub,data2)
					out<-.estimableFixedRandom(fixedPara,contrast)
	
					## any error for out, just NA
					if(is.null(out)){
						out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
					}else{
						out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
					}
				}
	
				allout<-rbind(allout, out)
	
			} # end loop for contrast

		} ## restricted, fixed subject	


		if(scopeOfBioReplication=="expanded"){

			# case-control
			if (!repeated){
				if(!TechReplicate){
					message("error : can't analyze with expanded scope of BioReplicates. Use the restricted scope of BioReplicates")
				}else{
					fit.full<-lmer(ABUNDANCE ~ GROUP+(1|SUBJECT) , data = data2)
					df.full<-lm(ABUNDANCE ~ GROUP+SUBJECT , data = data2)$df.residual
				}

			}else{ ## time-course
				if(!TechReplicate){
					fit.full<-lmer(ABUNDANCE ~ GROUP+(1|SUBJECT) , data = data2)
					df.full<-lm(ABUNDANCE ~ GROUP+SUBJECT , data = data2)$df.residual
				}else{
					fit.full<-lmer(ABUNDANCE ~ GROUP+(1|SUBJECT)+(1|GROUP:SUBJECT), data = data2) ## SUBJECT==SUBJECT_NESTED here
					df.full<-lm(ABUNDANCE ~ GROUP+SUBJECT+GROUP:SUBJECT , data = data2)$df.residual
				}
			} ## time-course

	
			## get parameter from model
			randomPara<-.getParameterRandom(fit.full,df.full)
	
			#### each comparison
			allout<-NULL
			
			for(k in 1:nrow(contrast.matrix)){
	
				# choose each comparison
				contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
				row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
				GroupComparisonAgreement<-.chechGroupComparisonAgreement(data2,contrast.matrix.sub)

				if(GroupComparisonAgreement$sign==TRUE){
					message("error : results of Protein ", levels(data2$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

					out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
				}else{
					contrast<-.make.contrast.free.single(fit.full,contrast.matrix.sub,data2)
					out<-.estimableFixedRandom(randomPara,contrast)
	
					## any error for out, just NA
					if(is.null(out)){
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
					}else{
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
					}
				}
	
				allout<-rbind(allout, out)
	
			} # end loop for comparion

		}	## expanded, random subject

	} ## label-free

	if(class(fit.full)=="lm"){  ### lm model
		finalresid<-fit.full$residuals
		finalfitted<-fit.full$fitted.values
	}else{   ### lmer model
		finalresid<-resid(fit.full)
		finalfitted<-fitted(fit.full)
	}
		
	finalout<-list(result=allout,valueresid=finalresid, valuefitted=finalfitted)	
	return(finalout)
			
} ## .fit.model.single



########################################################
.fit.model<-function(contrast.matrix,data,labeled,scopeOfBioReplication,scopeOfTechReplication,interference,repeated,MissGroupByFeature,MissRunByFeature,UnequalSubject,singleSubject,TechReplicate,featureVar,origGroup){
	

	nrepeats=3
	labeled=labeled
	i<-NULL
	
	## label-free : use only light intensity, need to change contrast.matrix without Group0	
	if(!labeled){
		## use only light
		data2<-subset(data,LABEL=="L")
		data2$GROUP<-factor(data2$GROUP)
		data2$SUBJECT<-factor(data2$SUBJECT)

		# when single subject, there is no need to distinguish case-control and time course; fixed subject or random subject, and no need to have GXF since there is not enough degree of freedom, hence no separatation of interference	
		if (singleSubject){
			if(!TechReplicate){
				fit.full<-lm(ABUNDANCE ~ FEATURE + GROUP , data = data2)
			}else{
				fit.full<-lm(ABUNDANCE ~ FEATURE + GROUP + GROUP:FEATURE , data = data2)
			}
			
			## make equal variance for feature
			if(featureVar){
				fit.full<-.iter.wls.fit.model(data=data2,fit=fit.full,nrepeats)
			}

			## get parameter from model
			fixedPara<-.getParameterFixed(fit.full)
	
			#### each comparison
			allout<-NULL
				
			for(k in 1:nrow(contrast.matrix)){
	
				# choose each comparison
				contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
				row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
				GroupComparisonAgreement<-.chechGroupComparisonAgreement(data2,contrast.matrix.sub)

				if(GroupComparisonAgreement$sign==TRUE){
					message("error : results of Protein ", levels(data2$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

					out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
				}else{
					contrast<-.make.contrast.free(fit.full,contrast.matrix.sub,data2)
					out<-.estimableFixedRandom(fixedPara,contrast)
	
					## any error for out, just NA
					if(is.null(out)){
						out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
					}else{
						out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
					}
				}
	
				allout<-rbind(allout, out)
	
			} # end loop for contrast
	
		} ## singleSubject is TRUE

		if (!singleSubject){
	
		# when subject is fixed, it is ok using lm function. When we use lm, model with nested or crossed samples have the same results. In the future, if we want to have sample quantification. We need to separate the cases of nested or crossed samples.

		#===========
		# (1) Subject:F; case-control/time course
		# since it's fixed effect model, the results of nested model and crossed model are the same
		# but when there is unequal subject, the SE is slightly different
		# also no missing patterns (GXF, G, S) will stop the test
		#==========	
			if (scopeOfBioReplication=="restricted"){
	
				# case-control
				if (!repeated){	
					if(interference){
     					fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + FEATURE:GROUP , data = data2)
					}else{
	  					fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP, data = data2)
					}
				}else{	## time-course
					if(interference){
      					fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + FEATURE:GROUP , data = data2)
      					#fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT + SUBJECT:GROUP + GROUP + FEATURE:GROUP , data = data2)
					}else{
	  					fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP, data = data2)
	  					#fit.full<-lm(ABUNDANCE ~ FEATURE + SUBJECT + SUBJECT:GROUP + GROUP, data = data2)

					}
				}
	
	     		## make equal variance for feature
	     		if(featureVar){
	     			fit.full<-.iter.wls.fit.model(data=data2,fit=fit.full,nrepeats)
	     		}
	     		
	     		
				## get parameter from model before all comparison
				fixedPara<-.getParameterFixed(fit.full)
	
				#### each comparison
				allout<-NULL
				
				for(k in 1:nrow(contrast.matrix)){
	
					# choose each comparison
					contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
					row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
					GroupComparisonAgreement<-.chechGroupComparisonAgreement(data2,contrast.matrix.sub)

					if(GroupComparisonAgreement$sign==TRUE){
						message("error : results of Protein ", levels(data2$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

						out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
					}else{
						contrast<-.make.contrast.free(fit.full,contrast.matrix.sub,data2)
						out<-.estimableFixedRandom(fixedPara,contrast)
	
						## any error for out, just NA
						if(is.null(out)){
							out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
						}else{
							out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
						}
					}
	
					allout<-rbind(allout, out)
	
				} # end loop for contrast

			} ## bioRep is restricted


	#===========
	# (2) Subject:R; 
	# since it's mixed effect model, the results of nested model and crossed model are not the same
	# df of nested and crossed are not the same as well.
	# also missing patterns (GXF) will stop the test, but (G, S) will not
	#==========	
			if (scopeOfBioReplication=="expanded"){

				# case-control
				if (!repeated){
					if(interference){
						if(!MissGroupByFeature){
							
							# SE is slightly different	
      						fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + FEATURE:GROUP , data = data2)
      						temp<-unique(data2[,c("GROUP","SUBJECT")])
	  						temp1<-xtabs(~GROUP,data=temp)
      						df.full<-sum(temp1-1)
      					}else{  ## MissGroupByFeature==TRUE
      						
							# SE is slightly different	
      						fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP, data = data2)
      						temp<-unique(data2[,c("GROUP","SUBJECT")])
	  						temp1<-xtabs(~GROUP,data=temp)
      						df.full<-sum(temp1-1)
      					}
					}else{ ## interference==FALSE
	  					fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP, data = data2)
	  					temp<-unique(data2[,c("GROUP","SUBJECT")])
	  					temp1<-xtabs(~GROUP,data=temp)
      					df.full<-sum(temp1-1)
					}
				}else{  ##time-course
					if(interference){
						if(!MissGroupByFeature){
									
							# SE is slightly different	
      						fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP)  + GROUP + FEATURE:GROUP , data = data2)
      						temp<-unique(data2[,c("GROUP","SUBJECT")])
	  						temp1<-xtabs(~GROUP,data=temp)
      						df.full<-sum(temp1[-1]-1)
      					}else{ ## MissGroupByFeature==TRUE
      						
							# SE is slightly different	
      						fit.full<-lmer(ABUNDANCE ~ FEATURE +   (1|SUBJECT) + (1|SUBJECT:GROUP)  + GROUP, data = data2)
     						temp<-unique(data2[,c("GROUP","SUBJECT")])
	  						temp1<-xtabs(~GROUP,data=temp)
      						df.full<-sum(temp1[-1]-1)
      					}
					}else{  ## interference==FALSE
	  					fit.full<-lmer(ABUNDANCE ~ FEATURE +   (1|SUBJECT) + (1|SUBJECT:GROUP)  + GROUP, data = data2)
	  					temp<-unique(data2[,c("GROUP","SUBJECT")])
	  					temp1<-xtabs(~GROUP,data=temp)
      					df.full<-sum(temp1[-1]-1)
					}
				}
	
	     		## make equal variance for feature
	     		if(featureVar){
	     			fit.full<-.iter.wls.fit.model(data=data2,fit=fit.full,nrepeats)
	     		}
	     		
	     		
				## get parameter from model
				randomPara<-.getParameterRandom(fit.full,df.full)
	
				#### each comparison
				allout<-NULL
					
				for(k in 1:nrow(contrast.matrix)){
	
					# choose each comparison
					contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
					row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
					GroupComparisonAgreement<-.chechGroupComparisonAgreement(data2,contrast.matrix.sub)

					if(GroupComparisonAgreement$sign==TRUE){
						message("error : results of Protein ", levels(data2$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

						out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
					}else{
						contrast<-.make.contrast.free(fit.full,contrast.matrix.sub,data2)
						out<-.estimableFixedRandom(randomPara,contrast)
	
						## any error for out, just NA
						if(is.null(out)){
							out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
						}else{
							out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
						}
					}
					allout<-rbind(allout, out)
	
				} # end loop for contrast
	
			}## bio rep is expanded

		}## singleSub=FALSE

	} ## label -free	



	if(labeled){	
	
	# when single subject, there is no need to distinguish case-control and time course; fixed subject or random subject, and no need to have GXF RXF since there is not enough degree of freedom, hence no separatation of interference	
	
		if (singleSubject){
			
			message(paste("*** ",unique(data$PROTEIN)," has single subject per condition. Only restricted scope of conclusion is possible."))
			
			if(TechReplicate){
				fit.full<-lm(ABUNDANCE ~ FEATURE + GROUP + RUN + GROUP:FEATURE + RUN:FEATURE, data = data)
			}else{
				fit.full<-lm(ABUNDANCE ~ FEATURE + GROUP + RUN, data = data)
			}
			
	     	## make equal variance for feature
	     	if(featureVar){
	   			fit.full<-.iter.wls.fit.model(data=data,fit=fit.full,nrepeats)
	   		}
	     		
	     		
			## get parameter from model
			fixedPara<-.getParameterFixed(fit.full)
	
			#### each comparison
			allout<-NULL
			
			for(k in 1:nrow(contrast.matrix)){
	
				# choose each comparison
				contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
				row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
				GroupComparisonAgreement<-.chechGroupComparisonAgreement(data,contrast.matrix.sub)

				if(GroupComparisonAgreement$sign==TRUE){
					message("error : results of Protein ", levels(data$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

					out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
				}else{
					contrast<-.make.contrast.based(fit.full,contrast.matrix.sub,data)
					out<-.estimableFixedRandom(fixedPara,contrast)
	
					## any error for out, just NA
					if(is.null(out)){
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
					}else{
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
					}
				}
	
				allout<-rbind(allout, out)
	
			} # end loop for contrast

		} ## singleSub is True

		if (!singleSubject){
			
		##########################	
		# (1) Subject: F, Run: F
		##########################	

			if(scopeOfBioReplication=="restricted" & scopeOfTechReplication=="restricted"){
	
				## (1.1) case-control
				if (!repeated){
					if(interference){
						fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = data)
					}else{
						fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + RUN, data = data)
					}
				}else{  ## (1.2) time-course
					if(interference){
						fit.full<-lm(ABUNDANCE ~ FEATURE  +  SUBJECT_NESTED + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = data)
					}else{
						fit.full<-lm(ABUNDANCE ~ FEATURE  +  SUBJECT_NESTED  + GROUP + RUN, data = data)
					}
				}	
	
		     	## make equal variance for feature
	     		if(featureVar){
	     			fit.full<-.iter.wls.fit.model(data=data,fit=fit.full,nrepeats)
	     		}
	     		
	     		
				## get parameter from model
				fixedPara<-.getParameterFixed(fit.full)
	
				#### each comparison
				allout<-NULL
				for(k in 1:nrow(contrast.matrix)){
	
					# choose each comparison
					contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
					row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
					GroupComparisonAgreement<-.chechGroupComparisonAgreement(data,contrast.matrix.sub)

					if(GroupComparisonAgreement$sign==TRUE){
						message("error : results of Protein ", levels(data$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
					}else{
						contrast<-.make.contrast.based(fit.full,contrast.matrix.sub,data)
						out<-.estimableFixedRandom(fixedPara,contrast)
	
						## any error for out, just NA
						if(is.null(out)){
							out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
						}else{
							out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
						}
					}
	
					allout<-rbind(allout, out)	
				} # end loop for comparison
			} ## fSfR
	
	
		##########################	
		# (2) Subject: R, Run: F
		##########################	

			if(scopeOfBioReplication=="expanded"&scopeOfTechReplication=="restricted"){
	
				## (2.1) case-control
				if (!repeated){
					if(interference){	
						if(!MissGroupByFeature & !MissRunByFeature){
							
							fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = data)
							# calculate DF of SUBJECT_NESTED
							temp<-unique(data[,c("GROUP","SUBJECT")])
							temp1<-xtabs(~GROUP,data=temp)
							df.full<-sum(temp1[-1]-1)
						}else{
								
							fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + RUN, data = data)
							# calculate DF of SUBJECT_NESTED
							temp<-unique(data[,c("GROUP","SUBJECT")])
							temp1<-xtabs(~GROUP,data=temp)
							df.full<-sum(temp1[-1]-1)		
						}
					}else{ ##interference==FALSE
						
						fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + RUN, data = data)
						# calculate DF of SUBJECT_NESTED
						temp<-unique(data[,c("GROUP","SUBJECT")])
						temp1<-xtabs(~GROUP,data=temp)
						df.full<-sum(temp1[-1]-1)
					}
				}else{  ## (2.2) time-course
					if(interference){	
						if(!MissGroupByFeature & !MissRunByFeature){
							
							fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = data)
							# calculate DF of SUBJECT:GROUP
							temp<-unique(data[,c("GROUP","SUBJECT")])
							temp1<-xtabs(~GROUP,data=temp)
							df.full<-sum(temp1[-c(1,2)]-1)
						}else{
							fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + RUN, data = data)
							# calculate DF of SUBJECT:GROUP
							temp<-unique(data[,c("GROUP","SUBJECT")])
							temp1<-xtabs(~GROUP,data=temp)
							df.full<-sum(temp1[-c(1,2)]-1)	
						}
					}else{  ##interference==FALSE
						fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP)  + GROUP + RUN, data = data)
						# calculate DF of SUBJECT:GROUP
						temp<-unique(data[,c("GROUP","SUBJECT")])
						temp1<-xtabs(~GROUP,data=temp)
						df.full<-sum(temp1[-c(1,2)]-1)
					}	
				}
	
		     	## make equal variance for feature
	     		if(featureVar){
	     			fit.full<-.iter.wls.fit.model(data=data,fit=fit.full,nrepeats)
	     		}
	     		
	     		
				## get parameter from model
				randomPara<-.getParameterRandom(fit.full,df.full)
	
				#### each comparison
				allout<-NULL
				
				for(k in 1:nrow(contrast.matrix)){
	
					# choose each comparison
					contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
					row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
					GroupComparisonAgreement<-.chechGroupComparisonAgreement(data,contrast.matrix.sub)

					if(GroupComparisonAgreement$sign==TRUE){
						message("error : results of Protein ", levels(data$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
					}else{
						contrast<-.make.contrast.based(fit.full,contrast.matrix.sub,data)
						out<-.estimableFixedRandom(randomPara,contrast)
	
						## any error for out, just NA
						if(is.null(out)){
							out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
						}else{
						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
						}	
					}
	
					allout<-rbind(allout, out)	
					
				} # end loop for comparison
	
			}	#rSfR
	
		##########################	
		# (3) Subject: F, Run: R
		# only SUBJECT_ORIGINAL + SUBJECT_ORIGINAL:GROUP can work for model fitting
		# therefore, there is no need to separate case-control and time course
		# however, when there is unequal subject, we need to fit additive model of removing SUBJECT by GROUP
		##########################	

			if(scopeOfBioReplication=="restricted"&scopeOfTechReplication=="expanded"){
	
				## (3.1) equal subject per group
				if (!UnequalSubject){
					if(interference){
						if(!MissGroupByFeature & !MissRunByFeature){
							
							fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL + SUBJECT_ORIGINAL:GROUP + GROUP + (1|RUN) + FEATURE:GROUP +  (1|FEATURE:RUN), data = data)
							df.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT + GROUP + GROUP:SUBJECT + RUN + FEATURE:GROUP +  FEATURE:RUN, data = data)$df.residual
						}else{
			
							fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL + SUBJECT_ORIGINAL:GROUP + GROUP + (1|RUN), data = data)
							df.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT + GROUP + GROUP:SUBJECT + RUN, data = data)$df.residual
						}
					}else{  ## interference
						fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL + SUBJECT_ORIGINAL:GROUP + GROUP + (1|RUN), data = data)
						df.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT + GROUP + GROUP:SUBJECT + RUN, data = data)$df.residual
					}
				}else{ ## (3.2) unequal subject per group
					if(interference){		
						if(!MissGroupByFeature & !MissRunByFeature){

							fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL + GROUP + (1|RUN) + FEATURE:GROUP +  (1|FEATURE:RUN), data = data)
							df.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT + GROUP + GROUP:SUBJECT + RUN + FEATURE:GROUP +  FEATURE:RUN, data = data)$df.residual
						}else{
							fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL  + GROUP + (1|RUN), data = data)
							df.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT + GROUP + GROUP:SUBJECT  + RUN, data = data)$df.residual
						}			
					}else{ ## interference==FALSE
						fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL  + GROUP + (1|RUN), data = data)
						df.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT + GROUP + GROUP:SUBJECT + RUN, data = data)$df.residual
					}
				}	
	
		     	## make equal variance for feature
	     		if(featureVar){
	     			fit.full<-.iter.wls.fit.model(data=data,fit=fit.full,nrepeats)
	     		}
	     		
	     		
				## get parameter from model
				randomPara<-.getParameterRandom(fit.full,df.full)
	
				#### each comparison
				allout<-NULL
				
				for(k in 1:nrow(contrast.matrix)){
	
					# choose each comparison
					contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
					row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
					GroupComparisonAgreement<-.chechGroupComparisonAgreement(data,contrast.matrix.sub)

					if(GroupComparisonAgreement$sign==TRUE){
						message("error : results of Protein ", levels(data$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
					}else{
						contrast<-.make.contrast.based(fit.full,contrast.matrix.sub,data)
						out<-.estimableFixedRandom(randomPara,contrast)
	
						## any error for out, just NA
						if(is.null(out)){
							out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
						}else{
							out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
						}
					}
					allout<-rbind(allout, out)	
				} # end loop for comparison
			} #fSrR	
	
		##########################	
		# (4) Subject: R, Run: R
		##########################	

			if(scopeOfBioReplication=="expanded" & scopeOfTechReplication=="expanded"){
	
				## (4.1) case-control
				if (!repeated){
					if(interference){
						if(!MissGroupByFeature & !MissRunByFeature){
							
							fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + (1|RUN) + FEATURE:GROUP +  (1|FEATURE:RUN), data = data)
							# calculate DF of SUBJECT_NESTED
							temp<-unique(data[,c("GROUP","SUBJECT")])
							temp1<-xtabs(~GROUP,data=temp)
							df.full<-sum(temp1[-1]-1)
						}else{
							
							fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED)  + GROUP + (1|RUN), data = data)
							# calculate DF of SUBJECT_NESTED
							temp<-unique(data[,c("GROUP","SUBJECT")])
							temp1<-xtabs(~GROUP,data=temp)
							df.full<-sum(temp1[-1]-1)
						}
					}else{ ## interference==FALSE
						fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + (1|RUN), data = data)
						# calculate DF of SUBJECT_NESTED
						temp<-unique(data[,c("GROUP","SUBJECT")])
						temp1<-xtabs(~GROUP,data=temp)
						df.full<-sum(temp1[-1]-1)
					}
				}else{  ## (4.2) time-course
					if(interference){	
						if(!MissGroupByFeature & !MissRunByFeature){
							
							fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + (1|RUN) + FEATURE:GROUP +  (1|FEATURE:RUN), data = data)
							# calculate DF of SUBJECT:GROUP
							temp<-unique(data[,c("GROUP","SUBJECT")])
							temp1<-xtabs(~GROUP,data=temp)
							df.full<-sum(temp1[-c(1,2)]-1)
						}else{
							fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + (1|RUN), data = data)
							# calculate DF of SUBJECT:GROUP
							temp<-unique(data[,c("GROUP","SUBJECT")])
							temp1<-xtabs(~GROUP,data=temp)
							df.full<-sum(temp1[-c(1,2)]-1)
						}
					}else{  ## interference==FALSE
						fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP)  + GROUP + (1|RUN), data = data)
						# calculate DF of SUBJECT:GROUP
						temp<-unique(data[,c("GROUP","SUBJECT")])
						temp1<-xtabs(~GROUP,data=temp)
						df.full<-sum(temp1[-c(1,2)]-1)
					}
				}	
	
		     	## make equal variance for feature
	     		if(featureVar){
	     			fit.full<-.iter.wls.fit.model(data=data,fit=fit.full,nrepeats)
	     		}
	     		
	     		
				## get parameter from model
				randomPara<-.getParameterRandom(fit.full,df.full)
	
				#### each comparison
				allout<-NULL
				
				for(k in 1:nrow(contrast.matrix)){
	
					# choose each comparison
					contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
					row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
	
					GroupComparisonAgreement<-.chechGroupComparisonAgreement(data,contrast.matrix.sub)

					if(GroupComparisonAgreement$sign==TRUE){
						message("error : results of Protein ", levels(data$PROTEIN)[i], " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")

						out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
					}else{
						contrast<-.make.contrast.based(fit.full,contrast.matrix.sub,data)
						out<-.estimableFixedRandom(randomPara,contrast)
	
						## any error for out, just NA
						if(is.null(out)){
							out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
						}else{
							out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
						}
					}
	
					allout<-rbind(allout, out)	
				} # end loop for comparison	
			} #rSrR		
		} ## singleSub is FALSE
	}## label-based
	
	if(class(fit.full)=="lm"){  ### lm model
		finalresid<-fit.full$residuals
		finalfitted<-fit.full$fitted.values
	}else{   ### lmer model
		finalresid<-resid(fit.full)
		finalfitted<-fitted(fit.full)
	}
		
	finalout<-list(result=allout,valueresid=finalresid, valuefitted=finalfitted)	

	return(finalout)
}	



#############################################
#############################################
# Part 3 groupComparisonPlots
#############################################
#############################################


groupComparisonPlots<-function(data=data,type=type,sig=0.05,FCcutoff=FALSE,ylimUp=FALSE,ylimDown=FALSE,xlimUp=FALSE,axis.size=10,text.size=4,ProteinName=TRUE,ProteinNameLoc=1,width=10, height=10, which.Comparison="all",address=""){


		#### choose comparison to draw plots
		
		if(which.Comparison!="all"){
			## check which.comparison is name of comparison
			if(is.character(which.Comparison)){
				
				temp.name<-which.Comparison
		
				## message if name of comparison is wrong.
				if(length(setdiff(temp.name,unique(data$Label)))>0)
					stop(message(paste("error : Please check labels of comparions. Result does not have this comparison. -", paste(temp.name, collapse=", "),sep=" ")))
			}
					
					
			## check which.comparison is order number of comparison
			if(is.numeric(which.Comparison)){
				
				temp.name<-levels(data$Label)[which.Comparison]
				
				## message if name of comparison is wrong.
				if(length(levels(data$Label))<max(which.Comparison))
					stop(message(paste("error : Please check your selection of comparisons. There are ", length(levels(data$Label))," comparisons in this result.",sep=" ")))
			}

			
			## use only assigned proteins
			data<-data[which(data$Label %in% temp.name),]
			
			data$Protein<-factor(data$Protein)
			data$Label<-factor(data$Label)
		}else{
			
			data$Protein<-factor(data$Protein)
			data$Label<-factor(data$Label)
		}
		
		
#######################
## Heatmap
#######################

	if (type=="Heatmap"){
		
		if(address!=FALSE) pdf(paste(address,"Heatmap.pdf",sep=""),width=width,height=height)

		y.limUp <-30
		if(is.numeric(ylimUp)) y.limUp<-ylimUp 

		## when NA, change it
		#data$adj.pvalue[is.na(data$adj.pvalue)]<-1 ## missing will be grey
		
		data$adj.pvalue[data$adj.pvalue< 2^(-y.limUp)]<-2^(-y.limUp)

		## if FCcutoff is assigned, make p-value insignificant.
		if(is.numeric(FCcutoff)){
			if(colnames(data)[3]=="log2FC"){
				data$adj.pvalue[data[,3]<log2(FCcutoff) & data[,3]>(-log2(FCcutoff)) ]<-1
			}
			if(colnames(data)[3]=="log10FC"){
				data$adj.pvalue[data[,3]<log10(FCcutoff) & data[,3]>(-log10(FCcutoff))]<-1
			}
		}

		final<-NULL

		### based on p-value
		for (i in 1:nlevels(data$Label)){
			
			sub<-data[data$Label==levels(data$Label)[i],]
			temp<- -log2(sub$adj.pvalue)*sign(sub[,3])
			
			final<-	data.frame(cbind(final,temp))
		}
		
		obj<-final
		data$Protein<-factor(data$Protein)
		rownames(obj)<-levels(data$Protein)
		colnames(obj)<-levels(data$Label)

		## remove if whole rows or columns are NA
		obj<-obj[rowSums(!is.na(obj))!=0, colSums(!is.na(obj))!=0]

		### remove any NA
		#obj<-obj[!is.na(obj[,1]),]

		blue.red.18 <- maPalette(low = "blue", high = "red", mid = "white", k = 12)
		my.colors <-blue.red.18
		my.colors[my.colors=="#FFFFFF"]<-"gold"

		up<-ceiling(-log10(2^-y.limUp ))
		temp<-10^(-sort(ceiling(seq(3,up,length=10)[c(1,2,3,5,10)]),decreasing = TRUE))
		breaks<-c(temp,sig)
		neg.breaks <- log(breaks, 2)
		my.breaks  <- c(neg.breaks,0,-neg.breaks[6:1])

		### draw color key
		blocks<-c(-breaks,1, breaks[6:1])
		x.at<-seq(-0.05,1.05,length.out=13)

		par(mar=c(3,3,3,3), mfrow=c(3,1),oma=c(3,0,3,0))
		plot.new()
		image(z = matrix(seq(1:length(my.colors)+1), ncol = 1), col = my.colors, xaxt = "n", yaxt = "n")
		mtext("Color Key", side=3,line=1, cex=3)
		mtext("(sign) Adjusted p-value", side=1,line=3,at=0.5, cex=1.7)
		mtext(blocks,side=1,line=1,at=x.at, cex=1)

		### draw heatmap
		par(oma=c(3,0,0,4))
		heatmap.2(as.matrix(obj),
			col=my.colors,
			Colv=FALSE,
			dendrogram="none",
			breaks=my.breaks,
			trace="none",
			na.color="grey", ## missing value will be grey
			cexCol=1.2,cexRow=1.2, # assign text.size as option
			key=FALSE,
			lhei=c(0.1,0.9),lwid=c(0.1,0.9)
		)
		
		if(address!=FALSE) dev.off()
	}
	
	
#######################
## VolcanoPlot
#######################
	if (type=="VolcanoPlot"){
		
		if(address!=FALSE) pdf(paste(address,"VolcanoPlot.pdf",sep=""))

		y.limUp <-30
		if(is.numeric(ylimUp)) y.limUp<-ylimUp 

		## remove the result, NA
		data<-data[!is.na(data$adj.pvalue),]
			
		### group for coloring dots
		if(!FCcutoff){  
			data[data$adj.pvalue>=sig,"colgroup"]<-"black"
			data[data$adj.pvalue<sig & data[,3]>0,"colgroup"]<-"red"
			data[data$adj.pvalue<sig & data[,3]<0,"colgroup"]<-"blue" 
		}

		if(is.numeric(FCcutoff)){
			data$colgroup<-"black"
			
			if(colnames(data)[3]=="log2FC"){
				data[data$adj.pvalue<sig & data[,3]>log2(FCcutoff),"colgroup"]<-"red"
				data[data$adj.pvalue<sig & data[,3]<(-log2(FCcutoff)),"colgroup"]<-"blue"
			}
	
			if(colnames(data)[3]=="log10FC"){
				data[data$adj.pvalue<sig & data[,3]>log10(FCcutoff),"colgroup"]<-"red"
				data[data$adj.pvalue<sig & data[,3]<(-log10(FCcutoff)),"colgroup"]<-"blue"
			}
		}

		data$colgroup<-factor(data$colgroup, levels=c("black","blue","red"))

		## for multiple volcano plots, 
		for(i in 1:nlevels(data$Label)){
				
			sub<-data[data$Label==levels(data$Label)[i],]
			sub$adj.pvalue[sub$adj.pvalue<2^(-y.limUp)]<-2^(-y.limUp)

			## ylimUp
			y.limup<-ceiling(max(-log2(sub[!is.na(sub$adj.pvalue),"adj.pvalue"])))
			if(y.limup < (-log2(sig))) y.limup<- (-log2(sig)+1) ## for too small y.lim
	
			## ylimDown
			y.limdown<-0 ## default is zero
			if(is.numeric(ylimDown)) y.limdown<-ylimDown

			## x.lim
			x.lim<-ceiling(max(abs(sub[!is.na(sub[,3]),3]))) ## log2FC or log10FC
			if(x.lim<3) x.lim<-3
			if(is.numeric(xlimUp)) x.lim<-xlimUp

			### for assigning x in ggplot2
			subtemp<-sub
			colnames(subtemp)[3]<-"logFC"
			subtemp$log2adjp<-(-log2(subtemp$adj.pvalue))

			### Plotting
			ptemp<-ggplot()+geom_point(aes_string(x='logFC', y='log2adjp',color='colgroup',label='Protein'),size=4, data=subtemp)+scale_colour_manual(values=c("black","blue","red"),limits=c("black","blue","red"),breaks=c("black","blue","red"),labels=c("No regulation","Down-regulated","Up-regulated"))+scale_y_continuous('-Log2 (adjusted p-value)', limit=c(y.limdown, y.limup))+labs(title=unique(sub$Label))+theme(
			panel.background=element_rect(fill='white', colour="black"),
			panel.grid.minor = element_blank(),
			axis.text.x=element_text(size=axis.size,colour="black"),
			axis.text.y=element_text(size=axis.size,colour="black"),
			axis.ticks=element_line(colour="black"),
			axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
			axis.title.y=element_text(size=axis.size+5,vjust=0.3),
			title=element_text(size=axis.size+8,vjust=1.5),
			legend.position="bottom",
			legend.key=element_rect(fill='white',colour='white'),
			legend.title=element_blank())

			## x-axis labeling
			if(colnames(sub)[3]=="log2FC") ptemp<-ptemp+scale_x_continuous('Log2 fold change',limit=c(-x.lim,x.lim))
			if(colnames(sub)[3]=="log10FC") ptemp<-ptemp+scale_x_continuous('Log10 fold change',limit=c(-x.lim,x.lim))

			## add protein name
			if(ProteinName){
				ptemp<-ptemp+geom_text(aes_string(x='logFC', y='log2adjp',label='Protein'),data=subtemp,color="black",guide="none",hjust=0.5-sign(sub[,3])*ProteinNameLoc, vjust=0.5,size=text.size)
			} 

			## For legend of linetype for cutoffs
			## first assign line type
			ltypes<-c("type1"="twodash","type2"="dotted")

			## cutoff lines, FDR only
			if(!FCcutoff){ 
	
				sigcut<-data.frame(logFC=seq(-x.lim, x.lim, length.out=20), log2adjp=(-log2(sig)))
	 
				pfinal<-ptemp+geom_line(data=sigcut,aes_string(x='logFC',y='log2adjp'), linetype="twodash", colour="darkgrey",size=0.6)+
                                    scale_linetype_manual(values=ltypes,labels=paste("Adj p-value cutoff(",sig,")",sep=""))
			}

			# cutoff lines, FDR and Fold change cutoff
			if(is.numeric(FCcutoff)){
				if(colnames(sub)[3]=="log2FC"){
		
					## three different lines
					sigcut<-data.frame(logFC=seq(-x.lim, x.lim, length.out=10), log2adjp=(-log2(sig)))
					FCcutpos<-data.frame(logFC=log2(FCcutoff), log2adjp=seq(y.limdown, y.limup, length.out=10))
					FCcutneg<-data.frame(logFC=(-log2(FCcutoff)), log2adjp=seq(y.limdown, y.limup, length.out=10))
		
					## three lines, with order color first and then assign linetype manual
					pfinal<-ptemp+geom_line(data=sigcut,aes_string(x='logFC',y='log2adjp'),linetype="twodash",colour="darkgrey",size=0.6)+geom_line(data=FCcutpos,aes_string(x='logFC',y='log2adjp'),linetype='dotted',colour="darkgrey",size=0.6)+geom_line(data=FCcutneg,aes_string(x='logFC',y='log2adjp'), linetype='dotted',colour="darkgrey",size=0.6)+guides(colour=guide_legend(order=1),linetype=guide_legend(order=2))+scale_linetype_manual(values=ltypes,labels=c(paste("Adj p-value cutoff(",sig,")",sep=""),paste("Fold change cutoff(",FCcutoff,")",sep="")))
				}

				if(colnames(sub)[3]=="log10FC"){
					
					## three different lines
					sigcut<-data.frame(logFC=seq(-x.lim, x.lim, length.out=10), log2adjp=(-log2(sig)))
					FCcutpos<-data.frame(logFC=log10(FCcutoff), log2adjp=seq(y.limdown, y.limup, length.out=10))
					FCcutneg<-data.frame(logFC=(-log10(FCcutoff)), log2adjp=seq(y.limdown, y.limup, length.out=10))
		
					## three lines, with order color first and then assign linetype manual
					pfinal<-ptemp+geom_line(data=sigcut,aes_string(x='logFC',y='log2adjp') ,linetype="twodash",colour="darkgrey",size=0.6)+geom_line(data=FCcutpos,aes_string(x='logFC',y='log2adjp'),linetype='dotted',colour="darkgrey",size=0.6)+geom_line(data=FCcutneg,aes_string(x='logFC',y='log2adjp'),linetype="dotted",colour="darkgrey",size=0.6)+guides(colour=guide_legend(order=1),linetype=guide_legend(order=2))+scale_linetype_manual(values=ltypes,labels=c(paste("Adj p-value cutoff(",sig,")",sep=""),paste("Fold change cutoff(",FCcutoff,")",sep="")))
				}
			}
			print(pfinal)
		} ## end-loop

		if(address!=FALSE) dev.off()
	}	

#######################
## Comparison Plot
#######################
	if (type=="ComparisonPlot"){
		
		
		datatemp<-data[!is.na(data$adj.pvalue),]
		datatemp$Protein<-factor(datatemp$Protein)
		colnames(datatemp)[3]<-"logFC"
				
		if(address!=FALSE) pdf(paste(address,"ComparisonPlot.pdf",sep=""))


		for (i in 1:nlevels(datatemp$Protein)){
			
  			sub<-datatemp[datatemp$Protein==levels(datatemp$Protein)[i],]
  			sub$ciw<-qt(1-sig/2,sub$DF)*sub$SE
  
  			## ylimUp
			y.limup<-ceiling(max(sub$logFC+sub$ciw))
			if(is.numeric(ylimUp)) y.limup<-ylimUp 

			## ylimDown
			y.limdown<-floor(min(sub$logFC-sub$ciw))
			if(is.numeric(ylimDown)) y.limdown<-ylimDown 
  
			ptemp<-ggplot(aes_string(x='Label', y='logFC'), data=sub)+geom_errorbar(aes(ymax = logFC + ciw, ymin=logFC - ciw),data=sub, width=0.1,colour="red")+geom_point(size=3,colour="darkred")+scale_x_discrete('Comparison')+geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+labs(title=levels(datatemp$Protein)[i])+theme(
			panel.background=element_rect(fill='white', colour="black"),
			panel.grid.major.y = element_line(colour="grey95"),
			panel.grid.minor.y = element_blank(),
			axis.text.x=element_text(size=axis.size,colour="black"),
			axis.text.y=element_text(size=axis.size,colour="black"),
			axis.ticks=element_line(colour="black"),
			axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
			axis.title.y=element_text(size=axis.size+5,vjust=0.3),
			title=element_text(size=axis.size+8,vjust=1.5))
	
			if(colnames(data)[3]=="log2FC") ptemp<-ptemp+scale_y_continuous("Log2-Fold Change",limit=c(y.limdown, y.limup))
			if(colnames(data)[3]=="log10FC") ptemp<-ptemp+scale_y_continuous("Log10-Fold Change",limit=c(y.limdown, y.limup))

			print(ptemp)
		} ## end-loop

	if(address!=FALSE) dev.off()
	} ## end Comparison plot

}



#############################################
#############################################
# Part : designSampleSize
#############################################
#############################################


designSampleSize<-function(data=data,numSample=numSample,numPep=numPep,numTran=numTran,desiredFC=desiredFC,FDR=FDR,power=power,scopeOfBioReplication="restricted",interference=TRUE){

	## save process output in each step
	processout<-NULL

	## data format
	rawinput<-c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Run","Intensity")
	
	if(length(setdiff(toupper(rawinput),toupper(colnames(data))))==0){
		processout<-rbind(processout,c(paste("The required input - data : did not process from dataProcess function. - stop")))
		colnames(processout)<-c("Process-designSampleSize function")

		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)

		stop(message("Please use 'dataProcess' first. Then use output of dataProcess function as input in groupComparison."))
	}

	
## check one TRUE or not
if( sum(isTRUE(numSample),isTRUE(numPep),isTRUE(numTran),isTRUE(power))!=1 ){
	processout<-rbind(processout,c(paste("The required input - number of sample or features : Only one value should be TRUE. - stop")))
	write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)

	 stop(message("One of (numSample, numPep, numTran, power) needs to be TRUE"))
}


## check other options
	if(!(scopeOfBioReplication=="restricted" | scopeOfBioReplication=="expanded")){
		processout<-rbind(processout,c(paste("The required input - scopeOfBioReplication : 'scopeOfBioReplication' value is wrong. - stop")))
		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)
		
		stop(message("'scopeOfBioReplication' must be one of \"restricted\" or \"expanded\".
	"))
	} 
	
	if(!(interference==TRUE | interference==FALSE) | !is.logical(interference)){
		processout<-rbind(processout,c(paste("The required input - interference : 'interference' value is wrong. - stop")))
		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)

		 stop(message("'interference' must be one of TRUE or FALSE as a logical value."))
	}
	
	## all input
		processout<-rbind(processout,c(paste("number of sample = ",numSample,sep="")))
		processout<-rbind(processout,c(paste("number of peptide per protein = ",numPep,sep="")))
		processout<-rbind(processout,c(paste("number of transition per peptide = ", numTran,sep="")))
		processout<-rbind(processout,c(paste("Desired fold change = ",paste(desiredFC,collapse=" - "),sep="")))
		processout<-rbind(processout,c(paste("FDR = ",FDR,sep="")))
		processout<-rbind(processout,c(paste("Power = ", power,sep="")))
		processout<-rbind(processout,c(paste("scopeOfBioReplication = ",scopeOfBioReplication,sep="")))
		processout<-rbind(processout,c(paste("interference = ", interference,sep="")))
		colnames(processout)<-c("Process-designSampleSize function")

		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)



## check whether case-control(FALSE) or time-course(TRUE)
	repeated<-.checkRepeated(data)
	
	if(repeated){ 
		processout<-rbind(processout,c(paste("Time course design of experiment - okay")))
	}else{
		processout<-rbind(processout,c(paste("Case control design of experiment - okay")))

	}
	write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)

	
	if (nlevels(data$LABEL)==1){
		
		sigma.error<-NULL	
		VarComponent<-data.frame(Protein=levels(data$PROTEIN),Error=NA,Subject=NA,GroupBySubject=NA)

		for (i in 1:nlevels(data$PROTEIN)){
			
			sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]

			sub<-sub[!is.na(sub$ABUNDANCE),]
			sub$GROUP<-factor(sub$GROUP)
			sub$SUBJECT<-factor(sub$SUBJECT)
			sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
			sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)
			sub$SUBJECT_NESTED<-factor(sub$SUBJECT_NESTED)
			sub$FEATURE<-factor(sub$FEATURE)	
			sub$RUN<-factor(sub$RUN)

			singleFeature<-.checkSingleFeature(sub)
			MissGroupByFeature<-.checkMissGroupByFeature(sub)
			MissRunByFeature<-.checkMissRunByFeature(sub)

			# note: when run is fixed, we can obtain the same variance of error for both case-control and time course studies.
			if(singleFeature){
				fit.full<-lm(ABUNDANCE ~ GROUP , data = sub)
			}else{
				if (nlevels(data$SUBJECT_ORIGINAL)==1){
					fit.full<-lm(ABUNDANCE ~ FEATURE + GROUP, data = sub)
				}
				
				if (nlevels(data$SUBJECT_ORIGINAL)>1){
		
					# (1) Subject:F
					if(scopeOfBioReplication=="restricted"){
							
						##case-control
						if(!repeated & interference){
							fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP  + FEATURE:GROUP , data = sub)
						}
						if(!repeated & !interference){
							fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP  , data = sub)
						}
			
						## time-course
						if(repeated & interference){
							fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT + SUBJECT:GROUP + GROUP  + FEATURE:GROUP , data = sub)
						}
						if(repeated & !interference){
							fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT + GROUP , data = sub)
						}
					}
		
					# (2) Subject:R
					if(scopeOfBioReplication=="expanded"){
			
						## case-control
						if(!repeated){
							if(!MissGroupByFeature & interference){
								fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + FEATURE:GROUP , data = sub)
							}
							if(MissGroupByFeature | !interference){
								fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP, data = sub)
							}
						}else{  ## time-course
							if(!MissGroupByFeature & interference){
								fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) +(1|SUBJECT:GROUP) + GROUP + FEATURE:GROUP , data = sub)
							}
							if(MissGroupByFeature| !interference){
								fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) +(1|SUBJECT:GROUP) + GROUP, data = sub)
							}
						}			
					}
				}
			} ## singleSub is FALSE	

			if(class(fit.full)!="mer"){
				VarComponent[i,"Error"]<-summary.lm(fit.full)$sigma^2
			}else{
				stddev <- c(sapply(VarCorr(fit.full), function(el) attr(el, "stddev")),attr(VarCorr(fit.full), "sc"))
				VarComponent[i,"Error"]<-stddev[names(stddev)==""]^2
				if(sum(names(stddev)%in%"SUBJECT_NESTED.(Intercept)")>0){
					VarComponent[i,"Subject"]<-stddev[names(stddev)=="SUBJECT_NESTED.(Intercept)"]^2
				}
				if(sum(names(stddev)%in%"SUBJECT.(Intercept)")>0){
					VarComponent[i,"Subject"]<-stddev[names(stddev)=="SUBJECT.(Intercept)"]^2
				}
				if(sum(names(stddev)%in%"SUBJECT:GROUP.(Intercept)")>0){
					VarComponent[i,"GroupBySubject"]<-stddev[names(stddev)=="SUBJECT:GROUP.(Intercept)"]^2
				}
			}
		} ## end-loop

		VarComponent[is.na(VarComponent)]<-0	
		median.sigma.error<-median(VarComponent[,"Error"],na.rm=TRUE)
		if(sum(VarComponent[,"GroupBySubject"])>0){
			median.sigma.subject<-median(VarComponent[,"GroupBySubject"],na.rm=TRUE)
		}else{
			median.sigma.subject<-median(VarComponent[,"Subject"],na.rm=TRUE)
		}
		
		###
		processout<-rbind(processout,c("Calculated variance component. - okay"))
		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)


		### power calculation
		if(isTRUE(power)){
			delta<-log2(seq(desiredFC[1],desiredFC[2],0.025))
			desiredFC<-2^delta
			m0_m1=99
			t<-delta/sqrt(2*median.sigma.error/numPep/numTran/numSample+median.sigma.subject/numSample)
   			 powerTemp<-seq(0,1,0.01)
   			 
    		power<-NULL
   			for(i in 1:length(t)){
    			diff<-qnorm(powerTemp)+qnorm(1-powerTemp*FDR/(1+(1-FDR)*m0_m1)/2)-t[i]
    			min(abs(diff),na.rm=TRUE)
    			power[i]<-powerTemp[order(abs(diff))][1]
    		}
    		
    		###
		processout<-rbind(processout,c("Power is calculated. - okay"))
		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)

    		out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power=power)
    		return(out)
    		
    	}

		if(is.numeric(power)){

			# Large portion of proteins are not changing
			m0_m1=99 ## it means m0/m1=99, m0/(m0+m1)=0.99
			alpha<-power*FDR/(1+(1-FDR)*m0_m1)

			### Num Sample calculation
			if (isTRUE(numSample)){
				delta<-log2(seq(desiredFC[1],desiredFC[2],0.025))
				desiredFC<-2^delta
				z_alpha<-qnorm(1-alpha/2)
				z_beta<-qnorm(power)
				aa<-(delta/(z_alpha+z_beta))^2
				numSample<-round((2*median.sigma.error/numPep/numTran+median.sigma.subject)/aa,0)
				CV<-round((2*median.sigma.error/(numSample*numPep*numTran)+median.sigma.subject/numSample)/desiredFC,3)
				
				###
		processout<-rbind(processout,c("The number of sample is calculated. - okay"))
		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)

				out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power,CV)
				return(out)
			}

			### Num Peptide calculation
			if (isTRUE(numPep)){
				delta<-log2(seq(desiredFC[1],desiredFC[2],0.025))
				desiredFC<-2^delta
				z_alpha<-qnorm(1-alpha/2)
				z_beta<-qnorm(power)
				aa<-(delta/(z_alpha+z_beta))^2
				numPep<-round((2*median.sigma.error/numSample/numTran+median.sigma.subject/numSample)/aa,0)
				CV<-round((2*median.sigma.error/(numSample*numPep*numTran)+median.sigma.subject/numSample)/desiredFC,3)
				
				###
		processout<-rbind(processout,c("The number of peptide per protein is calculated. - okay"))
		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)

				out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power,CV)
				return(out)
			}

			### Num Transition calculation
			if (isTRUE(numTran)){
				delta<-log2(seq(desiredFC[1],desiredFC[2],0.025))
				desiredFC<-2^delta
				z_alpha<-qnorm(1-alpha/2)
				z_beta<-qnorm(power)
				aa<-(delta/(z_alpha+z_beta))^2
				numTran<-round((2*median.sigma.error/numSample/numPep+median.sigma.subject/numSample)/aa,0)
				CV<-round((2*median.sigma.error/(numSample*numPep*numTran)+median.sigma.subject/numSample)/desiredFC,3)
				
				###
		processout<-rbind(processout,c("The number of transition per peptide is calculated. - okay"))
		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)

				out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power,CV)
				return(out)
			}
		} # when power is numeric
	} ## label-free


	if (nlevels(data$LABEL)==2){
		
		sigma.error<-NULL	
		VarComponent<-data.frame(Protein=levels(data$PROTEIN),Error=NA,Subject=NA,GroupBySubject=NA)

		for (i in 1:nlevels(data$PROTEIN)){
			
			sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]

			sub<-sub[!is.na(sub$ABUNDANCE),]
			sub$GROUP<-factor(sub$GROUP)
			sub$SUBJECT<-factor(sub$SUBJECT)
			sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
			sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)
			sub$SUBJECT_NESTED<-factor(sub$SUBJECT_NESTED)
			sub$FEATURE<-factor(sub$FEATURE)	
			sub$RUN<-factor(sub$RUN)

			singleFeature<-.checkSingleFeature(sub)
			MissGroupByFeature<-.checkMissGroupByFeature(sub)
			MissRunByFeature<-.checkMissRunByFeature(sub)

			# note: when run is fixed, we can obtain the same variance of error for both case-control and time course studies.

			if(singleFeature){
				fit.full<-lm(ABUNDANCE ~ GROUP , data = sub)
			}else{
				if (nlevels(sub$SUBJECT_ORIGINAL)==1){
					fit.full<-lm(ABUNDANCE ~ FEATURE + GROUP + RUN, data = sub)
				}
	
				if (nlevels(sub$SUBJECT_ORIGINAL)>1){
		
					# (1) Subject: F
					if(scopeOfBioReplication=="restricted"){
			
						## case-control
						if(!repeated){
							if(interference){
							fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = sub)
							}else{
							fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + RUN, data = sub)		
							}
						}else{ ## time-course
							if(interference){
								fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT + SUBJECT:GROUP + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = sub)
							}else{
							fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT + SUBJECT:GROUP + GROUP + RUN, data = sub)
							}	
						}
					}
		
					# (2) Subject: R
					if(scopeOfBioReplication=="expanded"){
			
						## case-control
						if(!repeated){
							if(!MissGroupByFeature & !MissRunByFeature & interference){
								fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = sub)
							}else{
								fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + RUN, data = sub)
							}
						}else{  # time-course
							if(!MissGroupByFeature & !MissRunByFeature & interference){
								fit.full<-lmer(ABUNDANCE ~ FEATURE + (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = sub)
							}else{
								fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + RUN, data = sub)
							}
						}
					}
				}
			}

			if(class(fit.full)!="mer"){
				VarComponent[i,"Error"]<-summary.lm(fit.full)$sigma^2
			}else{
				stddev <- c(sapply(VarCorr(fit.full), function(el) attr(el, "stddev")),attr(VarCorr(fit.full), "sc"))
				VarComponent[i,"Error"]<-stddev[names(stddev)==""]^2
				if(sum(names(stddev)%in%"SUBJECT_NESTED.(Intercept)")>0){
					VarComponent[i,"Subject"]<-stddev[names(stddev)=="SUBJECT_NESTED.(Intercept)"]^2
				}
				if(sum(names(stddev)%in%"SUBJECT.(Intercept)")>0){
					VarComponent[i,"Subject"]<-stddev[names(stddev)=="SUBJECT.(Intercept)"]^2
				}
				if(sum(names(stddev)%in%"SUBJECT:GROUP.(Intercept)")>0){
					VarComponent[i,"GroupBySubject"]<-stddev[names(stddev)=="SUBJECT:GROUP.(Intercept)"]^2
				}
			}
		} ## end-loop

		VarComponent[is.na(VarComponent)]<-0	
		median.sigma.error<-median(VarComponent[,"Error"],na.rm=TRUE)
		if(sum(VarComponent[,"GroupBySubject"])>0){
			median.sigma.subject<-median(VarComponent[,"GroupBySubject"],na.rm=TRUE)
		}else{
			median.sigma.subject<-median(VarComponent[,"Subject"],na.rm=TRUE)
		}
		
		###
		processout<-rbind(processout,c("Calculated variance component. - okay"))
		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)


		### power calculation
		if(isTRUE(power)){
			delta<-log2(seq(desiredFC[1],desiredFC[2],0.025))
			desiredFC<-2^delta
			m0_m1=99
			t<-delta/sqrt((4*median.sigma.error/numPep/numTran/numSample)+(2*median.sigma.subject/numSample))
			
			powerTemp<-seq(0,1,0.01)
			power<-NULL
			
			for(i in 1:length(t)){
				diff<-qnorm(powerTemp)+qnorm(1-powerTemp*FDR/(1+(1-FDR)*m0_m1)/2)-t[i]
				min(abs(diff),na.rm=TRUE)
				power[i]<-powerTemp[order(abs(diff))][1]
			}
    
        	###
		processout<-rbind(processout,c("Power is calculated. - okay"))
		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)

  			out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power=power)
    		return(out)
    	}

		if(is.numeric(power)){
	
			# Large portion of proteins are not changing
			m0_m1=99 ## it means m0/m1=99, m0/(m0+m1)=0.99
			alpha<-power*FDR/(1+(1-FDR)*m0_m1)

			### Num Sample calculation
			if(isTRUE(numSample)){
				delta<-log2(seq(desiredFC[1],desiredFC[2],0.025))
				desiredFC<-2^delta
				z_alpha<-qnorm(1-alpha/2)
				z_beta<-qnorm(power)
				aa<-(delta/(z_alpha+z_beta))^2
				numSample<-round(((4*median.sigma.error/numPep/numTran)+(2*median.sigma.subject))/aa,0)
				CV<-round((4*median.sigma.error/(numSample*numPep*numTran)+(2*median.sigma.subject/numSample))/desiredFC,3)
				
				###
		processout<-rbind(processout,c("The number of sample is calculated. - okay"))
		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)

				out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power,CV)
				return(out)
			}

			### Num Peptide calculation
			if(isTRUE(numPep)){
				delta<-log2(seq(desiredFC[1],desiredFC[2],0.025))
				desiredFC<-2^delta
				z_alpha<-qnorm(1-alpha/2)
				z_beta<-qnorm(power)
				aa<-(delta/(z_alpha+z_beta))^2
				numPep<-round((4*median.sigma.error/numSample/numTran+2*median.sigma.subject/numSample)/aa,0)
				CV<-round((4*median.sigma.error/(numSample*numPep*numTran)+2*median.sigma.subject/numSample)/desiredFC,3)
				
				###
		processout<-rbind(processout,c("The number of peptide per protein is calculated. - okay"))
		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)

				out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power,CV)
				return(out)
			}

			### Num Transition calculation
			if(isTRUE(numTran)){
				delta<-log2(seq(desiredFC[1],desiredFC[2],0.025))
				desiredFC<-2^delta
				z_alpha<-qnorm(1-alpha/2)
				z_beta<-qnorm(power)
				aa<-(delta/(z_alpha+z_beta))^2
				numTran<-round((4*median.sigma.error/numSample/numPep+2*median.sigma.subject/numSample)/aa,0)
				CV<-round((4*median.sigma.error/(numSample*numPep*numTran)+2*median.sigma.subject/numSample)/desiredFC,3)
				
				###
		processout<-rbind(processout,c("The number of transition per peptide is calculated. - okay"))
		write.table(processout, file="Process_designSampleSize.txt",row.names=FALSE)

				out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power,CV)
				return(out)
			}
		} # power is numeric
	} ## label-based
}



#############################################
#############################################
# Part designSampleSizePlots
#############################################
#############################################

designSampleSizePlots<-function(data=data){
	
	if (length(unique(data$numSample))>1) index<-"numSample"
	if (length(unique(data$numPep))>1) index<-"numPep"
	if (length(unique(data$numTran))>1) index<-"numTran"
	if (length(unique(data$power))>1) index<-"power"

	text.size<-1.2	
	axis.size<-1.3	
	lab.size<-1.7

	if (index=="numSample"){

	with(data, {
plot(desiredFC,numSample,lwd=2,xlab="",ylab="",cex.axis=axis.size,type="l",xaxt="n")})
	axis(1,at=seq(min(data$desiredFC),max(data$desiredFC),0.05),labels=seq(min(data$desiredFC),max(data$desiredFC),0.05),cex.axis=axis.size)
	axis(3,at=seq(min(data$desiredFC),max(data$desiredFC),0.05),labels=data$CV[which(data$desiredFC%in%seq(min(data$desiredFC),max(data$desiredFC),0.05))],cex.axis=axis.size)
	mtext("Coefficient of variation, CV",3,line=2.5,cex=lab.size)
	mtext("Desired fold change",1,line=3.5,cex=lab.size)
	mtext("Minimal number of biological replicates",2,line=2.5,cex=lab.size)
	legend("topright",c(paste("Number of peptides is",unique(data$numPep),sep=" "),paste("Number of transitions is",unique(data$numTran),sep=" "),paste("FDR is",unique(data$FDR),sep=" "),paste("Statistical power is",unique(data$power),sep=" ")),bty="n",cex=text.size)
	}

	if (index=="numPep"){

	with(data, {
  plot(desiredFC,numPep,lwd=2,xlab="",ylab="",cex.axis=axis.size,type="l",xaxt="n")})
	axis(1,at=seq(min(data$desiredFC),max(data$desiredFC),0.05),labels=seq(min(data$desiredFC),max(data$desiredFC),0.05),cex.axis=axis.size)
	axis(3,at=seq(min(data$desiredFC),max(data$desiredFC),0.05),labels=data$CV[which(data$desiredFC%in%seq(min(data$desiredFC),max(data$desiredFC),0.05))],cex.axis=axis.size)
	mtext("Coefficient of variation, CV",3,line=2.5,cex=lab.size)
	mtext("Desired fold change",1,line=3.5,cex=lab.size)
	mtext("Minimal number of peptides",2,line=2.5,cex=lab.size)
	legend("topright",c(paste("Number of replicates is",unique(data$numSample),sep=" "),paste("Number of transitions is",unique(data$numTran),sep=" "),paste("FDR is",unique(data$FDR),sep=" "),paste("Statistical power is",unique(data$power),sep=" ")),bty="n",cex=text.size)
	}

	if (index=="numTran"){

	with(data, {
  plot(desiredFC,numTran,lwd=2,xlab="",ylab="",cex.axis=axis.size,type="l",xaxt="n")})
	axis(1,at=seq(min(data$desiredFC),max(data$desiredFC),0.05),labels=seq(min(data$desiredFC),max(data$desiredFC),0.05),cex.axis=axis.size)
	axis(3,at=seq(min(data$desiredFC),max(data$desiredFC),0.05),labels=data$CV[which(data$desiredFC%in%seq(min(data$desiredFC),max(data$desiredFC),0.05))],cex.axis=axis.size)
	mtext("Coefficient of variation, CV",3,line=2.5,cex=lab.size)
	mtext("Desired fold change",1,line=3.5,cex=lab.size)
	mtext("Minimal number of transitions",2,line=2.5,cex=lab.size)
	legend("topright",c(paste("Number of replicates is",unique(data$numSample),sep=" "),paste("Number of peptides is",unique(data$numPep),sep=" "),paste("FDR is",unique(data$FDR),sep=" "),paste("Statistical power is",unique(data$power),sep=" ")),bty="n",cex=text.size)
	}

	if (index=="power"){

with(data, {
  plot(desiredFC,power,lwd=2,xlab="",ylab="",cex.axis=axis.size,type="l",xaxt="n")})
	axis(1,at=seq(min(data$desiredFC),max(data$desiredFC),0.05),labels=seq(min(data$desiredFC),max(data$desiredFC),0.05),cex.axis=axis.size)
	mtext("Desired fold change",1,line=3.5,cex=lab.size)
	mtext("Power",2,line=2.5,cex=lab.size)
	legend("bottomright",c(paste("Number of replicates is",unique(data$numSample),sep=" "),paste("Number of peptides is",unique(data$numPep),sep=" "),paste("Number of transitions is",unique(data$numTran),sep=" "),paste("FDR is",unique(data$FDR),sep=" ")),bty="n",cex=text.size)
	}

#	dev.off()
}


#############################################
#############################################
# Part 5 quantification
#############################################
#############################################


quantification<-function(data,type="Sample",format="matrix"){

	## save process output in each step
	processout<-NULL

	## data format
	rawinput<-c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Run","Intensity")
	
	if(length(setdiff(toupper(rawinput),toupper(colnames(data))))==0){
		processout<-rbind(processout,c(paste("The required input - data : did not process from dataProcess function. - stop")))
		colnames(processout)<-c("Process-quantification function")

		write.table(processout, file="Process_quantification.txt",row.names=FALSE)

		stop(message("Please use 'dataProcess' first. Then use output of dataProcess function as input in groupComparison."))
	}

## check other options
	if(!(type=="Sample" | type=="Group")){
		processout<-rbind(processout,c(paste("The required input - type : 'type' value is wrong. - stop")))
		write.table(processout, file="Process_quantification.txt",row.names=FALSE)
		
		stop(message("'type' must be one of \"Sample\" or \"Group\".
	"))
	} 
	
	if(!(format=="matrix" | format=="long")){
		processout<-rbind(processout,c(paste("The required input - format : 'format' value is wrong. - stop")))
		write.table(processout, file="Process_quantification.txt",row.names=FALSE)
		
		stop(message("'format' must be one of \"matrix\" or \"long\".
	"))
	} 

		## all input
		processout<-rbind(processout,c(paste("type of quantification = ",type,sep="")))
		processout<-rbind(processout,c(paste("format of output = ",format,sep="")))
		colnames(processout)<-c("Process-quantification function")

		write.table(processout, file="Process_quantification.txt",row.names=FALSE)


	## Group quantification
	if (type=="Group"){

		if(nlevels(data$LABEL)==2){	
		
			## set output format
			finalresult<-data.frame(Protein=rep(levels(data$PROTEIN),each=nlevels(data$GROUP)),Condition=rep(c(levels(data$GROUP_ORIGINAL),"Ref"),nlevels(data$PROTEIN)), LogIntensities=NA)
	
			result<-NULL

			for (i in 1: nlevels(data$PROTEIN)){
			
				sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
				sub<-sub[!is.na(sub$ABUNDANCE),]
				sub$GROUP<-factor(sub$GROUP)
				sub$SUBJECT<-factor(sub$SUBJECT)
				sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
				sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)
				sub$SUBJECT_NESTED<-factor(sub$SUBJECT_NESTED)
				sub$FEATURE<-factor(sub$FEATURE)	
				sub$RUN<-factor(sub$RUN)

				##  make result matrix per protein in case some subjects can miss in certain protein	
				sub.result<-data.frame(Protein=rep(levels(data$PROTEIN)[i],each=nlevels(sub$GROUP)),Condition=rep(c(levels(sub$GROUP_ORIGINAL),"Ref"),1),LogIntensities=NA)

				singleFeature<-.checkSingleFeature(sub)
				singleSubject<-.checkSingleSubject(sub)

				if(singleFeature==TRUE){
					fit<-lm(ABUNDANCE ~ GROUP + RUN , data = sub)
				}	
	
				if(singleFeature==FALSE){
	
					if (singleSubject==TRUE){
						fit<-lm(ABUNDANCE ~ FEATURE + GROUP + RUN, data = sub)
					}
					if (singleSubject==FALSE){
						fit<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = sub)
					}
				}

## extrac coefficient : 
## if we add lmer for random
# if(class(fit)=="lm"){ cf <- summary.lm(fit)$coefficients
#	}else{ cf <- as.matrix(fixef(fit)) }
# but now only fixed,

				cf <- summary.lm(fit)$coefficients

				# calculate group quantification for all levels of group
				a=1
		
				for (j in 1:nlevels(sub$GROUP_ORIGINAL)){
					contrast.matrix<-rep(0,nlevels(sub$GROUP_ORIGINAL))
					contrast.matrix[j]<-1
					contrast<-.make.contrast.group.quantification(fit,contrast.matrix,sub)
		
					## instead of result, use sub.result
					sub.result[a,3]<-.estimableQuantification(cf,contrast)
					a=a+1
				}
				contrast<-.make.contrast.group.quantification.reference(fit,contrast.matrix,sub)
				sub.result[a,3]<-.estimableQuantification(cf,contrast)
				a=a+1
		
				result<-rbind(result, sub.result)
				
				##
	processout<-rbind(processout,c(paste("Finished a comparison for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")")))
	write.table(processout, file="Process_quantification.txt",row.names=FALSE)

			}	## end loop for each protein
		}	## label-based
	

		if(nlevels(data$LABEL)==1){	

			finalresult<-data.frame(Protein=rep(levels(data$PROTEIN),each=nlevels(data$GROUP_ORIGINAL)),Condition=rep(c(levels(data$GROUP_ORIGINAL)),nlevels(data$PROTEIN)),LogIntensities=NA)
	
			result<-NULL

			for (i in 1: nlevels(data$PROTEIN)){
				sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
				sub<-sub[!is.na(sub$ABUNDANCE),]
				sub$GROUP<-factor(sub$GROUP)
				sub$SUBJECT<-factor(sub$SUBJECT)
				sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
				sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)
				sub$SUBJECT_NESTED<-factor(sub$SUBJECT_NESTED)
				sub$FEATURE<-factor(sub$FEATURE)	
				sub$RUN<-factor(sub$RUN)
				
				sub.result<-data.frame(Protein=rep(levels(data$PROTEIN)[i],each=nlevels(sub$GROUP_ORIGINAL)),Condition=rep(c(levels(sub$GROUP_ORIGINAL)),1),LogIntensities=NA)

				singleFeature<-.checkSingleFeature(sub)
				singleSubject<-.checkSingleSubject(sub)

				if(singleFeature==TRUE){
					fit<-lm(ABUNDANCE ~ GROUP , data = sub)
				}	
				if(singleFeature==FALSE){
	
					if(singleSubject==TRUE){
						fit<-lm(ABUNDANCE ~ FEATURE + GROUP, data = sub)
					}
					if (singleSubject==FALSE){
						fit<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + FEATURE:GROUP, data = sub)
					}
				}

## extrac coefficient : 
## if we add lmer for random
# if(class(fit)=="lm"){ cf <- summary.lm(fit)$coefficients
#	}else{ cf <- as.matrix(fixef(fit)) }
# but now only fixed,

				cf <- summary.lm(fit)$coefficients

				# calculate group quantification for all levels of group
				a=1
				
				for(j in 1:nlevels(sub$GROUP_ORIGINAL)){
					contrast.matrix<-rep(0,nlevels(sub$GROUP_ORIGINAL))
					contrast.matrix[j]<-1
					contrast<-.make.contrast.group.quantification(fit,contrast.matrix,sub)
					sub.result[a,3]<-.estimableQuantification(cf,contrast)
					a=a+1
				}
		
				result<-rbind(result, sub.result)
				
				##
	processout<-rbind(processout,c(paste("Finished a comparison for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")")))
	write.table(processout, file="Process_quantification.txt",row.names=FALSE)
	
			} ## end-loop for each protein	
		}	## label-free

		## make new data.frame in case there are missing group
		for(k in 1:nrow(result)){	
			finalresult[finalresult$Protein==result[k,"Protein"] & finalresult$Condition==result[k,"Condition"],"LogIntensities"]<-result[k,"LogIntensities"]
		}


		### output1 : long format	
		#write.csv(finalresult,file=paste(address,"GroupQuantification_longFormat.csv",sep=""))	

		### output2 : dataMatrix
		finalresult$Condition<-factor(finalresult$Condition, levels=rep(c(levels(data$GROUP_ORIGINAL),"Ref")))
		finalmatrix<-as.data.frame.matrix(xtabs(LogIntensities~Protein+Condition, data=finalresult))
		
		#write.csv(finalmatrix,file=paste(address,"GroupQuantification_dataMatrix.csv",sep=""))	

				##
	processout<-rbind(processout,c("Finish group quantificiation - okay."))
	write.table(processout, file="Process_quantification.txt",row.names=FALSE)
	
	
		if(format=="long") return(finalresult)
		if(format=="matrix") return(finalmatrix)
	}	
	
	
###### sample quantification		
	if (type=="Sample"){
	
		if(nlevels(data$LABEL)==2){	

			data$PROTEIN<-factor(data$PROTEIN)

			finalresult<-data.frame(Protein=rep(levels(data$PROTEIN),each=nlevels(data$SUBJECT_NESTED)),Subject=rep(c(levels(data$SUBJECT_NESTED)[-1],"Ref"),nlevels(data$PROTEIN)),Condition=NA, BioReplicate=NA,LogIntensities=NA,NumFeature=NA,NumPeaks=NA)
	
			result<-NULL

			for (i in 1: nlevels(data$PROTEIN)){
				
				sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
				sub<-sub[!is.na(sub$ABUNDANCE),]
				sub$GROUP<-factor(sub$GROUP)
				sub$SUBJECT<-factor(sub$SUBJECT)
				sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
				sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)
				sub$SUBJECT_NESTED<-factor(sub$SUBJECT_NESTED)
				sub$FEATURE<-factor(sub$FEATURE)	
				sub$RUN<-factor(sub$RUN)
	
				temp<-data.frame(xtabs(~SUBJECT_NESTED, data=sub))
	
				sub.result<-data.frame(Protein=rep(levels(data$PROTEIN)[i],each=nlevels(sub$SUBJECT_NESTED)),Subject=rep(c(levels(sub$SUBJECT_NESTED)[-1],"Ref"),1),LogIntensities=NA, NumFeature=length(unique(sub$FEATURE)),NumPeaks=c(temp[-1,"Freq"],temp[1,"Freq"]))

				singleFeature<-.checkSingleFeature(sub)
				singleSubject<-.checkSingleSubject(sub)

				if(singleFeature==TRUE){
					#fit<-lm(ABUNDANCE ~ GROUP + RUN , data = sub)
					fit<-lm(ABUNDANCE ~ SUBJECT_NESTED , data = sub)

					message("** Protein (",levels(data$PROTEIN)[i],") has only single transition. Results may be unreliable.", "\n")
				}	
	
				if(singleFeature==FALSE){
					if (singleSubject==TRUE){
						fit<-lm(ABUNDANCE ~ FEATURE + GROUP + RUN, data = sub)
					}
					if (singleSubject==FALSE){
						fit<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = sub)
					}
				}

## extrac coefficient : 
## if we add lmer for random
# if(class(fit)=="lm"){ cf <- summary.lm(fit)$coefficients
#	}else{ cf <- as.matrix(fixef(fit)) }
# but now only fixed,

				cf <- summary.lm(fit)$coefficients

				# calculate sample quantification for all levels of sample
				a=1	
					
				for(j in 1:(nlevels(sub$SUBJECT_NESTED)-1)){
					contrast.matrix<-rep(0,(nlevels(sub$SUBJECT_NESTED)-1))
					contrast.matrix[j]<-1
					
					if(singleFeature==TRUE){
						contrast<-.make.contrast.subject.quantification.single(fit,contrast.matrix,sub)
					}else{
						contrast<-.make.contrast.subject.quantification(fit,contrast.matrix,sub)
					}
					
					sub.result[a,3]<-.estimableQuantification(cf,contrast)
					a=a+1
				} ## end subject_nested
				
				contrast<-.make.contrast.subject.quantification.reference(fit,contrast.matrix,sub)
				sub.result[a,3]<-.estimableQuantification(cf,contrast)
				
				result<-rbind(result, sub.result)
				
				##
	processout<-rbind(processout,c(paste("Finished a comparison for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")")))
	write.table(processout, file="Process_quantification.txt",row.names=FALSE)

			} ## end-loop
		} ##label-based	

		if(nlevels(data$LABEL)==1){
			
			finalresult<-data.frame(Protein=rep(levels(data$PROTEIN),each=nlevels(data$SUBJECT_NESTED)),Subject=rep(c(levels(data$SUBJECT_NESTED)),nlevels(data$PROTEIN)),Condition=NA, BioReplicate=NA, LogIntensities=NA,NumFeature=NA,NumPeaks=NA)
			
			result<-NULL
			
			for(i in 1: nlevels(data$PROTEIN)){
				
				sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
				sub<-sub[!is.na(sub$ABUNDANCE),]
				sub$GROUP<-factor(sub$GROUP)
				sub$SUBJECT<-factor(sub$SUBJECT)
				sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
				sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)
				sub$SUBJECT_NESTED<-factor(sub$SUBJECT_NESTED)
				sub$FEATURE<-factor(sub$FEATURE)	
				sub$RUN<-factor(sub$RUN)
				
				temp<-data.frame(xtabs(~SUBJECT_NESTED, data=sub))
				
				sub.result<-data.frame(Protein=rep(levels(data$PROTEIN)[i],each=nlevels(sub$SUBJECT_NESTED)),Subject=rep(c(levels(sub$SUBJECT_NESTED)),1),LogIntensities=NA, NumFeature=length(unique(sub$FEATURE)),NumPeaks=temp$Freq)
				
				singleFeature<-.checkSingleFeature(sub)
				singleSubject<-.checkSingleSubject(sub)

				if(singleFeature==TRUE){
					fit<-lm(ABUNDANCE~ SUBJECT_NESTED , data = sub)
					
					message("** Protein (",unique(sub$PROTEIN),") has only single transition. Results may be unreliable.", "\n")
				}	
				if(singleFeature==FALSE){
					if (singleSubject==TRUE){
						fit<-lm(ABUNDANCE ~ FEATURE + GROUP, data = sub)
					}
					if (singleSubject==FALSE){
						fit<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + FEATURE:GROUP, data = sub)
					}
				}

## extrac coefficient : 
## if we add lmer for random
# if(class(fit)=="lm"){ cf <- summary.lm(fit)$coefficients
#	}else{ cf <- as.matrix(fixef(fit)) }
# but now only fixed,

				cf <- summary.lm(fit)$coefficients

				# calculate sample quantification for all levels of sample
				a=1	
					
				for(j in 1:nlevels(sub$SUBJECT_NESTED)){
					contrast.matrix<-rep(0,nlevels(sub$SUBJECT_NESTED))
					contrast.matrix[j]<-1
					
					if(singleFeature==TRUE){
						contrast<-.make.contrast.subject.quantification.single(fit,contrast.matrix,sub)
					}else{
						contrast<-.make.contrast.subject.quantification(fit,contrast.matrix,sub)
					}
					
					sub.result[a,3]<-.estimableQuantification(cf,contrast)
					a=a+1
				}
				result<-rbind(result, sub.result)
				
				##
	processout<-rbind(processout,c(paste("Finished a comparison for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")")))
	write.table(processout, file="Process_quantification.txt",row.names=FALSE)

			} ## end-loop for each protein	
		} ## label-free

		## make new data.frame in case there are missing group
		name<-unique(data[,c("SUBJECT_NESTED","GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
		name<-name[name$SUBJECT_NESTED!="0.0",]
		name<-name[order(name$SUBJECT_ORIGINAL),]  ## sort by bioReplicates
		name$SUBJECT_NESTED<-as.character(name$SUBJECT_NESTED)
		name$SUBJECT_ORIGINAL<-as.character(name$SUBJECT_ORIGINAL)
		name$GROUP_ORIGINAL<-as.character(name$GROUP_ORIGINAL)

		for(i in 1:nrow(name)){
			finalresult[finalresult$Subject==name$SUBJECT_NESTED[i],"Condition"]<-name$GROUP_ORIGINAL[i]
		}

		for(i in 1:nrow(name)){
			finalresult[finalresult$Subject==name$SUBJECT_NESTED[i],"BioReplicate"]<-name$SUBJECT_ORIGINAL[i]
		}

		for(k in 1:nrow(result)){
			finalresult[finalresult$Protein==result[k,"Protein"] & finalresult$Subject==result[k,"Subject"],c("LogIntensities","NumFeature","NumPeaks")]<-result[k,c("LogIntensities","NumFeature","NumPeaks")]
		}

		## use the level of original SUBJECT_ORIGINAL factoring
		finalresult$BioReplicate<-factor(finalresult$BioReplicate, levels=levels(data$SUBJECT_ORIGINAL))

		## sort by protein and BioReplicate
		finalresult<-finalresult[order(finalresult$Protein,finalresult$BioReplicate),]

		finalresult$Subject<-as.character(finalresult$Subject)
		finalresult$BioReplicate<-as.character(finalresult$BioReplicate)

		finalresult[finalresult$Subject=="Ref","Condition"]<-"Ref"
		finalresult[finalresult$Subject=="Ref","BioReplicate"]<-"Ref"

		## remove Subject (Subject_Nested column)
		finalresultout<-finalresult[,-2] 

		### output1 : long format		
		#write.csv(finalresultout,file=paste(address,"SampleQuantification_longformat.csv",sep=""))	

		### output2 : dataMatrix
		finalresult<-data.frame(finalresult, uniqueID=paste(finalresult$BioReplicate,finalresult$Condition,sep="_"))

		finalresult$uniqueID<-factor(finalresult$uniqueID, levels=unique(finalresult$uniqueID))
		finalmatrix<-as.data.frame.matrix(xtabs(LogIntensities~Protein+uniqueID, data=finalresult))

		if(sum(colnames(finalmatrix)=="Ref_Ref")!=0){
			colnames(finalmatrix)<-c(colnames(finalmatrix)[-length(colnames(finalmatrix))],"Ref")
		}

		#write.csv(finalmatrix,file=paste(address,"SampleQuantification_dataMatrix.csv",sep=""))		##
	processout<-rbind(processout,c("Finish sample quantificiation - okay."))
	write.table(processout, file="Process_quantification.txt",row.names=FALSE)


		if(format=="long") return(finalresultout)
		if(format=="matrix") return(finalmatrix)
	}## end sample quantification		
}
