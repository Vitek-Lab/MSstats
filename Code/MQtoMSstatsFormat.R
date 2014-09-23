library(reshape2)

################################################
################################################
################################################
## evidence : evidence.txt
## annotation : annotation.txt - Raw.file, Condition, BioReplicate, Run, (IsotopeLabelType)
## proteinGroups : proteinGroups.txt . if proteinGroups=NULL, use 'Proteins'. if not, use proteinGroups information for matching Protein group ID

## proteinID : Proteins or proteinGroupID
## useUniquePeptide : remove peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
## summary : max or sum - when there are multiple measurements for certain feature and certain fun, use highest or sum of all.
## fewMeasurements : if 1 or 2 measurements across runs per feature, 'remove' will remove those featuares. It can affected for unequal variance analysis.
## experiment : "DDA" or "SILAC"


MaxQtoMSstatsFormat<-function(evidence, annotation,proteinGroups, proteinID="Proteins", useUniquePeptide=TRUE, summary=max, fewMeasurements="remove", experiment="DDA"){
	
	### evidence.txt file
	infile<-evidence
	
	### annotation.txt : Raw.file, Condition, BioReplicate, Run, (IsotopeLabelType)
	annot<-annotation
	
	################################################
	### 1. remove contaminant, reverse proteinID 
	### Contaminant, Reverse column in evidence
	infile<-infile[-which(infile$Contaminant %in% "+" | infile$Reverse %in% "+") ,]

	### ? Only.identified.by.site column in proteinGroupID? : sometimes, it is not in evidence.txt
	if(is.element("Only.identified.by.site", colnames(infile))){
		infile<-infile[-which(infile$Only.identified.by.site %in% "+") ,]
	}
	
	### need to check proteinGroupID in evidence and proteinGroup.txt the same
	### 'Proteins' vs 'Protein.group.IDs' in evidence file and proteinGroup.txt files ????

	### get useful information
	### ? can remove Retention.time column later
	if(experiment=="SILAC"){
		
		infile <- infile[c("Proteins","Protein.group.IDs","Sequence","Modified.sequence","Charge","Raw.file","Intensity.L","Intensity.H","Retention.time","id")]
		infile.l<-infile[,!(colnames(infile) %in% "Intensity.H")]
		infile.h<-infile[,!(colnames(infile) %in% "Intensity.L")]

		colnames(infile.l)[colnames(infile.l)=="Intensity.L"]<-"Intensity"
		colnames(infile.h)[colnames(infile.h)=="Intensity.H"]<-"Intensity"
		
		### new IsotopeLabelType column
		infile.l$IsotopeLabelType<-"L"
		infile.h$IsotopeLabelType<-"H"

		infile<-rbind(infile.l, infile.h)
		
		rm(infile.l)
		rm(infile.h)
		
	}else{
		
		infile <- infile[c("Proteins","Protein.group.IDs","Sequence","Modified.sequence","Charge","Raw.file","Intensity","Retention.time","id")]
		
	}
	
	## remove "_" at the beginning and end
	infile$Modified.sequence<-gsub("_","",infile$Modified.sequence)
	
	################################################
	## 2. remove peptides which are used in more than one protein
	## we assume to use unique peptide
	################################################
	if(useUniquePeptide){
		
		pepcount<-unique(infile[,c("Proteins","Modified.sequence")]) ## Protein.group.IDs or Sequence
		pepcount$Modified.sequence<-factor(pepcount$Modified.sequence)
		
		## count how many proteins are assigned for each peptide
		structure<-aggregate(Proteins~.,data=pepcount, length)
		removePep<-structure[structure$Proteins!=1,]
		
		## remove the peptides which are used in more than one protein
		if(length(removePep$Proteins!=1)!=0){
			infile<-infile[-which(infile$Modified.sequence %in% removePep$Modified.sequence),]
		}
		
	}
	
	################################################
	## 3. duplicated rows for certain feature and certain runs
	## 	3.1) take highest intensity
	##  3.2) take sum of intensities
	################################################

	## Let's find duplicates
	# first remove NA intensity
	infile<-infile[!is.na(infile$Intensity),]
	### length(unique(infile$Proteins)) ## 5424
	### length(unique(infile$Protein.group.IDs)) ## 4477

	#########################
	### 2.1) general Label-free : one measurement for a feature and a run
	if(experiment=="DDA"){
		## count the number of intensities for feature by runs
		#infile$Feature<-paste(infile$Modified.sequence, infile$Charge, sep="_")
		#structure<-dcast(Feature ~ Raw.file, data=infile, value.var='Intensity')
		#flagduplicate = sum(structure>1)>0	
	
		## take the highest intensity among duplicated or sum of intensities 
		## summary="max" or "sum
		infile_w<-castMaxQToWide.GLF(infile, aggregateFun=summary)
	
		## *** remove features which has less than 2 measurements across runs
		## good to remove before reformatting to long-format
	
		if(fewMeasurements=="remove"){
			infile_w<-removeFeatureWithfew.GLF(infile_w)
		}
	
		## then, go back to long-format
		# good to fill rows with NAs, then now can have balanced data-structure.
		infile_l<-meltMaxQToLong.GLF(infile_w)
		
		## need to set 'IsotopeLabelType' because SILAC already has it.
		infile_l$IsotopeLabelType <- "L"

	}
	
	#########################
	### 2.2) label-free : however, several runs for a sample.
	
	
	
	#########################
	### 2.3) SILAC : two measurements for a feature and a run -> one measurements for a feature and a run and condition
	
	if(experiment=="SILAC"){
		## count the number of intensities for feature by runs
		#infile$Feature<-paste(infile$Modified.sequence, infile$Charge, sep="_")
		#structure<-dcast(Feature ~ Raw.file, data=infile, value.var='Intensity')
		#flagduplicate = sum(structure>1)>0	
	
		## take the highest intensity among duplicated or sum of intensities 
		## summary="max" or "sum
		infile_w<-castMaxQToWide.SILAC(infile, aggregateFun=summary)
	
		## *** remove features which has less than 2 measurements across runs
		## good to remove before reformatting to long-format
	
		if(fewMeasurements=="remove"){
			
			## it is the same across experiments. # measurement per feature. 
			infile_w<-removeFeatureWithfew(infile_w)
		}
	
		## then, go back to long-format
		# good to fill rows with NAs, then now can have balanced data-structure.
		infile_l<-meltMaxQToLong.SILAC(infile_w)
	}


	
	################################################	
	### merge all information
	colnames(infile_l)[1] <- "ProteinName"
	colnames(infile_l)[2] <- "PeptideSequence"
	colnames(infile_l)[3] <- "PrecursorCharge"

	#Add in columns for FramentIon & ProductCharge (all values are NA)
	#Add column for IsotopeLabelType (all "L")
	infile_l$FragmentIon <- NA
	infile_l$ProductCharge <- NA

	#Create Condition & Bioreplicate columns; TODO: fill in with correct values
	if(length(unique(infile_l$IsotopeLabelType))>1){
		infile_l<-merge(infile_l, annot, by=c("Raw.file","IsotopeLabelType"))

	}else{
		infile_l<-merge(infile_l, annot, by="Raw.file")
	}
	
	
	infile_l <- infile_l[,c(c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Raw.file","Intensity"))]
	colnames(infile_l)[9] <- "Run"

	infile_l$PeptideSequence<-factor(infile_l$PeptideSequence)
	infile_l$ProteinName<-factor(infile_l$ProteinName)

	return(infile_l)
}

castMaxQToWide.GLF= function(d_long, aggregateFun=aggregateFun){
  data_w = dcast( Proteins + Modified.sequence + Charge ~ Raw.file, data=d_long, value.var='Intensity', fun.aggregate=aggregateFun, keep=TRUE) 
  ## keep=TRUE : will keep the data.frame value as 1 even though there is no values for certain feature and certain run.
  ## when there is completely missing in certain feature and certain run, '1' will be filled. Therefore put NA instead of 1.
  data_w[data_w==1]<-NA
  return(data_w)
}

castMaxQToWide.SILAC= function(d_long, aggregateFun=aggregateFun){
	
	## check any cell has more than 1
	#data_w = dcast( Proteins + Modified.sequence + Charge + IsotopeLabelType ~ Raw.file, data=d_long, value.var='Intensity') 
	#temp<-data_w[,c(5:ncol(data_w))]
	#head(temp)
	#sum(temp>1)
	#which(temp>1, arr.ind=TRUE)
	#data_w[16300,]
	#d_long[d_long$Modified.sequence=="HIILVLSGK" & d_long$Charge=="2" & d_long$IsotopeLabelType=="L",]
	
  data_w = dcast( Proteins + Modified.sequence + Charge + IsotopeLabelType ~ Raw.file, data=d_long, value.var='Intensity', fun.aggregate=aggregateFun, keep=TRUE) 
  ## keep=TRUE : will keep the data.frame value as 1 even though there is no values for certain feature and certain run.
  
  ## when there is completely missing in certain feature and certain run, '1' will be filled. Therefore put NA instead of 1.
  data_w[data_w==1]<-NA
  return(data_w)
}



meltMaxQToLong.GLF = function(d_wide){
  data_l = melt(d_wide, id.vars=c('Proteins','Modified.sequence','Charge'))
  colnames(data_l)[colnames(data_l)==c("variable","value")]<-c('Raw.file','Intensity')
  return(data_l)
}


meltMaxQToLong.SILAC = function(d_wide){
  data_l = melt(d_wide, id.vars=c('Proteins','Modified.sequence','Charge',"IsotopeLabelType"))
  colnames(data_l)[colnames(data_l)==c("variable","value")]<-c('Raw.file','Intensity')
  return(data_l)
}


removeFeatureWithfew=function(x){
	count.measure = apply (x[,!(colnames(x) %in% c("Proteins","Modified.sequence","Charge"))], 1, function ( x ) length ( x[!is.na(x)] ) ) 
	removeFeatureName<-x[count.measure<3,c("Proteins","Modified.sequence","Charge")]

	x$Feature<-paste(x$Proteins, x$Modified.sequence,x$Charge, sep="_")
	removeFeatureName$Feature<-paste(removeFeatureName$Proteins, removeFeatureName$Modified.sequence, removeFeatureName$Charge, sep="_")

	x<-x[-which(x$Feature %in% removeFeatureName$Feature),]
	x<-x[,-ncol(x)]
	
	return(x)
}