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
    stop("Only MSnSet class can be converted to input format for MSstats.")
  
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
    stop(paste("Please check the variable name. The provided variable name",paste(diff.name,collapse=","), "is not present in the data set.",sep=" "))
  
  
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

dataProcess<-function(raw,logTrans=2, normalization="equalizeMedians",nameStandards=NULL, betweenRunInterferenceScore=FALSE,address="", fillIncompleteRows=TRUE, featureSubset="all", summaryMethod="TMP",equalFeatureVar=TRUE, filterLogOfSum=TRUE, censoredInt="NA",cutoffCensored="minFeatureNRun", MBimpute=TRUE,remove50missing=FALSE,skylineReport=FALSE){
  
  ## save process output in each step
  allfiles<-list.files()
  
  num<-0
  filenaming<-"msstats"
  finalfile<-"msstats.log"
  
  while(is.element(finalfile,allfiles)){
    num<-num+1
    finalfile<-paste(paste(filenaming,num,sep="-"),".log",sep="")
  }
  
  session<-sessionInfo()
  sink("sessionInfo.txt")
  print(session)
  sink()
  
  processout<-as.matrix(read.table("sessionInfo.txt", header=T, sep="\t"))
  write.table(processout, file=finalfile, row.names=FALSE)
  
  processout<-rbind(processout,as.matrix(c(" "," ","MSstats - dataProcess function"," "),ncol=1))
  
  #######################################	
  #### make case-insensitive for function options
  #######################################	
  
  normalization<-toupper(normalization)
  #	nameStandards<-toupper(nameStandards)
  #	betweenRunInterferenceScore<-toupper(betweenRunInterferenceScore)
  #	fillIncompleteRows<-toupper(fillIncompleteRows)
  
  
  
  ####################################
  #### Check correct option or input
  
  #### check right column in input
 
  	requiredinput<-c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Run","Intensity")

  
  ## PeptideModifiedSequence is also allowed.
  tempinput<-setdiff(toupper(requiredinput),toupper(colnames(raw)))
  
  if(sum(length(tempinput)!=0 & tempinput!="PEPTIDESEQUENCE")>0){
    processout<-rbind(processout,c(paste("ERROR : The required input : ",paste(tempinput[tempinput!="PEPTIDESEQUENCE"],collapse=", ")," are not provided in input - stop")))
    write.table(processout, file=finalfile,row.names=FALSE)
    
    stop("Please check the required input. The required input needs (ProteinName, PeptideSequence (or PeptideModifiedSequence), PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity)")
  }else{
    processout<-rbind(processout,c("The required input : provided - okay"))
    write.table(processout, file=finalfile,row.names=FALSE)
  }
  
  ## if "PeptideModifiedSequence", change the name as PeptideSequence
  if(is.element(toupper("PeptideModifiedSequence"), toupper(colnames(raw)))){
    colnames(raw)[toupper(colnames(raw))==toupper("PeptideModifiedSequence")]<-"PeptideSequence"	
  }
  
  ####  check whether class of intensity is factor or chaterer, if yes, neec to chage as numeric
  if(is.factor(raw$Intensity) | is.character(raw$Intensity)){	
    suppressWarnings(raw$Intensity<-as.numeric(as.character(raw$Intensity)))
  }
  
  #### check whether the intensity has 0 value or negative value
  if(length(which(raw$Intensity<=0))>0 & !skylineReport){
  	
  	if(is.null(censoredInt)){
  		processout<-rbind(processout,c("ERROR : There are some intensities which are zero or negative values. need to change them. - stop"))
    write.table(processout, file=finalfile,row.names=FALSE)
    
    stop("Intensity has 0 or negative values. Please check these intensities and change them. \n")
  		
  	}else{
  		if(censoredInt=="NA"){
  			processout<-rbind(processout,c("ERROR : There are some intensities which are zero or negative values. need to change them. - stop"))
    write.table(processout, file=finalfile,row.names=FALSE)
    
    stop("Intensity has 0 or negative values. Please check these intensities and change them. \n")
  		}	
  	}
    
  }
  
  #### however skyline report, keep zero and replace with 1,then relace with NA for truncated
  if(skylineReport){
  	
  	## if there is 'Truncated' column,
  	if(is.element(toupper("Truncated"), toupper(colnames(raw)))){
  		## remove truncated peaks
  		raw[!is.na(raw$Truncated) & raw$Truncated==TRUE,"Intensity"]<-NA
  	
    	processout<-rbind(processout,c("There are some truncated peaks. They replaced with NA."))
    	write.table(processout, file=finalfile,row.names=FALSE)
    }
    
  }
  
  #### check unexpected token(":")
  if(length(grep(":",raw$PeptideSequence))!=0){
    processout<-rbind(processout,c("ERROR : colon(:) is invalid in PeptideSequence column. Please replace with other entry. - stop"))
    write.table(processout, file=finalfile,row.names=FALSE)
    
    stop("Colon(:) is invalid in PeptideSequence column. Please replace with other entry. \n")
  }
  
  if(length(grep(":",raw$FragmentIon))!=0){
    processout<-rbind(processout,c("ERROR : colon(:) is invalid in FragmentIon column. Please replace with other entry. - stop"))
    write.table(processout, file=finalfile,row.names=FALSE)
    
    stop("Colon(:) is invalid in FragmentIon column. Please replace with other entry. \n")
  }
  
  ##### check logTrans is 2,10 or not
  if(logTrans!=2 & logTrans!=10){
    processout<-rbind(processout,c("ERROR : Logarithm transformation : log2 or log10 only - stop"))
    write.table(processout, file=finalfile,row.names=FALSE)
    
    stop("Only log2 or log10 are posssible.\n")
  }
  
  ##### check no row for some feature : balanced structure or not
  
  if(!(fillIncompleteRows==TRUE | fillIncompleteRows==FALSE) | !is.logical(fillIncompleteRows)){
    processout<-rbind(processout,c(paste("The required input - fillIncompleteRows : 'fillIncompleteRows' value is wrong. It should be either TRUE or FALSE. - stop")))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    stop("'fillIncompleteRows' must be one of TRUE or FALSE as a logical value.")
  }
  
  ##### check input for summaryMethod
  
  if(sum(summaryMethod==c("linear","TMP","logOfSum"))==0){
    processout<-rbind(processout,c("The required input - summaryMethod : 'summaryMethod' value is wrong. It should be one of 'linear','TMP','logOfSum'. - stop"))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    stop("'summaryMethod' value is wrong. It should be one of 'linear','TMP','logOfSum'.")
  }else{
  	processout<-rbind(processout,c(paste("summaryMethod : ",as.character(summaryMethod), sep="")))
    write.table(processout, file=finalfile, row.names=FALSE)
  }
  
  ##### check input for cutoffCensored
  
  if(sum(cutoffCensored==c("minFeature","minRun","minFeatureNRun"))==0){
    processout<-rbind(processout,c("The required input - cutoffCensored : 'cutoffCensored' value is wrong. It should be one of 'minFeature','minRun','minFeatureNRun'. - stop"))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    stop("'cutoffCensored' value is wrong. It should be one of 'minFeature','minRun','minFeatureNRun'.")
  }else{
  	processout<-rbind(processout,c(paste("cutoffCensored : ",as.character(cutoffCensored), sep="")))
    write.table(processout, file=finalfile, row.names=FALSE)
  }
  
  ##### check input for censoredInt
  
  if(sum(censoredInt==c("0","NA"))==0 & !is.null(censoredInt)){
    processout<-rbind(processout,c("The required input - censoredInt : 'censoredInt' value is wrong. It should be one of '0','NA', NULL. - stop"))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    stop("'censoredInt' value is wrong. It should be one of '0','NA',NULL.")
  }else{
  	processout<-rbind(processout,c(paste("censoredInt : ",as.character(censoredInt), sep="")))
    write.table(processout, file=finalfile, row.names=FALSE)
  }
  
  ##### need the names of global standards
  if(!is.element("NONE",normalization) & !is.element("FALSE",normalization) & is.element("GLOBALSTANDARDS",normalization) & is.null(nameStandards)){
    
    processout<-rbind(processout,c("ERROR : For normalization with global standards, the names of global standards are needed. Please add 'nameStandards' input."))
    write.table(processout, file=finalfile,row.names=FALSE)
    
    stop ("For normalization with global standards, the names of global standards are needed. Please add 'nameStandards' input." )
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
  
  ## before remove, get PeptideSequence and combination of PeptideSequence and precursorcharge for global standard normalization
  tempPeptide<-unique(raw[,c("PEPTIDESEQUENCE", "PRECURSORCHARGE")])
  tempPeptide$PEPTIDE<-paste(tempPeptide$PEPTIDESEQUENCE,tempPeptide$PRECURSORCHARGE,sep="_")
  
  rm(raw)
  
  
  ## assign peptide, transition
  raw.temp<-data.frame(raw.temp,PEPTIDE=paste(raw.temp$PEPTIDESEQUENCE,raw.temp$PRECURSORCHARGE,sep="_"), TRANSITION=paste(raw.temp$FRAGMENTION, raw.temp$PRODUCTCHARGE,sep="_"))
  
  if(length(unique(raw.temp$ISOTOPELABELTYPE))>2){
    processout<-rbind(processout,c("ERROR : There are more than two levels of labeling. So far, only label-free or reference-labeled experiment are supported. - stop"))
    write.table(processout, file=finalfile,row.names=FALSE)
    
    stop("Statistical tools in MSstats are only proper for label-free or with reference peptide experiments.")
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
  write.table(processout, file=finalfile,row.names=FALSE)
  
  
    
  work<-data.frame(work,INTENSITY=raw.temp$Intensity)
  
  
  ############## log transformation
  work$ABUNDANCE<-work$INTENSITY
  
  ## now, INTENSITY keeps original values.
    
  ## change zero with 1 for log transformation
  ## NA means no observation. assume that spectral tools are not report if no observation. zero means detected but zero. 
  work[!is.na(work$ABUNDANCE) & work$ABUNDANCE==0,"ABUNDANCE"]<-1
  	
  processout<-rbind(processout,c("There are some intensities which are zero. Intensities with zero are replaced with 1 in order to do log transformation."))
  write.table(processout, file=finalfile,row.names=FALSE)
    
  ## based on logTrans option, assign log transformation
  ## remove log2 or log10 intensity
  if(logTrans==2){
    work$ABUNDANCE<-log2(work$ABUNDANCE)
  } else if(logTrans==10){
    work$ABUNDANCE<-log10(work$ABUNDANCE)
  } 	
  
  processout<-rbind(processout,c(paste("Logarithm transformation: log",logTrans," transformation is done - okay", sep="")))
  write.table(processout, file=finalfile,row.names=FALSE)
  
  
  ############## Check multi-method or not : multiple run for a replicate
  
  #standardFeature<-unique(work[work$RUN=="1","FEATURE"]) ## if some feature are missing for this spedific run, it could be error. that is why we need balanced design.
  
  #countdiff = tapply (work$FEATURE, work$RUN, function ( x ) length(setdiff(unique(x),standardFeature)) ) 
  
  ## whether multirun or not : we assume that different method has completely different feature
  work$RUN<-factor(work$RUN)
  multirun<-.countMultiRun(work)
  checkMultirun<-any(multirun==0)
  
  ##############################
  #### if multirun, make new column 'method'
  
  if(checkMultirun){ ## when checkMultirun is TRUE, it means there are more than 1 method.
    
    ## make column 'method'
    work$METHOD<-1
    numincreasing=1
    
    ## get run which has different feature names from run1
    while(length(multirun[multirun==0])!=0){ ## until there is no more unique feature per run
      nextmethod<-names(multirun[multirun==0])
      numincreasing<-numincreasing+1
      work[which(work$RUN %in% nextmethod),"METHOD"]<-numincreasing
      
      worktemp<-work[which(work$RUN %in% nextmethod),]
      worktemp$RUN<-factor(worktemp$RUN)
      
      multirun<-.countMultiRun(worktemp)	
    }
    
    processout<-rbind(processout,c(paste("Multiple methods are existed : ", length(unique(work$METHOD))," methods per MS replicate.")))
    write.table(processout, file=finalfile,row.names=FALSE)
    
  }else{
    work$METHOD<-1
  }
  
  ### check messingness for multirun
  
  
  
  ##############################
  ## check no value for some feature : balanced structure or not
  ## need to separate label-free or label-based
  
  processout<-rbind(processout,c(paste("fillIncompleteRows = ",fillIncompleteRows,sep="")))
  write.table(processout, file=finalfile, row.names=FALSE)
  
  
  ## only 1 method
  
  if(!checkMultirun){
    
    ## label-free experiments
    if(nlevels(work$LABEL)==1){
      
      ## get feature by Run count of data
      structure = tapply ( work$ABUNDANCE, list ( work$FEATURE, work$RUN ) , function ( x ) length ( x ) ) 
      
      ## structure value should be 1 for label-free, if not there are missingness. if more there are duplicates.
      
      flagmissing = sum(is.na(structure))>0
      flagduplicate = sum(structure[!is.na(structure)]>1)>0
      
      ### if there is missing rows
      if ( flagmissing ){
        processout<-rbind(processout,c("CAUTION: the input dataset has incomplete rows. If missing peaks occur they should be included in the dataset as separate rows, and the missing intensity values should be indicated with ’NA’. The incomplete rows are listed below."))
        write.table(processout, file=finalfile,row.names=FALSE)
        
        message("CAUTION : the input dataset has incomplete rows. If missing peaks occur they should be included in the dataset as separate rows, and the missing intensity values should be indicated with ’NA’. The incomplete rows are listed below.")
        
        ## first, which run has missing	
        runstructure<-apply ( structure, 2, function ( x ) sum ( is.na ( x ) ) ) > 0
        
        ## get the name of Run
        runID<-names(runstructure[runstructure==TRUE])
        
        ## for missign row, need to assign before looping
        missingwork<-NULL
        
        ## then for each run, which features are missing,
        for(j in 1:length(runID)){
          
          # get subject, group information for this run
          nameID<-unique(work[work$RUN==runID[j],c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
          
          # get feature ID
          featureID<-structure[,colnames(structure)==runID[j]]
          
          # get feature ID which has no measuremnt.
          finalfeatureID<-featureID[is.na(featureID)]
          
          # print features ID	 	
          message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
          
          ## save in process file.
          processout<-rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
          write.table(processout, file=finalfile,row.names=FALSE)
          
          ### add missing rows if option is TRUE
          if(fillIncompleteRows){
            
            tempTogetfeature<-work[which(work$FEATURE %in% names(finalfeatureID)),]
            
            # get PROTEIN and FEATURE infomation
            tempfeatureID<-unique(tempTogetfeature[,c("PROTEIN","PEPTIDE","TRANSITION","FEATURE")])
            
            # merge feature info and run info as 'work' format
            tempmissingwork<-data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, RUN=nameID$RUN, GROUP=nameID$GROUP, SUBJECT=nameID$SUBJECT, SUBJECT_NESTED=nameID$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD)	
            
            # merge with tempary space, missingwork
            missingwork<-rbind(missingwork,tempmissingwork)
          } ## end fillIncompleteRows options
        } ## end loop for run ID
        
        if(fillIncompleteRows){
          
          # merge with work
          ## in future, use rbindlist?? rbindlist(list(work, missingwork))
          work<-rbind(work,missingwork)
          
          ## print message
          message("\n DONE : Incomplete rows for missing peaks are added with intensity values=NA. \n")
          
          ## save in process file.
          processout<-rbind(processout,"Incomplete rows for missing peaks are added with intensity values=NA. - done, Okay")
          write.table(processout, file=finalfile,row.names=FALSE)
          
        }else{
          
          ## save in process file.
          processout<-rbind(processout,"Please check whether features in the list are generated from spectral processing tool. Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
          write.table(processout, file=finalfile,row.names=FALSE)
          
          stop("Please check whether features in the list are generated from spectral processing tool or not. Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
          
        }
      } ## end for flag missing
      
      ########### if there are duplicates measurements
      if(flagduplicate){
        
        ## first, which run has duplicates
        runstructure<-apply ( structure, 2, function ( x ) sum ( x[!is.na(x)] > 1 )>0 )
        
        runID<-names(runstructure[runstructure==TRUE])
        
        ## then for each run, which features have duplicates,
        for(j in 1:length(runID)){
          
          nameID<-unique(work[work$RUN==runID[j],c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
          
          featureID<-structure[,colnames(structure)==runID[j]]
          finalfeatureID<-featureID[!is.na(featureID) & featureID>1]
          
          message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has multiple rows (duplicate rows) for some features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
          
          ## save in process file.
          processout<-rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has multiple rows (duplicate rows) for some features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
          write.table(processout, file=finalfile,row.names=FALSE)
        }
        
        ## save in process file.
        processout<-rbind(processout,"Please remove duplicate rows in the list above. ")
        write.table(processout, file=finalfile,row.names=FALSE)
        
        stop("Please remove duplicate rows in the list above.\n")		
      }## end flag duplicate
      
      ## no missing and no duplicates
      if(!flagmissing & !flagduplicate){
        processout<-rbind(processout,c("Balanced data format with NA for missing feature intensities - okay"))
        write.table(processout, file=finalfile,row.names=FALSE)
      } 
      
      ## end label-free
    }else{ 
      
      ########### label-based experiment
      
      ### count the reference and endobenous separately
      work.l<-work[work$LABEL=="L",]
      work.h<-work[work$LABEL=="H",]
      
      ## get feature by Run count of data
      structure.l = tapply ( work.l$ABUNDANCE, list ( work.l$FEATURE, work.l$RUN ) , function ( x ) length ( x ) ) 
      structure.h = tapply ( work.h$ABUNDANCE, list ( work.h$FEATURE, work.h$RUN ) , function ( x ) length ( x ) ) 
      
      
      ###### first, check some features which completely missing across run
      missingcomplete.l<-NULL	
      missingcomplete.h<-NULL	
      
      
      ### 1. reference peptides
      featurestructure.h<-apply(structure.h, 1, function (x) sum(is.na(x)))
      
      # get feature ID of reference which are completely missing across run
      featureID.h<-names(featurestructure.h[featurestructure.h==ncol(structure.h)])
      
      if(length(featureID.h)>0){
        # print message
        message(paste("CAUTION : some REFERENCE features have missing intensities in all the runs. The completely missing REFERENCE features are ", paste(featureID.h, collapse=", "),". Please check whether features in the list are correctly generated from spectral processing tool. \n",sep=""))
        
        ## save in process file.
        processout<-rbind(processout,c(paste("CAUTION : some REFERENCE features have missing intensities in all the runs. The completely missing REFERENCE features are ", paste(featureID.h, collapse=", "),". Please check whether features in the list are correctly generated from spectral processing tool.", sep="")))
        write.table(processout, file=finalfile,row.names=FALSE)
        
        ### add missing rows if option is TRUE
        if(fillIncompleteRows){
          
          # get unique Run information
          nameID<-unique(work.h[,c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
          
          # get PROTEIN and FEATURE information
          # here use whole work dataset
          tempTogetfeature<-work[which(work$FEATURE %in% featureID.h),]
          tempfeatureID<-unique(tempTogetfeature[,c("PROTEIN","PEPTIDE","TRANSITION","FEATURE")])
          
          ## then generate data.frame for missingness,
          for(j in 1:nrow(nameID)){
            
            # merge feature info and run info as 'work' format
            tempmissingwork<-data.frame(tempfeatureID, LABEL="H",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL[j], SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL[j], RUN=nameID$RUN[j], GROUP=nameID$GROUP[j], SUBJECT=nameID$SUBJECT[j], SUBJECT_NESTED=nameID$SUBJECT_NESTED[j], INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD[j])	
            
            # merge with tempary space, missingwork
            missingcomplete.h<-rbind(missingcomplete.h,tempmissingwork)
          }
        }	## end fillIncompleteRows option
        
      } ## end for reference peptides
      
      ### 2. endogenous peptides
      featurestructure.l<-apply(structure.l, 1, function (x) sum(is.na(x)))
      
      # get feature ID of reference which are completely missing across run
      featureID.l<-names(featurestructure.l[featurestructure.l==ncol(structure.l)])
      
      if(length(featureID.l)>0){
        # print message
        message(paste("CAUTION : some ENDOGENOUS features have missing intensities in all the runs. The completely missing ENDOGENOUS features are ", paste(featureID.l, collapse=", "), ". Please check whether features in the list are correctly generated from spectral processing tool. \n",sep=""))
        
        ## save in process file.
        processout<-rbind(processout,c(paste("CAUTION : some ENDOGENOUS features have missing intensities in all the runs. The completely missing ENDOGENOUS features are ", paste(featureID.l, collapse=", "),". Please check whether features in the list are correctly generated from spectral processing tool. \n",sep="")))
        write.table(processout, file=finalfile,row.names=FALSE)
        
        ### add missing rows if option is TRUE
        if(fillIncompleteRows){
          
          # get unique Run information
          nameID<-unique(work.l[,c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
          
          # get PROTEIN and FEATURE information
          # here use whole work dataset
          tempTogetfeature<-work[which(work$FEATURE %in% featureID.l),]
          tempfeatureID<-unique(tempTogetfeature[,c("PROTEIN","PEPTIDE","TRANSITION","FEATURE")])
          
          ## then generate data.frame for missingness,
          for(j in 1:nrow(nameID)){
            
            # merge feature info and run info as 'work' format
            tempmissingwork<-data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL[j], SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL[j], RUN=nameID$RUN[j], GROUP=nameID$GROUP[j], SUBJECT=nameID$SUBJECT[j], SUBJECT_NESTED=nameID$SUBJECT_NESTED[j], INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD[j])	
            
            # merge with tempary space, missingwork
            missingcomplete.l<-rbind(missingcomplete.l,tempmissingwork)
          }
        }	## end fillIncompleteRows option
      } ## end endogenous peptides
      
      
      ###### second, check other some missingness
      
      ## for missign row, need to assign before looping. need to assign at the beginning because it need either cases, with missingness or not
      missingwork.l<-NULL
      missingwork.h<-NULL
      
      ## structure value should be 1 for reference and endogenous separately, if not there are missingness. if more there are duplicates.
      
      ## if count of NA is not zero and not number of run (excluding complete missingness across runs)
      
      missing.l<-names(featurestructure.l[featurestructure.l!=ncol(structure.l) & featurestructure.l!=0])
      missing.h<-names(featurestructure.h[featurestructure.h!=ncol(structure.h) & featurestructure.h!=0])
      
      flagmissing.l = length(missing.l)>0
      flagmissing.h = length(missing.h)>0
      
      ## structure value is greater than 1, there are duplicates
      flagduplicate.l = sum(structure.l[!is.na(structure.l)]>1)>0
      flagduplicate.h = sum(structure.h[!is.na(structure.h)]>1)>0
      
      ### if there is missing rows for endogenous
      if ( flagmissing.l | flagmissing.h){
        processout<-rbind(processout,c("CAUTION: the input dataset has incomplete rows. If missing peaks occur they should be included in the dataset as separate rows, and the missing intensity values should be indicated with ’NA’. The incomplete rows are listed below."))
        write.table(processout, file=finalfile,row.names=FALSE)
        
        message("CAUTION : the input dataset has incomplete rows. If missing peaks occur they should be included in the dataset as separate rows, and the missing intensity values should be indicated with ’NA’. The incomplete rows are listed below.")
        
        #### endogenous intensities
        if(flagmissing.l){
          
          ## first, which run has missing	
          runstructure<-apply ( structure.l[which(rownames(structure.l) %in% missing.l),], 2, function ( x ) sum ( is.na ( x ) ) ) > 0
          
          ## get the name of Run
          runID<-names(runstructure[runstructure==TRUE])
          
          
          ## then for each run, which features are missing,
          for(j in 1:length(runID)){
            
            # get subject, group information for this run
            nameID<-unique(work.l[work.l$RUN==runID[j],c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
            
            # get feature ID
            featureID<-structure.l[which(rownames(structure.l) %in% missing.l),colnames(structure.l)==runID[j]]
            
            # get feature ID which has no measuremnt.
            finalfeatureID<-featureID[is.na(featureID)]
            
            # print features ID	 	
            message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some ENDOGENOUS features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
            
            ## save in process file.
            processout<-rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some ENDOGENOUS features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
            write.table(processout, file=finalfile,row.names=FALSE)
            
            ### add missing rows if option is TRUE
            if(fillIncompleteRows){
              
              tempTogetfeature<-work.l[which(work.l$FEATURE %in% names(finalfeatureID)),]
              
              # get PROTEIN and FEATURE infomation
              tempfeatureID<-unique(tempTogetfeature[,c("PROTEIN","PEPTIDE","TRANSITION","FEATURE")])
              
              # merge feature info and run info as 'work' format
              tempmissingwork<-data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, RUN=nameID$RUN, GROUP=nameID$GROUP, SUBJECT=nameID$SUBJECT, SUBJECT_NESTED=nameID$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD)	
              
              # merge with tempary space, missingwork
              missingwork.l<-rbind(missingwork.l,tempmissingwork)
            } ## end fillIncompleteRows options
          } ## end loop for run ID
        } ## end for endogenous
        
        #### reference intensities
        if(flagmissing.h){
          
          ## first, which run has missing	
          runstructure<-apply ( structure.h[which(rownames(structure.h) %in% missing.h),], 2, function ( x ) sum ( is.na ( x ) ) ) > 0
          
          ## get the name of Run
          runID<-names(runstructure[runstructure==TRUE])
          
          ## then for each run, which features are missing,
          for(j in 1:length(runID)){
            
            # get subject, group information for this run
            nameID<-unique(work.h[work.h$RUN==runID[j],c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
            
            # get feature ID
            featureID<-structure.h[which(rownames(structure.h) %in% missing.h),colnames(structure.h)==runID[j]]
            
            # get feature ID which has no measuremnt.
            finalfeatureID<-featureID[is.na(featureID)]
            
            # print features ID	 	
            message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some REFERENCE features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
            
            ## save in process file.
            processout<-rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some REFERENCE features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
            write.table(processout, file=finalfile,row.names=FALSE)
            
            ### add missing rows if option is TRUE
            if(fillIncompleteRows){
              
              tempTogetfeature<-work.h[which(work.h$FEATURE %in% names(finalfeatureID)),]
              
              # get PROTEIN and FEATURE infomation
              tempfeatureID<-unique(tempTogetfeature[,c("PROTEIN","PEPTIDE","TRANSITION","FEATURE")])
              
              # merge feature info and run info as 'work' format
              tempmissingwork<-data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, RUN=nameID$RUN, GROUP=nameID$GROUP, SUBJECT=nameID$SUBJECT, SUBJECT_NESTED=nameID$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD)	
              
              # merge with tempary space, missingwork
              missingwork.h<-rbind(missingwork.h,tempmissingwork)
            } ## end fillIncompleteRows options
          } ## end loop for run ID
        } ## end for endogenous
        
      } ## end for flag missing
      
      ## merge missing rows if fillIncompleteRows=TRUE or message.
      if(fillIncompleteRows){
        
        # merge with work
        ## in future, use rbindlist?? rbindlist(list(work, missingwork))
        work<-rbind(work,missingcomplete.l, missingcomplete.h, missingwork.l,missingwork.h)
        
        ## print message
        message("\n DONE : Incomplete rows for missing peaks are added with intensity values=NA. \n")
        
        ## save in process file.
        processout<-rbind(processout,"Incomplete rows for missing peaks are added with intensity values=NA. - done, Okay")
        write.table(processout, file=finalfile,row.names=FALSE)
        
      }else if(!is.null(missingcomplete.l) | !is.null(missingcomplete.h) | !is.null(missingwork.l) | !is.null(missingwork.l) ){
        
        ## save in process file.
        processout<-rbind(processout,"Please check whether features in the list are generated from spectral processing tool. Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
        write.table(processout, file=finalfile,row.names=FALSE)
        
        stop("Please check whether features in the list are generated from spectral processing tool or not. Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
        
      }
      
      ########### if there are duplicates measurements
      if(flagduplicate.h){
        
        ## first, which run has duplicates
        runstructure<-apply ( structure.h, 2, function ( x ) sum ( x[!is.na(x)] > 1 )>0 )
        
        runID<-names(runstructure[runstructure==TRUE])
        
        ## then for each run, which features have duplicates,
        for(j in 1:length(runID)){
          
          nameID<-unique(work[work$RUN==runID[j],c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
          
          featureID<-structure.h[,colnames(structure.h)==runID[j]]
          finalfeatureID<-featureID[!is.na(featureID) & featureID>1]
          
          message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has multiple rows (duplicate rows) for some REFERENCE features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
          
          ## save in process file.
          processout<-rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has multiple rows (duplicate rows) for some REFERENCE features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
          write.table(processout, file=finalfile,row.names=FALSE)
        }
        
        ## save in process file.
        processout<-rbind(processout,"Please remove duplicate rows in the list above. ")
        write.table(processout, file=finalfile,row.names=FALSE)
        
        stop("Please remove duplicate rows in the list above.\n")		
      }## end flag duplicate for reference
      
      if(flagduplicate.l){
        
        ## first, which run has duplicates
        runstructure<-apply ( structure.l, 2, function ( x ) sum ( x[!is.na(x)] > 1 )>0 )
        
        runID<-names(runstructure[runstructure==TRUE])
        
        ## then for each run, which features have duplicates,
        for(j in 1:length(runID)){
          
          nameID<-unique(work[work$RUN==runID[j],c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
          
          featureID<-structure.l[,colnames(structure.l)==runID[j]]
          finalfeatureID<-featureID[!is.na(featureID) & featureID>1]
          
          message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has multiple rows (duplicate rows) for some ENDOGENOUS features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
          
          ## save in process file.
          processout<-rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has multiple rows (duplicate rows) for some ENDOGENOUS features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
          write.table(processout, file=finalfile,row.names=FALSE)
        }
        
        ## save in process file.
        processout<-rbind(processout,"ERROR : Please remove duplicate rows in the list above. ")
        write.table(processout, file=finalfile,row.names=FALSE)
        
        stop("ERROR : Please remove duplicate rows in the list above.\n")		
      }## end flag duplicate for endogenous
      
      ## no missing and no duplicates
      if(!flagmissing.h & !flagmissing.l & !flagduplicate.h & !flagduplicate.l){
        processout<-rbind(processout,c("Balanced data format with NA for missing feature intensities - okay"))
        write.table(processout, file=finalfile,row.names=FALSE)
      } 		 
    }
    
    ## end 1 method
    
  }else{ ## multiple methods
    
    allflagmissing<-NULL
    allflagduplicate<-NULL
    
    ## check each method
    for(k in 1:length(unique(work$METHOD))){
      
      worktemp<-work[work$METHOD==k,]
      worktemp$RUN<-factor(worktemp$RUN)
      worktemp$FEATURE<-factor(worktemp$FEATURE)
      
      structure = tapply ( worktemp$ABUNDANCE, list ( worktemp$FEATURE, worktemp$RUN ) , function ( x ) length ( x ) ) 
      
      ## structure value should be 2 for labeled, 1 for label-free, if not there are missingness
      if(nlevels(worktemp$LABEL)==2){ ## label-based
        flag = sum(is.na(structure))>0 | sum(structure[!is.na(structure)]<2)>0
      }else{  ## label-free
        flag = sum(is.na(structure))>0
      }
      
      allflagmissing<-c(allflagmissing,flag)
      
      ## for duplicate
      if(nlevels(worktemp$LABEL)==2){ ## label-based
        worktemp.h<-worktemp[worktemp$LABEL=="H",]
        worktemp.l<-worktemp[worktemp$LABEL=="L",]
        
        structure.h = tapply ( worktemp.h$ABUNDANCE, list ( worktemp.h$FEATURE, worktemp.h$RUN ) , function ( x ) length ( x ) ) 
        structure.l = tapply ( worktemp.l$ABUNDANCE, list ( worktemp.l$FEATURE, worktemp.l$RUN ) , function ( x ) length ( x ) ) 
        
        flagduplicate = sum(structure.h[!is.na(structure.h)]>1)>0 | sum(structure.l[!is.na(structure.l)]>1)>0
        
      }else{  ## label-free
        flagduplicate = sum(structure[!is.na(structure)]>1)>0
      }
      
      allflagduplicate<-c(allflagduplicate,flag)
      
    } ## end to check any flag among methods
    
    if ( sum(allflagmissing)!=0 ){
      processout<-rbind(processout,c("CAUTION: the input dataset has incomplete rows. Missing feature intensities should be present in the dataset, and their intensities should be indicated with 'NA'. The incomplete rows are listed below."))
      write.table(processout, file=finalfile,row.names=FALSE)
      
      message("CAUTION : the input dataset has incomplete rows. Missing feature intensities should be present in the dataset, and their intensities should be indicated with 'NA'. The incomplete rows are listed below.")
      
      ## for missign row, need to assign before looping
      missingwork<-NULL
      
      missingcomplete.h<-NULL
      missingcomplete.l<-NULL
      missingwork.h<-NULL
      missingwork.l<-NULL
      
      for(k in 1:length(unique(work$METHOD))){
        
        ## see which method has missing rows
        if(allflagmissing[k]){
          worktemp<-work[work$METHOD==k,]
          worktemp$RUN<-factor(worktemp$RUN)
          worktemp$FEATURE<-factor(worktemp$FEATURE)
          
          if(nlevels(worktemp$LABEL)==1){ ## label-free
            
            structure = tapply ( worktemp$ABUNDANCE, list ( worktemp$FEATURE, worktemp$RUN ) , function ( x ) length ( x ) ) 
            
            ## first, which run has missing	
            runstructure<-apply ( structure, 2, function ( x ) sum ( is.na ( x ) ) ) > 0
            
            ## get the name of Run
            runID<-names(runstructure[runstructure==TRUE])
            
            ## then for each run, which features are missing,
            for(j in 1:length(runID)){
              
              nameID<-unique(worktemp[worktemp$RUN==runID[j],c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
              
              # get feature ID
              featureID<-structure[,colnames(structure)==runID[j]]
              
              # get feature ID which has no measuremnt.
              finalfeatureID<-featureID[is.na(featureID)]
              
              # print features ID	 	
              message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
              
              ## save in process file.
              processout<-rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
              write.table(processout, file=finalfile,row.names=FALSE)
              
              ### add missing rows if option is TRUE
              if(fillIncompleteRows){
                
                tempTogetfeature<-work[which(work$FEATURE %in% names(finalfeatureID)),]
                
                # get PROTEIN and FEATURE infomation
                tempfeatureID<-unique(tempTogetfeature[,c("PROTEIN","PEPTIDE","TRANSITION","FEATURE")])
                
                # merge feature info and run info as 'work' format
                tempmissingwork<-data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, RUN=nameID$RUN, GROUP=nameID$GROUP, SUBJECT=nameID$SUBJECT, SUBJECT_NESTED=nameID$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD)	
                
                # merge with tempary space, missingwork
                missingwork<-rbind(missingwork,tempmissingwork)
              } ## end fillIncompleteRows options
            } ## end loop for run
            
          }else{## end label-free
            
            ########### label-based
            ### count the reference and endobenous separately
            work.l<-worktemp[worktemp$LABEL=="L",]
            work.h<-worktemp[worktemp$LABEL=="H",]
            
            ## get feature by Run count of data
            structure.l = tapply ( work.l$ABUNDANCE, list ( work.l$FEATURE, work.l$RUN ) , function ( x ) length ( x ) ) 
            structure.h = tapply ( work.h$ABUNDANCE, list ( work.h$FEATURE, work.h$RUN ) , function ( x ) length ( x ) ) 
            
            ### 1. reference peptides
            featurestructure.h<-apply(structure.h, 1, function (x) sum(is.na(x)))
            
            # get feature ID of reference which are completely missing across run
            featureID.h<-names(featurestructure.h[featurestructure.h==ncol(structure.h)])
            
            if(length(featureID.h)>0){
              # print message
              message(paste("CAUTION : some REFERENCE features have missing intensities in all the runs. The completely missing REFERENCE features are ", paste(featureID.h, collapse=", "),". Please check whether features in the list are correctly generated from spectral processing tool. \n",sep=""))
              
              ## save in process file.
              processout<-rbind(processout,c(paste("CAUTION : some REFERENCE features have missing intensities in all the runs. The completely missing REFERENCE features are ", paste(featureID.h, collapse=", "),". Please check whether features in the list are correctly generated from spectral processing tool.", sep="")))
              write.table(processout, file=finalfile,row.names=FALSE)
              
              ### add missing rows if option is TRUE
              if(fillIncompleteRows){
                
                # get unique Run information
                nameID<-unique(work.h[,c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
                
                # get PROTEIN and FEATURE information
                # here use whole worktemp dataset
                tempTogetfeature<-worktemp[which(worktemp$FEATURE %in% featureID.h),]
                tempfeatureID<-unique(tempTogetfeature[,c("PROTEIN","PEPTIDE","TRANSITION","FEATURE")])
                
                ## then generate data.frame for missingness,
                for(j in 1:nrow(nameID)){
                  
                  # merge feature info and run info as 'work' format
                  tempmissingwork<-data.frame(tempfeatureID, LABEL="H",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL[j], SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL[j], RUN=nameID$RUN[j], GROUP=nameID$GROUP[j], SUBJECT=nameID$SUBJECT[j], SUBJECT_NESTED=nameID$SUBJECT_NESTED[j], INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD[j])	
                  
                  # merge with tempary space, missingwork
                  missingcomplete.h<-rbind(missingcomplete.h,tempmissingwork)
                }
              }	## end fillIncompleteRows option
              
            } ## end for reference peptides
            
            ### 2. endogenous peptides
            featurestructure.l<-apply(structure.l, 1, function (x) sum(is.na(x)))
            
            # get feature ID of reference which are completely missing across run
            featureID.l<-names(featurestructure.l[featurestructure.l==ncol(structure.l)])
            
            if(length(featureID.l)>0){
              # print message
              message(paste("CAUTION : some ENDOGENOUS features have missing intensities in all the runs. The completely missing ENDOGENOUS features are ", paste(featureID.l, collapse=", "), ". Please check whether features in the list are correctly generated from spectral processing tool. \n",sep=""))
              
              ## save in process file.
              processout<-rbind(processout,c(paste("CAUTION : some ENDOGENOUS features have missing intensities in all the runs. The completely missing ENCOGENOUS features are ", paste(featureID, collapse=", "),". Please check whether features in the list are correctly generated from spectral processing tool. \n",sep="")))
              write.table(processout, file=finalfile,row.names=FALSE)
              
              ### add missing rows if option is TRUE
              if(fillIncompleteRows){
                
                # get unique Run information
                nameID<-unique(work.l[,c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
                
                # get PROTEIN and FEATURE information
                # here use whole worktemp dataset
                tempTogetfeature<-worktemp[which(worktemp$FEATURE %in% featureID.l),]
                tempfeatureID<-unique(tempTogetfeature[,c("PROTEIN","PEPTIDE","TRANSITION","FEATURE")])
                
                ## then generate data.frame for missingness,
                for(j in 1:nrow(nameID)){
                  
                  # merge feature info and run info as 'work' format
                  tempmissingwork<-data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL[j], SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL[j], RUN=nameID$RUN[j], GROUP=nameID$GROUP[j], SUBJECT=nameID$SUBJECT[j], SUBJECT_NESTED=nameID$SUBJECT_NESTED[j], INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD[j])	
                  
                  # merge with tempary space, missingwork
                  missingcomplete.l<-rbind(missingcomplete.l,tempmissingwork)
                }
              }	## end fillIncompleteRows option
            } ## end endogenous peptides
            
            
            ###### second, check other some missingness
            
            ## structure value should be 1 for reference and endogenous separately, if not there are missingness. if more there are duplicates.
            
            ## if count of NA is not zero and not number of run (excluding complete missingness across runs)
            missing.l<-names(featurestructure.l[featurestructure.l!=ncol(structure.l) & featurestructure.l!=0])
            missing.h<-names(featurestructure.h[featurestructure.h!=ncol(structure.h) & featurestructure.h!=0])
            
            flagmissing.l = length(missing.l)>0
            flagmissing.h = length(missing.h)>0
            
            ## structure value is greater than 1, there are duplicates
            flagduplicate.l = sum(structure.l[!is.na(structure.l)]>1)>0
            flagduplicate.h = sum(structure.h[!is.na(structure.h)]>1)>0
            
            ### if there is missing rows for endogenous
            if ( flagmissing.l | flagmissing.h){
              processout<-rbind(processout,c("CAUTION: the input dataset has incomplete rows. If missing peaks occur they should be included in the dataset as separate rows, and the missing intensity values should be indicated with ’NA’. The incomplete rows are listed below."))
              write.table(processout, file=finalfile,row.names=FALSE)
              
              message("CAUTION : the input dataset has incomplete rows. If missing peaks occur they should be included in the dataset as separate rows, and the missing intensity values should be indicated with ’NA’. The incomplete rows are listed below.")
              
              #### endogenous intensities
              if(flagmissing.l){
                
                ## first, which run has missing	
                runstructure<-apply ( structure.l[-which(rownames(structure.l) %in% featureID.l),], 2, function ( x ) sum ( is.na ( x ) ) ) > 0
                
                ## get the name of Run
                runID<-names(runstructure[runstructure==TRUE])
                
                ## then for each run, which features are missing,
                for(j in 1:length(runID)){
                  
                  # get subject, group information for this run
                  nameID<-unique(work.l[work.l$RUN==runID[j],c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
                  
                  # get feature ID
                  featureID<-structure.l[-which(rownames(structure.l) %in% featureID.l),colnames(structure.l)==runID[j]]
                  
                  # get feature ID which has no measuremnt.
                  finalfeatureID<-featureID[is.na(featureID)]
                  
                  # print features ID	 	
                  message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some ENDOGENOUS features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
                  
                  ## save in process file.
                  processout<-rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some ENDOGENOUS features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
                  write.table(processout, file=finalfile,row.names=FALSE)
                  
                  ### add missing rows if option is TRUE
                  if(fillIncompleteRows){
                    
                    tempTogetfeature<-work.l[which(work.l$FEATURE %in% names(finalfeatureID)),]
                    
                    # get PROTEIN and FEATURE infomation
                    tempfeatureID<-unique(tempTogetfeature[,c("PROTEIN","PEPTIDE","TRANSITION","FEATURE")])
                    
                    # merge feature info and run info as 'work' format
                    tempmissingwork<-data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, RUN=nameID$RUN, GROUP=nameID$GROUP, SUBJECT=nameID$SUBJECT, SUBJECT_NESTED=nameID$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD)	
                    
                    # merge with tempary space, missingwork
                    missingwork.l<-rbind(missingwork.l,tempmissingwork)
                  } ## end fillIncompleteRows options
                } ## end loop for run ID
              } ## end for endogenous
              
              #### reference intensities
              if(flagmissing.h){
                
                ## first, which run has missing	
                runstructure<-apply ( structure.h[-which(rownames(structure.h) %in% featureID.h),], 2, function ( x ) sum ( is.na ( x ) ) ) > 0
                
                ## get the name of Run
                runID<-names(runstructure[runstructure==TRUE])
                
                ## then for each run, which features are missing,
                for(j in 1:length(runID)){
                  
                  # get subject, group information for this run
                  nameID<-unique(work.h[work.h$RUN==runID[j],c("SUBJECT_ORIGINAL","GROUP_ORIGINAL","GROUP","SUBJECT","SUBJECT_NESTED","RUN","METHOD")])
                  
                  # get feature ID
                  featureID<-structure.h[-which(rownames(structure.h) %in% featureID.h),colnames(structure.h)==runID[j]]
                  
                  # get feature ID which has no measuremnt.
                  finalfeatureID<-featureID[is.na(featureID)]
                  
                  # print features ID	 	
                  message(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some REFERENCE features (", paste(names(finalfeatureID), collapse=", "),")", sep="" ))
                  
                  ## save in process file.
                  processout<-rbind(processout,c(paste("*** Subject : ", as.character(nameID[,"SUBJECT_ORIGINAL"]) ,", Condition : ", as.character(nameID[,"GROUP_ORIGINAL"]), " has incomplete rows for some REFERENCE features (", paste(names(featureID[is.na(featureID)]), collapse=", "),")", sep="" )))
                  write.table(processout, file=finalfile,row.names=FALSE)
                  
                  ### add missing rows if option is TRUE
                  if(fillIncompleteRows){
                    
                    tempTogetfeature<-work.h[which(work.h$FEATURE %in% names(finalfeatureID)),]
                    
                    # get PROTEIN and FEATURE infomation
                    tempfeatureID<-unique(tempTogetfeature[,c("PROTEIN","PEPTIDE","TRANSITION","FEATURE")])
                    
                    # merge feature info and run info as 'work' format
                    tempmissingwork<-data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL, SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL, RUN=nameID$RUN, GROUP=nameID$GROUP, SUBJECT=nameID$SUBJECT, SUBJECT_NESTED=nameID$SUBJECT_NESTED, INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD)	
                    
                    # merge with tempary space, missingwork
                    missingwork.h<-rbind(missingwork.h,tempmissingwork)
                  } ## end fillIncompleteRows options
                } ## end loop for run ID
              } ## end for endogenous
            } ## end any missingness
            
          } ## end label-based
        } ## if only any flag for method
      } ## end loop for methods
      
      if(fillIncompleteRows){
        
        # merge with work
        ## in future, use rbindlist?? rbindlist(list(work, missingwork))
        if(nlevels(worktemp$LABEL)==1){
          work<-rbind(work,missingwork)
        }else{
          work<-rbind(work,missingcomplete.l, missingcomplete.h, missingwork.l,missingwork.h)
        }
        
        ## print message
        message("\n DONE : Incomplete rows for missing peaks are added with intensity values=NA. \n")
        
        ## save in process file.
        processout<-rbind(processout,"Incomplete rows for missing peaks are added with intensity values=NA. - done, Okay")
        write.table(processout, file=finalfile,row.names=FALSE)
        
      }else if(!is.null(missingcomplete.l) | !is.null(missingcomplete.h) | !is.null(missingwork.l) | !is.null(missingwork.l) | !is.null(missingwork)){
        
        ## save in process file.
        processout<-rbind(processout,"Please check whether features in the list are generated from spectral processing tool. Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
        write.table(processout, file=finalfile,row.names=FALSE)
        
        stop("Please check whether features in the list are generated from spectral processing tool. Or the option, fillIncompleteRows=TRUE, will add incomplete rows for missing peaks with intensity=NA.")
        
      }
      
      
    }else{
      processout<-rbind(processout,c("Balanced data format with NA for missing feature intensities - okay"))
      write.table(processout, file=finalfile,row.names=FALSE)
    }
    
    ### for duplicate, in future
    
  } ## end multiple method
  
  
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
  write.table(processout, file=finalfile,row.names=FALSE)
  
  
  ############################################################
  ## Normalization : 
  
  if(!(normalization=="NONE" | normalization=="FALSE" | normalization=="EQUALIZEMEDIANS" | normalization=="QUANTILE" | normalization=="GLOBALSTANDARDS")){
    processout<-rbind(processout,c(paste("The required input - normalization : 'normalization' value is wrong. - stop")))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    message("'normalization' must be one of \"None\", \"FALSE\", \"equalizeMedians\", \"quantile\", or \"globalStandards\". Here \"equalizeMedian\" will be used. If you want different normalization, please assign 'normalization' again.")
  } 
  
  ############################################################
  ## Normalization : option 0. none
  if(is.element("NONE",normalization) | is.element("FALSE",normalization)){ ## after 'toupper', FALSE becomes character.
    processout<-rbind(processout,c("Normalization : no normalization - okay"))
    write.table(processout, file=finalfile,row.names=FALSE)
  }
  
  ############################################################
  ## Normalization : option 1. constant normalization , equalize medians
  if(!is.element("NONE",normalization) & !is.element("FALSE",normalization) & is.element("EQUALIZEMEDIANS",normalization)){
    
    #	if(!checkMultirun){
    #		if(nlevels(work$LABEL)==1){
    #			### Constant normalization by endogenous
    #			median.run<-tapply(work$ABUNDANCE,work[,"RUN"], function(x) median(x,na.rm=TRUE))
    #			median.all<-median(work$ABUNDANCE, na.rm=TRUE)
    
    #			for (i in 1:length(unique(work$RUN))){
    #				## ABUNDANCE is normalized
    #				work[work$RUN==i,"ABUNDANCE"]<-work[work$RUN==i,"ABUNDANCE"]-median.run[i]+median.all
    #			}
    #		}
    
    #		if(nlevels(work$LABEL)==2){
    ### Constant normalization by heavy standard
    #			h<-work[work$LABEL=="H",]
    #			median.run<-tapply(h$ABUNDANCE,h[,"RUN"], function(x) median(x,na.rm=TRUE))
    #			median.all<-median(h$ABUNDANCE, na.rm=TRUE)
    
    #			for (i in 1:length(unique(work$RUN))){
    #				## ABUNDANCE is normalized
    #				work[work$RUN==i,"ABUNDANCE"]<-work[work$RUN==i,"ABUNDANCE"]-median.run[i]+median.all
    #			}		
    #		}
    #	}else{ ## for multi-method case
    if(nlevels(work$LABEL)==1){
      ### Constant normalization by endogenous per method
      
      median.run<-tapply(work$ABUNDANCE,work$RUN, function(x) median(x,na.rm=TRUE))
      
      median.method<-tapply(work$ABUNDANCE,work$METHOD, function(x) median(x,na.rm=TRUE))
      
      nmethod<-unique(work$METHOD)
      
      for(j in 1:length(nmethod)){
        namerun<-unique(work[work$METHOD==nmethod[j],"RUN"])
        
        for (i in 1:length(namerun)){
          ## ABUNDANCE is normalized
          work[work$RUN==namerun[i],"ABUNDANCE"]<-work[work$RUN==namerun[i],"ABUNDANCE"]-median.run[names(median.run)==namerun[i]]+median.method[j]
        }
      }
    }
    
    
    if(nlevels(work$LABEL)==2){
      
      ### Constant normalization by heavy standard per method
      h<-work[work$LABEL=="H",]
      median.run<-tapply(h$ABUNDANCE,h[,"RUN"], function(x) median(x,na.rm=TRUE))
      median.method<-tapply(h$ABUNDANCE,h$METHOD, function(x) median(x,na.rm=TRUE))
      
      nmethod<-unique(work$METHOD)
      
      for(j in 1:length(nmethod)){
        namerun<-unique(work[work$METHOD==nmethod[j],"RUN"])
        
        for (i in 1:length(namerun)){
          ## ABUNDANCE is normalized
          work[work$RUN==namerun[i],"ABUNDANCE"]<-work[work$RUN==namerun[i],"ABUNDANCE"]-median.run[names(median.run)==namerun[i]]+median.method[j]
        }
      } ## end loop method		
    } # for labe-based 
    
    
    processout<-rbind(processout,c("Normalization : Constant normalization (equalize medians) - okay"))
    write.table(processout, file=finalfile,row.names=FALSE)
  } ## end equaliemedian normalization	
  
  ############################################################
  ## Normalization : option 2. quantile normalization
  if(!is.element("NONE",normalization) & !is.element("FALSE",normalization) & is.element("QUANTILE",normalization)){
    
    
    if(nlevels(work$LABEL)==1){
      
      ## for label-free, just use endogenous
      
      nmethod<-unique(work$METHOD)
      quantileall<-NULL
      
      for(j in 1:length(nmethod)){
        namerun<-unique(work[work$METHOD==nmethod[j],"RUN"])
        
        worktemp<-work[which(work$RUN %in% namerun & !is.na(work$INTENSITY)),]
        worktemp$RUN<-factor(worktemp$RUN)
        worktemp$FEATURE<-factor(worktemp$FEATURE)
        
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
        rm(quantiledone)
        
        #			quantileall<-rbindlist(list(quantileall,quantilelong))
        quantileall<-rbind(quantileall,quantilelong)
        
        rm(quantilelong)
      }
      
      work<-merge(work,quantileall, by=c("FEATURE","RUN"))
      rm(quantileall)
      
      ## reorder
      work<-data.frame("PROTEIN"=work$PROTEIN, "PEPTIDE"=work$PEPTIDE, "TRANSITION"=work$TRANSITION, "FEATURE"=work$FEATURE, "LABEL"=work$LABEL, "GROUP_ORIGINAL"=work$GROUP_ORIGINAL, "SUBJECT_ORIGINAL"=work$SUBJECT_ORIGINAL, "RUN"=work$RUN, "GROUP"=work$GROUP, "SUBJECT"=work$SUBJECT, "SUBJECT_NESTED"=work$SUBJECT_NESTED, "INTENSITY"=work$INTENSITY, "ABUNDANCE"=work$ABUNDANCE_quantile, "METHOD"=work$METHOD)
      
      work<-work[with(work,order(LABEL,GROUP_ORIGINAL,SUBJECT_ORIGINAL,RUN,PROTEIN,PEPTIDE,TRANSITION)),]
      
    }
    
    if(nlevels(work$LABEL)==2){
      
      nmethod<-unique(work$METHOD)
      quantileall<-NULL
      
      for(j in 1:length(nmethod)){
        namerun<-unique(work[work$METHOD==nmethod[j],"RUN"])
        
        ##### for label-based, make quantile normalization for reference
        #worktemp<-work[which(work$RUN %in% namerun & work$LABEL=="H" & !is.na(work$INTENSITY)),] ## because for sparse of reference
        worktemp<-work[which(work$RUN %in% namerun & work$LABEL=="H"),]
        worktemp$RUN<-factor(worktemp$RUN)
        worktemp$FEATURE<-factor(worktemp$FEATURE)
        
        quantiletemp<-as.matrix(xtabs(ABUNDANCE~FEATURE+RUN, data=worktemp))
        rm(worktemp)
        
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
        #worktemp.l<-work[which(work$RUN %in% namerun & work$LABEL=="L" & !is.na(work$INTENSITY)),] ## because for sparse of reference
        worktemp.l<-work[which(work$RUN %in% namerun & work$LABEL=="L"),]
        worktemp.l$RUN<-factor(worktemp.l$RUN)
        worktemp.l$FEATURE<-factor(worktemp.l$FEATURE)
        
        quantiletemp.l<-as.matrix(xtabs(ABUNDANCE~FEATURE+RUN, data=worktemp.l))
        rm(worktemp.l)
        
        ## need to put NA for missing value in endogenous
        quantiletemp.l[quantiletemp.l==0]<-NA
        
        ## apply the difference from reference
        quantiledone.l<-quantiletemp.l-(quantiletemp-quantiledone)
        
        ## get quantiled to long format for apply difference endogenous
        quantilelong.l<-melt(quantiledone.l, id=rownames(quantiledone.l))
        colnames(quantilelong.l)<-c("FEATURE","RUN","ABUNDANCE_quantile")
        quantilelong.l<-data.frame(quantilelong.l, LABEL="L")
        
        rm(quantiletemp)
        rm(quantiledone)
        rm(quantiletemp.l)
        rm(quantiledone.l)
        
        #			quantileall<-rbindlist(list(quantileall,quantilelong.h, quantilelong.l))
        quantileall<-rbind(quantileall,quantilelong.h, quantilelong.l)
        
      }
      
      ##### merge with original data
      work<-merge(work,quantileall, by=c("FEATURE","RUN","LABEL"))
      
      ## reorder
      work<-data.frame("PROTEIN"=work$PROTEIN, "PEPTIDE"=work$PEPTIDE, "TRANSITION"=work$TRANSITION, "FEATURE"=work$FEATURE, "LABEL"=work$LABEL, "GROUP_ORIGINAL"=work$GROUP_ORIGINAL, "SUBJECT_ORIGINAL"=work$SUBJECT_ORIGINAL, "RUN"=work$RUN, "GROUP"=work$GROUP, "SUBJECT"=work$SUBJECT, "SUBJECT_NESTED"=work$SUBJECT_NESTED, "INTENSITY"=work$INTENSITY, "ABUNDANCE"=work$ABUNDANCE_quantile,"METHOD"=work$METHOD)
      
      work<-work[with(work,order(LABEL,GROUP_ORIGINAL,SUBJECT_ORIGINAL,RUN,PROTEIN,PEPTIDE,TRANSITION)),]
      
    }
    
    processout<-rbind(processout,c("Normalization : Quantile normalization - okay"))
    write.table(processout, file=finalfile,row.names=FALSE)
  }
  
  
  
  ############################################################
  ## Normalization : option 3. global standards - for endogenous
  if(!is.element("NONE",normalization) & !is.element("FALSE",normalization) & is.element("GLOBALSTANDARDS",normalization)){
    
    work$RUN<-factor(work$RUN)
    combine<-data.frame(RUN=levels(work$RUN))
    allPeptide<-unique(work$PEPTIDE)
    allProtein<-unique(work$PROTEIN)
    
    for(i in 1:length(nameStandards)){
      
      ## if Peptides
      #namePeptide<-allPeptide[grep(nameStandards[i],allPeptide)] ## cannot grep for modified peptide sequence, [,],+ sign
      namePeptide<-tempPeptide[tempPeptide$PEPTIDESEQUENCE==nameStandards[i],"PEPTIDE"]
      
      if(length(namePeptide)!=0){
        tempStandard<-work[work$PEPTIDE==namePeptide,]
      }else{
        
        ## if Proteins
        nameProtein<-allProtein[grep(nameStandards[i],allProtein)]
        
        if(length(nameProtein)!=0){
          tempStandard<-work[work$PROTEIN==nameProtein,]
        }else{
          processout<-rbind(processout,c(paste("global standard peptides or proteins, ",nameStandards[i] ,", is not in dataset. Please check whether 'nameStandards' input is correct or not.")))
          write.table(processout, file=finalfile,row.names=FALSE)
          
          stop(paste("global standard peptides or proteins, ",nameStandards[i] ,", is not in dataset. Please check whether 'nameStandards' input is correct or not."))
        }	
        
      }
      
      ### here, by RUN, but need to check !!!
      tempStandard<-tempStandard[tempStandard$GROUP!="0",]
      tempStandard$RUN<-factor(tempStandard$RUN)
      
      tempStandard<-tempStandard[!is.na(tempStandard$ABUNDANCE),]
      meanStandard<-tapply(tempStandard$ABUNDANCE, tempStandard$RUN, function(x) mean(x, na.rm=TRUE))
      
      meanStandard<-data.frame(RUN=names(meanStandard),meanStandard)
      combine<-merge(combine, meanStandard, by="RUN", all=TRUE)
      colnames(combine)[i+1]<-paste("meanStandard",i,sep="")
    }
    
    rownames(combine)<-combine$RUN
    combine<-subset(combine, select=-c(RUN))
    
    ## get mean among global standards
    allmean<-apply(combine,1,mean)
    #allmean[is.na(allmean)]<-0
    
    allmeantemp<-data.frame(RUN=names(allmean),allmean)
    allrun<-unique(work[,c("RUN","METHOD")])
    
    allmeantemp<-merge(allmeantemp, allrun,by="RUN")
    median.all<-tapply(allmeantemp$allmean, allmeantemp$METHOD, function(x) median(x,na.rm=TRUE))
    
    #### adjust
    nmethod<-unique(work$METHOD)
    
    for(j in 1:length(nmethod)){
      namerun<-unique(work[work$METHOD==nmethod[j],"RUN"])
      
      for (i in 1:length(namerun)){
        ## ABUNDANCE is normalized			
        if(!is.na(allmean[names(allmean)==namerun[i]])) work[work$RUN==namerun[i] & work$LABEL=="L","ABUNDANCE"]<-work[work$RUN==namerun[i] & work$LABEL=="L","ABUNDANCE"]-allmean[names(allmean)==namerun[i]]+median.all[j]
      }
    } ## end loop method
    
    
    processout<-rbind(processout,c("Normalization : normalization with global standards protein - okay"))
    write.table(processout, file=finalfile,row.names=FALSE)
    
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
    write.table(processout, file=finalfile,row.names=FALSE)
    
  }else{
    processout<-rbind(processout,c("Between Run Interference Score is not calculated."))
    write.table(processout, file=finalfile,row.names=FALSE)
  }
  
  
  
  
  ##############################
  ## featureSubset
  ##  !! need to decide how to present : keep original all data and make new column to mark, or just present selected subset    
  
  if(featureSubset=="all"){
 	 message("* Use all features that the dataset origianally has.")
 	 
 	 processout<-rbind(processout,c("* Use all features that the dataset origianally has."))
     write.table(processout, file=finalfile, row.names=FALSE)
  } 
  
  
  if(featureSubset=="top3"){
  		message("* Use top3 features that have highest average of log2(intensity) across runs.")
  		
  		processout<-rbind(processout,c("* Use top3 features that have highest average of log2(intensity) across runs."))
        write.table(processout, file=finalfile, row.names=FALSE)

  	 
  	 	## INTENSITY vs ABUNDANCE?
  		## how to decide top3 for DIA?
	
		### 2015.09.05
		temp1<-aggregate(INTENSITY~PROTEIN+FEATURE,data=work, function(x) mean(x, na.rm=TRUE))

		temp2<-split(temp1, temp1$PROTEIN)

		temp3<-lapply(temp2, function(x){ 
			x<-x[order(x$INTENSITY,decreasing=T),]
			x<-x$FEATURE[1:3]
			})
	
		selectfeature<-unlist(temp3,use.names=FALSE)
		selectfeature<-selectfeature[!is.na(selectfeature)]
		
		# get subset
		work<-work[which(work$FEATURE %in% selectfeature),]	

	 #### old version with looping
  	 #selectfeature<-NULL
		
	#	for(i in 1:length(unique(work$PROTEIN))){
	
	#		sub<-work[work$PROTEIN==unique(work$PROTEIN)[i],]
	#		sub$FEATURE<-factor(sub$FEATURE)
			
	#		message(paste("Getting top3 featuares per protein for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(work$PROTEIN)),")"))
			
	#		subtemp<-sub[!is.na(sub$INTENSITY),]
	#		if(nrow(subtemp)==0) next
			
	#		## from raw scale 
	#		maxvalue<-aggregate(INTENSITY~FEATURE,data=sub, function(x) mean(x, na.rm=TRUE))
	#		maxvalue<-maxvalue[order(maxvalue$INTENSITY, decreasing=T),]

	#		## choose top n
	#		maxfeature<-maxvalue$FEATURE[1:3]
	
	#		selectfeature<-c(selectfeature, as.character(maxfeature))		
	#	}
		
	#	selectfeature<-selectfeature[!is.na(selectfeature)]
	#	work<-work[which(work$FEATURE %in% selectfeature),]		
  }
  
  
  if(featureSubset=="highQuality"){
  	  message("* Use feature selection algorithm in order to get high quality features.")
  	  
  	  processout<-rbind(processout,c("* Use feature selection algorithm in order to get high quality features."))
      write.table(processout, file=finalfile, row.names=FALSE)
  	 
      tempfeature<-try(.FeatureSelection1(work,lambda,eta, address),silent=TRUE)
    
      if(class(tempfeature)=="try-error") {
        message("*** error : can't select the best features. Now use all features.")
      
        processout<-rbind(processout,c(paste("error : can't select the best features. Now use all features.")))
        write.table(processout, file=finalfile, row.names=FALSE)
      
        work<-work
      
      }else{
        work<-tempfeature
      }
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
  write.table(processout, file=finalfile,row.names=FALSE)
  
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
  
  write.table(processout, file=finalfile,row.names=FALSE)
  
  #### protein list with 1 feature
  temp<-unique(work[,c("PROTEIN","FEATURE")])
  temp1<-xtabs(~PROTEIN,data=temp)
  temp2<-as.data.frame(temp1[temp1==1])
  if(nrow(temp2)>0) message("\n","** Protein (",paste(rownames(temp2),collapse = ", "),") has only single transition : Consider excluding this protein from the dataset.", "\n")
  
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
  c.tech<-round(summary.s[1,]/(summary.s[2,]* length(unique(work$METHOD))))
  #summary.s[3,]<-ifelse(c.tech==1,0,c.tech)
  summary.s[3,]<-c.tech
  
  colnames(summary.s)<-unique(work$GROUP_ORIGINAL)
  rownames(summary.s)<-c("# of MS runs","# of Biological Replicates", "# of Technical Replicates")
  
  print(summary.s)
  
  
  
  message("\n Summary of Missingness :\n" )
  message("  # transitions are completely missing in one condition: ", sum(final.decision!=0), "\n")
  if(sum(final.decision!=0)!=0) message("    -> ", paste(names(final.decision[final.decision!=0]),collapse = ", "))
  
  without<-xtabs(~RUN,work)
  withall<-xtabs(~RUN,all.work)
  run.missing<-without/withall
  message("\n  # run with 75% missing observations: ", sum(run.missing<0.25), "\n")
  if(sum(run.missing<0.25)!=0) message("    -> ", paste("RUN",names(without[run.missing<0.25]),sep=" "))
  
  ## output process
  processout<-rbind(processout,c("Summary of Missingness :"))
  processout<-rbind(processout,c(paste("  # transitions are completely missing in one condition: ", sum(final.decision!=0), sep="")))
  if(sum(final.decision!=0)!=0) processout<-rbind(processout,"    -> ", paste(names(final.decision[final.decision!=0]),collapse = ", "))
  
  processout<-rbind(processout,c(paste("  # run with 75% missing observations: ", sum(run.missing<0.25), sep="")))
  if(sum(run.missing<0.25)!=0) processout<-rbind(processout,"    -> ", paste("RUN",names(without[run.missing<0.25]),sep=" "))
  
  write.table(processout, file=finalfile,row.names=FALSE)
  
  
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
  write.table(processout, file=finalfile,row.names=FALSE)
  
   ##### after normalization, zero intensity could be negative
    work[!is.na(work$ABUNDANCE) & work$ABUNDANCE<0,"ABUNDANCE"]<-0
    
    work[!is.na(work$INTENSITY) & work$INTENSITY==0,"ABUNDANCE"]<-0

	
	###==================================================
 	### get the summarization per subplot (per RUN)   

	message("\n == Start the summarization per subplot...")

	rqresult<-try(.runQuantification(work,summaryMethod,equalFeatureVar,filterLogOfSum,cutoffCensored,censoredInt,remove50missing,MBimpute),silent=TRUE)

	if(class(rqresult)=="try-error") {
      message("*** error : can't summarize per subplot with ", summary, ".")
     
      processout<-rbind(processout,c(paste("error : can't summarize per subplot with ", summary, ".", sep = "")))
      write.table(processout, file=finalfile, row.names=FALSE)
	
	  rqall<-NULL
	  rqmodelqc<-NULL
	  workpred<-NULL
	  
	 }else{
      
		label<-nlevels(work$LABEL)==2
	
		if(sum(is.element(colnames(rqresult$rqdata),"RUN"))==0){
			## logsum is summarization per subject
			lab<-unique(work[,c("GROUP","GROUP_ORIGINAL","SUBJECT_ORIGINAL","SUBJECT_NESTED","SUBJECT")])
	
			if(label) lab<-lab[lab$GROUP!=0,]

			rqall<-merge(rqresult$rqdata, lab, by="SUBJECT_ORIGINAL")
			
		}else{
			lab<-unique(work[,c("RUN","GROUP","GROUP_ORIGINAL","SUBJECT_ORIGINAL","SUBJECT_NESTED","SUBJECT")])
	
			if(label) lab<-lab[lab$GROUP!=0,]

			rqall<-merge(rqresult$rqdata, lab, by="RUN")
		}
		
		rqall$GROUP<-factor(rqall$GROUP)
		rqall$Protein<-factor(rqall$Protein)
		
		rqmodelqc<-rqresult$ModelQC
		
		workpred<-rqresult$PredictedBySurvival
		
		message("\n == the summarization per subplot is done.")
		
		processout<-rbind(processout,c(paste("the summarization per subplot is done.- okay : ",summaryMethod, sep="")))
		write.table(processout, file=finalfile, row.names=FALSE)

	 }
	 
  ## return work data.frame	and run quantification
	
	processedquant<-list(ProcessedData=work, RunlevelData=rqall, SummaryMethod=summaryMethod, ModelQC=rqmodelqc, PredictBySurvival=workpred)
    return(processedquant)
  
}




#############################################
#############################################
# Part 2 dataProcessPlots
#############################################
#############################################


dataProcessPlots<-function(data=data,type=type,featureName="Transition",ylimUp=FALSE,ylimDown=FALSE,scale=FALSE,interval="CI",x.axis.size=10,y.axis.size=10,text.size=4,text.angle=0,legend.size=7,dot.size.profile=2,dot.size.condition=3,width=10,height=10, which.Protein="all", originalPlot=TRUE, summaryPlot=TRUE, address=""){
	
  datafeature<-data$ProcessedData
  datarun<-data$RunlevelData
  
  datafeature$PROTEIN<-factor(datafeature$PROTEIN)	
  datarun$Protein<-factor(datarun$Protein)	
  
  if(!is.element("SUBJECT_NESTED",colnames(datafeature))){
    stop("Input for dataProcessPlots function should be processed by dataProcess function previously. Please use 'dataProcess' function first.")
  }
  
  if(length(setdiff(toupper(type),c(toupper("ProfilePlot"),toupper("QCPlot"),toupper("ConditionPlot"))))!=0){
    stop(paste("Input for type=",type,". However,'type' should be one of \"ProfilePlot\", \"QCPlot\",\"ConditionPlot\".",sep=""))
  }
  

  #################
  ## Profile plot
  #################
  if (toupper(type)=="PROFILEPLOT"){
    
    #### choose Proteins or not
    if(which.Protein!="all"){
      ## check which.Protein is name of Protein
      if(is.character(which.Protein)){
        temp.name<-which.Protein
        
        ## message if name of Protein is wrong.
        if(length(setdiff(temp.name,unique(datafeature$PROTEIN)))>0)
          stop(paste("Please check protein name. Data set does not have this protein. -", paste(temp.name, collapse=", "),sep=" "))
      }
      
      ## check which.Protein is order number of Protein
      if(is.numeric(which.Protein)){
        temp.name<-levels(datafeature$PROTEIN)[which.Protein]
        
        ## message if name of Protein is wrong.
        if(length(levels(datafeature$PROTEIN))<max(which.Protein))
          stop(paste("Please check your selection of proteins. There are ", length(levels(datafeature$PROTEIN))," proteins in this dataset.",sep=" "))
      }
      
      ## use only assigned proteins
      datafeature<-datafeature[which(datafeature$PROTEIN %in% temp.name),]
      datafeature$PROTEIN<-factor(datafeature$PROTEIN)
      
      datarun<-datarun[which(datarun$PROTEIN %in% temp.name),]
      datarun$PROTEIN<-factor(datarun$PROTEIN)
    }
    
    # assign upper or lower limit
    ## ylimUp
    y.limup<-30
    if(is.numeric(ylimUp)) y.limup<-ylimUp 
    
    ## ylimDown
    y.limdown=-1
    if(is.numeric(ylimDown)) y.limdown<-ylimDown 
    
    datafeature<-datafeature[with(datafeature,order(GROUP_ORIGINAL,SUBJECT_ORIGINAL,LABEL)),]
    datafeature$RUN<-factor(datafeature$RUN,levels=unique(datafeature$RUN),labels=seq(1,length(unique(datafeature$RUN))))
    datafeature$RUN<-as.numeric(datafeature$RUN)
    tempGroupName<-unique(datafeature[,c("GROUP_ORIGINAL","RUN")])
    
    groupAxis<-as.numeric(xtabs(~GROUP_ORIGINAL,tempGroupName))
    cumGroupAxis<-cumsum(groupAxis)
    lineNameAxis<-cumGroupAxis[-nlevels(datafeature$GROUP_ORIGINAL)]
    
    groupName<-data.frame(RUN=c(0,lineNameAxis)+groupAxis/2+0.5,y=rep(y.limup-1,length(groupAxis)),Name=levels(datafeature$GROUP_ORIGINAL))
    
    if(length(unique(datafeature$LABEL))==2){
      datafeature$LABEL<-factor(datafeature$LABEL,labels=c("Reference","Endogenous"))	
    }else{
      if(unique(datafeature$LABEL)=="L"){
        datafeature$LABEL<-factor(datafeature$LABEL,labels=c("Endogenous"))	
      }
      if(unique(datafeature$LABEL)=="H"){
        datafeature$LABEL<-factor(datafeature$LABEL,labels=c("Reference"))
      }
    }
    
    ## need to fill in incomplete rows for Runlevel data
    haverun<-FALSE
    
    if(sum(is.element(colnames(datarun),"RUN"))!=0){
    		datamat = dcast( Protein ~ RUN, data=datarun, value.var='LogIntensities', keep=TRUE) 
    
    		datarun = melt(datamat, id.vars=c('Protein'))
		colnames(datarun)[colnames(datarun) %in% c("variable","value")]<-c('RUN','ABUNDANCE')
		
		haverun<-TRUE
    }
    
    ## remove the column called 'SuggestToFilter' if there.
    if(any(is.element(colnames(datafeature),"SuggestToFilter"))) datafeature$SuggestToFilter<-NULL
    
    
    #### save the plots as pdf or not
    # If there are the file with the same name, add next numbering at the end of file name		
    if(originalPlot){
    if(address!=FALSE){
      allfiles<-list.files()
      
      num<-0
      filenaming<-paste(address,"ProfilePlot",sep="")
      finalfile<-paste(address,"ProfilePlot.pdf",sep="")
      
      while(is.element(finalfile,allfiles)){
        num<-num+1
        finalfile<-paste(paste(filenaming,num,sep="-"),".pdf",sep="")
      }	
      
      pdf(finalfile, width=width, height=height)
    }
    
    for (i in 1:nlevels(datafeature$PROTEIN)){	
      sub<-datafeature[datafeature$PROTEIN==levels(datafeature$PROTEIN)[i],]
      sub$FEATURE<-factor(as.character(sub$FEATURE))	
      sub$SUBJECT<-factor(sub$SUBJECT)	
      sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
      sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)
      sub$PEPTIDE<-factor(as.character(sub$PEPTIDE))
      
      
      #sub<- sub[with(sub, order(LABEL,RUN,FEATURE)), ]
      
      ## if all measurements are NA,
      if(nrow(sub)==sum(is.na(sub$ABUNDANCE))){
        message(paste("Can't the Profile plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),") because all measurements are NAs."))
        next()
      }
      
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
      
      ## if there is 'imputed' column,
#      if(any(is.element(colnames(sub),"imputed"))){
      	## need to keep 'sub' to present imputed values in the second plot
#      	subtemp<-sub
#      	subtemp[subtemp$imputed==TRUE,"ABUNDANCE"]<-NA
#      }else{
#      	subtemp<-sub
#      }

   
      #options(show.error.messages = FALSE)
      
      if(toupper(featureName)=="TRANSITION"){
        
        ## 1st plot for original plot
        ptemp<-ggplot(aes_string(x='RUN', y='ABUNDANCE', color='FEATURE',linetype='FEATURE'), data=sub)+facet_grid(~LABEL)+geom_point(size=dot.size.profile)+geom_line(size=0.5)+scale_colour_manual(values=s)+scale_linetype_manual(values=ss,guide="none")+scale_x_continuous('MS runs',breaks=cumGroupAxis)+annotate("text",x=groupName$RUN,y=groupName$y,label=groupName$Name,size=text.size,angle=text.angle)+geom_vline(xintercept=lineNameAxis+0.5,colour="grey",linetype="longdash")+labs(title=unique(sub$PROTEIN))+theme(
          panel.background=element_rect(fill='white', colour="black"),
          legend.key=element_rect(fill='white',colour='white'),
          panel.grid.minor = element_blank(),
          strip.background=element_rect(fill='gray95'),
          strip.text.x=element_text(colour=c("#00B0F6"),size=14),
          axis.text.x=element_text(size=x.axis.size,colour="black"),
          axis.text.y=element_text(size=y.axis.size,colour="black"),
          axis.ticks=element_line(colour="black"),
          axis.title.x=element_text(size=x.axis.size+5,vjust=-0.4),
          axis.title.y=element_text(size=y.axis.size+5,vjust=0.3),
          title=element_text(size=x.axis.size+8,vjust=1.5),
          legend.position="top",
          legend.text=element_text(size=legend.size))
        
        ## y-axis labeling
        temp<-sub[!is.na(sub[,"ABUNDANCE"]) & !is.na(sub[,"INTENSITY"]),]
        temptest<-abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])
        
        if(temptest){
          ptemp<-ptemp+scale_y_continuous('Log2-intensities',limit=c(y.limdown, y.limup))
        }else{
          ptemp<-ptemp+scale_y_continuous('Log10-intensities',limit=c(y.limdown, y.limup))
        }	
        
        ## draw point here,
        ## if there are imputed data
        if(any(is.element(colnames(sub),"imputed"))){
			ptemp<-ptemp+geom_point(data=sub,aes(x=RUN, y=ABUNDANCE,shape=imputed, size=imputed))+scale_shape_manual(values=c(20,4),labels=c("Detected data","Imputed data"))+scale_size_manual(values=c(1.5,2), guide="none")+guides(color=guide_legend(title=paste("# peptide:",nlevels(sub$PEPTIDE)),ncol=3), shape=guide_legend(title=NULL))
		}else{
			ptemp<-ptemp+geom_point(size=1.5)+guides(color=guide_legend(title=paste("# peptide:",nlevels(sub$PEPTIDE)),ncol=3))
		}
		
		print(ptemp)
       
        message(paste("Drew the Profile plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),")"))
        
      }
      
      if(toupper(featureName)=="PEPTIDE"){
        
        ptemp<-ggplot(aes_string(x='RUN', y='ABUNDANCE', color='PEPTIDE',linetype='FEATURE'), data=sub)+facet_grid(~LABEL)+geom_point(size=dot.size.profile)+geom_line(size=0.5)+scale_colour_manual(values=unique(s))+scale_linetype_manual(values=ss,guide="none")+scale_x_continuous('MS runs',breaks=cumGroupAxis)+annotate("text",x=groupName$RUN,y=groupName$y,label=groupName$Name,size=text.size,angle=text.angle)+geom_vline(xintercept=lineNameAxis+0.5,colour="grey",linetype="longdash")+labs(title=unique(sub$PROTEIN))+theme(
          panel.background=element_rect(fill='white', colour="black"),
          legend.key=element_rect(fill='white',colour='white'),
          panel.grid.minor = element_blank(),
          strip.background=element_rect(fill='gray95'),	
          strip.text.x=element_text(colour=c("#00B0F6"),size=14),
          axis.text.x=element_text(size=x.axis.size,colour="black"),
          axis.text.y=element_text(size=y.axis.size,colour="black"),
          axis.ticks=element_line(colour="black"),
          axis.title.x=element_text(size=x.axis.size+5,vjust=-0.4),
          axis.title.y=element_text(size=y.axis.size+5,vjust=0.3),
          title=element_text(size=x.axis.size+8,vjust=1.5),
          legend.position="top",
          legend.text=element_text(size=legend.size))
        
        ## y-axis labeling
        temp<-sub[!is.na(sub[,"ABUNDANCE"]) & !is.na(sub[,"INTENSITY"]),]
        temptest<-abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])
        
        if(temptest){
          ptemp<-ptemp+scale_y_continuous('Log2-intensities',limit=c(y.limdown, y.limup))	
        }else{
          ptemp<-ptemp+scale_y_continuous('Log10-intensities',limit=c(y.limdown, y.limup))
        }
        
        ## draw point here,
        ## if there are imputed data
        if(any(is.element(colnames(sub),"imputed"))){
			ptemp<-ptemp+geom_point(data=sub,aes(x=RUN, y=ABUNDANCE,shape=imputed, size=imputed))+scale_shape_manual(values=c(20,4),labels=c("Detected data","Imputed data"))+scale_size_manual(values=c(1.5,2), guide="none")+guides(color=guide_legend(title=paste("# peptide:",nlevels(sub$PEPTIDE)),ncol=3), shape=guide_legend(title=NULL))
		}else{
			ptemp<-ptemp+guides(color=guide_legend(title=paste("# peptide:",nlevels(sub$PEPTIDE)),ncol=3))
		}
        
		print(ptemp)

        message(paste("Drew the Profile plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),")"))
        
      }
      
      if(toupper(featureName)=="NA"){
        
        ptemp<-ggplot(aes_string(x='RUN', y='ABUNDANCE', color='PEPTIDE',linetype='FEATURE'), data=sub)+facet_grid(~LABEL)+geom_point(size=dot.size.profile)+geom_line(size=0.5)+scale_colour_manual(values=unique(s),guide="none")+scale_linetype_manual(values=ss,guide="none")+scale_x_continuous('MS runs',breaks=cumGroupAxis)+annotate("text",x=groupName$RUN,y=groupName$y,label=groupName$Name,size=text.size,angle=text.angle)+geom_vline(xintercept=lineNameAxis+0.5,colour="grey",linetype="longdash")+labs(title=unique(sub$PROTEIN))+theme(
          panel.background=element_rect(fill='white', colour="black"),
          legend.key=element_rect(fill='white',colour='white'),
          panel.grid.minor = element_blank(),
          strip.background=element_rect(fill='gray95'),	
          strip.text.x=element_text(colour=c("#00B0F6"),size=14),
          axis.text.x=element_text(size=x.axis.size,colour="black"),
          axis.text.y=element_text(size=y.axis.size,colour="black"),
          axis.ticks=element_line(colour="black"),
          axis.title.x=element_text(size=x.axis.size+5,vjust=-0.4),
          axis.title.y=element_text(size=y.axis.size+5,vjust=0.3),
          title=element_text(size=x.axis.size+8,vjust=1.5),
          legend.position="top",
          legend.text=element_text(size=legend.size))
        
        ## y-axis labeling
        temp<-sub[!is.na(sub[,"ABUNDANCE"]) & !is.na(sub[,"INTENSITY"]),]
        temptest<-abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])
        
        if(temptest){
          ptemp<-ptemp+scale_y_continuous('Log2-intensities',limit=c(y.limdown, y.limup))
        }else{
          ptemp<-ptemp+scale_y_continuous('Log10-intensities',limit=c(y.limdown, y.limup))
        }	
        
        ## draw point here,
        ## if there are imputed data
        if(any(is.element(colnames(sub),"imputed"))){
			ptemp<-ptemp+geom_point(data=sub,aes(x=RUN, y=ABUNDANCE,shape=imputed, size=imputed))+scale_shape_manual(values=c(20,4),labels=c("Detected data","Imputed data"))+scale_size_manual(values=c(1.5,2), guide="none")+guides(shape=guide_legend(title=NULL))
		}else{
			ptemp<-ptemp
		}
        
       	print(ptemp)
        	
        message(paste("Drew the Profile plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),")"))
        
      }
    } ## end-loop for each protein
    if(address!=FALSE) dev.off() 
    } # end original plot
    
    ###############
    ## 2st plot for original plot : summary
    
    if(summaryPlot){
		if(address!=FALSE){
     	 	allfiles<-list.files()
      
     	 	num<-0
      		filenaming<-paste(address,"ProfilePlot_wSummarization",sep="")
      		finalfile<-paste(address,"ProfilePlot_wSummarization.pdf",sep="")
      
      		while(is.element(finalfile,allfiles)){
        		num<-num+1
        		finalfile<-paste(paste(filenaming,num,sep="-"),".pdf",sep="")
      		}	
      
     	 	pdf(finalfile, width=width, height=height)
    	}

    	for (i in 1:nlevels(datafeature$PROTEIN)){	
     		sub<-datafeature[datafeature$PROTEIN==levels(datafeature$PROTEIN)[i],]
      		sub$FEATURE<-factor(as.character(sub$FEATURE))	
      		sub$SUBJECT<-factor(sub$SUBJECT)	
      		sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
      		sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)
      		sub$PEPTIDE<-factor(as.character(sub$PEPTIDE))
      
      		#sub<- sub[with(sub, order(LABEL,RUN,FEATURE)), ]
      
      		## if all measurements are NA,
     		 if(nrow(sub)==sum(is.na(sub$ABUNDANCE))){
       			 message(paste("Can't the Profile plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),") because all measurements are NAs."))
        		next()
     		 }
      
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
      
      		if(haverun){
        		subrun<-datarun[datarun$Protein==levels(datafeature$PROTEIN)[i],]
			
				if(nrow(subrun)!=0){
					quantrun<-data.frame(PROTEIN=subrun$Protein, PEPTIDE="Run summary",TRANSITION="Run summary",FEATURE="Run summary",LABEL="Endogenous",GROUP_ORIGINAL=NA,SUBJECT_ORIGINAL=NA,RUN=subrun$RUN,GROUP=NA,SUBJECT=NA,SUBJECT_NESTED=NA, INTENSITY=NA, ABUNDANCE=subrun$ABUNDANCE,METHOD=1)
				}else{ ## if there is only one Run measured across all runs, no Run information for linear with censored
					quantrun<-data.frame(PROTEIN=levels(datafeature$PROTEIN)[i], PEPTIDE="Run summary",TRANSITION="Run summary",FEATURE="Run summary",LABEL="Endogenous",GROUP_ORIGINAL=NA,SUBJECT_ORIGINAL=NA,RUN=unique(datafeature$RUN)[1],GROUP=NA,SUBJECT=NA,SUBJECT_NESTED=NA, INTENSITY=NA, ABUNDANCE=NA,METHOD=1)
				}
			

				if(any(is.element(colnames(sub),"imputed"))){
					quantrun$imputed<-FALSE
				}

				quantrun$analysis<-"Run summary"
				sub$analysis<-"Processed feature-level data"
		
				final<-rbind(sub, quantrun)
				final$analysis<-factor(final$analysis)
				final$FEATURE<-factor(final$FEATURE)
				final$RUN<-as.numeric(final$RUN)

				ptempall<-ggplot(aes_string(x='RUN', y='ABUNDANCE', color='analysis',linetype='FEATURE',size='analysis',shape='analysis'), data=final)+facet_grid(~LABEL)+geom_point(size=dot.size.profile)+geom_line(size=0.5)+scale_colour_manual(values=c("lightgray","darkred"))+scale_shape_manual(values=c(20,20))+scale_size_manual(values=c(1,2))+scale_linetype_manual(values=c(rep(1,times=length(unique(final$FEATURE))-1),2),guide="none")+scale_x_continuous('MS runs',breaks=cumGroupAxis)+annotate("text",x=groupName$RUN,y=groupName$y,label=groupName$Name,size=text.size,angle=text.angle)+geom_vline(xintercept=lineNameAxis+0.5,colour="grey",linetype="longdash")+labs(title=unique(final$PROTEIN))+theme(
					panel.background=element_rect(fill='white', colour="black"),
					legend.key=element_rect(fill='white',colour='white'),
					panel.grid.minor = element_blank(),
					strip.background=element_rect(fill='gray95'),
					strip.text.x=element_text(colour=c("#00B0F6"),size=14),
					axis.text.x=element_text(size=x.axis.size,colour="black"),
					axis.text.y=element_text(size=y.axis.size,colour="black"),
					axis.ticks=element_line(colour="black"),
					axis.title.x=element_text(size=x.axis.size+5,vjust=-0.4),
					axis.title.y=element_text(size=y.axis.size+5,vjust=0.3),
					title=element_text(size=x.axis.size+8,vjust=1.5),
					legend.position="top",
					legend.text=element_text(size=legend.size),
					legend.title=element_blank())

					ptempall<-ptempall+scale_y_continuous('Log2-intensities',limit=c(y.limdown, y.limup))
				
					## draw point again because some red summary dots could be hiden
					ptempall<-ptempall+geom_point(data=final,aes(x=RUN, y=ABUNDANCE,shape=analysis, size=analysis, color=analysis))
			
				print(ptempall)
			
				message(paste("Drew the Profile plot with summarization for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),")"))
			
			}

    	} ## end-loop for each protein
    	if(address!=FALSE) dev.off() 
  	}  
  } ## end Profile plot	
  
  
  #################### 
  ## QC plot (Quality control plot)
  #################### 
  if (toupper(type)=="QCPLOT"){
    
    #### save the plots as pdf or not
    # If there are the file with the same name, add next numbering at the end of file name		
    if(address!=FALSE){
      allfiles<-list.files()
      
      num<-0
      filenaming<-paste(address,"QCPlot",sep="")
      finalfile<-paste(address,"QCPlot.pdf",sep="")
      
      while(is.element(finalfile,allfiles)){
        num<-num+1
        finalfile<-paste(paste(filenaming,num,sep="-"),".pdf",sep="")
      }	
      
      pdf(finalfile, width=width, height=height)
    }
    
    ## options(warn=-1)
    
    # assign upper or lower limit
    ## ylimUp
    y.limup<-30
    if(is.numeric(ylimUp)) y.limup<-ylimUp 
    
    ## ylimDown
    y.limdown=-1
    if(is.numeric(ylimDown)) y.limdown<-ylimDown 
    
    # relabel the Run (make it sorted by group first)
    datafeature<-datafeature[with(datafeature,order(GROUP_ORIGINAL,SUBJECT_ORIGINAL)),]
    datafeature$RUN<-factor(datafeature$RUN,levels=unique(datafeature$RUN),labels=seq(1,length(unique(datafeature$RUN))))
    
    if(length(unique(datafeature$LABEL))==2){
      datafeature$LABEL<-factor(datafeature$LABEL,labels=c("Reference","Endogenous"))	
      label.color<-c("darkseagreen1","lightblue")
    }else{
      if(unique(datafeature$LABEL)=="L"){
        datafeature$LABEL<-factor(datafeature$LABEL,labels=c("Endogenous"))
        label.color<-c("lightblue")	
      }
      if(unique(datafeature$LABEL)=="H"){
        datafeature$LABEL<-factor(datafeature$LABEL,labels=c("Reference"))
        label.color<-c("darkseagreen1")
      }
    }
    
    tempGroupName<-unique(datafeature[,c("GROUP_ORIGINAL","RUN")])
    datafeature<-datafeature[with(datafeature,order(LABEL,GROUP_ORIGINAL,SUBJECT_ORIGINAL)),]
    
    groupAxis<-as.numeric(xtabs(~GROUP_ORIGINAL,tempGroupName))
    cumGroupAxis<-cumsum(groupAxis)
    lineNameAxis<-cumGroupAxis[-nlevels(datafeature$GROUP_ORIGINAL)]
    
    groupName<-data.frame(RUN=c(0,lineNameAxis)+groupAxis/2+0.5,y=rep(y.limup-1,length(groupAxis)),Name=levels(datafeature$GROUP_ORIGINAL))
    
    #### all protein
    ptemp<-ggplot(aes_string(x='RUN', y='ABUNDANCE'), data=datafeature)+facet_grid(~LABEL)+geom_boxplot(aes_string(fill='LABEL'),outlier.shape=1,outlier.size=1.5)+scale_fill_manual(values=label.color, guide="none")+scale_x_discrete('MS runs',breaks=cumGroupAxis)+annotate("text",x=groupName$RUN,y=groupName$y,label=groupName$Name,size=text.size,angle=text.angle)+geom_vline(xintercept=lineNameAxis+0.5,colour="grey",linetype="longdash")+labs(title="All proteins")+theme(
      panel.background=element_rect(fill='white', colour="black"),
      legend.key=element_rect(fill='white',colour='white'),
      panel.grid.minor = element_blank(),
      strip.background=element_rect(fill='gray95'),	
      strip.text.x=element_text(colour=c("#00B0F6"),size=14),
      axis.text.x=element_text(size=x.axis.size,colour="black"),
      axis.text.y=element_text(size=y.axis.size,colour="black"),
      axis.ticks=element_line(colour="black"),
      axis.title.x=element_text(size=x.axis.size+5,vjust=-0.4),
      axis.title.y=element_text(size=y.axis.size+5,vjust=0.3),
      title=element_text(size=x.axis.size+8,vjust=1.5))
    
    ## y-axis labeling
    temp<-datafeature[!is.na(datafeature[,"ABUNDANCE"]) & !is.na(datafeature[,"INTENSITY"]),]
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
        if(length(setdiff(temp.name,unique(datafeature$PROTEIN)))>0){
          dev.off()
          stop(paste("Please check protein name. Data set does not have this protein. -", paste(temp.name, collapse=", "),sep=" "))
        }
      }
      
      ## check which.Protein is order number of Protein
      if(is.numeric(which.Protein)){
        temp.name<-levels(datafeature$PROTEIN)[which.Protein]
        
        ## message if name of Protein is wrong.
        if(length(levels(datafeature$PROTEIN))<max(which.Protein)){
          dev.off()
          stop(paste("Please check your selection of proteins. There are ", length(levels(datafeature$PROTEIN))," proteins in this dataset.",sep=" "))
        }
      }
      
      ## use only assigned proteins
      datafeature<-datafeature[which(datafeature$PROTEIN %in% temp.name),]
      datafeature$PROTEIN<-factor(datafeature$PROTEIN)
    }
    
    for (i in 1:nlevels(datafeature$PROTEIN)){	
      sub<-datafeature[datafeature$PROTEIN==levels(datafeature$PROTEIN)[i],]
      subTemp<-sub[!is.na(sub$ABUNDANCE),]
      #subTemp$LABEL<-factor(subTemp$LABEL)	
      sub<-sub[with(sub, order(LABEL,RUN)),]
      
      ## if all measurements are NA,
      if(nrow(sub)==sum(is.na(sub$ABUNDANCE))){
        message(paste("Can't the Quality Control plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),") because all measurements are NAs."))
        next()
      }
      
      
      #			if(length(unique(subTemp$LABEL))==2){
      #				label.color<-c("darkseagreen1","lightblue")
      #			}else{
      #				if(unique(subTemp$LABEL)=="Endogenous"){
      #					label.color<-c("lightblue")	
      #				}
      #				if(unique(subTemp$LABEL)=="Reference"){
      #					label.color<-c("darkseagreen1")
      #				}
      #			}
      
      ##options(show.error.messages = FALSE)
      
      ptemp<-ggplot(aes_string(x='RUN', y='ABUNDANCE'), data=sub)+facet_grid(~LABEL)+geom_boxplot(aes_string(fill='LABEL'),outlier.shape=1,outlier.size=1.5)+scale_fill_manual(values=label.color, guide="none")+scale_x_discrete('MS runs',breaks=cumGroupAxis)+annotate("text",x=groupName$RUN,y=groupName$y,label=groupName$Name,size=text.size,angle=text.angle)+geom_vline(xintercept=lineNameAxis+0.5,colour="grey",linetype="longdash")+labs(title=unique(sub$PROTEIN))+theme(
        panel.background=element_rect(fill='white', colour="black"),
        legend.key=element_rect(fill='white',colour='white'),
        panel.grid.minor = element_blank(),
        strip.background=element_rect(fill='gray95'),	
        strip.text.x=element_text(colour=c("#00B0F6"),size=14),
        axis.text.x=element_text(size=x.axis.size,colour="black"),
        axis.text.y=element_text(size=y.axis.size,colour="black"),
        axis.ticks=element_line(colour="black"),
        axis.title.x=element_text(size=x.axis.size+5,vjust=-0.4),
        axis.title.y=element_text(size=y.axis.size+5,vjust=0.3),
        title=element_text(size=x.axis.size+8,vjust=1.5))
      
      ## y-axis labeling
      if(temptest){
        ptemp<-ptemp+scale_y_continuous('Log2-intensities',limit=c(y.limdown, y.limup))
      }else{
        ptemp<-ptemp+scale_y_continuous('Log10-intensities',limit=c(y.limdown, y.limup))
      }	
      
      print(ptemp)
      
      message(paste("Drew the Quality Contol plot(boxplot) for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),")"))
      
    } ## end-loop
    if(address!=FALSE) dev.off()
  } # end QC plot	
  
  
  #################### 
  ### Condition plot
  #################### 
  if(toupper(type)=="CONDITIONPLOT"){
    
    #### choose Proteins or not
    if(which.Protein!="all"){
      ## check which.Protein is name of Protein
      if(is.character(which.Protein)){
        
        temp.name<-which.Protein
        
        ## message if name of Protein is wrong.
        if(length(setdiff(temp.name,unique(datafeature$PROTEIN)))>0)
          stop(paste("Please check protein name. Dataset does not have this protein. -", paste(temp.name, collapse=", "),sep=" "))
      }
      
      ## check which.Protein is order number of Protein
      if(is.numeric(which.Protein)){
        
        temp.name<-levels(datafeature$PROTEIN)[which.Protein]
        
        ## message if name of Protein is wrong.
        if(length(levels(datafeature$PROTEIN))<max(which.Protein))
          stop(paste("Please check your selection of proteins. There are ", length(levels(datafeature$PROTEIN))," proteins in this dataset.",sep=" "))
      }
      
      ## use only assigned proteins
      datafeature<-datafeature[which(datafeature$PROTEIN %in% temp.name),]
      datafeature$PROTEIN<-factor(datafeature$PROTEIN)
    }
    
    #### save the plots as pdf or not
    # If there are the file with the same name, add next numbering at the end of file name		
    if(address!=FALSE){
      allfiles<-list.files()
      
      num<-0
      filenaming<-paste(address,"ConditionPlot",sep="")
      finalfile<-paste(address,"ConditionPlot.pdf",sep="")
      
      while(is.element(finalfile,allfiles)){
        num<-num+1
        finalfile<-paste(paste(filenaming,num,sep="-"),".pdf",sep="")
      }	
      
      pdf(finalfile, width=width, height=height)
    }
    
    if(nlevels(datafeature$LABEL)==1){
      for (i in 1:nlevels(datafeature$PROTEIN)){	
        sub<-datafeature[datafeature$PROTEIN==levels(datafeature$PROTEIN)[i],]
        sub<-na.omit(sub)	
        sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
        sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)	
        sub$FEATURE<-factor(sub$FEATURE)	
        
        ## if all measurements are NA,
        if(nrow(sub)==sum(is.na(sub$ABUNDANCE))){
          message(paste("Can't the Condition plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),") because all measurements are NAs."))
          next()
        }
        
        ######## statistics
        sub.mean<-by(sub$ABUNDANCE, sub$GROUP_ORIGINAL, function(x) mean(x,na.rm=TRUE))
        sub.sd<-by(sub$ABUNDANCE, sub$GROUP_ORIGINAL, sd)
        sub.len<-by(sub$ABUNDANCE, sub$GROUP_ORIGINAL, length)
        if(interval=="CI") ciw<-qt(0.975,sub.len)*sub.sd/sqrt(sub.len)
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
          
          ptemp<-ggplot(aes_string(x='Label', y='mean'), data=tempsummary)+geom_errorbar(aes(ymax = mean + ciw, ymin=mean - ciw),data=tempsummary, width=0.1,colour="red")+geom_point(size=dot.size.condition,colour="darkred")+scale_x_discrete('Condition')+geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+labs(title=unique(sub$PROTEIN))+theme(
            panel.background=element_rect(fill='white', colour="black"),
            panel.grid.major.y = element_line(colour="grey95"),
            panel.grid.minor.y = element_blank(),
            axis.text.x=element_text(size=x.axis.size,colour="black",angle=text.angle),
            axis.text.y=element_text(size=y.axis.size,colour="black"),
            axis.ticks=element_line(colour="black"),
            axis.title.x=element_text(size=x.axis.size+5,vjust=-0.4),
            axis.title.y=element_text(size=y.axis.size+5,vjust=0.3),
            title=element_text(size=x.axis.size+8,vjust=1.5))
          
        }else{
          #extract numeric value, don't use levels (because T1,T10,T3,...)
          # ?? sort??
          
          ## reformat as data.frame
          tempsummary<-data.frame(Label=as.numeric(gsub("\\D","",unique(sub$GROUP_ORIGINAL))),mean=as.vector(sub.mean),ciw=as.vector(ciw))
          
          ptemp<-ggplot(aes_string(x='Label', y='mean'), data=tempsummary)+geom_errorbar(aes(ymax = mean + ciw, ymin=mean - ciw),data=tempsummary, width=0.1,colour="red")+geom_point(size=dot.size.condition,colour="darkred")+scale_x_continuous('Condition')+geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+labs(title=unique(sub$PROTEIN))+theme(
            panel.background=element_rect(fill='white', colour="black"),
            panel.grid.major.y = element_line(colour="grey95"),
            panel.grid.minor.y = element_blank(),
            axis.text.x=element_text(size=x.axis.size,colour="black",angle=text.angle),
            axis.text.y=element_text(size=y.axis.size,colour="black"),
            axis.ticks=element_line(colour="black"),
            axis.title.x=element_text(size=x.axis.size+5,vjust=-0.4),
            axis.title.y=element_text(size=y.axis.size+5,vjust=0.3),
            title=element_text(size=x.axis.size+8,vjust=1.5))
          
        } ## scale : true
        
        ## y-axis labeling, find log 2 or log 10
        temp<-sub[!is.na(sub[,"ABUNDANCE"]) & !is.na(sub[,"INTENSITY"]),]
        temptest<-abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])
        
        if(temptest){
          ptemp<-ptemp+scale_y_continuous("Log2-intensities",limit=c(y.limdown, y.limup))
        }else{
          ptemp<-ptemp+scale_y_continuous("Log10-intensities",limit=c(y.limdown, y.limup))
        }
        
        print(ptemp)
        
        message(paste("Drew the condition plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),")"))
        
      } ## end-loop
    } ## label-free
    
    if(nlevels(datafeature$LABEL)==2){
      for (i in 1:nlevels(datafeature$PROTEIN)){	
        
        sub<-datafeature[datafeature$PROTEIN==levels(datafeature$PROTEIN)[i],]
        sub<-sub[!is.na(sub$ABUNDANCE),]	
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
        if(interval=="CI") ciw<-qt(0.975,sub.len)*sub.sd/sqrt(sub.len)
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
          
          ptemp<-ggplot(aes_string(x='Label', y='mean'), data=tempsummary)+geom_errorbar(aes(ymax = mean + ciw, ymin=mean - ciw),data=tempsummary, width=0.1,colour="red")+geom_point(size=dot.size.condition,colour="darkred")+scale_x_discrete('Condition')+geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+labs(title=unique(sub$PROTEIN))+theme(
            panel.background=element_rect(fill='white', colour="black"),
            panel.grid.major.y = element_line(colour="grey95"),
            panel.grid.minor.y = element_blank(),
            axis.text.x=element_text(size=x.axis.size,colour="black",angle=text.angle),
            axis.text.y=element_text(size=y.axis.size,colour="black"),
            axis.ticks=element_line(colour="black"),
            axis.title.x=element_text(size=x.axis.size+5,vjust=-0.4),
            axis.title.y=element_text(size=y.axis.size+5,vjust=0.3),
            title=element_text(size=x.axis.size+8,vjust=1.5))
        }else{
          ## reformat as data.frame
          tempsummary<-data.frame(Label=as.numeric(gsub("\\D","",unique(sub.l$GROUP_ORIGINAL))),mean=as.vector(sub.mean),ciw=as.vector(ciw))
          
          ptemp<-ggplot(aes_string(x='Label', y='mean'), data=tempsummary)+geom_errorbar(aes(ymax = mean + ciw, ymin=mean - ciw),data=tempsummary, width=0.1,colour="red")+geom_point(size=dot.size,colour="darkred")+scale_x_continuous('Condition')+geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+labs(title=unique(sub$PROTEIN))+theme(
            panel.background=element_rect(fill='white', colour="black"),
            panel.grid.major.y = element_line(colour="grey95"),
            panel.grid.minor.y = element_blank(),
            axis.text.x=element_text(size=x.axis.size,colour="black",angle=text.angle),
            axis.text.y=element_text(size=y.axis.size,colour="black"),
            axis.ticks=element_line(colour="black"),
            axis.title.x=element_text(size=x.axis.size+5,vjust=-0.4),
            axis.title.y=element_text(size=y.axis.size+5,vjust=0.3),
            title=element_text(size=x.axis.size+8,vjust=1.5))
        } ## scale : true
        
        ## y-axis labeling, find log 2 or log 10
        temp<-sub[!is.na(sub[,"ABUNDANCE"]) & !is.na(sub[,"INTENSITY"]),]
        temptest<-abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])
        
        if(temptest){
          ptemp<-ptemp+scale_y_continuous("Log2-Ratio(L/H)",limit=c(y.limdown, y.limup))
        }else{
          ptemp<-ptemp+scale_y_continuous("Log10-Ratio(L/H)",limit=c(y.limdown, y.limup))
        }
        
        print(ptemp)
        message(paste("Drew the condition plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),")"))
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

groupComparison<-function(contrast.matrix=contrast.matrix,data=data){
  
  scopeOfBioReplication<-"expanded"
  scopeOfTechReplication<-"restricted"
  
  ## save process output in each step
  allfiles<-list.files()
  filenaming<-"msstats"
  
  
  if(length(grep(filenaming,allfiles))==0){
    
    finalfile<-"msstats.log"
    processout<-NULL
    
  }else{
    
    num<-0
    finalfile<-"msstats.log"
    
    while(is.element(finalfile,allfiles)){
      num<-num+1
      lastfilename<-finalfile ## in order to rea
      finalfile<-paste(paste(filenaming,num,sep="-"),".log",sep="")
    }
    
    finalfile<-lastfilename
    processout<-as.matrix(read.table(finalfile, header=T, sep="\t"))
  }
  
  processout<-rbind(processout,as.matrix(c(" "," ","MSstats - groupComparison function"," "),ncol=1))
  
  
  ## check input is correct
  ## data format
  rawinput<-c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Run","Intensity")
  
  if(length(setdiff(toupper(rawinput),toupper(colnames(data))))==0){
    processout<-rbind(processout,c(paste("The required input - data : did not process from dataProcess function. - stop")))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    stop("Please use 'dataProcess' first. Then use output of dataProcess function as input in groupComparison.")
  }
  
  
  ## contrast. matrix
  if(ncol(contrast.matrix)!=length(unique(data$ProcessedData$GROUP_ORIGINAL))){
    processout<-rbind(processout,c(paste("The required input - contrast.matrix: the number of column and the number of group are not the same. - stop")))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    stop("Please check contrast matrix. The number of group in data set is different with columns of contrast.matrix.")
  }
  
  # check whether row.names of contrast.matrix.sub exists or not
  if(sum(is.null(row.names(contrast.matrix)))>0){
    processout<-rbind(processout,c(paste("The required input - contrast.matrix: need row names of contrast.matrix . - stop")))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    stop("No row.names of comparison exist.\n")
  }
  
  
  ## other option value
#  if(!(scopeOfTechReplication=="restricted" | scopeOfTechReplication=="expanded")){
#    processout<-rbind(processout,c(paste("The required input - scopeOfTechReplication : 'scopeOfTechReplication' value is wrong. - stop")))
#    write.table(processout, file=finalfile, row.names=FALSE)
    
#    stop("'scopeOfTechReplication' must be one of \"restricted\" or \"expanded\".
#	")
#  }
  
  
  if(!(scopeOfBioReplication=="restricted" | scopeOfBioReplication=="expanded")){
    processout<-rbind(processout,c(paste("The required input - scopeOfBioReplication : 'scopeOfBioReplication' value is wrong. - stop")))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    stop("'scopeOfBioReplication' must be one of \"restricted\" or \"expanded\".
	")
  } 
  
  labeled<-ifelse(length(unique(data$ProcessedData$LABEL))==1, FALSE, TRUE)
  
  ## all input
  processout<-rbind(processout,c(paste("labeled = ",labeled,sep="")))
  processout<-rbind(processout,c(paste("scopeOfBioReplication = ",scopeOfBioReplication,sep="")))
#  processout<-rbind(processout,c(paste("scopeOfTechReplication = ", scopeOfTechReplication,sep="")))
  
  write.table(processout, file=finalfile, row.names=FALSE)
  
  
  ## check whether case-control(FALSE) or time-course(TRUE)
  repeated<-.checkRepeated(data$ProcessedData)
  
  if(repeated){ 
    processout<-rbind(processout,c(paste("Time course design of experiment - okay")))
  }else{
    processout<-rbind(processout,c(paste("Case control design of experiment - okay")))
    
  }
  write.table(processout, file=finalfile, row.names=FALSE)
  
  
  #data$PROTEIN<-factor(data$PROTEIN)	
  
  ## for final result report
  out<-NULL
  outsummary<-NULL
  outfitted<-NULL
  dataafterfit<-NULL
  
  ### check input for labeled
  
  ###==================================================
  ### start to analyze by protein ID
 	message(paste("\n Start to test and get inference in whole plot..."))
 	
 	if(data$SummaryMethod=="skyline"){
 		message(paste("\n * Use t-test for log sum summarization per subject."))
 		processout<-rbind(processout,c(paste("Use t-test for log sum summarization per subject.")))
 		write.table(processout, file=finalfile, row.names=FALSE)

 		if(repeated){
 			
 			processout<-rbind(processout,c(paste("Only group comparsion is available for t-test.- stop")))
 			write.table(processout, file=finalfile, row.names=FALSE)

 			stop("\n * Only group comparsion is available for t-test.")
 		}
 	}
 
  ### need original group information
  rqall<-data$RunlevelData

  origGroup<-unique(rqall$GROUP_ORIGINAL)
  
  for (i in 1:nlevels(rqall$Protein)){
    
    sub<-rqall[rqall$Protein==levels(rqall$Protein)[i],]
    colnames(sub)[colnames(sub)=="LogIntensities"]<-"ABUNDANCE"
    colnames(sub)[colnames(sub)=="Protein"]<-"PROTEIN"

    # it is important to remove NA first, before we have the correct structure for each factor
    sub<-sub[!is.na(sub$ABUNDANCE),]
    
    ### 1. for logsum, t-test
    if(data$SummaryMethod=="skyline"){
    	
    		## subject level
    		sub$GROUP<-factor(sub$GROUP)
   		sub$SUBJECT<-factor(sub$SUBJECT)
    		sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	

    
    		########### testing and inference in whole plot
    		message(paste("Testing a comparison for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(rqall$Protein)),")"))
    	
     	temp<-try(.ttest.logsum(contrast.matrix,sub,origGroup),silent=TRUE)
     	
     }else{ ### linear model
     	
    		sub$GROUP<-factor(sub$GROUP)
   	 	sub$SUBJECT<-factor(sub$SUBJECT)
    		sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
    		sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)
    		sub$SUBJECT_NESTED<-factor(sub$SUBJECT_NESTED)
    		# sub$FEATURE<-factor(sub$FEATURE)	
    		sub$RUN<-factor(sub$RUN)
    
   	 	# singleFeature<-.checkSingleFeature(sub)
    		singleSubject<-.checkSingleSubject(sub)
    		TechReplicate<-.checkTechReplicate(sub) ## use for label-free model
    
    		# MissGroupByFeature<-.checkMissGroupByFeature(sub)
    		# MissRunByFeature<-.checkMissRunByFeature(sub)
    		MissSubjectByGroup<-.checkRunbyFeature(sub)
    		UnequalSubject<-.checkUnequalSubject(sub)
        
    
    		########### testing and inference in whole plot
    		message(paste("Testing a comparison for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(rqall$Protein)),")"))
    
     	## fit the model 
     	temp<-try(.fit.model.single(contrast.matrix,sub,scopeOfTechReplication,scopeOfBioReplication,TechReplicate,singleSubject,repeated,origGroup),silent=TRUE)
     	
     }

    ## fix(apr 16)
    if(class(temp)=="try-error") {
      message("*** error : can't analyze ", levels(rqall$Protein)[i], " for comparison.")
      
      processout<-rbind(processout,c(paste("error : can't analyze ", levels(rqall$Protein)[i], " for comparison.", sep = "")))
      write.table(processout, file=finalfile, row.names=FALSE)
      
      tempresult<-list(result=NULL,valueresid=NULL, valuefitted=NULL, fittedmodel=NULL)
      
      for(k in 1:nrow(contrast.matrix)){	
        tempresult$result<-rbind(tempresult$result, data.frame(Protein=levels(rqall$Protein)[i],Label=row.names(contrast.matrix)[k], logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA))
      }
    }else{
      tempresult<-temp
    }
    
    ## comparison result table
    #	out<-rbindlist(list(out,tempresult$result))
    out<-rbind(out,tempresult$result)
    
    ## for checking model assumptions
    ## add residual and fitted after fitting the model
    if(class(temp)=="try-error") {
      if(nrow(sub)!=0){
        sub$residuals<-NA
        sub$fitted<-NA
      }
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
    
    #	dataafterfit<-rbindlist(list(dataafterfit,sub))
    dataafterfit<-rbind(dataafterfit,sub)
    
    ## save fitted model
    outfitted<-c(outfitted, list(tempresult$fittedmodel))
    
    ##
    processout<-rbind(processout,c(paste("Finished a comparison for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(rqall$Protein)),")")))
    write.table(processout, file=finalfile, row.names=FALSE)
    
  } ### end protein loop
  
  ##
  processout<-rbind(processout,c("Comparisons for all proteins are done.- okay"))
  write.table(processout, file=finalfile, row.names=FALSE)
  
  
  ##### finalize result
  ## need to FDR per comparison
  out.all<-NULL
  
  out$Label<-as.character(out$Label)
  
  for(i in 1:length(unique(out$Label))){
    outsub<-out[out$Label==unique(out$Label)[i],]
    outsub<-data.frame(outsub,adj.pvalue=p.adjust(outsub$pvalue,method="BH"))
    #	out.all<-rbindlist(list(out.all, outsub))
    out.all<-rbind(out.all, outsub)
    
  }
  out.all$Label<-factor(out.all$Label)
  
  
  ##
  processout<-rbind(processout,c("Adjust p-values per comparison are calculated - okay."))
  write.table(processout, file=finalfile, row.names=FALSE)
  
  
  temp<-data$ProcessedData[!is.na(data$ProcessedData[,"ABUNDANCE"]) & !is.na(data$ProcessedData[,"INTENSITY"]),]
  
  if(abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])<
       abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])){
    colnames(out.all)[3]<-"log2FC"
  }
  
  if(abs(log2(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])>
       abs(log10(temp[1,"INTENSITY"])-temp[1,"ABUNDANCE"])){
    colnames(out.all)[3]<-"log10FC"
  }
  
  ## change the format as data.frame
  out.all<-data.frame(out.all)
  
  ##
  processout<-rbind(processout,c("Group comparison is done. - okay"))
  write.table(processout, file=finalfile, row.names=FALSE)
  
  
  finalout<-list(ComparisonResult=out.all, ModelQC=dataafterfit, fittedmodel=outfitted)
  return(finalout)	
  
}


#############################################
#############################################
### Model-based quality control plot 
#############################################
#############################################

modelBasedQCPlots<-function(data,type,axis.size=10,dot.size=3,text.size=7,legend.size=7,width=10, height=10,featureName=TRUE,feature.QQPlot="all",which.Protein="all",address=""){
  
  if(length(setdiff(toupper(type),c("QQPLOTS","RESIDUALPLOTS")))!=0){
    stop(paste("Input for type=",type,". However,'type' should be one of \"QQPlots\", \"ResidualPlots\" .",sep=""))
  }
  
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
        stop(paste("Please check protein name. Dataset does not have this protein. -", paste(temp.name, collapse=", "),sep=" "))
    }
    
    ## check which.Protein is order number of Protein		
    if(is.numeric(which.Protein)){
      
      temp.name<-levels(data$PROTEIN)[which.Protein]
      
      ## message if name of Protein is wrong.
      if(length(levels(data$PROTEIN))<max(which.Protein))
        stop(paste("Please check your selection of proteins. There are ", length(levels(data$PROTEIN))," proteins in this dataset.",sep=" "))
    }
    
    ## use only assigned proteins
    data<-data[which(data$PROTEIN %in% temp.name),]
    data$PROTEIN<-factor(data$PROTEIN)
  }
  
  #############################################
  ### normality, QQ plot
  #############################################
  if(toupper(type)=="QQPLOTS"){
    
    ### test feature.QQPlot is something wrong.
    
    
    #### one QQplot per protein, all features together
    if(toupper(feature.QQPlot)=="ALL"){
      
      #### save the plots as pdf or not
      # If there are the file with the same name, add next numbering at the end of file name	
      if(address!=FALSE){
        allfiles<-list.files()
        
        num<-0
        filenaming<-paste(address,"QQPlot_allFeatures",sep="")
        finalfile<-paste(address,"QQPlot_allFeatures.pdf",sep="")
        
        while(is.element(finalfile,allfiles)){
          num<-num+1
          finalfile<-paste(paste(filenaming,num,sep="-"),".pdf",sep="")
        }	
        
        pdf(finalfile, width=width, height=height)
      }
      
      for (i in 1:nlevels(data$PROTEIN)){	
        
        sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
        
        ## get slope and intercept for qline
        y<-quantile(sub$residuals[!is.na(sub$residuals)], c(0.25, 0.75))
        x<-qnorm(c(0.25, 0.75))
        slope<-diff(y)/diff(x)
        int<-y[1L]-slope*x[1L]
        
        ptemp<-ggplot(sub, aes_string(sample='residuals'))+geom_point(stat="qq",alpha=0.8,shape=1,size=dot.size)+scale_shape(solid=FALSE)+ geom_abline(slope = slope, intercept = int,colour="red")+scale_y_continuous('Sample Quantiles')+scale_x_continuous('Theoretical Quantiles')+labs(title=paste("Normal Q-Q Plot (",unique(sub$PROTEIN),")"))+theme(
          panel.background=element_rect(fill='white', colour="black"),
          panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor =element_blank(),
          axis.text.x=element_text(size=axis.size,colour="black"),
          axis.text.y=element_text(size=axis.size,colour="black"),
          axis.ticks=element_line(colour="black"),
          axis.title.x=element_text(size=axis.size+5,vjust=-0.4),
          axis.title.y=element_text(size=axis.size+5,vjust=0.3),
          title=element_text(size=axis.size+8,vjust=1.5),
          legend.position="none")
        
        print(ptemp)
        
        message(paste("Drew the QQ plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
        
      } ## end loop
      
      if(address!=FALSE) dev.off()
    } ## end QQplot by all feature
    
    
    #### panel for each feature,
    if(toupper(feature.QQPlot)=="BYFEATURE"){
      
      #### save the plots as pdf or not
      # If there are the file with the same name, add next numbering at the end of file name	
      if(address!=FALSE){
        allfiles<-list.files()
        
        num<-0
        filenaming<-paste(address,"QQPlot_byFeatures",sep="")
        finalfile<-paste(address,"QQPlot_byFeatures.pdf",sep="")
        
        while(is.element(finalfile,allfiles)){
          num<-num+1
          finalfile<-paste(paste(filenaming,num,sep="-"),".pdf",sep="")
        }	
        
        pdf(finalfile, width=width, height=height)
      }
      
      for (i in 1:nlevels(data$PROTEIN)){	
        
        sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
        
        
        ## label-free
        if(length(unique(sub$LABEL))==1){
          
          ## need to update for qline per feature
          ptemp<-ggplot(sub, aes_string(sample='residuals',color='FEATURE'))+geom_point(stat="qq",alpha=0.8,size=dot.size)+facet_wrap(~FEATURE)+scale_y_continuous('Sample Quantiles')+scale_x_continuous('Theoretical Quantiles')+labs(title=paste("Normal Q-Q Plot (",unique(sub$PROTEIN),")"))+theme(
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
          
          message(paste("Drew the QQ plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
          
        }
        
        ## label-based : seperate endogenous and reference
        if(length(unique(sub$LABEL))==2){
          ## need to update for qline per feature
          
          ### endogenous intensities
          sub.l<-sub[sub$LABEL=="L",]
          ptemp<-ggplot(sub.l, aes_string(sample='residuals',color='FEATURE'))+geom_point(stat="qq",alpha=0.8,size=dot.size)+facet_wrap(~FEATURE)+scale_y_continuous('Sample Quantiles')+scale_x_continuous('Theoretical Quantiles')+labs(title=paste("Normal Q-Q Plot (",unique(sub$PROTEIN),") - Endogenous Intensities"))+theme(
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
          ptemp<-ggplot(sub.h, aes_string(sample='residuals',color='FEATURE'))+geom_point(stat="qq",alpha=0.8,size=dot.size)+facet_wrap(~FEATURE)+scale_y_continuous('Sample Quantiles')+scale_x_continuous('Theoretical Quantiles')+labs(title=paste("Normal Q-Q Plot (",unique(sub$PROTEIN),") - Reference Intensities"))+theme(
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
          
          message(paste("Drew the QQ plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
          
        }
        
      } ## end loop
      
      if(address!=FALSE) dev.off()	
    } ### end for QQ by FEATURE
    
  } ### end QQplots
  
  
  #############################################
  #### Residual plot
  #############################################
  if(toupper(type)=="RESIDUALPLOTS"){
    
    y.limdown<-min(data$residuals,na.rm=TRUE)
    y.limup<-max(data$residuals,na.rm=TRUE)
    
    x.limdown<-min(data$fitted,na.rm=TRUE)
    x.limup<-max(data$fitted,na.rm=TRUE)
    
    
    #### save the plots as pdf or not
    # If there are the file with the same name, add next numbering at the end of file name	
    if(address!=FALSE){
      allfiles<-list.files()
      
      num<-0
      filenaming<-paste(address,"ResidualPlot",sep="")
      finalfile<-paste(address,"ResidualPlot.pdf",sep="")
      
      while(is.element(finalfile,allfiles)){
        num<-num+1
        finalfile<-paste(paste(filenaming,num,sep="-"),".pdf",sep="")
      }	
      
      pdf(finalfile, width=width, height=height)
    }
    
    for (i in 1:nlevels(data$PROTEIN)){	
      
      sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
      sub$PEPTIDE<-factor(sub$PEPTIDE)
      sub$FEATURE<-factor(sub$FEATURE)
      
      ptemp<-ggplot(aes_string(x='fitted', y='residuals', color='FEATURE', shape='LABEL'), data=sub)+geom_point(size=dot.size,alpha=0.5)+geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+scale_y_continuous('Residuals',limit=c(y.limdown,y.limup))+scale_x_continuous('Predicted Abundance',limit=c(x.limdown,x.limup))+labs(title=levels(data$PROTEIN)[i])
      
      if(length(unique(sub$LABEL))==2){
        ptemp<-ptemp+scale_shape_manual(values=c(2,19),name="",labels=c("Reference","Endogenous"))
      }else{
        ptemp<-ptemp+scale_shape_manual(values=c(19),name="",labels=c("Endogenous"))
      }
      
      
      if(featureName){
        ptemp<-ptemp+theme(
          panel.background=element_rect(fill='white', colour="black"),
          legend.key=element_rect(fill='white',colour='white'),
          legend.text=element_text(size=legend.size),
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
      
      message(paste("Drew the residual plot for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
      
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
# check whether there are multiple runs for a replicate
# if yes, normalization should be different way.
#############################################

.countMultiRun<-function(data){
  
 standardFeature<-unique(data[data$RUN==unique(data$RUN[1]),"FEATURE"]) ## if some feature are missing for this spedific run, it could be error. that is why we need balanced design.
  
  ## get overlapped feature ID
  countdiff = tapply (data$FEATURE, data$RUN, function ( x ) length(intersect(unique(x),standardFeature)) ) 
  
  return(countdiff)
}

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
  
  temp<-unique(data[,c("SUBJECT_NESTED","RUN")])
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
.runQuantification<-function(data,summaryMethod,equalFeatureVar,filterLogOfSum,cutoffCensored,censoredInt,remove50missing,MBimpute){
	
    data$LABEL<-factor(data$LABEL)
    label<-nlevels(data$LABEL)==2
    
   	# set ref which is distinguish reference and endogenous. any reference=0. endogenous is the same as RUN
	if(label){
		data$ref<-0
		data$ref[data$LABEL!="H"]<-data$RUN[data$LABEL!="H"]
		data$ref<-factor(data$ref)
#		unique(data[,c("RUN","LABEL","GROUP","ref")])
	}
	      
#    finalresult<-data.frame(Protein=rep(levels(data$PROTEIN),each=nlevels(data$RUN)),RUN=rep(c(levels(data$RUN)),nlevels(data$PROTEIN)),Condition=NA, BioReplicate=NA,LogIntensities=NA,NumFeature=NA,NumPeaks=NA)

	# for saving predicting value for impute option
	predAbundance<-NULL
	
	###################################
	## method 1 : model based summarization
	if(summaryMethod=="linear"  & is.null(censoredInt)){
		
		data<-data[!is.na(data$ABUNDANCE),]
    	data$PROTEIN<-factor(data$PROTEIN)
    	data$RUN<-factor(data$RUN)
    
		result<-NULL
	    dataafterfit<-NULL
	    
		for(i in 1: nlevels(data$PROTEIN)){
        
     		sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]

#      sub$GROUP<-factor(sub$GROUP)
#      sub$SUBJECT<-factor(sub$SUBJECT)
#      sub$GROUP_ORIGINAL<-factor(sub$GROUP_ORIGINAL)	
#      sub$SUBJECT_ORIGINAL<-factor(sub$SUBJECT_ORIGINAL)
     		sub$SUBJECT_NESTED<-factor(sub$SUBJECT_NESTED)
      		sub$FEATURE<-factor(sub$FEATURE)	
      		sub$RUN<-factor(sub$RUN)	        
      
      		if(!label){
      			temp<-data.frame(xtabs(~RUN, data=sub))
	
      			sub.result<-data.frame(Protein=rep(levels(data$PROTEIN)[i],each=nlevels(sub$RUN)),RUN=rep(c(levels(sub$RUN)),1),LogIntensities=NA, NumFeature=length(unique(sub$FEATURE)),NumPeaks=temp$Freq)
      	
      		}else{
      	
      			sub$ref<-factor(sub$ref)			

      			temp<-data.frame(xtabs(~ref, data=sub))

      			sub.result<-data.frame(Protein=rep(levels(data$PROTEIN)[i],each=nlevels(sub$ref)),RUN=rep(c(levels(sub$ref)[-1],"Ref"),1),LogIntensities=NA, NumFeature=length(unique(sub$FEATURE)),NumPeaks=c(temp[-1,"Freq"],temp[1,"Freq"]))
      	
      		}
      
      		singleFeature<-.checkSingleFeature(sub)
      		singleSubject<-.checkSingleSubject(sub)
      		TechReplicate<-.checkTechReplicate(sub) ## use for label-free model
        
      		##### fit the model
      		message(paste("Getting the summarization per subplot for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
        
      		fit<-try(.fit.quantification.run(sub,singleFeature,singleSubject, TechReplicate,labeled=label, equalFeatureVar), silent=TRUE)
             
      		if(class(fit)=="try-error") {
        		message("*** error : can't fit the model for ", levels(data$PROTEIN)[i])
          
        		result<-rbind(result, sub.result)
          
         		if(nrow(sub)!=0){
        			sub$residuals<-NA
        			sub$fitted<-NA
      			}
          
      		}else{
          
       			if(class(fit)=="lm"){
          			cf <- summary(fit)$coefficients
       	 		}else{
         			cf <- fixef(fit)
         		}
          
        		# calculate sample quantification for all levels of sample
        		a=1	
          
        		for(j in 1:nlevels(sub$RUN)){
         			contrast.matrix<-rep(0,nlevels(sub$RUN))
          			contrast.matrix[j]<-1
            
          			contrast<-.make.contrast.run.quantification(fit,contrast.matrix,sub, labeled=label)
            
         			if(class(fit)=="lm"){
           				sub.result[a,3]<-.estimableFixedQuantification(cf,contrast)
          			}else{
            			sub.result[a,3]<-.estimableRandomQuantification(cf,contrast)
          			}
          
          			a=a+1
        		}
          
          		## for label-based case, need reference quantification
        		if(label){
        			contrast<-.make.contrast.run.quantification.reference(fit,contrast.matrix,sub)
        			if(class(fit)=="lm"){
              			sub.result[a,3]<-.estimableFixedQuantification(cf,contrast)
          			}else{
              			sub.result[a,3]<-.estimableRandomQuantification(cf,contrast)
          			}
        		}
          
        		result<-rbind(result, sub.result)

 				if(class(fit)=="lm"){  ### lm model
 					sub$residuals<-fit$residuals
      				sub$fitted<-fit$fitted.values
  				}else{   ### lmer model
    				sub$residuals<-resid(fit)
    				sub$fitted<-fitted(fit)
  				}
  			
     			dataafterfit<-rbind(dataafterfit,sub)

          	}
        
      	} ## end-loop for each protein	
	} ## for linear model summary
	
	###################################
	## Method 2 : Tukey Median Polish	
	if(summaryMethod=="TMP"){
		
		#data<-data[!is.na(data$ABUNDANCE),]
   	 	data$PROTEIN<-factor(data$PROTEIN)
    	data$RUN<-factor(data$RUN)
    
		result<-NULL
	  
		for(i in 1: nlevels(data$PROTEIN)){
              		
     		sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
     		
     		message(paste("Getting the summarization by Tukey's median polish per subplot for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))

      		sub$FEATURE<-factor(sub$FEATURE)	
      		sub$feature.label<-paste(sub$FEATURE, sub$LABEL, sep="_")
      		sub$run.label<-paste(sub$RUN, sub$LABEL, sep="_")
      		
      		## if all measurements are NA,
      		if(nrow(sub)==sum(is.na(sub$ABUNDANCE))){
       			message(paste("Can't summarize for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),") because all measurements are NAs."))
        		next()
      		}
      		
      		## remove features which are completely NAs
			subtemp<-sub[sub$LABEL=="L" & !is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
			countfeature<-xtabs(~FEATURE, subtemp)
			namefeature<-names(countfeature)[countfeature==0]
				
			if(length(namefeature)!=0){
				sub<-sub[-which(sub$FEATURE %in% namefeature),]
				sub$FEATURE<-factor(sub$FEATURE)
			}
      		
      		## remove run which has no measurement at all 
			subtemp<-sub[sub$LABEL=="L" & !is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
			count<-aggregate(ABUNDANCE~RUN,data=subtemp, length)
			norun<-setdiff(unique(data$RUN),count$RUN)
				
			if(length(norun)!=0 & length(intersect(norun, as.character(unique(sub$RUN))))){ # removed NA rows already, if there is no overlapped run, error
				sub<-sub[-which(sub$RUN %in% norun),]
				sub$RUN<-factor(sub$RUN)
			}

			
			if(remove50missing){
				# count # feature per run
        			if(!is.null(censoredInt)){
						if(censoredInt=="NA"){
							subtemp<-sub[sub$LABEL=="L" & !is.na(sub$INTENSITY),]
						}
					
						if(censoredInt=="0"){
							subtemp<-sub[sub$LABEL=="L" & !is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
						}
					}
					
					numFea<-xtabs(~RUN, subtemp) ## RUN or run.label?
					numFea<-numFea/length(unique(subtemp$FEATURE))
					numFea<-numFea<=0.5
					removerunid<-names(numFea)[numFea]
					
					## if all measurements are NA,
      				if(length(removerunid)==length(numFea)){
       					message(paste("Can't summarize for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),") because all runs have more than 50% NAs and are removed with the option, remove50missing=TRUE."))
        				next()
      				}

			}
			
			##### how to decide censored or not
			if(!is.null(censoredInt)){
				## 1. censored 
				if(censoredInt=="0"){
					sub$cen<-ifelse(!is.na(sub$INTENSITY) & sub$INTENSITY==0,0,1)
				}
				
				### 2. all censored missing
				if(censoredInt=="NA"){
					sub$cen<-ifelse(is.na(sub$INTENSITY),0,1)
				}
				
				### check whether we need to impute or not.
				if(sum(sub$cen==0)>0){
			
					##### cutoffCensored
					## 1. put 0 to censored
					#if(cutoffCensored=="0"){
					#	if(censoredInt=="NA"){
					#		sub[is.na(sub$INTENSITY),"ABUNDANCE"]<-0
					#	}
				
					#	if(censoredInt=="0"){
					#		sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"]<-0
					#	}
					#}
			
					## 2. put minimum in feature level to NA
					if(cutoffCensored=="minFeature"){
						if(censoredInt=="NA"){
							cut<-aggregate(ABUNDANCE~feature.label,data=sub, function(x) min(x, na.rm=TRUE))
							## cutoff for each feature is little less than minimum abundance in a run.
							cut$ABUNDANCE<-0.99*cut$ABUNDANCE
						
							## remove runs which has more than 50% missing values
							if(remove50missing){
								if(length(removerunid)!=0){
									sub<-sub[-which(sub$RUN %in% removerunid),]
									sub$RUN<-factor(sub$RUN)
								}
							}

							for(j in 1:length(unique(sub$feature.label))){
								sub[is.na(sub$INTENSITY) & sub$feature.label==cut$feature.label[j],"ABUNDANCE"]<-cut$ABUNDANCE[j]
							}
						}
					
						if(censoredInt=="0"){
							subtemptemp<-sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
							cut<-aggregate(ABUNDANCE~feature.label,data=subtemptemp, FUN=min)
							## cutoff for each feature is little less than minimum abundance in a run.
							cut$ABUNDANCE<-0.99*cut$ABUNDANCE
						
						
							## remove runs which has more than 50% missing values
							if(remove50missing){
								if(length(removerunid)!=0){
									sub<-sub[-which(sub$RUN %in% removerunid),]
									sub$RUN<-factor(sub$RUN)
								}
							}
						
							for(j in 1:length(unique(sub$feature.label))){
								sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0  & sub$feature.label==cut$feature.label[j],"ABUNDANCE"]<-cut$ABUNDANCE[j]
							}
						}
					}
				
					## 3. put minimum in RUN to NA
					if(cutoffCensored=="minRun"){
				
						## remove runs which has more than 50% missing values
						if(remove50missing){
							if(length(removerunid)!=0){
								sub<-sub[-which(sub$RUN %in% removerunid),]
								sub$RUN<-factor(sub$RUN)
							}
						}
						
						if(censoredInt=="NA"){
							cut<-aggregate(ABUNDANCE~run.label,data=sub, function(x) min(x, na.rm=TRUE))
							## cutoff for each Run is little less than minimum abundance in a run.
							cut$ABUNDANCE<-0.99*cut$ABUNDANCE

							for(j in 1:length(unique(sub$run.label))){
								sub[is.na(sub$INTENSITY) & sub$run.label==cut$run.label[j],"ABUNDANCE"]<-cut$ABUNDANCE[j]
							}
						}
					
						if(censoredInt=="0"){
							subtemptemp<-sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
							cut<-aggregate(ABUNDANCE~run.label,data=subtemptemp, FUN=min)
							cut$ABUNDANCE<-0.99*cut$ABUNDANCE

							for(j in 1:length(unique(sub$run.label))){
								sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$run.label==cut$run.label[j],"ABUNDANCE"]<-cut$ABUNDANCE[j]
							}
						}
					}
				
					## 20150829 : 4. put minimum RUN and FEATURE
					if(cutoffCensored=="minFeatureNRun"){
						if(censoredInt=="NA"){
						
							## cutoff for each feature is little less than minimum abundance in a run.

							cut.fea<-aggregate(ABUNDANCE~feature.label,data=sub, function(x) min(x, na.rm=TRUE))
							cut.fea$ABUNDANCE<-0.99*cut.fea$ABUNDANCE
						
							## remove runs which has more than 50% missing values
							## before removing, need to contribute min feature calculation
							if(remove50missing){
								if(length(removerunid)!=0){
									sub<-sub[-which(sub$RUN %in% removerunid),]
									sub$RUN<-factor(sub$RUN)
								}
							}
						
							## cutoff for each Run is little less than minimum abundance in a run.

							cut.run<-aggregate(ABUNDANCE~run.label,data=sub, function(x) min(x, na.rm=TRUE))
							cut.run$ABUNDANCE<-0.99*cut.run$ABUNDANCE
						
						
							if(length(unique(sub$feature.label))>1){
								for(j in 1:length(unique(sub$feature.label))){
									for(k in 1:length(unique(sub$run.label))){
										# get smaller value for min Run and min Feature
										finalcut<-min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
								
										sub[is.na(sub$INTENSITY) & sub$feature.label==cut.fea$feature.label[j] & sub$run.label==cut.run$run.label[k],"ABUNDANCE"]<-finalcut
									}
								}
							}
							# if single feature, not impute
						}
					
						if(censoredInt=="0"){
							subtemptemp<-sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]

							cut.fea<-aggregate(ABUNDANCE~feature.label,data=subtemptemp, FUN=min)
							cut.fea$ABUNDANCE<-0.99*cut.fea$ABUNDANCE
												
							## remove runs which has more than 50% missing values
							## before removing, need to contribute min feature calculation
							if(remove50missing){
								if(length(removerunid)!=0){
									sub<-sub[-which(sub$RUN %in% removerunid),]
									sub$RUN<-factor(sub$RUN)
								}
							}

							cut.run<-aggregate(ABUNDANCE~run.label,data=subtemptemp, FUN=min)
							cut.run$ABUNDANCE<-0.99*cut.run$ABUNDANCE
	
							if(length(unique(sub$feature.label))>1){
								for(j in 1:length(unique(sub$feature.label))){
									for(k in 1:length(unique(sub$run.label))){
										# get smaller value for min Run and min Feature
										finalcut<-min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
								
										sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$feature.label==cut.fea$feature.label[j] & sub$run.label==cut.run$run.label[k],"ABUNDANCE"]<-finalcut
									}
								}
							}else{ # single feature
					
								sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"]<-cut.fea$ABUNDANCE
							
							}
						}
					}
				
					if(MBimpute){
					
						if(nrow(sub[sub$cen==0,])>0){
							## impute by survival model
							subtemp<-sub[!is.na(sub$ABUNDANCE),]
							countdf<-nrow(subtemp)<(length(unique(subtemp$FEATURE))+length(unique(subtemp$RUN))-1)
				
							### fit the model
							if(length(unique(sub$FEATURE))==1){
								fittest<-survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN,data=sub, dist='gaussian')
							}else{
								if(countdf){
									fittest<-survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN,data=sub, dist='gaussian')
								}else{
									fittest<-survreg(Surv(ABUNDANCE, cen, type='left') ~ FEATURE+RUN,data=sub, dist='gaussian')
								}
							}
					
							# get predicted value from survival
							sub<-data.frame(sub, pred=predict(fittest, newdata=sub, type="response"))
					
							# the replace censored value with predicted value
							sub[sub$cen==0,"ABUNDANCE"]<-sub[sub$cen==0,"pred"]	
					
							# save predicted value
							predAbundance<-c(predAbundance,predict(fittest, newdata=sub, type="response"))				
						}
					}	
				}
			}
			
			## then, finally remove NA in abundance
			sub<-sub[!is.na(sub$ABUNDANCE),]
					       
					        
      		if(nlevels(sub$FEATURE)>1){ ## for more than 1 features
      			if(!label){ ## label-free
      			
      				data_w = dcast(RUN ~ FEATURE, data=sub, value.var='ABUNDANCE', keep=TRUE)
  					rownames(data_w)<-data_w$RUN
  					data_w<-data_w[,-1]
  					data_w[data_w==1]<-NA
  					
  					meddata <- medpolish(data_w,na.rm=TRUE,trace.iter = FALSE)
					tmpresult<-meddata$overall + meddata$row
					
					# count # feature per run
					if(!is.null(censoredInt)){
						if(censoredInt=="NA"){
							subtemp<-sub[!is.na(sub$INTENSITY),]
						}
					
						if(censoredInt=="0"){
							subtemp<-sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
						}
					}else{
						subtemp<-sub[!is.na(sub$INTENSITY),]
					}
					
					numFea<-xtabs(~RUN, subtemp)
					numFea<-numFea/length(unique(subtemp$FEATURE))
					numFea<-numFea<=0.5
					
					sub.result<-data.frame(Protein=unique(sub$PROTEIN),LogIntensities=tmpresult, RUN=names(tmpresult), more50missing=numFea)
			
      				result<-rbind(result, sub.result)
      			
      			}else{ ## labeled
					
					data_w = dcast(run.label ~ FEATURE, data=sub, value.var='ABUNDANCE', keep=TRUE)
  					rownames(data_w)<-data_w$run.label
  					data_w<-data_w[,-1]
  					#data_w[data_w==1]<-NA
  					
  					meddata <- medpolish(data_w,na.rm=TRUE,trace.iter = FALSE)
					tmpresult<-meddata$overall + meddata$row
					
					reformresult<-data.frame(tmpresult)
					end<-nchar(rownames(reformresult))
					reformresult$LABEL<-substr(rownames(reformresult),end,end)
					reformresult$RUN<-substr(rownames(reformresult),1,end-2)
					colnames(reformresult)[1]<-"ABUNDANCE"
					
					## now single feature, adjust reference feature difference
      				h<-reformresult[reformresult$LABEL=="H",]
 					allmed<-median(h$ABUNDANCE, na.rm=TRUE)

        			for (i in 1:length(unique(reformresult$RUN))){
          				## ABUNDANCE is normalized
          				reformresult[reformresult$RUN==unique(reformresult$RUN)[i],"ABUNDANCE"]<-reformresult[reformresult$RUN==unique(reformresult$RUN)[i],"ABUNDANCE"]-reformresult[reformresult$RUN==unique(reformresult$RUN)[i] & reformresult$LABEL=="H","ABUNDANCE"]+allmed
        			}
        			
        			reformresult<-reformresult[reformresult$LABEL=="L",]
        			
        			subtemp<-reformresult[!is.na(reformresult$ABUNDANCE),]
        			
        			# count # feature per run
        			if(!is.null(censoredInt)){
						if(censoredInt=="NA"){
							subtemp<-sub[sub$LABEL=="L" & !is.na(sub$INTENSITY),]
						}
					
						if(censoredInt=="0"){
							subtemp<-sub[sub$LABEL=="L" & !is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
						}
					}
					
					numFea<-xtabs(~RUN, subtemp)
					numFea<-numFea/length(unique(subtemp$FEATURE))
					numFea<-numFea<=0.5

      				sub.result<-data.frame(Protein=unique(sub$PROTEIN),LogIntensities=reformresult$ABUNDANCE, RUN=reformresult$RUN, more50missing=numFea)
      				
      				result<-rbind(result, sub.result)
      			}
      		
      		}else{ ## single feature
      			if(label){ ## label-based
      				## single feature, adjust reference feature difference
      				h<-sub[sub$LABEL=="H",]
 					allmed<-median(h$ABUNDANCE, na.rm=TRUE)

        			for (i in 1:length(unique(sub$RUN))){
          				## ABUNDANCE is normalized
          				sub[sub$RUN==unique(sub$RUN)[i],"ABUNDANCE"]<-sub[sub$RUN==unique(sub$RUN)[i],"ABUNDANCE"]-sub[sub$RUN==unique(sub$RUN)[i] & sub$LABEL=="H","ABUNDANCE"]+allmed
        			}
        			
        			sub<-sub[sub$LABEL=="L",]
      			}
      			
				## single feature, use original values
				
				
      			subtemp<-sub[!is.na(sub$ABUNDANCE),]
      			
      			if(!is.null(censoredInt)){
      				if(censoredInt=="NA"){
						subtempcount<-sub[sub$LABEL=="L" & !is.na(sub$INTENSITY),]
					}
					
					if(censoredInt=="0"){
						subtempcount<-sub[sub$LABEL=="L" & !is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
					}
				}else{
					subtempcount<-subtemp
				}

				
				numFea<-xtabs(~RUN, subtempcount)
				numFea<-numFea/length(unique(subtempcount$FEATURE))
				numFea<-numFea<=0.5
					
      			sub.result<-data.frame(Protein=subtemp$PROTEIN,LogIntensities=subtemp$ABUNDANCE, RUN=subtemp$RUN,  more50missing=numFea)
      				
      			result<-rbind(result, sub.result)
      		}
      	}  ## loop for proteins
      	
      	dataafterfit<-NULL
	}
		
	###################################
	## Method 3 : log sum	
	if(summaryMethod=="logOfSum"){
		
		if(label){
			message("* For log2(sum of intensities) summarization with label-based experiment, ratio between endogenous intensity and reference intensity is used.")
		}
		
		data<-data[!is.na(data$ABUNDANCE),]
    	data$PROTEIN<-factor(data$PROTEIN)
    	data$RUN<-factor(data$RUN)
    
		result<-NULL
		
		single<-unique(data[,c("GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
	  	singlesubject<-any(xtabs(~GROUP_ORIGINAL,single)==1)
	  	
		for(i in 1: nlevels(data$PROTEIN)){
              		
     		sub<-data[data$PROTEIN==levels(data$PROTEIN)[i],]
     		
     		message(paste("Getting the summarization by log2 (sum of intensities) per subject for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))

 			## 1.remove runs which have any missing values : same as Skyline Filtering
 			if(filterLogOfSum){
 				sub$FEATURE<-factor(sub$FEATURE)
				countmissing<-xtabs(~RUN, data=sub) ## already removed NAs
			
				if(label){
					nomissingrun<-names(countmissing[countmissing==(2*length(levels(sub$FEATURE)))])

				}else{
					nomissingrun<-names(countmissing[countmissing==length(levels(sub$FEATURE))])
				}
			
			###### if there is any run which have any missing values, remove those runs

				sub<-sub[which(sub$RUN %in% nomissingrun),]
			}
			
			if(nrow(sub)==0){
				
				message(paste("* All replicates in ",levels(data$PROTEIN)[i], " have any missing values. Can't summarize by log2 (sum of intensities)."))
				
				next
			}
			
		
			## 2. sum of normalized intensities per run
		
			#### back to original scale from normalized abundance
			sub$INTENSITY<-2^(sub$ABUNDANCE)
		
			#### 2.1 label based
			if(label){
				sub.l<-sub[sub$LABEL=="L",]
				sub.h<-sub[sub$LABEL=="H",]
				
				data_w.l = dcast(SUBJECT_ORIGINAL+RUN ~ FEATURE, data=sub.l, value.var='ABUNDANCE', keep=TRUE)

				data_w.h = dcast(SUBJECT_ORIGINAL+RUN ~ FEATURE, data=sub.h, value.var='ABUNDANCE', keep=TRUE)
				
				if(length(unique(sub$FEATURE))==1){
					
					data_w.l$sumInt<-data_w.l[,3]/data_w.h[,3]
					
					subsum<-data_w.l[,c("SUBJECT_ORIGINAL","RUN","sumInt")]
				
				}else{
					
					data_w<-data_w.l[,c(3:ncol(data_w.l))]/data_w.h[,c(3:ncol(data_w.h))]
					
					data_w<-cbind(SUBJECT_ORIGINAL=data_w.l$SUBJECT_ORIGINAL,RUN=data_w.l$RUN,data_w)
				
					data_w$sumInt<-apply(data_w[,-c(1:2)],1,sum)
				
					subsum<-data_w[,c("SUBJECT_ORIGINAL","RUN","sumInt")]

				}
				
			}else{
			#### 2.2 label-free

				subsum<-aggregate(sub$INTENSITY,list(SUBJECT_ORIGINAL=sub$SUBJECT_ORIGINAL,RUN=sub$RUN),sum, na.rm=TRUE)
				colnames(subsum)[3]<-"sumInt"
			}
			
			## 3. log2
			subsum$logsum<-log2(subsum$sumInt)
		
			if(!singlesubject){
				
				## 4. average per subject
				sublogsum<-aggregate(subsum$logsum,list(SUBJECT_ORIGINAL=subsum$SUBJECT_ORIGINAL),mean, na.rm=TRUE)
				colnames(sublogsum)[2]<-"LogIntensities"
		
				## 5. make output
				sub.result<-data.frame(Protein=unique(sub$PROTEIN),LogIntensities=sublogsum$LogIntensities, SUBJECT_ORIGINAL=sublogsum$SUBJECT_ORIGINAL)
				
			}else{
				
				## make output per run
				sub.result<-data.frame(Protein=unique(sub$PROTEIN),LogIntensities=subsum$logsum, RUN=subsum$RUN)

			}
			
      		result<-rbind(result, sub.result)
		} ## end loop for each protein
		
		dataafterfit<-NULL
	}
	
	###################################
	## method 4 : survival model for censored missing values
	if(summaryMethod=="linear" & !is.null(censoredInt)){
			
		#data<-data[!is.na(data$ABUNDANCE),]
    	data$PROTEIN<-factor(data$PROTEIN)
    	data$RUN<-factor(data$RUN)
    	
		if(label){
			
			result<-NULL

			for(i in 1:length(unique(data$PROTEIN))){
	
				sub<-data[data$PROTEIN==unique(data$PROTEIN)[i],]
				
				message(paste("Getting the summarization for censored missing values per subplot for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
				
				sub$FEATURE<-factor(sub$FEATURE)
				sub$feature.label<-paste(sub$FEATURE, sub$LABEL, sep="_")
				sub$run.label<-paste(sub$RUN, sub$LABEL, sep="_")

				## if all measurements are NA,
      			if(nrow(sub)==sum(is.na(sub$ABUNDANCE))){
       				message(paste("Can't summarize for ",unique(sub$PROTEIN), "(",i," of ",length(unique(datafeature$PROTEIN)),") because all measurements are NAs."))
        			next()
      			}
      			
				## remove run which has no measurement at all
				subtemp<-sub[sub$LABEL=="L" & !is.na(sub$INTENSITY),]
				count<-aggregate(ABUNDANCE~RUN,data=subtemp, length)
				norun<-setdiff(unique(data$RUN),count$RUN)
				
				if(length(norun)!=0 & length(intersect(norun, as.character(unique(sub$RUN))))){ # removed NA rows already, if there is no overlapped run, error
					sub<-sub[-which(sub$RUN %in% norun),]
					sub$RUN<-factor(sub$RUN)
				}
				
				if(length(unique(sub$RUN))==1){
				
					message(paste("* Only 1 MS run in ",levels(data$PROTEIN)[i], " has measurement. Can't summarize with censored intensities."))
				
					next
				}	
				
				## remove features which are completely NAs or zero
				subtemp<-sub[sub$LABEL=="L" & !is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
				countfeature<-xtabs(~FEATURE, subtemp)
				namefeature<-names(countfeature)[countfeature==0]
				
				if(length(namefeature)!=0){
					sub<-sub[-which(sub$FEATURE %in% namefeature),]
					sub$FEATURE<-factor(sub$FEATURE)
				}		
				
				##### how to decide censored or not
				## 1. censored 
				if(censoredInt=="0"){
					sub$cen<-ifelse(!is.na(sub$INTENSITY) & sub$INTENSITY==0,0,1)
				}
				
				### 2. all censored missing
				if(censoredInt=="NA"){
					sub$cen<-ifelse(is.na(sub$INTENSITY),0,1)
				}
				
				##### cutoffCensored
				## 1. put minimum in protein level to NA
				#if(cutoffCensored=="minEachProtein"){
				#	if(censoredInt=="NA"){
				#		cut<-min(sub$ABUNDANCE, na.rm=TRUE) 
				#		sub[is.na(sub$INTENSITY),"ABUNDANCE"]<-cut
				#	} 
					
				#	if(censoredInt=="0"){
				#		cut<-min(sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,"ABUNDANCE"]) 
				#		sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"]<-cut
				#	}
				#}
				
				## 2. put minimum in feature level to NA
				if(cutoffCensored=="minFeature"){
					if(censoredInt=="NA"){
						cut<-aggregate(ABUNDANCE~feature.label,data=sub, function(x) min(x, na.rm=TRUE))
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE<-0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$feature.label))){
							sub[is.na(sub$INTENSITY) & sub$feature.label==cut$feature.label[j],"ABUNDANCE"]<-cut$ABUNDANCE[j]
						}
					}
					
					if(censoredInt=="0"){
						subtemptemp<-sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
						cut<-aggregate(ABUNDANCE~feature.label,data=subtemptemp, FUN=min)
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE<-0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$feature.label))){
							sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0  & sub$feature.label==cut$feature.label[j],"ABUNDANCE"]<-cut$ABUNDANCE[j]
						}
					}
				}
				
				## 3. put minimum in RUN to NA
				if(cutoffCensored=="minRun"){
					if(censoredInt=="NA"){
						cut<-aggregate(ABUNDANCE~run.label,data=sub, function(x) min(x, na.rm=TRUE))
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE<-0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$run.label))){
							sub[is.na(sub$INTENSITY) & sub$run.label==cut$run.label[j],"ABUNDANCE"]<-cut$ABUNDANCE[j]
						}
					}
					
					if(censoredInt=="0"){
						subtemptemp<-sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
						cut<-aggregate(ABUNDANCE~run.label,data=subtemptemp, FUN=min)
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE<-0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$run.label))){
							sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$run.label==cut$run.label[j],"ABUNDANCE"]<-cut$ABUNDANCE[j]
						}
					}
				}
				
				## 20150829 : 4. put minimum RUN and FEATURE
				if(cutoffCensored=="minFeatureNRun"){
					if(censoredInt=="NA"){
						
						## cutoff for each feature is little less than minimum abundance in a run.

						cut.fea<-aggregate(ABUNDANCE~feature.label,data=sub, function(x) min(x, na.rm=TRUE))
						cut.fea$ABUNDANCE<-0.99*cut.fea$ABUNDANCE
						
						## cutoff for each Run is little less than minimum abundance in a run.

						cut.run<-aggregate(ABUNDANCE~run.label,data=sub, function(x) min(x, na.rm=TRUE))
						cut.run$ABUNDANCE<-0.99*cut.run$ABUNDANCE
						
						
						if(length(unique(sub$feature.label))>1){
							for(j in 1:length(unique(sub$feature.label))){
								for(k in 1:length(unique(sub$run.label))){
									# get smaller value for min Run and min Feature
									finalcut<-min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
								
									sub[is.na(sub$INTENSITY) & sub$feature.label==cut.fea$feature.label[j] & sub$run.label==cut.run$run.label[k],"ABUNDANCE"]<-finalcut
								}
							}
						}
							# if single feature, not impute
					}
					
					if(censoredInt=="0"){
						subtemptemp<-sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]

						cut.fea<-aggregate(ABUNDANCE~feature.label,data=subtemptemp, FUN=min)
						cut.fea$ABUNDANCE<-0.99*cut.fea$ABUNDANCE
												
						## remove runs which has more than 50% missing values
						## before removing, need to contribute min feature calculation
						if(remove50missing){
							if(length(removerunid)!=0){
								sub<-sub[-which(sub$RUN %in% removerunid),]
								sub$RUN<-factor(sub$RUN)
							}
						}

						cut.run<-aggregate(ABUNDANCE~run.label,data=subtemptemp, FUN=min)
						cut.run$ABUNDANCE<-0.99*cut.run$ABUNDANCE
	
						if(length(unique(sub$feature.label))>1){
							for(j in 1:length(unique(sub$feature.label))){
								for(k in 1:length(unique(sub$run.label))){
									# get smaller value for min Run and min Feature
									finalcut<-min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
								
									sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$feature.label==cut.fea$feature.label[j] & sub$run.label==cut.run$run.label[k],"ABUNDANCE"]<-finalcut
								}
							}
						}else{ # single feature
					
							sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"]<-cut.fea$ABUNDANCE
							
						}
					}
				}	
	
				## when number of measurement is less than df, error for fitting
				subtemp<-sub[!is.na(sub$ABUNDANCE),]
				countdf<-nrow(subtemp)<(length(unique(subtemp$FEATURE))+length(unique(subtemp$RUN))-1)
				
				
				### fit the model
				if(length(unique(sub$FEATURE))==1){
					
					# with single feature, not converge, wrong intercept
					# need to check
					fittest<-survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN+ref,data=sub, dist='gaussian')
				}else{
					if(countdf){
						fittest<-survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN+ref,data=sub, dist='gaussian')
					}else{
						fittest<-survreg(Surv(ABUNDANCE, cen, type='left') ~ FEATURE+RUN+ref,data=sub, dist='gaussian')
					}
				}
				
		
				sub.result<-data.frame(Protein=unique(sub$PROTEIN),RUN=rep(c(levels(sub$RUN)),1),LogIntensities=NA)

				# get the parameters
				cf <- summary(fittest)$coefficients

				# calculate sample quantification for all levels of sample
      			a=1	
          
       			for(j in 1:nlevels(sub$RUN)){
       				contrast.matrix<-rep(0,nlevels(sub$RUN))
        			contrast.matrix[j]<-1
          					contrast<-.make.contrast.run.quantification.Survival(fittest,contrast.matrix,sub, labeled=TRUE)
	
         			 sub.result[a,3]<-.estimableFixedQuantificationSurvival(cf,contrast)
         			 a=a+1
        		}

				result<-rbind(result, sub.result)
			}

			datamat = dcast( Protein ~ RUN, data=result, value.var='LogIntensities', keep=TRUE)
			datamat = melt(datamat, id.vars=c('Protein'))
			colnames(datamat)<-c('Protein','RUN','LogIntensities')
			result<-datamat
			
		}else{
			
			result<-NULL

			for(i in 1:length(unique(data$PROTEIN))){
	
				sub<-data[data$PROTEIN==unique(data$PROTEIN)[i],]
				
				message(paste("Getting the summarization for censored missing values per subplot for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
				
				sub$FEATURE<-factor(sub$FEATURE)
	
				## if all measurements are NA,
      			if(nrow(sub)==sum(is.na(sub$ABUNDANCE))){
       				message(paste("Can't summarize for ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),") because all measurements are NAs."))
        			next
      			}
      			
				## remove run which has no measurement at all
				subtemp<-sub[!is.na(sub$INTENSITY),]
				count<-aggregate(ABUNDANCE~RUN,data=subtemp, length)
				norun<-setdiff(unique(data$RUN),count$RUN)
				
				if(length(norun)!=0 & length(intersect(norun, as.character(unique(sub$RUN))))!=0){ # removed NA rows already, if there is no overlapped run, error
					sub<-sub[-which(sub$RUN %in% norun),]
					sub$RUN<-factor(sub$RUN)
				}
				
				if(length(unique(sub$RUN))==1){
				
					message(paste("* Only 1 MS run in ",levels(data$PROTEIN)[i], " has measurement. Can't summarize with censored intensities."))
				
					next
				}	
						
				
				## remove features which are (completely NAs or zero) 
				subtemp<-sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0 ,]
				countfeature<-xtabs(~FEATURE, subtemp)
				namefeature<-names(countfeature)[countfeature==0]
				
				if(length(namefeature)!=0){
					sub<-sub[-which(sub$FEATURE %in% namefeature),]
					sub$FEATURE<-factor(sub$FEATURE)
				}
				
				if(nrow(sub)==0){
				
					message(paste("* All measurements are NAs or only one measurement per feature in ",levels(data$PROTEIN)[i], ". Can't summarize with censored intensities."))
				
					next
				}	

				##### how to decide censored or not
				## 1. censored 
				if(censoredInt=="0"){
					sub$cen<-ifelse(!is.na(sub$INTENSITY) & sub$INTENSITY==0,0,1)
				}
				
				### 2. all censored missing
				if(censoredInt=="NA"){
					sub$cen<-ifelse(is.na(sub$INTENSITY),0,1)
				}

				## 3. decide random above some point
	
				## get runs which has all features
				#subtemp<-sub[!is.na(sub$INTENSITY),]
				#count<-aggregate(ABUNDANCE~RUN,data=subtemp, length)
	
				#completerun<-count[count$ABUNDANCE==length(unique(sub$FEATURE)),"RUN"]
	
				#if(length(completerun)!=0){
				#	subtemp<-sub[which(sub$RUN %in% completerun),]
				#}else{
				#	subtemp<-sub[which(sub$RUN %in% count[which.max(count$ABUNDANCE),"RUN"]),]
				#}
				
				# get feature mean and make order of feature
				# mean or median?
				#featureorder<-aggregate(ABUNDANCE~FEATURE,data=subtemp, mean)
				#featureorder<-featureorder[with(featureorder, order(ABUNDANCE, decreasing=T)),]
	
				# runs which has any missing
				#if(length(completerun)!=0){
				#	incompleterun<-count[count$ABUNDANCE!=length(unique(sub$FEATURE)),"RUN"]
				#}else{
				#	incompleterun<-count[-which.max(count$ABUNDANCE),"RUN"]
				#}
	
				#if(length(incompleterun)!=0){
				#	for(j in 1:length(incompleterun)){
				#		temp<-sub[sub$RUN==incompleterun[j],]
				#		temptemp<-temp[!is.na(temp$INTENSITY),]
		
				#		minfeature<-temptemp[which.min(temptemp$ABUNDANCE),"FEATURE"]
				#		abovefeature<-featureorder[1:which(featureorder$FEATURE==minfeature),"FEATURE"]
	
				#		sub[which(sub$RUN==incompleterun[j] & sub$FEATURE %in% abovefeature & is.na(sub$INTENSITY)),"ABUNDANCE"]<-NA
				#		sub[which(sub$RUN==incompleterun[j] & sub$FEATURE %in% abovefeature & is.na(sub$INTENSITY)),"cen"]<-1
				#	}
				#}
				
				##### cutoffCensored
				## 1. put minimum in protein level to NA
				#if(cutoffCensored=="minEachProtein"){
				#	if(censoredInt=="NA"){
				#		cut<-min(sub$ABUNDANCE, na.rm=TRUE) 
				#		sub[is.na(sub$INTENSITY),"ABUNDANCE"]<-cut
				#	} 
					
				#	if(censoredInt=="0"){
				#		cut<-min(sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,"ABUNDANCE"]) 
				#		sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"]<-cut
				#	}
				#}
				
				## 2. put minimum in feature level to NA
				if(cutoffCensored=="minFeature"){
					if(censoredInt=="NA"){
						cut<-aggregate(ABUNDANCE~FEATURE,data=sub, function(x) min(x, na.rm=TRUE))
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE<-0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$FEATURE))){
							sub[is.na(sub$INTENSITY) & sub$FEATURE==cut$FEATURE[j],"ABUNDANCE"]<-cut$ABUNDANCE[j]
						}
					}
					
					if(censoredInt=="0"){
						subtemptemp<-sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0 ,]
						cut<-aggregate(ABUNDANCE~FEATURE,data=subtemptemp, FUN=min)
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE<-0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$FEATURE))){
							sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$FEATURE==cut$FEATURE[j],"ABUNDANCE"]<-cut$ABUNDANCE[j]
						}
					}
				}
				
				## 3. put minimum in RUN to NA
				if(cutoffCensored=="minRun"){
					if(censoredInt=="NA"){
						cut<-aggregate(ABUNDANCE~RUN,data=sub, function(x) min(x, na.rm=TRUE))
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE<-0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$RUN))){
							sub[is.na(sub$INTENSITY) & sub$RUN==cut$RUN[j],"ABUNDANCE"]<-cut$ABUNDANCE[j]
						}
					}
					
					if(censoredInt=="0"){
						subtemptemp<-sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]
						cut<-aggregate(ABUNDANCE~RUN,data=subtemptemp, FUN=min)
						
						## cutoff for each Run is little less than minimum abundance in a run.
						cut$ABUNDANCE<-0.99*cut$ABUNDANCE

						for(j in 1:length(unique(cut$RUN))){
							sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$RUN==cut$RUN[j],"ABUNDANCE"]<-cut$ABUNDANCE[j]
						}
					}
				}	
				
				## 20150829 : 4. put minimum RUN and FEATURE
				if(cutoffCensored=="minFeatureNRun"){
					if(censoredInt=="NA"){
						
						## cutoff for each feature is little less than minimum abundance in a run.
						cut.fea<-aggregate(ABUNDANCE~FEATURE,data=sub, function(x) min(x, na.rm=TRUE))
						cut.fea$ABUNDANCE<-0.99*cut.fea$ABUNDANCE
												
						## cutoff for each Run is little less than minimum abundance in a run.

						cut.run<-aggregate(ABUNDANCE~RUN,data=sub, function(x) min(x, na.rm=TRUE))
						cut.run$ABUNDANCE<-0.99*cut.run$ABUNDANCE
						
						if(length(unique(sub$FEATURE))>1){
							for(j in 1:length(unique(sub$FEATURE))){
								for(k in 1:length(unique(sub$RUN))){
									# get smaller value for min Run and min Feature
									finalcut<-min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
								
									sub[is.na(sub$INTENSITY) & sub$FEATURE==cut.fea$FEATURE[j] & sub$RUN==cut.run$RUN[k],"ABUNDANCE"]<-finalcut
								}
							}
						}
							# if single feature, not impute
					}
					
					if(censoredInt=="0"){
						subtemptemp<-sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,]

						cut.fea<-aggregate(ABUNDANCE~FEATURE,data=subtemptemp, FUN=min)
						cut.fea$ABUNDANCE<-0.99*cut.fea$ABUNDANCE
						
						cut.run<-aggregate(ABUNDANCE~RUN,data=subtemptemp, FUN=min)
						cut.run$ABUNDANCE<-0.99*cut.run$ABUNDANCE
	
						if(length(unique(sub$FEATURE))>1){
							for(j in 1:length(unique(sub$FEATURE))){
								for(k in 1:length(unique(sub$RUN))){
									# get smaller value for min Run and min Feature
									finalcut<-min(cut.fea$ABUNDANCE[j],cut.run$ABUNDANCE[k])
								
									sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0 & sub$FEATURE==cut.fea$FEATURE[j] & sub$RUN==cut.run$RUN[k],"ABUNDANCE"]<-finalcut
								}
							}
						}else{ # single feature
					
							sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"]<-cut.fea$ABUNDANCE
							
						}
					}
				}
					
				
				## when number of measurement is less than df, error for fitting
				subtemp<-sub[!is.na(sub$ABUNDANCE),]
				countdf<-nrow(subtemp)<(length(unique(subtemp$FEATURE))+length(unique(subtemp$RUN))-1)
				
				### fit the model
				if(length(unique(sub$FEATURE))==1){
					fittest<-survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN,data=sub, dist='gaussian')
				}else{
					if(countdf){
						fittest<-survreg(Surv(ABUNDANCE, cen, type='left') ~ RUN,data=sub, dist='gaussian')
					}else{
						fittest<-survreg(Surv(ABUNDANCE, cen, type='left') ~ FEATURE+RUN,data=sub, dist='gaussian')
					}
				}
				
		
				sub.result<-data.frame(Protein=unique(sub$PROTEIN),RUN=rep(c(levels(sub$RUN)),1),LogIntensities=NA)

				# get the parameters
				cf <- summary(fittest)$coefficients

				# calculate sample quantification for all levels of sample
      			a=1	
          
       			for(j in 1:nlevels(sub$RUN)){
       				contrast.matrix<-rep(0,nlevels(sub$RUN))
        			contrast.matrix[j]<-1
          					contrast<-.make.contrast.run.quantification.Survival(fittest,contrast.matrix,sub, labeled=FALSE)
	
         			 sub.result[a,3]<-.estimableFixedQuantificationSurvival(cf,contrast)
         			 a=a+1
        		}

				result<-rbind(result, sub.result)
			}

			datamat = dcast( Protein ~ RUN, data=result, value.var='LogIntensities', keep=TRUE)
			datamat = melt(datamat, id.vars=c('Protein'))
			colnames(datamat)<-c('Protein','RUN','LogIntensities')
			result<-datamat
		}
		dataafterfit<-NULL	
	}
	
	###################################
	## final result
	finalout<-list(rqdata=result,ModelQC=dataafterfit, PredictedBySurvival=predAbundance)
	return(finalout)
}



##########################################################################################
## updated v3
.fit.quantification.run<-function(sub,singleFeature,singleSubject, TechReplicate,labeled,equalFeatureVar){
	
	if(!labeled){ ### label-free case
		## for single Feature, original value is the run quantification
		if(singleFeature){
			fit.full<-lm(ABUNDANCE ~ RUN , data = sub)
		}else{
			fit.full<-lm(ABUNDANCE ~ FEATURE + RUN , data = sub)
		}
		
	}else{ ### labeled-based case
		### update v3
		if(singleFeature){
			fit.full<-lm(ABUNDANCE ~ RUN+ref , data = sub)
		}else{ ## several subjects
			fit.full<-lm(ABUNDANCE ~ FEATURE+RUN+ref , data = sub)
		}
	}
	
	## make equal variance for feature : need to be updated
    if(!equalFeatureVar){
       fit.full<-.iter.wls.fit.model(data=sub,fit=fit.full,nrepeats=1)
    }
	
	return(fit.full)
}



########################################################
.fit.model.single<-function(contrast.matrix,data,scopeOfTechReplication,scopeOfBioReplication,TechReplicate,singleSubject,repeated,origGroup){
  
  ## input : output of run quantification
  	
	data2<-data
    data2$GROUP<-factor(data2$GROUP)
    data2$SUBJECT<-factor(data2$SUBJECT)
    
    # when subject is fixed, it is ok using lm function.
    
    # when single feature, consider technical replicates for time-course.
    
    if (scopeOfBioReplication=="restricted"){
      
      # case-control: impossible for GROUP+SUBJECT for case-control (parameter for the last SUBJECT is NA even though with technical replicates)
      if (!repeated){ ## case-control
        if(!TechReplicate | singleSubject){
          fit.full<-lm(ABUNDANCE ~ GROUP , data = data2)
        }else{
          fit.full<-lm(ABUNDANCE ~ GROUP + SUBJECT , data = data2)
          
        }
      }else{ ### repeated==TRUE, time-course
      	if(singleSubject){
          fit.full<-lm(ABUNDANCE ~ GROUP , data = data2)
      	}else{ ## no single subject
        	if(!TechReplicate){
          		fit.full<-lm(ABUNDANCE ~ GROUP+SUBJECT , data = data2)
        	}else{
          		fit.full<-lm(ABUNDANCE ~ GROUP+SUBJECT+GROUP:SUBJECT, data = data2) ## SUBJECT==SUBJECT_NESTED here
        	}
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
          message("*** error : results of Protein ", unique(data2$PROTEIN), " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")
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
         if(!TechReplicate | singleSubject){
           fit.full<-lm(ABUNDANCE ~ GROUP , data = data2)
         }else{
          fit.full<-lmer(ABUNDANCE ~ GROUP+(1|SUBJECT) , data = data2)
          df.full<-lm(ABUNDANCE ~ GROUP+SUBJECT , data = data2)$df.residual
         }
        
      }else{ ## time-course
      	if(singleSubject){
          fit.full<-lm(ABUNDANCE ~ GROUP , data = data2)
      	}else{ ## no single subject
        	if(!TechReplicate){
          	fit.full<-lmer(ABUNDANCE ~ GROUP+(1|SUBJECT) , data = data2)
          	df.full<-lm(ABUNDANCE ~ GROUP+SUBJECT , data = data2)$df.residual
        	}else{
          	fit.full<-lmer(ABUNDANCE ~ GROUP+(1|SUBJECT)+(1|GROUP:SUBJECT), data = data2) ## SUBJECT==SUBJECT_NESTED here
          	df.full<-lm(ABUNDANCE ~ GROUP+SUBJECT+GROUP:SUBJECT , data = data2)$df.residual
        	}
        }
      } ## time-course
      
      ## get parameter from model
      	if(class(fit.full)=="lm"){
			Para<-.getParameterFixed(fit.full)
		}else{
			Para<-.getParameterRandom(fit.full,df.full)
		}
      
      
      #### each comparison
      allout<-NULL
      
      for(k in 1:nrow(contrast.matrix)){
        
        # choose each comparison
        contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
        row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
        
        GroupComparisonAgreement<-.chechGroupComparisonAgreement(data2,contrast.matrix.sub)
        
        if(GroupComparisonAgreement$sign==TRUE){
          message("*** error : results of Protein ", unique(data2$PROTEIN), " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")
          
          out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
        }else{
          contrast<-.make.contrast.free.single(fit.full,contrast.matrix.sub,data2)
          out<-.estimableFixedRandom(Para,contrast)
          
          ## any error for out, just NA
          if(is.null(out)){
            out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
          }else{
            out<-data.frame(Protein=unique(data2$PROTEIN),Label=row.names(contrast.matrix.sub),out)	
          }
        }
        
        allout<-rbind(allout, out)
        
      } # end loop for comparion
      
    }	## expanded, random subject
    
  
  if(class(fit.full)=="lm"){  ### lm model
    finalresid<-fit.full$residuals
    finalfitted<-fit.full$fitted.values
  }else{   ### lmer model
    finalresid<-resid(fit.full)
    finalfitted<-fitted(fit.full)
  }
  
  finalout<-list(result=allout,valueresid=finalresid, valuefitted=finalfitted, fittedmodel=fit.full)	
  return(finalout)
  
} ## .fit.model.single

#############################################
.ttest.logsum<-function(contrast.matrix,data,origGroup){
	
	#### each comparison
      allout<-NULL
      
      for(k in 1:nrow(contrast.matrix)){
        
        # choose each comparison
        contrast.matrix.sub<-matrix(contrast.matrix[k,], nrow=1)
        row.names(contrast.matrix.sub)<-row.names(contrast.matrix)[k]
        
        #GroupComparisonAgreement<-.chechGroupComparisonAgreement(data,contrast.matrix.sub)
        
 #       if(GroupComparisonAgreement$sign==TRUE){
 #         message("*** error : results of Protein ", unique(data$PROTEIN), " for comparison ",row.names(contrast.matrix.sub), " are NA because measurements in Group ", origGroup[GroupComparisonAgreement$positionMiss], " are missing completely.")
          
 #         out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)		
 #       }else{
        	
        	## get two groups from contrast.matrix
        	datasub<-data[which(data$GROUP_ORIGINAL %in% origGroup[contrast.matrix.sub!=0]),]
          
          ## t test
          sumresult<-try(t.test(datasub$ABUNDANCE~datasub$GROUP_ORIGINAL, var.equal=TRUE),silent=TRUE)

			if(class(sumresult)=="try-error") {
				out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
			}else{
				out<-data.frame(Protein=unique(data$PROTEIN),Label=row.names(contrast.matrix.sub), logFC=NA,SE=NA,Tvalue=NA,DF=NA,pvalue=NA)
				
				out<-data.frame(Protein=unique(data$PROTEIN),Label=paste(names(sumresult$estimate)[1]," - ",names(sumresult$estimate)[2],sep=""), logFC=sumresult$estimate[1]-sumresult$estimate[2],SE=(sumresult$estimate[1]-sumresult$estimate[2])/sumresult$statistic, Tvalue=sumresult$statistic, DF=sumresult$parameter,pvalue=sumresult$p.value)
				rownames(out)<-NULL
			}

 #       }
        
        allout<-rbind(allout, out)
        
      } # end loop for comparion
      
      finalout<-list(result=allout,valueresid=NULL, valuefitted=NULL, fittedmodel=NULL)	
  return(finalout)

}




#############################################
#############################################
# Part 3 groupComparisonPlots
#############################################
#############################################


groupComparisonPlots<-function(data=data,type=type,sig=0.05,FCcutoff=FALSE,logBase.pvalue=10,ylimUp=FALSE,ylimDown=FALSE,xlimUp=FALSE,x.axis.size=10,y.axis.size=10,dot.size=3,text.size=4, legend.size=7,ProteinName=TRUE,ProteinNameLoc=1, numProtein=100, clustering="both", width=10, height=10, which.Comparison="all", address=""){
  
  ## save process output in each step
  allfiles<-list.files()
  filenaming<-"msstats"
  
  if(length(grep(filenaming,allfiles))==0){
    
    finalfile<-"msstats.log"
    processout<-NULL
    
  }else{
    
    num<-0
    finalfile<-"msstats.log"
    
    while(is.element(finalfile,allfiles)){
      num<-num+1
      lastfilename<-finalfile ## in order to rea
      finalfile<-paste(paste(filenaming,num,sep="-"),".log",sep="")
    }
    
    finalfile<-lastfilename
    processout<-as.matrix(read.table(finalfile, header=T, sep="\t"))
  }
  
  processout<-rbind(processout,as.matrix(c(" "," ","MSstats - groupComparisonPlots function"," "),ncol=1))
  
  
  ### make upper letter
  type<-toupper(type)
  
  if(length(setdiff(type,c("HEATMAP","VOLCANOPLOT","COMPARISONPLOT")))!=0){
    
    processout<-rbind(processout,c(paste("Input for type=",type,". However,'type' should be one of \"Heatmap\", \"VolcanoPlot\",\"ComparisonPlot\".",sep="")))
    write.table(processout, file=finalfile,row.names=FALSE)
    
    stop(paste("Input for type=",type,". However,'type' should be one of \"Heatmap\", \"VolcanoPlot\",\"ComparisonPlot\".",sep=""))
  }
  
  ##### check logBase.pvalue is 2,10 or not
  if(logBase.pvalue!=2 & logBase.pvalue!=10){
    processout<-rbind(processout,c("ERROR : (-) Logarithm transformation for adjusted p-values : log2 or log10 only - stop"))
    write.table(processout, file=finalfile,row.names=FALSE)
    
    stop("Only -log2 or -log10 for logarithm transformation for adjusted p-values are posssible.\n")
  }
  
  
  #### choose comparison to draw plots
  
  if(which.Comparison!="all"){
    ## check which.comparison is name of comparison
    if(is.character(which.Comparison)){
      
      temp.name<-which.Comparison
      
      ## message if name of comparison is wrong.
      if(length(setdiff(temp.name,unique(data$Label)))>0)
        
        processout<-rbind(processout,paste("Please check labels of comparions. Result does not have this comparison. -", paste(temp.name, collapse=", "),sep=" "))
      write.table(processout, file=finalfile,row.names=FALSE)
      
      stop(paste("Please check labels of comparions. Result does not have this comparison. -", paste(temp.name, collapse=", "),sep=" "))
    }
    
    
    ## check which.comparison is order number of comparison
    if(is.numeric(which.Comparison)){
      
      temp.name<-levels(data$Label)[which.Comparison]
      
      ## message if name of comparison is wrong.
      if(length(levels(data$Label))<max(which.Comparison))
        stop(paste("Please check your selection of comparisons. There are ", length(levels(data$Label))," comparisons in this result.",sep=" "))
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
  
  if (type=="HEATMAP"){
    
    #### check whether there is only one protein or not,
    if(length(unique(data$Protein))<=1){
      stop("At least two proteins are needed for heatmaps.")
    }
    
    #### check whether there is only one comparison or not,
    if(length(unique(data$Label))<=1){
      stop("At least two comparisons are needed for heatmaps.")
    }
    
    if(logBase.pvalue==2){
      y.limUp <-30
    }
    
    if(logBase.pvalue==10){
      y.limUp <-10
    }
    
    if(is.numeric(ylimUp)) y.limUp<-ylimUp 
    
    ## when NA, change it
    #data$adj.pvalue[is.na(data$adj.pvalue)]<-1 ## missing will be grey
    
    if(logBase.pvalue==2){
      data$adj.pvalue[data$adj.pvalue<2^(-y.limUp)]<-2^(-y.limUp)
    }
    
    if(logBase.pvalue==10){
      data$adj.pvalue[data$adj.pvalue<10^(-y.limUp)]<-10^(-y.limUp)
    }
    
    
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
      
      #if(logBase.pvalue==2){
      #  temp<- -log2(sub$adj.pvalue)*sign(sub[,3])
      #}
      
      #if(logBase.pvalue==10){
        temp<- -log10(sub$adj.pvalue)*sign(sub[,3])
      #}
      
      final<-	data.frame(cbind(final,temp))
    }
    
    obj<-final
    data$Protein<-factor(data$Protein)
    rownames(obj)<-levels(data$Protein)
    colnames(obj)<-levels(data$Label)
    
    ## remove if whole rows or columns are NA
    obj<-obj[rowSums(!is.na(obj))!=0, colSums(!is.na(obj))!=0]
    
    ### clustering for order
    tempobj<-obj
    tempobj[is.na(tempobj)]<-50
    
    if(toupper(clustering)=='PROTEIN'){
      obj<-obj[hclust(dist(tempobj),method="ward.D")$order,]
    }
    if(toupper(clustering)=='COMPARISON'){
      obj<-obj[hclust(dist(t(tempobj)),method="ward.D")$order,]
    }
    if(toupper(clustering)=='BOTH'){
      obj<-obj[hclust(dist(tempobj),method="ward.D")$order,hclust(dist(t(tempobj)),method="ward.D")$order]
    }
    if(toupper(clustering)=='NONE'){
      obj<-obj
    }
    
    rm(tempobj)
    
    ## change the order
    #obj$id<-seq(1:nrow(obj))
    #obj<-obj[order(obj$id,decreasing=TRUE),]
    #obj<-subset(obj, select=-c(id))
    
    ## color scale
    blue.red.18 <- maPalette(low = "blue", high = "red", mid = "black", k = 12)
    my.colors <-blue.red.18
    #my.colors[my.colors=="#FFFFFF"]<-"gold"
    my.colors<-c(my.colors,"grey") ## for NA
    
    if(logBase.pvalue==2){ up<-ceiling(-log10(2^-y.limUp )) }
    if(logBase.pvalue==10){ up<-ceiling(-log10(10^-y.limUp )) }
    
    temp<-10^(-sort(ceiling(seq(2,up,length=10)[c(1,2,3,5,10)]),decreasing = TRUE))
    breaks<-c(temp,sig)
    neg.breaks <- log(breaks, 10)
    my.breaks  <- c(neg.breaks,0,-neg.breaks[6:1],101)			
    
    ### draw color key
    blocks<-c(-breaks,1, breaks[6:1])
    x.at<-seq(-0.05,1.05,length.out=13)
    
    
    ### maximum number of proteins per heatmap
    namepro<-rownames(obj)
    totalpro<-length(namepro)
    numheatmap<-totalpro %/% numProtein +1
    
    
    # If there are the file with the same name, add next numbering at the end of file name
    if(address!=FALSE){
      allfiles<-list.files()
      
      num<-0
      filenaming<-paste(address,"Heatmap",sep="")
      finalfile<-paste(address,"Heatmap.pdf",sep="")
      
      while(is.element(finalfile,allfiles)){
        num<-num+1
        finalfile<-paste(paste(filenaming,num,sep="-"),".pdf",sep="")
      }	
      
      pdf(finalfile, width=width, height=height)
    }
    
    par(mar=c(3,3,3,3), mfrow=c(3,1),oma=c(3,0,3,0))
    plot.new()
    image(z = matrix(seq(1:(length(my.colors)-1)), ncol = 1), col = my.colors[-length(my.colors)], xaxt = "n", yaxt = "n")
    mtext("Color Key", side=3,line=1, cex=3)
    mtext("(sign) Adjusted p-value", side=1,line=3,at=0.5, cex=1.7)
    mtext(blocks,side=1,line=1,at=x.at, cex=1)
    
    ### draw heatmap
    
    ### loop for numProtein
    for(j in 1:numheatmap){
      
      ## get the number proteins needed
      if(j!=numheatmap){
        tempobj<-obj[((j-1)*numProtein+1):(j*numProtein),]
      }else{
        tempobj<-obj[((j-1)*numProtein+1):nrow(obj),]
      }
      
      par(oma=c(3,0,0,4))
      heatmap.2(as.matrix(tempobj),
                col=my.colors,
                Rowv=FALSE,
                Colv=FALSE,
                dendrogram="none",
                breaks=my.breaks,
                trace="none",
                na.color="grey", ## missing value will be grey
                cexCol=(x.axis.size/10),cexRow=(y.axis.size/10), # assign text.size as option
                key=FALSE,
                lhei=c(0.1,0.9),lwid=c(0.1,0.9)
      )
      
      
    } ## end loop for heatmap
    
    if(address!=FALSE) dev.off()
  }
  
  
  #######################
  ## VolcanoPlot
  #######################
  if (type=="VOLCANOPLOT"){
    
    # If there are the file with the same name, add next numbering at the end of file name		
    if(address!=FALSE){
      allfiles<-list.files()
      
      num<-0
      filenaming<-paste(address,"VolcanoPlot",sep="")
      finalfile<-paste(address,"VolcanoPlot.pdf",sep="")
      
      while(is.element(finalfile,allfiles)){
        num<-num+1
        finalfile<-paste(paste(filenaming,num,sep="-"),".pdf",sep="")
      }	
      
      pdf(finalfile, width=width, height=height)
    }
    
    if(logBase.pvalue==2){
      y.limUp <-30
    }
    
    if(logBase.pvalue==10){
      y.limUp <-10
    }
    
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
      
      if(logBase.pvalue==2){
        sub$adj.pvalue[sub$adj.pvalue<2^(-y.limUp)]<-2^(-y.limUp)
      }
      
      if(logBase.pvalue==10){
        sub$adj.pvalue[sub$adj.pvalue<10^(-y.limUp)]<-10^(-y.limUp)
      }
      
      sub<-as.data.frame(sub)
      
      ## ylimUp
      if(logBase.pvalue==2){
        y.limup<-ceiling(max(-log2(sub[!is.na(sub$adj.pvalue),"adj.pvalue"])))
        if(y.limup < (-log2(sig))) y.limup<- (-log2(sig)+1) ## for too small y.lim
      }
      
      if(logBase.pvalue==10){
        y.limup<-ceiling(max(-log10(sub[!is.na(sub$adj.pvalue),"adj.pvalue"])))
        if(y.limup < (-log10(sig))) y.limup<- (-log10(sig)+1) ## for too small y.lim
      }
      
      
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
      #setnames(subtemp,3,"logFC")	
      
      if(logBase.pvalue==2){
        subtemp$log2adjp<-(-log2(subtemp$adj.pvalue))
      }
      
      if(logBase.pvalue==10){
        subtemp$log10adjp<-(-log10(subtemp$adj.pvalue))
        
      }
      
      ### Plotting
      if(logBase.pvalue==2){
        ptemp<-ggplot()+geom_point(aes_string(x='logFC', y='log2adjp',color='colgroup',label='Protein'),size=dot.size, data=subtemp)+scale_colour_manual(values=c("black","blue","red"),limits=c("black","blue","red"),breaks=c("black","blue","red"),labels=c("No regulation","Down-regulated","Up-regulated"))+scale_y_continuous('-Log2 (adjusted p-value)', limit=c(y.limdown, y.limup))+labs(title=unique(sub$Label))
      }
      
      if(logBase.pvalue==10){
        ptemp<-ggplot()+geom_point(aes_string(x='logFC', y='log10adjp',color='colgroup',label='Protein'),size=dot.size, data=subtemp)+scale_colour_manual(values=c("black","blue","red"),limits=c("black","blue","red"),breaks=c("black","blue","red"),labels=c("No regulation","Down-regulated","Up-regulated"))+scale_y_continuous('-Log10 (adjusted p-value)', limit=c(y.limdown, y.limup))+labs(title=unique(sub$Label))
      }
      
      
      ## x-axis labeling
      if(colnames(sub)[3]=="log2FC") ptemp<-ptemp+scale_x_continuous('Log2 fold change',limit=c(-x.lim,x.lim))
      if(colnames(sub)[3]=="log10FC") ptemp<-ptemp+scale_x_continuous('Log10 fold change',limit=c(-x.lim,x.lim))
      
      ## add protein name
      if(ProteinName){
        if(logBase.pvalue==2){
          ptemp<-ptemp+geom_text(aes_string(x='logFC', y='log2adjp',label='Protein'),data=subtemp,color="black",guide="none",hjust=0.5-sign(sub[,3])*ProteinNameLoc, vjust=0.5,size=text.size)
        }
        
        if(logBase.pvalue==10){
          ptemp<-ptemp+geom_text(aes_string(x='logFC', y='log10adjp',label='Protein'),data=subtemp,color="black",guide="none",hjust=0.5-sign(sub[,3])*ProteinNameLoc, vjust=0.5,size=text.size)
        }
      } 
      
      ## For legend of linetype for cutoffs
      ## first assign line type
      ltypes<-c("type1"="twodash","type2"="dotted")
      
      ## cutoff lines, FDR only
      if(!FCcutoff){ 
        if(logBase.pvalue==2){
          sigcut<-data.frame(logFC=seq(-x.lim, x.lim, length.out=20), log2adjp=(-log2(sig)))
          
          pfinal<-ptemp+geom_line(data=sigcut,aes_string(x='logFC',y='log2adjp'), linetype="twodash", colour="darkgrey",size=0.6)+scale_linetype_manual(values=ltypes,labels=paste("Adj p-value cutoff(",sig,")",sep=""))
        }
        
        if(logBase.pvalue==10){
          sigcut<-data.frame(logFC=seq(-x.lim, x.lim, length.out=20), log10adjp=(-log10(sig)))
          
          pfinal<-ptemp+geom_line(data=sigcut,aes_string(x='logFC',y='log10adjp'), linetype="twodash", colour="darkgrey",size=0.6)+scale_linetype_manual(values=ltypes,labels=paste("Adj p-value cutoff(",sig,")",sep=""))
        }				
      }
      
      # cutoff lines, FDR and Fold change cutoff
      if(is.numeric(FCcutoff)){
        if(colnames(sub)[3]=="log2FC"){
          if(logBase.pvalue==2){
            
            ## three different lines
            sigcut<-data.frame(logFC=seq(-x.lim, x.lim, length.out=10), log2adjp=(-log2(sig)))
            FCcutpos<-data.frame(logFC=log2(FCcutoff), log2adjp=seq(y.limdown, y.limup, length.out=10))
            FCcutneg<-data.frame(logFC=(-log2(FCcutoff)), log2adjp=seq(y.limdown, y.limup, length.out=10))
            
            ## three lines, with order color first and then assign linetype manual
            pfinal<-ptemp+geom_line(data=sigcut,aes_string(x='logFC',y='log2adjp'),linetype="twodash",colour="darkgrey",size=0.6)+geom_line(data=FCcutpos,aes_string(x='logFC',y='log2adjp'),linetype='dotted',colour="darkgrey",size=0.6)+geom_line(data=FCcutneg,aes_string(x='logFC',y='log2adjp'), linetype='dotted',colour="darkgrey",size=0.6)+guides(colour=guide_legend(order=1),linetype=guide_legend(order=2))+scale_linetype_manual(values=ltypes,labels=c(paste("Adj p-value cutoff(",sig,")",sep=""),paste("Fold change cutoff(",FCcutoff,")",sep="")))
          }
          
          if(logBase.pvalue==10){
            
            ## three different lines
            sigcut<-data.frame(logFC=seq(-x.lim, x.lim, length.out=10), log10adjp=(-log10(sig)))
            FCcutpos<-data.frame(logFC=log2(FCcutoff), log10adjp=seq(y.limdown, y.limup, length.out=10))
            FCcutneg<-data.frame(logFC=(-log2(FCcutoff)), log10adjp=seq(y.limdown, y.limup, length.out=10))
            
            ## three lines, with order color first and then assign linetype manual
            pfinal<-ptemp+geom_line(data=sigcut,aes_string(x='logFC',y='log10adjp'),linetype="twodash",colour="darkgrey",size=0.6)+geom_line(data=FCcutpos,aes_string(x='logFC',y='log10adjp'),linetype='dotted',colour="darkgrey",size=0.6)+geom_line(data=FCcutneg,aes_string(x='logFC',y='log10adjp'), linetype='dotted',colour="darkgrey",size=0.6)+guides(colour=guide_legend(order=1),linetype=guide_legend(order=2))+scale_linetype_manual(values=ltypes,labels=c(paste("Adj p-value cutoff(",sig,")",sep=""),paste("Fold change cutoff(",FCcutoff,")",sep="")))
          }
          
        }
        
        if(colnames(sub)[3]=="log10FC"){
          if(logBase.pvalue==2){
            
            ## three different lines
            sigcut<-data.frame(logFC=seq(-x.lim, x.lim, length.out=10), log2adjp=(-log2(sig)))
            FCcutpos<-data.frame(logFC=log10(FCcutoff), log2adjp=seq(y.limdown, y.limup, length.out=10))
            FCcutneg<-data.frame(logFC=(-log10(FCcutoff)), log2adjp=seq(y.limdown, y.limup, length.out=10))
            
            ## three lines, with order color first and then assign linetype manual
            pfinal<-ptemp+geom_line(data=sigcut,aes_string(x='logFC',y='log2adjp') ,linetype="twodash",colour="darkgrey",size=0.6)+geom_line(data=FCcutpos,aes_string(x='logFC',y='log2adjp'),linetype='dotted',colour="darkgrey",size=0.6)+geom_line(data=FCcutneg,aes_string(x='logFC',y='log2adjp'),linetype="dotted",colour="darkgrey",size=0.6)+guides(colour=guide_legend(order=1),linetype=guide_legend(order=2))+scale_linetype_manual(values=ltypes,labels=c(paste("Adj p-value cutoff(",sig,")",sep=""),paste("Fold change cutoff(",FCcutoff,")",sep="")))
          }
          
          if(logBase.pvalue==10){
            
            ## three different lines
            sigcut<-data.frame(logFC=seq(-x.lim, x.lim, length.out=10), log10adjp=(-log10(sig)))
            FCcutpos<-data.frame(logFC=log10(FCcutoff), log10adjp=seq(y.limdown, y.limup, length.out=10))
            FCcutneg<-data.frame(logFC=(-log10(FCcutoff)), log10adjp=seq(y.limdown, y.limup, length.out=10))
            
            ## three lines, with order color first and then assign linetype manual
            pfinal<-ptemp+geom_line(data=sigcut,aes_string(x='logFC',y='log10adjp') ,linetype="twodash",colour="darkgrey",size=0.6)+geom_line(data=FCcutpos,aes_string(x='logFC',y='log10adjp'),linetype='dotted',colour="darkgrey",size=0.6)+geom_line(data=FCcutneg,aes_string(x='logFC',y='log10adjp'),linetype="dotted",colour="darkgrey",size=0.6)+guides(colour=guide_legend(order=1),linetype=guide_legend(order=2))+scale_linetype_manual(values=ltypes,labels=c(paste("Adj p-value cutoff(",sig,")",sep=""),paste("Fold change cutoff(",FCcutoff,")",sep="")))
          }
          
        }
      }
      
      pfinal<-pfinal+theme(
        panel.background=element_rect(fill='white', colour="black"),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=x.axis.size,colour="black"),
        axis.text.y=element_text(size=y.axis.size,colour="black"),
        axis.ticks=element_line(colour="black"),
        axis.title.x=element_text(size=x.axis.size+5,vjust=-0.4),
        axis.title.y=element_text(size=y.axis.size+5,vjust=0.3),
        title=element_text(size=x.axis.size+8,vjust=1.5),
        legend.position="bottom",
        legend.key=element_rect(fill='white',colour='white'),
        legend.text=element_text(size=legend.size),
        legend.title=element_blank())
      
      print(pfinal)
    } ## end-loop
    
    if(address!=FALSE) dev.off()
  }	
  
  #######################
  ## Comparison Plot
  #######################
  if (type=="COMPARISONPLOT"){
    
    
    datatemp<-data[!is.na(data$adj.pvalue),]
    datatemp$Protein<-factor(datatemp$Protein)
    setnames(datatemp,3,"logFC")
    
    # If there are the file with the same name, add next numbering at the end of file name		
    if(address!=FALSE){
      allfiles<-list.files()
      
      num<-0
      filenaming<-paste(address,"ComparisonPlot",sep="")
      finalfile<-paste(address,"ComparisonPlot.pdf",sep="")
      
      while(is.element(finalfile,allfiles)){
        num<-num+1
        finalfile<-paste(paste(filenaming,num,sep="-"),".pdf",sep="")
      }	
      
      pdf(finalfile, width=width, height=height)
    }
    
    for (i in 1:nlevels(datatemp$Protein)){
      
      sub<-datatemp[datatemp$Protein==levels(datatemp$Protein)[i],] 		
      #sub$ciw<-qt(1-sig/2,sub$DF)*sub$SE
      ## adjust for multiple comparison within protein
      sub$ciw<-qt(1-sig/(2*nrow(sub)),sub$DF)*sub$SE
      
      sub<-as.data.frame(sub)
      
      ## ylimUp
      y.limup<-ceiling(max(sub$logFC+sub$ciw))
      if(is.numeric(ylimUp)) y.limup<-ylimUp 
      
      ## ylimDown
      y.limdown<-floor(min(sub$logFC-sub$ciw))
      if(is.numeric(ylimDown)) y.limdown<-ylimDown 
      
      ptemp<-ggplot(aes_string(x='Label', y='logFC'), data=sub)+geom_errorbar(aes(ymax = logFC + ciw, ymin=logFC - ciw),data=sub, width=0.1,colour="red")+geom_point(size=dot.size,colour="darkred")+scale_x_discrete('Comparison')+geom_hline(yintercept=0, linetype="twodash", colour="darkgrey", size=0.6)+labs(title=levels(datatemp$Protein)[i])+theme(
        panel.background=element_rect(fill='white', colour="black"),
        panel.grid.major.y = element_line(colour="grey95"),
        panel.grid.minor.y = element_blank(),
        axis.text.x=element_text(size=x.axis.size,colour="black"),
        axis.text.y=element_text(size=y.axis.size,colour="black"),
        axis.ticks=element_line(colour="black"),
        axis.title.x=element_text(size=x.axis.size+5,vjust=-0.4),
        axis.title.y=element_text(size=y.axis.size+5,vjust=0.3),
        title=element_text(size=x.axis.size+8,vjust=1.5))
      
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


designSampleSize<-function(data=data,desiredFC=desiredFC,FDR=0.05,numSample=TRUE,power=0.9){
  
  # numPep=numPep,numTran=numTran
 
  labeled=FALSE
  scopeOfBioReplication="expanded"
  interference=TRUE
  equalFeatureVar=TRUE
  
  #nrepeats=3
  
  ## save process output in each step
  allfiles<-list.files()
  filenaming<-"msstats"
  
  
  if(length(grep(filenaming,allfiles))==0){
    
    finalfile<-"msstats.log"
    processout<-NULL
    
  }else{
    
    num<-0
    finalfile<-"msstats.log"
    
    while(is.element(finalfile,allfiles)){
      num<-num+1
      lastfilename<-finalfile ## in order to rea
      finalfile<-paste(paste(filenaming,num,sep="-"),".log",sep="")
    }
    
    finalfile<-lastfilename
    processout<-as.matrix(read.table(finalfile,header=T, sep="\t"))
  }
  
  processout<-rbind(processout,as.matrix(c(" "," ","MSstats - designSampleSize function"," "),ncol=1))
  
  
  
  
  ## check one TRUE or not
  #if( sum(isTRUE(numSample),isTRUE(numPep),isTRUE(numTran),isTRUE(power))!=1 ){
#    processout<-rbind(processout,c(paste("The required input - number of sample or features : Only one value should be TRUE. - stop")))
#    write.table(processout, file=finalfile, row.names=FALSE)
    
#    stop("One of (numSample, numPep, numTran, power) needs to be TRUE")
#  }
  
  ## labeled value
#  if(!(labeled==TRUE | labeled==FALSE) | !is.logical(labeled)){
#    processout<-rbind(processout,c(paste("The required input - labeled : 'labeled' value is wrong. - stop")))
#    write.table(processout, file=finalfile, row.names=FALSE)
    
#    stop("'labeled' must be one of TRUE or FALSE as a logical value.")
#  }\
  
  
  ## all input
#  processout<-rbind(processout,c(paste("number of sample = ",numSample,sep="")))
#  processout<-rbind(processout,c(paste("number of peptide per protein = ",numPep,sep="")))
#  processout<-rbind(processout,c(paste("number of transition per peptide = ", numTran,sep="")))
  processout<-rbind(processout,c(paste("Desired fold change = ",paste(desiredFC,collapse=" - "),sep="")))
  processout<-rbind(processout,c(paste("FDR = ",FDR,sep="")))
  processout<-rbind(processout,c(paste("Power = ", power,sep="")))
  
  write.table(processout, file=finalfile, row.names=FALSE)
  
  ## for label-free experiment
  if (!labeled){
    
    sigma.error<-NULL	
    VarComponent<-data.frame(Protein=seq(1,length(data)),Error=NA,Subject=NA,GroupBySubject=NA)
    
    for (i in 1:length(data)){
      
      # note: when run is fixed, we can obtain the same variance of error for both case-control and time course studies.
      
      fit.full<-data[[i]]
      
      ## if fit.full==NA (class(fit.full)=="try-error)
      if(is.null(fit.full)) {
        ## !!!!! but if we have NULL for last protein?
        next
        
      }else{
        
        ## get variance component
        
        if(class(fit.full)!="mer"){
          VarComponent[i,"Error"]<-summary(fit.full)$sigma^2
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
      }
      
    } ## end-loop
    
    #		VarComponent[is.na(VarComponent)]<-0	
    ## for label-free DDA, there are lots of missingness and lots of zero SE. So, remove NA SE.
    median.sigma.error<-median(VarComponent[,"Error"],na.rm=TRUE)
    if(sum(!is.na(VarComponent[,"GroupBySubject"]))>0){
      median.sigma.subject<-median(VarComponent[,"GroupBySubject"],na.rm=TRUE)
    }else{
      if(sum(!is.na(VarComponent[,"Subject"]))>0){
        median.sigma.subject<-median(VarComponent[,"Subject"],na.rm=TRUE)
      }else{
        median.sigma.subject<-0
      }
    }
    
    ###
    processout<-rbind(processout,c("Calculated variance component. - okay"))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    
    ### power calculation
    if(isTRUE(power)){
      delta<-log2(seq(desiredFC[1],desiredFC[2],0.025))
      desiredFC<-2^delta
      m0_m1=99
      t<-delta/sqrt(2*median.sigma.error/numSample+median.sigma.subject/numSample)
      #t<-delta/sqrt(2*median.sigma.error/numPep/numTran/numSample+median.sigma.subject/numSample)
      powerTemp<-seq(0,1,0.01)
      
      power<-NULL
      for(i in 1:length(t)){
        diff<-qnorm(powerTemp)+qnorm(1-powerTemp*FDR/(1+(1-FDR)*m0_m1)/2)-t[i]
        min(abs(diff),na.rm=TRUE)
        power[i]<-powerTemp[order(abs(diff))][1]
      }
      
      CV<-round((2*median.sigma.error/(numSample)+median.sigma.subject/numSample)/desiredFC,3)

      #CV<-round((2*median.sigma.error/(numSample*numPep*numTran)+median.sigma.subject/numSample)/desiredFC,3)
      
      ###
      processout<-rbind(processout,c("Power is calculated. - okay"))
      write.table(processout, file=finalfile, row.names=FALSE)
      
      out<-data.frame(desiredFC,numSample,FDR,power=power,CV)
      #out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power=power,CV)

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
        numSample<-round((2*median.sigma.error+median.sigma.subject)/aa,0)
        CV<-round((2*median.sigma.error/(numSample)+median.sigma.subject/numSample)/desiredFC,3)
 #       numSample<-round((2*median.sigma.error/numPep/numTran+median.sigma.subject)/aa,0)
 #       CV<-round((2*median.sigma.error/(numSample*numPep*numTran)+median.sigma.subject/numSample)/desiredFC,3)
        
        ###
        processout<-rbind(processout,c("The number of sample is calculated. - okay"))
        write.table(processout, file=finalfile, row.names=FALSE)
        
        out<-data.frame(desiredFC,numSample,FDR,power,CV)
        #out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power,CV)
        return(out)
      }
      
      ### Num Peptide calculation
#      if (isTRUE(numPep)){
#        delta<-log2(seq(desiredFC[1],desiredFC[2],0.025))
#        desiredFC<-2^delta
#        z_alpha<-qnorm(1-alpha/2)
 #       z_beta<-qnorm(power)
#        aa<-(delta/(z_alpha+z_beta))^2
#        numPep<-round((2*median.sigma.error/numSample/numTran+median.sigma.subject/numSample)/aa,0)
#        CV<-round((2*median.sigma.error/(numSample*numPep*numTran)+median.sigma.subject/numSample)/desiredFC,3)
        
        ###
#        processout<-rbind(processout,c("The number of peptide per protein is calculated. - okay"))
#        write.table(processout, file=finalfile, row.names=FALSE)
        
#        out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power,CV)
#        return(out)
#      }
      
      ### Num Transition calculation
#      if (isTRUE(numTran)){
#        delta<-log2(seq(desiredFC[1],desiredFC[2],0.025))
#        desiredFC<-2^delta
#        z_alpha<-qnorm(1-alpha/2)
#        z_beta<-qnorm(power)
 #       aa<-(delta/(z_alpha+z_beta))^2
 #       numTran<-round((2*median.sigma.error/numSample/numPep+median.sigma.subject/numSample)/aa,0)
 #       CV<-round((2*median.sigma.error/(numSample*numPep*numTran)+median.sigma.subject/numSample)/desiredFC,3)
        
        ###
#        processout<-rbind(processout,c("The number of transition per peptide is calculated. - okay"))
#        write.table(processout, file=finalfile, row.names=FALSE)
        
#        out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power,CV)
#        return(out)
#      }
     } # when power is numeric
  } ## label-free
  
  
  ## isotope labeled experiment
  if (labeled){
    
    sigma.error<-NULL	
    VarComponent<-data.frame(Protein=seq(1,length(data)),Error=NA,Subject=NA,GroupBySubject=NA)
    
    for (i in 1:length(data)){
      
      # note: when run is fixed, we can obtain the same variance of error for both case-control and time course studies.
      
      fit.full<-data[[i]]
      
      ## if fit.full==NA (class(fit.full)=="try-error)
      if(is.null(fit.full)) {
        next
        
      }else{
        
        if(class(fit.full)!="mer"){
          VarComponent[i,"Error"]<-summary(fit.full)$sigma^2
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
      }
    } ## end-loop
    
    ## label-based case, few of NA SE.
    VarComponent[is.na(VarComponent)]<-0	
    median.sigma.error<-median(VarComponent[,"Error"],na.rm=TRUE)
    if(sum(VarComponent[,"GroupBySubject"])>0){
      median.sigma.subject<-median(VarComponent[,"GroupBySubject"],na.rm=TRUE)
    }else{
      median.sigma.subject<-median(VarComponent[,"Subject"],na.rm=TRUE)
    }
    
    ###
    processout<-rbind(processout,c("Calculated variance component. - okay"))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    
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
      
      CV<-round((2*median.sigma.error/(numSample*numPep*numTran)+median.sigma.subject/numSample)/desiredFC,3)
      
      ###
      processout<-rbind(processout,c("Power is calculated. - okay"))
      write.table(processout, file=finalfile, row.names=FALSE)
      
      out<-data.frame(desiredFC,numSample,numPep,numTran,FDR,power=power,CV)
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
        write.table(processout, file=finalfile, row.names=FALSE)
        
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
        write.table(processout, file=finalfile, row.names=FALSE)
        
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
        write.table(processout, file=finalfile, row.names=FALSE)
        
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
#  if (length(unique(data$numPep))>1) index<-"numPep"
#  if (length(unique(data$numTran))>1) index<-"numTran"
  if (length(unique(data$power))>1) index<-"power"
  
  if (length(unique(data$numSample))==1 & length(unique(data$power))==1) index<-"numSample"
  
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
    legend("topright",c(paste("FDR is",unique(data$FDR),sep=" "),paste("Statistical power is",unique(data$power),sep=" ")),bty="n",cex=text.size)
  }
  
  # if (index=="numPep"){
    
    # with(data, {
      # plot(desiredFC,numPep,lwd=2,xlab="",ylab="",cex.axis=axis.size,type="l",xaxt="n")})
    # axis(1,at=seq(min(data$desiredFC),max(data$desiredFC),0.05),labels=seq(min(data$desiredFC),max(data$desiredFC),0.05),cex.axis=axis.size)
    # axis(3,at=seq(min(data$desiredFC),max(data$desiredFC),0.05),labels=data$CV[which(data$desiredFC%in%seq(min(data$desiredFC),max(data$desiredFC),0.05))],cex.axis=axis.size)
    # mtext("Coefficient of variation, CV",3,line=2.5,cex=lab.size)
    # mtext("Desired fold change",1,line=3.5,cex=lab.size)
    # mtext("Minimal number of peptides",2,line=2.5,cex=lab.size)
    # legend("topright",c(paste("Number of replicates is",unique(data$numSample),sep=" "),paste("Number of transitions is",unique(data$numTran),sep=" "),paste("FDR is",unique(data$FDR),sep=" "),paste("Statistical power is",unique(data$power),sep=" ")),bty="n",cex=text.size)
  # }
  
  # if (index=="numTran"){
    
    # with(data, {
      # plot(desiredFC,numTran,lwd=2,xlab="",ylab="",cex.axis=axis.size,type="l",xaxt="n")})
    # axis(1,at=seq(min(data$desiredFC),max(data$desiredFC),0.05),labels=seq(min(data$desiredFC),max(data$desiredFC),0.05),cex.axis=axis.size)
    # axis(3,at=seq(min(data$desiredFC),max(data$desiredFC),0.05),labels=data$CV[which(data$desiredFC%in%seq(min(data$desiredFC),max(data$desiredFC),0.05))],cex.axis=axis.size)
    # mtext("Coefficient of variation, CV",3,line=2.5,cex=lab.size)
    # mtext("Desired fold change",1,line=3.5,cex=lab.size)
    # mtext("Minimal number of transitions",2,line=2.5,cex=lab.size)
    # legend("topright",c(paste("Number of replicates is",unique(data$numSample),sep=" "),paste("Number of peptides is",unique(data$numPep),sep=" "),paste("FDR is",unique(data$FDR),sep=" "),paste("Statistical power is",unique(data$power),sep=" ")),bty="n",cex=text.size)
  # }
  
  if (index=="power"){
    
    with(data, {
      plot(desiredFC,power,lwd=2,xlab="",ylab="",cex.axis=axis.size,type="l",xaxt="n")})
    axis(1,at=seq(min(data$desiredFC),max(data$desiredFC),0.05),labels=seq(min(data$desiredFC),max(data$desiredFC),0.05),cex.axis=axis.size)
    mtext("Desired fold change",1,line=3.5,cex=lab.size)
    mtext("Power",2,line=2.5,cex=lab.size)
    legend("bottomright",c(paste("Number of replicates is",unique(data$numSample),sep=" "),paste("FDR is",unique(data$FDR),sep=" ")),bty="n",cex=text.size)
  }
  
  #	dev.off()
}


#############################################
#############################################
# Part 5 quantification
#############################################
#############################################


quantification<-function(data,type="Sample",format="matrix", scopeOfTechReplication="restricted", scopeOfBioReplication="restricted", interference=TRUE, missing.action="nointeraction",equalFeatureVar=TRUE){
  
  ## save process output in each step
  allfiles<-list.files()
  filenaming<-"msstats"
  
  
  if(length(grep(filenaming,allfiles))==0){
    
    finalfile<-"msstats.log"
    processout<-NULL
    
  }else{
    
    num<-0
    finalfile<-"msstats.log"
    
    while(is.element(finalfile,allfiles)){
      num<-num+1
      lastfilename<-finalfile ## in order to rea
      finalfile<-paste(paste(filenaming,num,sep="-"),".log",sep="")
    }
    
    finalfile<-lastfilename
    processout<-as.matrix(read.table(finalfile, header=T,sep="\t"))
  }
  
  processout<-rbind(processout,as.matrix(c(" "," ","MSstats - quantification function"," "),ncol=1))
  
  
  ## data format
  rawinput<-c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Run","Intensity")
  
  if(length(setdiff(toupper(rawinput),toupper(colnames(data))))==0){
    processout<-rbind(processout,c(paste("The required input - data : did not process from dataProcess function. - stop")))
    
    write.table(processout, file=finalfile, row.names=FALSE)
    
    stop("Please use 'dataProcess' first. Then use output of dataProcess function as input in groupComparison.")
  }
  
  ## check other options
  if(!(toupper(type)=="SAMPLE" | toupper(type)=="GROUP")){
    processout<-rbind(processout,c(paste("The required input - type : 'type' value is wrong. - stop")))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    stop("'type' must be one of \"Sample\" or \"Group\". ")
  } 
  
  if(!(toupper(format)=="MATRIX" | toupper(format)=="LONG")){
    processout<-rbind(processout,c(paste("The required input - format : 'format' value is wrong. - stop")))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    stop("'format' must be one of \"matrix\" or \"long\". ")
  } 
  
  ## all input
  processout<-rbind(processout,c(paste("type of quantification = ",type,sep="")))
  processout<-rbind(processout,c(paste("format of output = ",format,sep="")))
  
  write.table(processout, file=finalfile, row.names=FALSE)
  
  
  #################################
  ### how to handle missingness for endogenous
  
  data$PROTEIN<-factor(data$PROTEIN)
  
  data.l<-data[data$LABEL=="L",]
  data.h<-data[data$LABEL=="H",]
  
  missingPeptides<-.checkMissFeature(data.l)
  
  protein.list = tapply ( data.l$FEATURE, data.l$PROTEIN, function ( x ) unique ( as.character ( x ) ) )
  
  missing.results = sapply ( protein.list, function ( x ) any ( x %in% missingPeptides ) )
  
  
  
  ## Impute for missing endogenous intensity
  if ( toupper(missing.action) == "IMPUTE" ) { 
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
    
    #	data<-rbindlist(list(data.l,data.h))
    data<-rbind(data.l,data.h)
    
  }
  
  
  ## even though it is not the case, user can do with no interaction ( no else command)
  
  
  ## 
  if ( toupper(missing.action) == "REMOVE" ){
    data<-data[-which(data$FEATURE %in% missingPeptides),]
    
    message("The features that are missing intensities for an entire condition in Protein will be removed for fitting the model.")
    
  }
  
  
  ## for assigning interference
  missing.results = sapply ( protein.list, function ( x ) any ( x %in% missingPeptides ) )
  
  ##
  processout<-rbind(processout,c(paste("missing.action : ",missing.action," - okay",sep="")))
  write.table(processout, file=finalfile, row.names=FALSE)
  
  
  #################################
  ## Group quantification
  
  if (toupper(type)=="GROUP"){
    
    data<-data[!is.na(data$ABUNDANCE),]
    data$PROTEIN<-factor(data$PROTEIN)
    
    if(nlevels(data$LABEL)==2){	
      
      ## set output format
      data$GROUP<-factor(data$GROUP)
      data$GROUP_ORIGINAL<-factor(data$GROUP_ORIGINAL)
      
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
        
        repeated<-.checkRepeated(sub)
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
        
        
        ##### fit the model
        message(paste("Fitting the model for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
        
        
        fit<-try(.fit.quantification.group.labeled(sub,repeated,singleFeature, singleSubject, TechReplicate, MissGroupByFeature,MissRunByFeature,MissSubjectByGroup,UnequalSubject,scopeOfTechReplication,scopeOfBioReplication,remove.interaction,equalFeatureVar), silent=TRUE)
        
        ## extrac coefficient : 
        ## if we add lmer for random
        # if(class(fit)=="lm"){ cf <- summary.lm(fit)$coefficients
        #	}else{ cf <- as.matrix(fixef(fit)) }
        # but now only fixed,
        
        if(class(fit)=="try-error") {
          message("*** error : can't fit the model for ", levels(data$PROTEIN)[i])
          
          result<-rbind(result, sub.result)
          
          processout<-rbind(processout,c(paste("*** error : can't fit the model for ", levels(data$PROTEIN)[i],sep="")))
          write.table(processout, file=finalfile, row.names=FALSE)
          
        }else{
          if(class(fit)=="lm"){
            cf <- summary(fit)$coefficients
          }else{
            cf <- fixef(fit)
          }
          
          # calculate group quantification for all levels of group
          a=1
          
          for (j in 1:nlevels(sub$GROUP_ORIGINAL)){
            contrast.matrix<-rep(0,nlevels(sub$GROUP_ORIGINAL))
            contrast.matrix[j]<-1
            contrast<-.make.contrast.group.quantification(fit,contrast.matrix,sub)
            
            ## instead of result, use sub.result
            if(class(fit)=="lm"){
              sub.result[a,3]<-.estimableFixedQuantification(cf,contrast)
            }else{
              sub.result[a,3]<-.estimableRandomQuantification(cf,contrast)
            }
            a=a+1
          }
          contrast<-.make.contrast.group.quantification.reference(fit,contrast.matrix,sub)
          if(class(fit)=="lm"){
            sub.result[a,3]<-.estimableFixedQuantification(cf,contrast)
          }else{
            sub.result[a,3]<-.estimableRandomQuantification(cf,contrast)
          }
          a=a+1
          
          result<-rbind(result, sub.result)
          
          ##
          processout<-rbind(processout,c(paste("Finished to quantify for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")")))
          write.table(processout, file=finalfile, row.names=FALSE)
        }
        
      }	## end loop for each protein
    }	## label-based
    
    
    if(nlevels(data$LABEL)==1){	
      
      data$GROUP<-factor(data$GROUP)
      data$GROUP_ORIGINAL<-factor(data$GROUP_ORIGINAL)
      
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
        
        repeated<-.checkRepeated(sub)
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
        
        ##### fit the model
        message(paste("Fitting the model for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
        
        fit<-try(.fit.quantification.group.label.free(sub,repeated,singleFeature, singleSubject, TechReplicate, MissGroupByFeature,MissRunByFeature,MissSubjectByGroup,UnequalSubject,scopeOfTechReplication,scopeOfBioReplication,remove.interaction,equalFeatureVar), silent=TRUE)
        
        ## extrac coefficient : 
        ## if we add lmer for random
        # if(class(fit)=="lm"){ cf <- summary.lm(fit)$coefficients
        #	}else{ cf <- as.matrix(fixef(fit)) }
        # but now only fixed,
        
        
        if(class(fit)=="try-error") {
          message("*** error : can't fit the model for ", levels(data$PROTEIN)[i])
          
          result<-rbind(result, sub.result)
          
          processout<-rbind(processout,c(paste("*** error : can't fit the model for ", levels(data$PROTEIN)[i],sep="")))
          write.table(processout, file=finalfile, row.names=FALSE)
        }else{
          
          if(class(fit)=="lm"){
            cf <- summary(fit)$coefficients
          }else{
            cf <- fixef(fit)
          }
          
          # calculate group quantification for all levels of group
          a=1
          
          for(j in 1:nlevels(sub$GROUP_ORIGINAL)){
            contrast.matrix<-rep(0,nlevels(sub$GROUP_ORIGINAL))
            contrast.matrix[j]<-1
            contrast<-.make.contrast.group.quantification(fit,contrast.matrix,sub)
            
            if(class(fit)=="lm"){
              sub.result[a,3]<-.estimableFixedQuantification(cf,contrast)
            }else{
              sub.result[a,3]<-.estimableRandomQuantification(cf,contrast)
            }
            a=a+1
            
          }
          
          result<-rbind(result, sub.result)
          
          ##
          processout<-rbind(processout,c(paste("Finished to quantify for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")")))
          write.table(processout, file=finalfile, row.names=FALSE)
        }
        
      } ## end-loop for each protein	
    }	## label-free
    
    message("Preparing the output.")
    
    ## make new data.frame in case there are missing group
    for(k in 1:nrow(result)){	
      finalresult[finalresult$Protein==result[k,"Protein"] & finalresult$Condition==result[k,"Condition"],"LogIntensities"]<-result[k,"LogIntensities"]
    }
    
    
    ### output
    finalresult$Condition<-factor(finalresult$Condition, levels=rep(c(levels(data$GROUP_ORIGINAL))))
    
    finalmatrix<-as.data.frame.matrix(xtabs(LogIntensities~Protein+Condition, data=finalresult))
    
    ##
    message("Done.")
    
    processout<-rbind(processout,c("Finish group quantificiation - okay."))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    
    if(format=="long") return(finalresult)
    if(format=="matrix") return(finalmatrix)
  }	
  
  
  #################################	
  ###### sample quantification		
  if (toupper(type)=="SAMPLE"){
    
    data<-data[!is.na(data$ABUNDANCE),]
    data$PROTEIN<-factor(data$PROTEIN)
    
    if(nlevels(data$LABEL)==2){	
      
      data$SUBJECT<-factor(data$SUBJECT)
      data$SUBJECT_NESTED<-factor(data$SUBJECT_NESTED)
      
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
        
        repeated<-.checkRepeated(sub)
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
        
        ##### fit the model
        message(paste("Fitting the model for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
        
        
        fit<-try(.fit.quantification.sample.labeled(sub,repeated,singleFeature, singleSubject, TechReplicate, MissGroupByFeature,MissRunByFeature,MissSubjectByGroup,UnequalSubject,remove.interaction,scopeOfTechReplication,equalFeatureVar), silent=TRUE)
        
        ## extrac coefficient : 
        ## if we add lmer for random
        # if(class(fit)=="lm"){ cf <- summary.lm(fit)$coefficients
        #	}else{ cf <- as.matrix(fixef(fit)) }
        # but now only fixed,
        
        if(class(fit)=="try-error") {
          message("*** error : can't fit the model for ", levels(data$PROTEIN)[i])
          
          result<-rbind(result, sub.result)
          
          processout<-rbind(processout,c(paste("*** error : can't fit the model for ", levels(data$PROTEIN)[i],sep="")))
          write.table(processout, file=finalfile, row.names=FALSE)
          
        }else{
          
          if(class(fit)=="lm"){
            cf <- summary(fit)$coefficients
          }else{
            cf <- fixef(fit)
          }
          
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
            
            if(class(fit)=="lm"){
              sub.result[a,3]<-.estimableFixedQuantification(cf,contrast)
            }else{
              sub.result[a,3]<-.estimableRandomQuantification(cf,contrast)
            }
            a=a+1
          } ## end subject_nested
          
          contrast<-.make.contrast.subject.quantification.reference(fit,contrast.matrix,sub)
          if(class(fit)=="lm"){
            sub.result[a,3]<-.estimableFixedQuantification(cf,contrast)
          }else{
            sub.result[a,3]<-.estimableRandomQuantification(cf,contrast)
          }
          
          result<-rbind(result, sub.result)
          
          ##
          processout<-rbind(processout,c(paste("Finished a comparison for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")")))
          write.table(processout, file=finalfile, row.names=FALSE)
        }
        
      } ## end-loop
    } ##label-based	
    
    if(nlevels(data$LABEL)==1){
      
      data$SUBJECT<-factor(data$SUBJECT)
      data$SUBJECT_NESTED<-factor(data$SUBJECT_NESTED)
      
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
        
        repeated<-.checkRepeated(sub)
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
        
        ##### fit the model
        message(paste("Fitting the model for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
        
        fit<-try(.fit.quantification.sample.label.free(sub,singleFeature, singleSubject, TechReplicate, MissGroupByFeature,MissRunByFeature,MissSubjectByGroup,UnequalSubject,remove.interaction,equalFeatureVar), silent=TRUE)
        
        ## extrac coefficient : 
        ## if we add lmer for random
        # if(class(fit)=="lm"){ cf <- summary.lm(fit)$coefficients
        #	}else{ cf <- as.matrix(fixef(fit)) }
        # but now only fixed,
        
        if(class(fit)=="try-error") {
          message("*** error : can't fit the model for ", levels(data$PROTEIN)[i])
          
          result<-rbind(result, sub.result)
          
          processout<-rbind(processout,c(paste("*** error : can't fit the model for ", levels(data$PROTEIN)[i],sep="")))
          write.table(processout, file=finalfile, row.names=FALSE)
          
        }else{
          
          if(class(fit)=="lm"){
            cf <- summary(fit)$coefficients
          }else{
            cf <- fixef(fit)
          }
          
          # calculate sample quantification for all levels of sample
          a=1	
          
          for(j in 1:nlevels(sub$SUBJECT_NESTED)){
            contrast.matrix<-rep(0,nlevels(sub$SUBJECT_NESTED))
            contrast.matrix[j]<-1
            
            if(singleFeature){
              contrast<-.make.contrast.subject.quantification.single(fit,contrast.matrix,sub)
            }else{
              contrast<-.make.contrast.subject.quantification(fit,contrast.matrix,sub)
            }
            
            if(class(fit)=="lm"){
              sub.result[a,3]<-.estimableFixedQuantification(cf,contrast)
            }else{
              sub.result[a,3]<-.estimableRandomQuantification(cf,contrast)
            }
            a=a+1
          }
          result<-rbind(result, sub.result)
          
          ##
          processout<-rbind(processout,c(paste("Finished to quantify for protein ",unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")")))
          write.table(processout, file=finalfile, row.names=FALSE)
        }
        
      } ## end-loop for each protein	
    } ## label-free
    
    message("Preparing the output.")
    
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
    
    result$Subject<-as.character(result$Subject)
    finalresult$Subject<-as.character(finalresult$Subject)
    
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
    message("Done.")
    
    processout<-rbind(processout,c("Finish sample quantificiation - okay."))
    write.table(processout, file=finalfile, row.names=FALSE)
    
    
    if(format=="long") return(finalresultout)
    if(format=="matrix") return(finalmatrix)
  }## end sample quantification		
}


.fit.quantification.group.labeled<-function(sub,repeated,singleFeature, singleSubject, TechReplicate, MissGroupByFeature,MissRunByFeature,MissSubjectByGroup,UnequalSubject,scopeOfTechReplication,scopeOfBioReplication,remove.interaction,equalFeatureVar){
  
  if(singleFeature){
    #fit<-lm(ABUNDANCE ~ GROUP + RUN , data = sub)
    
    ### (1) fixed Subject, fixed Run
    if(scopeOfTechReplication=="restricted" & scopeOfBioReplication=="restricted"){
      
      # case-control
      if (!repeated){	
        if(TechReplicate){
          fit.full<-lm(ABUNDANCE ~ SUBJECT_NESTED + GROUP + RUN , data = sub)
        }else{
          fit.full<-lm(ABUNDANCE ~ GROUP + RUN , data = sub)
        }
      }else{ ### time-course
        if(TechReplicate){
          fit.full<-lm(ABUNDANCE ~ SUBJECT_NESTED + GROUP + RUN , data = sub)
        }else{
          fit.full<-lm(ABUNDANCE ~ SUBJECT + GROUP + RUN , data = sub)
        }
      }
    }
    
    ### (2) random Subject, fixed Run
    if(scopeOfTechReplication=="restricted" & scopeOfBioReplication=="expanded"){
      
      # case-control
      if (!repeated){	
        if(TechReplicate){
          fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT_NESTED) + GROUP + RUN , data = sub)	
        }else{
          message(paste(unique(data$PROTEIN),"This protein has single feature and no technical replicate. Therefore, we can't analysis with expanded scope of biolofical replication."))
        }
      }else{ ### time-course
        #if(TechReplicate){
        #	fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + RUN , data = sub)
        #}else{
        fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT) + GROUP + RUN , data = sub)
        #}
      }
    }
    
    ### (3) fixed subject, random Run
    if(scopeOfTechReplication=="expanded" & scopeOfBioReplication=="restricted"){
      
      # case-control
      if (!repeated){	
        if(TechReplicate){
          fit.full<-lmer(ABUNDANCE ~ SUBJECT_ORIGINAL + GROUP + (1|RUN) , data = sub)
        }else{
          fit.full<-lmer(ABUNDANCE ~ GROUP + (1|RUN) , data = sub)
        }
      }else{ ### time-course
        if(TechReplicate){
          fit.full<-lmer(ABUNDANCE ~ SUBJECT_ORIGINAL + SUBJECT_ORIGINAL:GROUP + GROUP + (1|RUN) , data = sub)
        }else{
          fit.full<-lmer(ABUNDANCE ~ SUBJECT_ORIGINAL + GROUP + (1|RUN) , data = sub)
        }
      }
    }
    
    ### (4) random subject, random Run
    if(scopeOfTechReplication=="expanded" & scopeOfBioReplication=="expanded"){
      
      # case-control
      if (!repeated){	
        if(TechReplicate){
          fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT_NESTED) + GROUP + (1|RUN) , data = sub)
        }else{
          fit.full<-lmer(ABUNDANCE ~ GROUP + (1|RUN) , data = sub)
        }
      }else{ ### time-course
        #if(TechReplicate){
        #	fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + (1|RUN) , data = sub)	
        #}else{
        fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT) + GROUP + (1|RUN) , data = sub)
        #}
      }
    }
    
  }else{ ## not single feature
    
    if (singleSubject){
      if(TechReplicate){
        fit.full<-lm(ABUNDANCE ~ FEATURE + GROUP + RUN + GROUP:FEATURE + RUN:FEATURE, data = sub)
      }else{
        fit.full<-lm(ABUNDANCE ~ FEATURE + GROUP + RUN, data = sub)
      }
    }else{
      if(scopeOfBioReplication=="restricted" & scopeOfTechReplication=="restricted"){
        
        ## (1.1) case-control
        if (!repeated){
          if(!remove.interaction){
            fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = sub)
          }else{
            fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + RUN, data = sub)
          }
        }else{  ## (1.2) time-course
          if(!remove.interaction){
            fit.full<-lm(ABUNDANCE ~ FEATURE  +  SUBJECT_NESTED + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = sub)
          }else{
            fit.full<-lm(ABUNDANCE ~ FEATURE  +  SUBJECT_NESTED  + GROUP + RUN, data = sub)
          }
        }	
      } ## fixed sub, run
      
      if(scopeOfBioReplication=="expanded"&scopeOfTechReplication=="restricted"){
        
        ## (2.1) case-control
        if (!repeated){
          if(!remove.interaction){	
            if(!MissGroupByFeature & !MissRunByFeature){
              
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = sub)
            }else{
              
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + RUN, data = sub)		
            }
          }else{ ##interference==FALSE
            
            fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + RUN, data = sub)
          }
        }else{  ## (2.2) time-course
          if(!remove.interaction){	
            if(!MissGroupByFeature & !MissRunByFeature){
              
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = sub)
            }else{
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + RUN, data = sub)
            }
          }else{  ##interference==FALSE
            fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP)  + GROUP + RUN, data = sub)
          }	
        }
      } # random sub, fix run
      
      if(scopeOfBioReplication=="restricted"&scopeOfTechReplication=="expanded"){
        
        ## (3.1) equal subject per group
        if (!UnequalSubject){
          if(!remove.interaction){
            if(!MissGroupByFeature & !MissRunByFeature){
              
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL + SUBJECT_ORIGINAL:GROUP + GROUP + (1|RUN) + FEATURE:GROUP +  (1|FEATURE:RUN), data = sub)
            }else{
              
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL + SUBJECT_ORIGINAL:GROUP + GROUP + (1|RUN), data = sub)
            }
          }else{  ## interference
            fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL + SUBJECT_ORIGINAL:GROUP + GROUP + (1|RUN), data = sub)
          }
        }else{ ## (3.2) unequal subject per group
          if(!remove.interaction){		
            if(!MissGroupByFeature & !MissRunByFeature){
              
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL + GROUP + (1|RUN) + FEATURE:GROUP +  (1|FEATURE:RUN), data = sub)
            }else{
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL  + GROUP + (1|RUN), data = sub)
            }			
          }else{ ## interference==FALSE
            fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL  + GROUP + (1|RUN), data = sub)
          }
        }
      } ## fix sub, random run
      
      if(scopeOfBioReplication=="expanded" & scopeOfTechReplication=="expanded"){
        
        ## (4.1) case-control
        if (!repeated){
          if(!remove.interaction){
            if(!MissGroupByFeature & !MissRunByFeature){
              
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + (1|RUN) + FEATURE:GROUP +  (1|FEATURE:RUN), data = sub)
            }else{
              
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED)  + GROUP + (1|RUN), data = sub)
            }
          }else{ ## interference==FALSE
            fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + (1|RUN), data = sub)
          }
        }else{  ## (4.2) time-course
          if(!remove.interaction){	
            if(!MissGroupByFeature & !MissRunByFeature){
              
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + (1|RUN) + FEATURE:GROUP +  (1|FEATURE:RUN), data = sub)
            }else{
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP) + GROUP + (1|RUN), data = sub)
            }
          }else{  ## interference==FALSE
            fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT)  + GROUP + (1|RUN), data = sub)
          }
        }	
      } ## random run, random sub
    } ## end single subject
  }
  
  ## make equal variance for feature
  if(!equalFeatureVar & scopeOfBioReplication=="restricted" & scopeOfTechReplication=="restricted"){
    nrepeats=3
    fit.full<-.iter.wls.fit.model(data=sub,fit=fit.full,nrepeats)
  }
  
  return(fit.full)
}


.fit.quantification.group.label.free<-function(sub,repeated,singleFeature, singleSubject, TechReplicate, MissGroupByFeature,MissRunByFeature,MissSubjectByGroup,UnequalSubject,scopeOfTechReplication,scopeOfBioReplication,remove.interaction,equalFeatureVar){
  
  ## for label-free : SUBJECT=SUBJECT_ORIGINAL
  if(singleFeature){
    #fit<-lm(ABUNDANCE ~ GROUP , data = sub)
    
    if (scopeOfBioReplication=="restricted"){
      if (!repeated){ ## case-control
        if(!TechReplicate){
          fit.full<-lm(ABUNDANCE ~ GROUP , data = sub)
        }else{
          fit.full<-lm(ABUNDANCE ~ SUBJECT_ORIGINAL + GROUP , data = sub)
        }
      }else{ ### repeated==TRUE, time-course
        if(!TechReplicate){
          fit.full<-lm(ABUNDANCE ~ SUBJECT_ORIGINAL + GROUP , data = sub)
        }else{
          fit.full<-lm(ABUNDANCE ~ SUBJECT_ORIGINAL + SUBJECT_ORIGINAL:GROUP + GROUP, data = sub) ## SUBJECT==SUBJECT_NESTED here
        }
      } ## time-course
    }
    
    if(scopeOfBioReplication=="expanded"){
      
      # case-control
      if (!repeated){
        #if(!TechReplicate){
        #	message("*** error : can't analyze with expanded scope of BioReplicates. Use the restricted scope of BioReplicates")
        #}else{
        fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT_ORIGINAL) + GROUP , data = sub)
        #}
      }else{ ## time-course
        #if(!TechReplicate){
        fit.full<-lmer(ABUNDANCE ~ (1|SUBJECT_ORIGINAL) + GROUP , data = sub)
        #}else{
        #	fit.full<-lmer(ABUNDANCE ~ GROUP+(1|SUBJECT)+(1|GROUP:SUBJECT), data = sub) ## SUBJECT==SUBJECT_NESTED here
        #}
      } ## time-course
    }
    
  }else{ ## not single feature
    
    if(singleSubject){
      if(!TechReplicate){
        fit.full<-lm(ABUNDANCE ~ FEATURE + GROUP , data = sub)
      }else{
        fit.full<-lm(ABUNDANCE ~ FEATURE + GROUP + GROUP:FEATURE , data = sub)
      }
      
    }else{ ## no random turn of run, can't do random subject
      if(scopeOfBioReplication=="restricted"){
        if(!remove.interaction){
          fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP, data = sub)
        }else{
          fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + FEATURE:GROUP, data = sub)
        }
      }
      
      if (scopeOfBioReplication=="expanded"){
        
        # case-control
        if (!repeated){
          if(!remove.interaction){
            if(!MissGroupByFeature){
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP + FEATURE:GROUP , data = sub)
            }else{  ## MissGroupByFeature==TRUE
              
              fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP, data = sub)
            }
          }else{ ## interference==FALSE
            fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT_NESTED) + GROUP, data = sub)
          }
        }else{  ##time-course
          #									if(interference){
          #										if(!MissGroupByFeature){
          #     										fit.full<-lmer(ABUNDANCE ~ FEATURE +  (1|SUBJECT) + (1|SUBJECT:GROUP)  + GROUP + FEATURE:GROUP , data = sub)
          #      									}else{ ## MissGroupByFeature==TRUE
          
          #     										fit.full<-lmer(ABUNDANCE ~ FEATURE +   (1|SUBJECT) + GROUP, data = sub)
          #   									}
          #								}else{  ## interference==FALSE
          fit.full<-lmer(ABUNDANCE ~ FEATURE +   (1|SUBJECT)  + GROUP, data = sub)
          #								}
        }
      } ## random sub
    }
  }
  
  if(!equalFeatureVar & scopeOfBioReplication=="restricted" & scopeOfTechReplication=="restricted"){
    nrepeats=3
    fit.full<-.iter.wls.fit.model(data=sub,fit=fit.full,nrepeats)
  }
  
  return(fit.full)
}

.fit.quantification.sample.labeled<-function(sub,repeated,singleFeature, singleSubject, TechReplicate, MissGroupByFeature,MissRunByFeature,MissSubjectByGroup,UnequalSubject,remove.interaction,scopeOfTechReplication,equalFeatureVar){
  
  if(singleFeature){
    
    ### (1) fixed Subject, fixed Run
    if(scopeOfTechReplication=="restricted" ){
      
      # case-control
      #				if (!repeated){	
      #					if(TechReplicate){
      #						fit.full<-lm(ABUNDANCE ~ SUBJECT_NESTED + RUN , data = sub)
      ## SUBJECT_NESTED+RUN = SUBJECT_NESTED for SUBJECT_NESTED value
      #					}else{
      fit.full<-lm(ABUNDANCE ~ SUBJECT_NESTED , data = sub)
      #					}
      #				}else{ ### time-course
      #					if(TechReplicate){
      #						fit.full<-lm(ABUNDANCE ~ SUBJECT_NESTED + RUN , data = sub)
      # same as GROUP+SUBJECT+SUBJECT:GROUP
      #					}else{
      #						fit.full<-lm(ABUNDANCE ~ SUBJECT_NESTED  , data = sub)
      #					}
      #				}
    }
    
    
    ### (3) fixed subject, random Run
    if(scopeOfTechReplication=="expanded" ){
      
      # case-control
      #				if (!repeated){	
      if(TechReplicate){
        fit.full<-lmer(ABUNDANCE ~ SUBJECT_NESTED + (1|RUN) , data = sub)
      }else{
        fit.full<-lm(ABUNDANCE ~ SUBJECT_NESTED , data = sub)
      }
      #				}else{ ### time-course
      #					if(TechReplicate){
      #						fit.full<-lmer(ABUNDANCE ~ SUBJECT_ORIGINAL + SUBJECT_ORIGINAL:GROUP + GROUP + (1|RUN) , data = sub)
      #					}else{
      #						fit.full<-lmer(ABUNDANCE ~ SUBJECT_ORIGINAL + GROUP + (1|RUN) , data = sub)
      #						fit.full<-lmer(ABUNDANCE ~ SUBJECT_NESTED + (1|RUN) , data = sub)
      #					}
      #				}
    }
    
    
  }else{ ## not single feature
    
    if (singleSubject){
      if(TechReplicate){
        if(scopeOfTechReplication=="restricted"){
          fit.full<-lm(ABUNDANCE ~ FEATURE + SUBJECT_NESTED + RUN + RUN:FEATURE, data = sub)
        }
        if(scopeOfTechReplication=="expanded"){
          fit.full<-lmer(ABUNDANCE ~ FEATURE + SUBJECT_NESTED + (1|RUN) + (1|RUN:FEATURE), data = sub)
        }
      }else{ 
        fit.full<-lm(ABUNDANCE ~ FEATURE + SUBJECT_NESTED, data = sub)
      }
    }else{
      if(scopeOfTechReplication=="restricted"){
        
        ## (1.1) case-control
        if (!repeated){
          if(!remove.interaction){
            fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + RUN + FEATURE:GROUP +  FEATURE:RUN, data = sub)
          }else{
            fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + RUN, data = sub)
          }
        }else{  ## (1.2) time-course
          if(!remove.interaction){
            fit.full<-lm(ABUNDANCE ~ FEATURE  +  SUBJECT_NESTED + RUN + FEATURE:GROUP +  FEATURE:RUN, data = sub)
          }else{
            fit.full<-lm(ABUNDANCE ~ FEATURE  +  SUBJECT_NESTED  + RUN, data = sub)
          }
        }	
      } ## fixed sub, run
      
      if(scopeOfTechReplication=="expanded"){
        
        ## (3.1) equal subject per group
        #					if (!UnequalSubject){
        #						if(!remove.interaction){
        #							if(!MissGroupByFeature & !MissRunByFeature){
        
        #								fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + (1|RUN) + FEATURE:GROUP +  (1|FEATURE:RUN), data = sub)
        #							}else{
        
        #								fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + (1|RUN), data = sub)
        #							}
        #						}else{  ## interference
        #							fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + (1|RUN), data = sub)
        #						}
        #					}else{ ## (3.2) unequal subject per group
        if(!remove.interaction){		
          if(!MissGroupByFeature & !MissRunByFeature){
            
            fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL + GROUP + (1|RUN) + FEATURE:GROUP +  (1|FEATURE:RUN), data = sub)
          }else{
            fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL  + GROUP + (1|RUN), data = sub)
          }			
        }else{ ## interference==FALSE
          fit.full<-lmer(ABUNDANCE ~ FEATURE +  SUBJECT_ORIGINAL  + GROUP + (1|RUN), data = sub)
        }
        #					}
      } ## fix sub, random run
      
    }
  }
  
  if(!equalFeatureVar & scopeOfTechReplication=="restricted"){
    nrepeats=3
    fit.full<-.iter.wls.fit.model(data=sub,fit=fit.full,nrepeats)
  }
  
  return(fit.full)
}


.fit.quantification.sample.label.free<-function(sub,singleFeature, singleSubject, TechReplicate, MissGroupByFeature,MissRunByFeature,MissSubjectByGroup,UnequalSubject,remove.interaction,equalFeatureVar){
  
  ## only for restricted subject (here is no run term.)
  if(singleFeature){
    
    #				if (!repeated){ ## case-control
    #					if(!TechReplicate){
    fit.full<-lm(ABUNDANCE ~ SUBJECT_NESTED , data = sub)
    #					}else{
    #						fit.full<-lm(ABUNDANCE ~ SUBJECT_NESTED , data = sub)
    #					}
    #				}else{ ### repeated==TRUE, time-course
    #					if(!TechReplicate){
    #						fit.full<-lm(ABUNDANCE ~ GROUP+SUBJECT , data = sub)
    #						fit.full<-lm(ABUNDANCE ~ SUBJECT_NESTED , data = sub)
    
    #					}else{
    #						fit.full<-lm(ABUNDANCE ~ GROUP+SUBJECT+GROUP:SUBJECT, data = sub) ## SUBJECT==SUBJECT_NESTED here
    #					}
    #				} ## time-course
  }else{ ## not single feature
    
    if(singleSubject){
      # when single subject, there is no need to distinguish case-control and time course; fixed subject or random subject, and no need to have GXF since there is not enough degree of freedom, hence no separatation of interference
      # GROUP=SUBJECT=SUBJECT_NESTED
      if(!TechReplicate){
        fit.full<-lm(ABUNDANCE ~ FEATURE + SUBJECT_NESTED , data = sub)
      }else{
        fit.full<-lm(ABUNDANCE ~ FEATURE + SUBJECT_NESTED + GROUP:FEATURE , data = sub)
      }
    }else{ ## no random turn of run, can't do random subject
      if(!remove.interaction){
        fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP, data = sub)
      }else{
        fit.full<-lm(ABUNDANCE ~ FEATURE +  SUBJECT_NESTED + GROUP + FEATURE:GROUP, data = sub)
      }			
    }
  }
  
  if(!equalFeatureVar){
    nrepeats=3
    fit.full<-.iter.wls.fit.model(data=sub,fit=fit.full,nrepeats)
  }
  
  return(fit.full)
  
}