
##################################
########### make contrast ########
##################################

##================================
## .make.contrast.free: 
## label-free: equal or unequal subjects per group
## fixed or random subject
##================================

.make.contrast.free<-function(fit,contrast.matrix,sub1){

	if(class(fit)=="lm"){
		coef.name<-names(coef(fit))
	}else{
		coef.name<-names(fixef(fit))
	}

#####
# intercept
#####
	temp<-coef.name[grep("Intercept",coef.name)]
	intercept.c<-rep(0,length(temp))
	names(intercept.c)<-temp
	if(length(temp)==0) intercept.c<-NULL

#####
# feature
#####
	temp<-coef.name[grep("FEATURE",coef.name)[!grep("FEATURE",coef.name)%in%grep(":",coef.name)]]
	feature.c<-rep(0,length(temp))
	names(feature.c)<-temp
	if(length(temp)==0) feature.c<-NULL

#####
# subject_nested 
#####
	temp<-coef.name[grep("SUBJECT_NESTED",coef.name)]
	if(length(temp)>0){
		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\.")),nrow=2))
		temp2<-as.vector(xtabs(~temp1[,1]))	
		
		tempdata<-fit$model
		levels<-unique(tempdata$GROUP)
		sub.contrast<-contrast.matrix[as.numeric(as.character(levels))]
		
		# the base is alway be the first SUBJECT_NESTED
		# (SUBJECT_NESTED1.1)
		temp3<-temp2
			if(length(temp2)==length(sub.contrast)){
				temp3[1]<-temp2[1]+1 ## this line first because if next is first,length of temp3 becomes >1
			}else{
				temp3<-c(1,temp3)
			}
		
		# subjectNested.c<-rep(contrast.matrix/(temp3),temp2) ## in case of unequal sample per group, wrong
		subjectNested.c<-rep(sub.contrast/(temp3),temp3)[-1]
		
		names(subjectNested.c)<-temp
	}
	if(length(temp)==0) subjectNested.c<-NULL

#####
# subject : for time-course
#####
#	temp<-coef.name[grep("SUBJECT",coef.name)[!grep("SUBJECT",coef.name)%in%c(grep(":",coef.name),grep("NESTED",coef.name))]]
#	if(length(temp)>0){
#		subject.c<-rep(0,length(temp))
#		names(subject.c)<-temp
#	}
#	if(length(temp)==0) subject.c<-NULL

	temp<-coef.name[grep("SUBJECT",coef.name)[!grep("SUBJECT",coef.name)%in%c(grep(":",coef.name),grep("NESTED",coef.name))]]
	
	if(length(temp)>0){
#		subject.c<-rep(0,length(temp))
#		names(subject.c)<-temp
	
	tempdata<-fit$model
	levels<-unique(tempdata$GROUP)
	labels<-paste("GROUP", levels, sep="")
	patients<-NULL
	for(i in 1:length(levels)){
		sub<-tempdata[tempdata$GROUP==levels[i],]
		sub.patients<-cbind(GROUP=paste("GROUP", levels[i], sep=""), SUBJECT=paste("SUBJECT", as.character(unique(sub$SUBJECT)), sep=""), Value=contrast.matrix[as.numeric(as.character(levels[i]))])
		patients<-data.frame(rbind(patients,sub.patients))
	}
	
	patient.count<-tapply(patients$SUBJECT, patients$GROUP, function(x) length(unique(x)))
	patient.seq<-rep(0, length(temp))
	for(i in 1:length(as.character(patients$SUBJECT))){
		match<-any(temp==as.character(patients$SUBJECT)[i])
		if(match & as.numeric(as.character(patients$Value[i]))!=0){
			res<-temp == as.character(patients$SUBJECT)[i]
			index<-which(res==TRUE)
			group<-as.character(patients[i,]$GROUP)
			count<-as.numeric(patient.count[names(patient.count)==group])
			value<-as.numeric(as.character(patients[i,]$Value))
			patient.value<-c(rep(0,index-1), value/count, rep(0, length(temp)-index))
		}else{
			patient.value<-rep(0,length(temp))
		}
		patient.seq<-patient.value+patient.seq
	}
	subject.c<-patient.seq
	names(subject.c)<-temp
	}
	
	if(length(temp)==0) subject.c<-NULL



#####
# subject by group : only for time-course - SUBJECT and GROUP (order) even GROUP:SUBJECT in model
#####
#	temp<-coef.name[intersect(grep("SUBJECT",coef.name),grep("GROUP",coef.name))]

#	tempSub<-unique(sub1[,c("GROUP","SUBJECT")])
#	tempSub1<-xtabs(~GROUP,data=tempSub)
#	tempSub2<-tempSub1[-1]
#	if(length(temp)>0){
#		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\:")),nrow=2))
#		temp2<-as.vector(xtabs(~temp1[,2]))	## count per GROUP
#		#gs.c<-rep(as.vector(contrast.matrix[-1]/(tempSub2)),temp2[1]) ## assume no missing for group and subject
		
#		# when Group completely  missing
#		sub.matrix<-contrast.matrix[unique(tempSub$GROUP)]
#		gs.c<-rep(as.vector(sub.matrix[-1]/(tempSub2)),each=temp2[1])
#		names(gs.c)<-temp
#	}
#	if(length(temp)==0) gs.c<-NULL

#####
# subject by group : only for time-course - SUBJECT and GROUP (order) even GROUP:SUBJECT in model
#####
	temp<-coef.name[intersect(grep("SUBJECT",coef.name),grep("GROUP",coef.name))]
	
	if(length(temp)>0){
#		subject.c<-rep(0,length(temp))
#		names(subject.c)<-temp
	
	tempdata<-fit$model
	levels<-unique(tempdata$GROUP)
	labels<-paste("GROUP", levels, sep="")
	patients<-NULL
	for(i in 1:length(levels)){
		sub<-tempdata[tempdata$GROUP==levels[i],]
		sub.patients<-cbind(GROUP=paste("GROUP", levels[i], sep=""), SUBJECT=paste("SUBJECT", as.character(unique(sub$SUBJECT)), sep=""), Value=contrast.matrix[as.numeric(as.character(levels[i]))])
		patients<-data.frame(rbind(patients,sub.patients))
	}
	
	patient.count<-tapply(patients$SUBJECT, patients$GROUP, function(x) length(unique(x)))
	interaction.seq<-rep(0, length(temp))
	interaction.labels<-paste(as.character(patients$GROUP),as.character(patients$SUBJECT), sep=":")
	
	for(i in 1:length(as.character(patients$SUBJECT))){
		match<-any(temp==interaction.labels[i])
		if(match & as.numeric(as.character(patients$Value[i]))!=0){
			res<-temp == interaction.labels[i]
			index<-which(res==TRUE)
			group<-as.character(patients[i,]$GROUP)
			count<-as.numeric(patient.count[names(patient.count)==group])
			value<-as.numeric(as.character(patients[i,]$Value))
			interaction.value<-c(rep(0,index-1), value/count, rep(0, length(temp)-index))
		}else{
			interaction.value<-rep(0,length(temp))
		}
		interaction.seq<-interaction.value+interaction.seq
	}
	gs.c<-interaction.seq
	names(gs.c)<-temp
	}
	
	if(length(temp)==0) gs.c<-NULL
	

#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]

### when there are some groups which are all missing
	tempSub<-as.numeric(as.character(unique(sub1[,c("GROUP")])))
	tempcontrast<-contrast.matrix[tempSub]

	group.c<-tempcontrast[-1] ## for label-free, need to remove first
	names(group.c)<-temp
	if(length(temp)==0) group.c<-NULL

####
# feature by group : different from labeled
#####

	temp<-coef.name[intersect(grep("GROUP",coef.name),grep("FEATURE",coef.name))]

	tempSub<-unique(sub1[,c("GROUP","FEATURE")])
	tempSub1<-xtabs(~GROUP,data=tempSub)
	tempSub2<-tempSub1[-1]
	if(length(temp)>0){
		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\:")),nrow=2))
		temp2<-as.vector(xtabs(~temp1[,2]))	
		gf.c<-rep(contrast.matrix[-1]/(tempSub2),temp2)
		names(gf.c)<-temp
	}

	if(length(temp)==0) gf.c<-NULL

	## be careful for order
	contrast<-c(intercept.c,feature.c,subjectNested.c,subject.c,group.c,gs.c,gf.c)

	if(class(fit)=="lm"){
		contrast1<-contrast[!is.na(coef(fit))]
	}else{
		contrast1<-contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}


##================================
## .make.contrast.based: 
## label-based: equal or unequal subjects per group
## fixed or random subject
##================================

.make.contrast.based<-function(fit,contrast.matrix,sub1){

	if(class(fit)=="lm"){
		coef.name<-names(coef(fit))
	}else{
		coef.name<-names(fixef(fit))
	}

#####
# intercept
#####
	temp<-coef.name[grep("Intercept",coef.name)]
	intercept.c<-rep(0,length(temp))
	names(intercept.c)<-temp
	if(length(temp)==0) intercept.c<-NULL

#####
# feature
#####
	temp<-coef.name[grep("FEATURE",coef.name)[!grep("FEATURE",coef.name)%in%grep(":",coef.name)]]
	feature.c<-rep(0,length(temp))
	names(feature.c)<-temp
	if(length(temp)==0) feature.c<-NULL

#####
# subject_nested
#####
	temp<-coef.name[grep("SUBJECT_NESTED",coef.name)]
	if(length(temp)>0){
		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\.")),nrow=2))
		temp2<-as.vector(xtabs(~temp1[,1]))	

		# the base is alway be the first SUBJECT_NESTED
		# (SUBJECT_NESTED0.0)
		temp3<-temp2
		# free:
		#temp3<-temp2
		#temp3[1]<-temp2[1]+1

		subjectNested.c<-rep(contrast.matrix/(temp3),temp2)
		names(subjectNested.c)<-temp
	}
	if(length(temp)==0) subjectNested.c<-NULL

#####
# subject
#####
	temp<-coef.name[grep("SUBJECT",coef.name)[!grep("SUBJECT",coef.name)%in%c(grep(":",coef.name),grep("NESTED",coef.name))]]
	if(length(temp)>0){
		subject.c<-rep(0,length(temp))
		names(subject.c)<-temp
	}
	if(length(temp)==0) subject.c<-NULL

#####
# subject by group
#####
	temp<-coef.name[intersect(grep("SUBJECT",coef.name),grep("GROUP",coef.name))]
	tempSub<-unique(sub1[,c("GROUP","SUBJECT")])
	tempSub1<-xtabs(~GROUP,data=tempSub)
	tempSub2<-tempSub1[-1]
	if(length(temp)>0){
		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\:")),nrow=2))
		temp2<-as.vector(xtabs(~temp1[,2]))	
		gs.c<-rep(as.numeric(contrast.matrix)/(tempSub2),temp2)
		names(gs.c)<-temp
	}
	if(length(temp)==0) gs.c<-NULL

#####
# subject_original_nested
#####
	temp<-coef.name[grep("SUBJECT_ORIGINAL_NESTED",coef.name)]
	if(length(temp)>0){
		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\.")),nrow=2))
		temp2<-as.vector(xtabs(~temp1[,1]))	

		# the base is alway be the first SUBJECT_ORIGINAL_NESTED
		# (SUBJECT_ORIGINAL_NESTED1.1)
		temp3<-temp2
		temp3[1]<-temp2[1]+1

		subjectOriginalNested.c<-rep(contrast.matrix/(temp3),temp2)
		names(subjectOriginalNested.c)<-temp
	}
	if(length(temp)==0) subjectOriginalNested.c<-NULL

#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]

	### when there are some groups which are all missing
	tempSub<-as.numeric(as.character(unique(sub1[,c("GROUP")])))
	group.c<-contrast.matrix[tempSub]

	# free
	#group.c<-contrast.matrix[-1]
	names(group.c)<-temp
	if(length(temp)==0) group.c<-NULL

#####
# run
#####
	temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	run.c<-rep(0,length(temp))
	names(run.c)<-temp
	if(length(temp)==0) run.c<-NULL

####
# feature by group : different from labeled
#####
	temp<-coef.name[intersect(grep("GROUP",coef.name),grep("FEATURE",coef.name))]

	tempSub<-unique(sub1[,c("GROUP","FEATURE")])
	tempSub1<-xtabs(~GROUP,data=tempSub)
	tempSub2<-tempSub1[-1]
	if(length(temp)>0){
		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\:")),nrow=2))
		temp2<-as.vector(xtabs(~temp1[,2]))	
		gf.c<-rep(as.numeric(contrast.matrix)/(tempSub2),temp2)
		# free
		#gf.c<-rep(contrast.matrix[-1]/(tempSub2),temp2)
		names(gf.c)<-temp
	}
	if(length(temp)==0) gf.c<-NULL

####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	if(length(temp)>0){
		rf.c<-rep(0,length(temp))
		names(rf.c)<-temp
	}
	if(length(temp)==0) rf.c<-NULL

	contrast<-c(intercept.c,feature.c,subjectNested.c,subject.c,subjectOriginalNested.c,group.c,run.c,gs.c,gf.c,rf.c)

	if(class(fit)=="lm"){
		contrast1<-contrast[!is.na(coef(fit))]
	}else{
		contrast1<-contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}



##================================
## .make.contrast.group.quantification: 
## label-based/label-free; single/multiple features
## all fixed subject and run
##================================

.make.contrast.group.quantification<-function(fit,contrast.matrix,sub1){

	if(class(fit)=="lm"){
		coef.name<-names(coef(fit))
	}else{
		coef.name<-names(fixef(fit))
	}

#####
# intercept
#####
	temp<-coef.name[grep("Intercept",coef.name)]
	intercept.c<-rep(1,length(temp))
	names(intercept.c)<-temp
	if(length(temp)==0) intercept.c<-NULL

#####
# feature
#####
	temp<-coef.name[grep("FEATURE",coef.name)[!grep("FEATURE",coef.name)%in%grep(":",coef.name)]]
	feature.c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
	names(feature.c)<-temp
	if(length(temp)==0) feature.c<-NULL

#####
# subject_nested
#####
	temp<-coef.name[grep("SUBJECT_NESTED",coef.name)]
	if(length(temp)>0){
		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\.")),nrow=2))
		temp2<-as.vector(xtabs(~temp1[,1]))	

		# the base is alway be the first SUBJECT_NESTED
		# (SUBJECT_NESTED0.0)
		if (nlevels(sub1$LABEL)==2){
			temp3<-temp2
			subjectNested.c<-rep(contrast.matrix/(temp3),temp2)
		}
		# free:
		if (nlevels(sub1$LABEL)==1){
			temp3<-temp2
			if(length(temp2)==length(contrast.matrix)){
				temp3[1]<-temp2[1]+1 ## this line first because if next is first,length of temp3 becomes >1
			}else{
				temp3<-c(1,temp3)
			}
			subjectNested.c<-rep(contrast.matrix/(temp3),temp3)[-1]
		}

		names(subjectNested.c)<-temp
	}
	if(length(temp)==0) subjectNested.c<-NULL

#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]

	if (nlevels(sub1$LABEL)==2){
		group.c<-contrast.matrix
	}
	# free
	if (nlevels(sub1$LABEL)==1){
		group.c<-contrast.matrix[-1]
	}
	names(group.c)<-temp
	if(length(temp)==0) group.c<-NULL

#####
# run
#####
	temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	run.c<-rep(1/nlevels(fit$model$RUN),length(temp))
	names(run.c)<-temp
	if(length(temp)==0) run.c<-NULL

####
# feature by group : different from labeled
#####
	temp<-coef.name[intersect(grep("GROUP",coef.name),grep("FEATURE",coef.name))]

	tempSub<-unique(sub1[,c("GROUP","FEATURE")])
	tempSub1<-xtabs(~GROUP,data=tempSub)
	tempSub2<-tempSub1[-1]
	if(length(temp)>0){
		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\:")),nrow=2))
		temp2<-as.vector(xtabs(~temp1[,2]))	

		if (nlevels(sub1$LABEL)==2){
			gf.c<-rep(as.numeric(contrast.matrix)/(tempSub2),temp2)
		}
		# free
		if (nlevels(sub1$LABEL)==1){
			gf.c<-rep(contrast.matrix[-1]/(tempSub2),temp2)
		}
		names(gf.c)<-temp
	}
	if(length(temp)==0) gf.c<-NULL

####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	tempSub<-dim(unique(sub1[,c("RUN","FEATURE")]))[1]
	if(length(temp)>0){
		rf.c<-rep(1/tempSub,length(temp))
		names(rf.c)<-temp
	}
	if(length(temp)==0) rf.c<-NULL

	contrast<-c(intercept.c,feature.c,subjectNested.c,group.c,run.c,gf.c,rf.c)

	if(class(fit)=="lm"){
		contrast1<-contrast[!is.na(coef(fit))]
	}else{
		contrast1<-contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}

##================================
## .make.contrast.group.quantification.reference: 
## label-based/label-free; single/multiple features
## all fixed subject and run
##================================

.make.contrast.group.quantification.reference<-function(fit,contrast.matrix,sub1){

	if(class(fit)=="lm"){
		coef.name<-names(coef(fit))
	}else{
		coef.name<-names(fixef(fit))
	}

#####
# intercept
#####
	temp<-coef.name[grep("Intercept",coef.name)]
	intercept.c<-rep(1,length(temp))
	names(intercept.c)<-temp
	if(length(temp)==0) intercept.c<-NULL

#####
# feature
#####
	temp<-coef.name[grep("FEATURE",coef.name)[!grep("FEATURE",coef.name)%in%grep(":",coef.name)]]
	feature.c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
	names(feature.c)<-temp
	if(length(temp)==0) feature.c<-NULL

#####
# subject_nested
#####
	temp<-coef.name[grep("SUBJECT_NESTED",coef.name)]
	if(length(temp)>0){
		subjectNested.c<-rep(0,length(temp))	
		names(subjectNested.c)<-temp
	}
	if(length(temp)==0) subjectNested.c<-NULL

#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]

	if(length(temp)>0){
		group.c<-rep(0,length(temp))	
		names(group.c)<-temp
	}
	if(length(temp)==0) group.c<-NULL

#####
# run
#####
	temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	run.c<-rep(1/nlevels(fit$model$RUN),length(temp))
	names(run.c)<-temp
	if(length(temp)==0) run.c<-NULL

####
# feature by group : different from labeled
#####
	temp<-coef.name[intersect(grep("GROUP",coef.name),grep("FEATURE",coef.name))]
	if(length(temp)>0){
		gf.c<-rep(0,length(temp))	
		names(gf.c)<-temp
	}
	if(length(temp)==0) gf.c<-NULL

####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	tempSub<-dim(unique(sub1[,c("RUN","FEATURE")]))[1]
	if(length(temp)>0){
		rf.c<-rep(1/tempSub,length(temp))
		names(rf.c)<-temp
	}
	if(length(temp)==0) rf.c<-NULL

	contrast<-c(intercept.c,feature.c,subjectNested.c,group.c,run.c,gf.c,rf.c)

	if(class(fit)=="lm"){
		contrast1<-contrast[!is.na(coef(fit))]
	}else{
		contrast1<-contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}


##================================
## .make.contrast.subject.quantification.single: 
## label-based/label-free; single features
## all fixed subject and run
##================================

.make.contrast.subject.quantification.single<-function(fit,contrast.matrix,sub1){

	if(class(fit)=="lm"){
		coef.name<-names(coef(fit))
	}else{
		coef.name<-names(fixef(fit))
	}

#contrastList<-unique(sub1[,c("GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
#contrastGroup<-rep(0,nlevels(sub1$GROUP_ORIGINAL))
#contrastGroup[as.numeric(contrastList[contrast.matrix==1,"GROUP_ORIGINAL"])]<-1

## for label-based	
	if (nlevels(sub1$LABEL)==2){
		contrastList<-unique(sub1[,c("GROUP","SUBJECT_ORIGINAL")])[-1,]## remove GROUP==0
		contrastList$GROUP<-factor(contrastList$GROUP) ## remove '0' group

		contrastGroup<-rep(0,(nlevels(sub1$GROUP)-1))
		contrastGroup[as.numeric(contrastList[contrast.matrix==1,"GROUP"])]<-1
	}else{ 
	## for label-free
		contrastList<-unique(sub1[,c("GROUP","SUBJECT_ORIGINAL")])

		contrastGroup<-rep(0,nlevels(sub1$GROUP))
		contrastGroup[as.numeric(contrastList[contrast.matrix==1,"GROUP"])]<-1
}


#####
# intercept
#####
	temp<-coef.name[grep("Intercept",coef.name)]
	intercept.c<-rep(1,length(temp))
	names(intercept.c)<-temp
	if(length(temp)==0) intercept.c<-NULL


#####
# feature
#####
	temp<-coef.name[grep("FEATURE",coef.name)[!grep("FEATURE",coef.name)%in%grep(":",coef.name)]]
	feature.c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
	names(feature.c)<-temp
	if(length(temp)==0) feature.c<-NULL

#####
# subject_nested
#####
	temp<-coef.name[grep("SUBJECT_NESTED",coef.name)]
	if(length(temp)>0){
		if (nlevels(sub1$LABEL)==2){
			subjectNested.c<-contrast.matrix
		}

		# free:
		if (nlevels(sub1$LABEL)==1){
			subjectNested.c<-contrast.matrix[-1]
		}
		names(subjectNested.c)<-temp
	}

	if(length(temp)==0) subjectNested.c<-NULL

#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]

	if (nlevels(sub1$LABEL)==2){
		group.c<-contrastGroup
	}
	# free
	if (nlevels(sub1$LABEL)==1){
		group.c<-contrastGroup[-1]
	}
	names(group.c)<-temp
	if(length(temp)==0) group.c<-NULL

#####
# run
#####
	temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]

#run.c<-rep(1/nlevels(fit$model$RUN),length(temp))
## when no technical replicate : subject_nested = run
	run.c<-contrast.matrix[-1]

## however, with technical replicate : subject_nested != run, need others

	names(run.c)<-temp
	if(length(temp)==0) run.c<-NULL

####
# feature by group : different from labeled
#####
	temp<-coef.name[intersect(grep("GROUP",coef.name),grep("FEATURE",coef.name))]

	tempSub<-unique(sub1[,c("GROUP","FEATURE")])
	tempSub1<-xtabs(~GROUP,data=tempSub)
	tempSub2<-tempSub1[-1]
	if(length(temp)>0){
		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\:")),nrow=2))
		temp2<-as.vector(xtabs(~temp1[,2]))	

		if (nlevels(sub1$LABEL)==2){
			gf.c<-rep(as.numeric(contrastGroup)/(tempSub2),temp2)
		}
		# free
		if (nlevels(sub1$LABEL)==1){
			gf.c<-rep(contrastGroup[-1]/(tempSub2),temp2)
		}
		names(gf.c)<-temp
	}
	if(length(temp)==0) gf.c<-NULL


####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	tempSub<-dim(unique(sub1[,c("RUN","FEATURE")]))[1]
	if(length(temp)>0){
		rf.c<-rep(1/tempSub,length(temp))
		names(rf.c)<-temp
	}
	if(length(temp)==0) rf.c<-NULL

	contrast<-c(intercept.c,feature.c,subjectNested.c,group.c,run.c,gf.c,rf.c)

	if(class(fit)=="lm"){
		contrast1<-contrast[!is.na(coef(fit))]
	}else{
		contrast1<-contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}


##================================
## .make.contrast.subject.quantification: 
## label-based/label-free; multiple features
## all fixed subject and run
##================================

.make.contrast.subject.quantification<-function(fit,contrast.matrix,sub1){

	if(class(fit)=="lm"){
		coef.name<-names(coef(fit))
	}else{
		coef.name<-names(fixef(fit))
	}

## when there are missing value in endogenous, there are error, because the number of group_original and fitted group are different
#contrastList<-unique(sub1[,c("GROUP_ORIGINAL","SUBJECT_ORIGINAL")])
#contrastGroup<-rep(0,nlevels(sub1$GROUP_ORIGINAL))
		#contrastGroup[as.numeric(contrastList[contrast.matrix==1,"GROUP_ORIGINAL"])]<-1

## for label-based	
	if (nlevels(sub1$LABEL)==2){
		contrastList<-unique(sub1[,c("GROUP","SUBJECT_ORIGINAL")])
		contrastList<-contrastList[contrastList$GROUP!="0",] ## remove GROUP==0
		contrastList$GROUP<-factor(contrastList$GROUP) ## remove '0' group

		contrastGroup<-rep(0,(nlevels(sub1$GROUP)-1))
		contrastGroup[as.numeric(contrastList[contrast.matrix==1,"GROUP"])]<-1
	}else{ 
	## for label-free
		contrastList<-unique(sub1[,c("GROUP","SUBJECT_ORIGINAL")])

		contrastGroup<-rep(0,nlevels(sub1$GROUP))
		contrastGroup[as.numeric(contrastList[contrast.matrix==1,"GROUP"])]<-1
	}

#####
# intercept
#####
	temp<-coef.name[grep("Intercept",coef.name)]
	intercept.c<-rep(1,length(temp))
	names(intercept.c)<-temp
	if(length(temp)==0) intercept.c<-NULL

#####
# feature
#####
	temp<-coef.name[grep("FEATURE",coef.name)[!grep("FEATURE",coef.name)%in%grep(":",coef.name)]]
	feature.c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
	names(feature.c)<-temp
	if(length(temp)==0) feature.c<-NULL

#####
# subject_nested
#####
	temp<-coef.name[grep("SUBJECT_NESTED",coef.name)]
	if(length(temp)>0){
		if (nlevels(sub1$LABEL)==2){
			subjectNested.c<-contrast.matrix
		}

		# free:
		if (nlevels(sub1$LABEL)==1){
			subjectNested.c<-contrast.matrix[-1]
		}
		names(subjectNested.c)<-temp
	}
	if(length(temp)==0) subjectNested.c<-NULL

#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]

	if (nlevels(sub1$LABEL)==2){
		group.c<-contrastGroup
	}
	# free
	if (nlevels(sub1$LABEL)==1){
		group.c<-contrastGroup[-1]
	}
	names(group.c)<-temp
	if(length(temp)==0) group.c<-NULL

#####
# run
#####
	temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	run.c<-rep(1/nlevels(fit$model$RUN),length(temp))
	names(run.c)<-temp
	if(length(temp)==0) run.c<-NULL

####
# feature by group : different from labeled
#####
	temp<-coef.name[intersect(grep("GROUP",coef.name),grep("FEATURE",coef.name))]

	tempSub<-unique(sub1[,c("GROUP","FEATURE")])
	tempSub1<-xtabs(~GROUP,data=tempSub)
	tempSub2<-tempSub1[-1]
	if(length(temp)>0){
		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\:")),nrow=2))
		temp2<-as.vector(xtabs(~temp1[,2]))	

		if (nlevels(sub1$LABEL)==2){
			gf.c<-rep(as.numeric(contrastGroup)/(tempSub2),temp2)
		}
		# free
		if (nlevels(sub1$LABEL)==1){
			gf.c<-rep(contrastGroup[-1]/(tempSub2),temp2)
		}
		names(gf.c)<-temp
	}
	if(length(temp)==0) gf.c<-NULL

####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	tempSub<-dim(unique(sub1[,c("RUN","FEATURE")]))[1]
	if(length(temp)>0){
		rf.c<-rep(1/tempSub,length(temp))
		names(rf.c)<-temp
	}
	if(length(temp)==0) rf.c<-NULL

	contrast<-c(intercept.c,feature.c,subjectNested.c,group.c,run.c,gf.c,rf.c)

	if(class(fit)=="lm"){
		contrast1<-contrast[!is.na(coef(fit))]
	}else{
		contrast1<-contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}



##================================
## .make.contrast.subject.quantification.reference: 
## label-based/label-free; single/multiple features
## all fixed subject and run
##================================

.make.contrast.subject.quantification.reference<-function(fit,contrast.matrix,sub1){

	if(class(fit)=="lm"){
		coef.name<-names(coef(fit))
	}else{
		coef.name<-names(fixef(fit))
	}

#####
# intercept
#####
	temp<-coef.name[grep("Intercept",coef.name)]
	intercept.c<-rep(1,length(temp))
	names(intercept.c)<-temp
	if(length(temp)==0) intercept.c<-NULL

#####
# feature
#####
	temp<-coef.name[grep("FEATURE",coef.name)[!grep("FEATURE",coef.name)%in%grep(":",coef.name)]]
	feature.c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
	names(feature.c)<-temp
	if(length(temp)==0) feature.c<-NULL

#####
# subject_nested
#####
	temp<-coef.name[grep("SUBJECT_NESTED",coef.name)]
	if(length(temp)>0){
		subjectNested.c<-rep(0,length(temp))	
		names(subjectNested.c)<-temp
	}
	if(length(temp)==0) subjectNested.c<-NULL

#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]

	if(length(temp)>0){
		group.c<-rep(0,length(temp))	
		names(group.c)<-temp
	}
	if(length(temp)==0) group.c<-NULL

#####
# run
#####
	temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	run.c<-rep(1/nlevels(fit$model$RUN),length(temp))
	names(run.c)<-temp
	if(length(temp)==0) run.c<-NULL

####
# feature by group : different from labeled
#####
	temp<-coef.name[intersect(grep("GROUP",coef.name),grep("FEATURE",coef.name))]
	if(length(temp)>0){
		gf.c<-rep(0,length(temp))	
		names(gf.c)<-temp
	}
	if(length(temp)==0) gf.c<-NULL

####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	tempSub<-dim(unique(sub1[,c("RUN","FEATURE")]))[1]
	if(length(temp)>0){
		rf.c<-rep(1/tempSub,length(temp))
		names(rf.c)<-temp
	}
	if(length(temp)==0) rf.c<-NULL

	contrast<-c(intercept.c,feature.c,subjectNested.c,group.c,run.c,gf.c,rf.c)

	if(class(fit)=="lm"){
		contrast1<-contrast[!is.na(coef(fit))]
	}else{
		contrast1<-contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}

##================================
## label-free, single
##================================
.make.contrast.free.single<-function(fit,contrast.matrix,sub1){

	if(class(fit)=="lm"){
		coef.name<-names(coef(fit))
	}else{
		coef.name<-names(fixef(fit))
	}

## change contrast.matrix without Group1
#	cons=1
#	if(contrast.matrix[1]==0) contrast.free<-contrast.matrix[-1]
#	if(contrast.matrix[1]<0){ 
#		contrast.free<-contrast.matrix[-1]/abs(contrast.matrix[1])
#		cons=abs(contrast.matrix[1])
#		}
#	if(contrast.matrix[1]>0){ 
#		contrast.free<-contrast.matrix[-1]*(-1)/contrast.matrix[1]
#		cons=contrast.matrix[1]*(-1)
#		}

#	if(class(fit)=="lm"){
#		coef.name<-names(coef(fit))
#	}else{
#		coef.name<-names(fixef(fit))
#	}

#####
# intercept
#####
	temp<-coef.name[grep("Intercept",coef.name)]
	intercept.c<-rep(0,length(temp))
	names(intercept.c)<-temp
	if(length(temp)==0) intercept.c<-NULL


#####
# subject
#####
	temp<-coef.name[grep("SUBJECT",coef.name)[!grep("SUBJECT",coef.name)%in%c(grep(":",coef.name),grep("NESTED",coef.name))]]
	
	if(length(temp)>0){
#		subject.c<-rep(0,length(temp))
#		names(subject.c)<-temp
	
	tempdata<-fit$model
	levels<-unique(tempdata$GROUP)
	labels<-paste("GROUP", levels, sep="")
	patients<-NULL
	for(i in 1:length(levels)){
		sub<-tempdata[tempdata$GROUP==levels[i],]
		sub.patients<-cbind(GROUP=paste("GROUP", levels[i], sep=""), SUBJECT=paste("SUBJECT", as.character(unique(sub$SUBJECT)), sep=""), Value=contrast.matrix[as.numeric(as.character(levels[i]))])
		patients<-data.frame(rbind(patients,sub.patients))
	}
	
	patient.count<-tapply(patients$SUBJECT, patients$GROUP, function(x) length(unique(x)))
	patient.seq<-rep(0, length(temp))
	for(i in 1:length(as.character(patients$SUBJECT))){
		match<-any(temp==as.character(patients$SUBJECT)[i])
		if(match & as.numeric(as.character(patients$Value[i]))!=0){
			res<-temp == as.character(patients$SUBJECT)[i]
			index<-which(res==TRUE)
			group<-as.character(patients[i,]$GROUP)
			count<-as.numeric(patient.count[names(patient.count)==group])
			value<-as.numeric(as.character(patients[i,]$Value))
			patient.value<-c(rep(0,index-1), value/count, rep(0, length(temp)-index))
		}else{
			patient.value<-rep(0,length(temp))
		}
		patient.seq<-patient.value+patient.seq
	}
	subject.c<-patient.seq
	names(subject.c)<-temp
	}
	
	if(length(temp)==0) subject.c<-NULL
	
	
#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]

### when there are some groups which are all missing
	tempSub<-as.numeric(as.character(unique(sub1[,c("GROUP")])))
	tempcontrast<-contrast.matrix[tempSub]

	group.c<-tempcontrast[-1] ## for label-free, need to remove first
	names(group.c)<-temp
	if(length(temp)==0) group.c<-NULL


#####
# subject by group : only for time-course - SUBJECT and GROUP (order) even GROUP:SUBJECT in model
#####
	temp<-coef.name[intersect(grep("SUBJECT",coef.name),grep("GROUP",coef.name))]
	
	if(length(temp)>0){
#		subject.c<-rep(0,length(temp))
#		names(subject.c)<-temp
	
	tempdata<-fit$model
	levels<-unique(tempdata$GROUP)
	labels<-paste("GROUP", levels, sep="")
	patients<-NULL
	for(i in 1:length(levels)){
		sub<-tempdata[tempdata$GROUP==levels[i],]
		sub.patients<-cbind(GROUP=paste("GROUP", levels[i], sep=""), SUBJECT=paste("SUBJECT", as.character(unique(sub$SUBJECT)), sep=""), Value=contrast.matrix[as.numeric(as.character(levels[i]))])
		patients<-data.frame(rbind(patients,sub.patients))
	}
	
	patient.count<-tapply(patients$SUBJECT, patients$GROUP, function(x) length(unique(x)))
	interaction.seq<-rep(0, length(temp))
	interaction.labels<-paste(as.character(patients$GROUP),as.character(patients$SUBJECT), sep=":")
	
	for(i in 1:length(as.character(patients$SUBJECT))){
		match<-any(temp==interaction.labels[i])
		if(match & as.numeric(as.character(patients$Value[i]))!=0){
			res<-temp == interaction.labels[i]
			index<-which(res==TRUE)
			group<-as.character(patients[i,]$GROUP)
			count<-as.numeric(patient.count[names(patient.count)==group])
			value<-as.numeric(as.character(patients[i,]$Value))
			interaction.value<-c(rep(0,index-1), value/count, rep(0, length(temp)-index))
		}else{
			interaction.value<-rep(0,length(temp))
		}
		interaction.seq<-interaction.value+interaction.seq
	}
	gs.c<-interaction.seq
	names(gs.c)<-temp
	}
	
	if(length(temp)==0) gs.c<-NULL
		
### combine all
	contrast<-c(intercept.c,group.c, subject.c, gs.c)

	if(class(fit)=="lm"){
		contrast1<-contrast[!is.na(coef(fit))]
	}else{
		contrast1<-contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}



##################################
##################################
########### estimate   ###########
##################################
##################################	
.getParameterFixed<-function(obj){
	temp1<-summary.lm(obj)
	cf <- temp1$coefficients
	vcv <- temp1$cov.unscaled * temp1$sigma^2	
	df <- obj$df.residual
	parameter<-list(cf=cf, vcv=vcv, df=df)
	
	return(parameter)
}

.getParameterRandom<-function(obj,df.full){
	cf <- as.matrix(fixef(obj))
	vcv <- as.matrix(vcov(obj))
	df <- df.full
	parameter<-list(cf=cf, vcv=vcv, df=df)
	
	return(parameter)
}

.estimableFixedRandom<-function(parameter,cm){
	cm <- matrix(cm, nrow=1)
	ct <- cm %*% parameter$cf[, 1]
	vc <- sqrt(diag(cm %*% parameter$vcv %*% t(cm)))
	prob <- 2 * (1 - pt(abs(ct/vc), parameter$df))
	result<-cbind(est=ct,stderr=vc,t=ct/vc,df=parameter$df,prob=prob)
	colnames(result)<-c("logFC","SE","Tvalue","DF","pvalue")

	return(result)
}

##########################################################################################

.estimableQuantification<-function(cf,cm){

	cm <- matrix(cm, nrow=1)
	ct <- cm %*% cf[, 1]
	result<-cbind(est=ct)
	colnames(result)<-c("log-intensities")
	
	return(result)
}


##########################################################################################

.iter.wls.fit.model<-function(data,fit,nrepeats){
	
	for(i in 1:nrepeats){	
		if(i==1){
			
			## lm or lmer
			if(class(fit)=="lm"){
				abs.resids<-data.frame(abs.resids=abs(fit$residuals))
				fitted<-data.frame(fitted=fit$fitted.values)
			}else{
				abs.resids<-data.frame(abs.resids=abs(resid(fit)))
				fitted<-data.frame(fitted=fitted(fit))
			}
	
			data<-data.frame(data,"abs.resids"=abs.resids,"fitted"=fitted)					
#			data<-merge(data,abs.resids,by="row.names",all=T)
#			rownames(data)<-data$Row.names
#			data<-merge(data, fitted, by="row.names",all=T)
#			rownames(data)<-data$Row.names
		}
		
		fit.loess<-loess(abs.resids~fitted, data=data)
		loess.fitted<-data.frame(loess.fitted=fitted(fit.loess))
		data<-data.frame(data,"loess.fitted"=loess.fitted)				

#		rownames(loess.fitted)<-names(resid(fit.loess))
#		data<-merge(data, loess.fitted, by="row.names",all=T)
#		rownames(data)<-data$Row.names
		
		## loess fitted valuaes are predicted sd
		data$weight<-1/(data$loess.fitted^2)
		
		data<-data[,-which(colnames(data) %in% "abs.resids")]
		
		## re-fit using weight
		
		if(class(fit)=="lm"){
			wls.fit<-lm(formula(fit),data=data,weights=weight)
		}else{
			wls.fit<-lmer(formula(fit),data=data,weights=weight)
		}

		## lm or lmer
#		if(class(fit)=="lm"){
#			residuals<-data.frame(residuals=wls.fit$residuals)
#		}else{
#			residuals<-data.frame(residuals=resid(wls.fit))
#		}
#		data<-merge(data,residuals,by="row.names",all=T)
#		rownames(data)<-data$Row.names
		
		## lm or lmer
		if(class(wls.fit)=="lm"){
			abs.resids<-data.frame(abs.resids=abs(wls.fit$residuals))
		}else{
			abs.resids<-data.frame(abs.resids=abs(resid(wls.fit)))
		}
		data<-data.frame(data,"abs.resids"=abs.resids)				

#		data<-merge(data,abs.resids,by="row.names",all=T)
#		rownames(data)<-data$Row.names
		
		data<-data[,-which(colnames(data) %in% c("loess.fitted","weight"))]
	}

	return(wls.fit)
}
