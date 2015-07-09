
##############################
## featureSelection
## (1) remove noisy feature
## (2) rank the feature
## (3) decision rule: minimize standard error



.featureSelection1<-function(work,lambda,eta, address){
	
	featureSelection<-data.frame(unique(work[,c("PROTEIN","PEPTIDE","FEATURE")]),perFeatureVariance=0,perProteinVariance=0,rank=NA,weights=NA,PRSS=NA,Decision="keep")
featureSelection$Decision<-as.vector(featureSelection$Decision)

	l<-subset(work,LABEL=="L")

	#=======================	
	##  generate per-feature weight

	for (i in 1:nlevels(l$PROTEIN)){
		sub<-subset(l,PROTEIN==levels(l$PROTEIN)[i])	
		sub<-within(sub,c(PEPTIDE<-factor(PEPTIDE),FEATURE<-factor(FEATURE)))

		# the sort of featue is based on sorting peptide and transition
		# it is the same as it is in feature selection and work
	
		# ABUNDANCE~FEATURE+SUBJECT_NESTED is the same as ABUNDANCE~GROUP+FEATURE+SUBJECT_NESTED
		# validated by SAS and R for per-feature variance

		# v10: with S(G) and not divide # of feature
		if (nlevels(sub$FEATURE)>1){
			contrasts(sub$FEATURE) = contr.sum(length(unique(sub$FEATURE)))
			fitSub<-lm(ABUNDANCE~GROUP+SUBJECT_NESTED+FEATURE,sub)
			featureSelection$perProteinVariance[featureSelection$PROTEIN==levels(l$PROTEIN)[i]]<-summary(fitSub)$sigma^2

			# since NA is removed from residual
			sub$res<-sub$ABUNDANCE
			sub$res[!is.na(sub$ABUNDANCE)]<-fitSub$residuals

			featureSelection$perFeatureVariance[featureSelection$PROTEIN==levels(l$PROTEIN)[i]]<-as.numeric(tapply(sub$res,sub[,c("FEATURE")], function(x) var(x,na.rm=TRUE)))

		}

		if (nlevels(sub$FEATURE)==1){
			contrasts(sub$FEATURE) = contr.sum(length(unique(sub$FEATURE)))      
			fitSub<-lm(ABUNDANCE~GROUP,sub)
			featureSelection$perProteinVariance[featureSelection$PROTEIN==levels(l$PROTEIN)[i]]<-summary(fitSub)$sigma^2
			featureSelection$perFeatureVariance[featureSelection$PROTEIN==levels(l$PROTEIN)[i]]<-summary(fitSub)$sigma^2
		}

		print(paste(i, 'out of ', nlevels(l$PROTEIN), 'proteins are done for computing the weight.', sep=' '))

	}    #End of the for loop over the proteins

	# up-to-here: run COV: add per-feature weights

	varCutoff<-min(c(quantile(featureSelection$perFeatureVariance,0.25,na.rm=TRUE),0.4))

	# add missingness per feature

	temp<-tapply(l$ABUNDANCE,l[,"FEATURE"],function(x) (1-mean(is.na(x),na.rm=TRUE)))
	temp1<-data.frame(FEATURE=names(temp),MissScore=as.numeric(temp))
	temp1$FEATURE<-factor(temp1$FEATURE,levels=unique(featureSelection$FEATURE))
	temp1<-temp1[with(temp1,order(FEATURE)),]
	featureSelection$MissScore<-temp1$MissScore

	# amount of penalty on missingness

	temp<-tapply(l[,"ABUNDANCE"],l[,"FEATURE"],function(x) mean(is.na(x),na.rm=TRUE))

	missPenalty<-ifelse(length(which(temp>0.5))>min(100,0.5*length(temp)),2,1)

	if(missPenalty==1){
		missPenalty<-ifelse(length(which(temp>0.5))>0.25*length(temp),1,0)
	}


## 20140514 : temporary ignore identificatino score.
#	if(IdentificationScore==TRUE){
	
#		temp<-unique(l[,c("PROTEIN","FEATURE","RUN","IDENTIFICATIONSCORE")])

#		# v13: select the worst score (the largest m_score) of a peptide across runs

#		temp1<-tapply(temp$IDENTIFICATIONSCORE,temp[,c("FEATURE")],function(x) min(x,na.rm=TRUE))
#		temp2<-data.frame(FEATURE=names(temp1),IdScore=as.numeric(temp1))

#		temp2$FEATURE<-factor(temp2$FEATURE,levels=unique(featureSelection$FEATURE))
#		temp2<-temp2[with(temp2,order(FEATURE)),]
#		featureSelection$IdScore<-temp2$IdScore

#		#featureSelection$weights<-ifelse((featureSelection$perProteinVariance/featureSelection$perFeatureVariance*featureSelection$IdScore)>3,3,(featureSelection$perProteinVariance/featureSelection$perFeatureVariance*featureSelection$IdScore))

#		featureSelection$weights<-featureSelection$perProteinVariance/ifelse(featureSelection$perFeatureVariance/(featureSelection$IdScore*featureSelection$MissScore^missPenalty)<varCutoff,varCutoff,featureSelection$perFeatureVariance/(featureSelection$IdScore*featureSelection$MissScore^missPenalty))

#	}else{

		featureSelection$weights<-featureSelection$perProteinVariance/ifelse(featureSelection$perFeatureVariance/featureSelection$MissScore^missPenalty<varCutoff,varCutoff,featureSelection$perFeatureVariance/featureSelection$MissScore^missPenalty)

#	}

	# force all the single feature protein has weights = 1
	temp1<-unique(l[,c("PROTEIN","FEATURE")])
	temp2<-xtabs(~PROTEIN,temp1)
	featureSelection[featureSelection$PROTEIN%in%names(temp2)[temp2==1],"weights"]<-1


	# #=======================	
	## (2.1) rank the feature 
	## based on perFeatureVariance
	
	featureSelection <- featureSelection[order(featureSelection$PROTEIN),]

	featureSelection[,"rank"]<-unlist(tapply(-featureSelection$weights,featureSelection[,c("PROTEIN")],function(x) rank(x,ties.method = c("min"))))

	# #=======================	
	## (2.2) save featureSelection results

	#write.csv(featureSelection,file=paste(address,"FeatureSelectionTable.csv",sep=""))
	#save(featureSelection,file=paste(address,"FeatureSelectionTable.RData",sep=""))


	# #=======================	
	# #=======================	
	# # Define PRSS: control based on Beta Feature and PRSS (without L1-penalty)
	# #=======================	
	# #=======================	
	# #=======================	

	# #=======================	
	## (2.1) different lambda control different Beta_Feature, PRSS is calculate based on the rank of weights

	temp1<-data.frame(Lambda=rep(lambda,each=length(featureSelection$FEATURE)),FEATURE=rep(featureSelection$FEATURE,length(lambda)))

	featureSelectionOut<-merge(featureSelection,temp1,by="FEATURE")	
	featureSelectionOut<-featureSelectionOut[with(featureSelectionOut,order(Lambda,FEATURE)),]
	#Mike: Don't want the duplicated rows
	featureSelectionOut <- featureSelectionOut[!duplicated(featureSelectionOut), ] 

	featureSelectionOut$PRSSprot<-NA


	# calculate PRSS

	temp<-featureSelection[,c("FEATURE","weights","rank")]
	l<-merge(l,temp,by="FEATURE",sort=FALSE)

	for (i in 1:nlevels(l$PROTEIN)){
		
		sub<-subset(l,PROTEIN==levels(l$PROTEIN)[i])	
		sub$FEATURE<-factor(sub$FEATURE)

		if(nlevels(sub$FEATURE)>1){
		# estimate Beta_Feature
		# since NA is removed from residuals
			sub$residual<-sub$ABUNDANCE
			sub$residual[!is.na(sub$ABUNDANCE)]<-lm(ABUNDANCE~GROUP+SUBJECT_NESTED,data=sub,weights=weights)$residuals	

			X<-model.matrix(~FEATURE,data=sub,contrasts=list(FEATURE="contr.sum"))

			X<-X[,-1]
			a<-lm(residual~X,data=sub)     #Same as lm(formula = residual ~ FEATURE, data = sub, contrasts = list(FEATURE = "contr.sum"))
			beta<-c(coef(a)[-1],0-sum(coef(a)[-1],na.rm=TRUE))    #Due to the zero-sum constraint
			feature.rank<-subset(featureSelection,PROTEIN==levels(l$PROTEIN)[i])$rank
		}

		if(nlevels(sub$FEATURE)==1){
			sub$residual<-sub$ABUNDANCE
			sub$residual[!is.na(sub$ABUNDANCE)]<-lm(ABUNDANCE~GROUP,data=sub,weights=weights)$residuals	

			feature.rank<-1
		}


		for (k in 1:length(lambda)){
			
			res.L<-lambda[k]

			if(nlevels(sub$FEATURE)>1){
				update.beta<-sign(beta)*ifelse(abs(beta)-res.L>0,abs(beta)-res.L,0)
				accept.rank<-sum(update.beta!=0,na.rm=TRUE)
			}

			if(nlevels(sub$FEATURE)==1){
				accept.rank<-1
			}

			# gurantee has one feature
			if(accept.rank==0){accept.rank<-1}

			# based on the feature.rank and select the number of non-0
			keep.feature<-levels(sub$FEATURE)[which(feature.rank<=accept.rank)]

			# new dataset with only keep.feature
			sub.new<-sub[sub$FEATURE%in%keep.feature,]
			sub.new$FEATURE<-factor(sub.new$FEATURE)

			if(nlevels(sub.new$FEATURE)>1){

				contrasts(sub.new$FEATURE) = contr.sum(length(unique(sub.new$FEATURE)))
				b<-lm(residual~FEATURE,data=sub.new)

				featureSelectionOut[featureSelectionOut$PROTEIN==levels(l$PROTEIN)[i]&featureSelectionOut$Lambda==res.L,"PRSS"]<-anova(b)$Mean[2]

			}

			if(nlevels(sub.new$FEATURE)=="1"){
		
				#featureSelectionOut[featureSelectionOut$PROTEIN==levels(l$PROTEIN)[i]&featureSelectionOut$Lambda==res.L,"PRSS"]<-sum((sub.new$residual-mean(sub.new$residual,na.rm=TRUE))^2,na.rm=TRUE)/(dim(sub.new[!is.na(sub.new$residual),])[1]-1)

				b2 <- lm(residual~1,data=sub.new)
				featureSelectionOut[featureSelectionOut$PROTEIN==levels(l$PROTEIN)[i]&featureSelectionOut$Lambda==res.L,"PRSS"]<-anova(b2)$Mean[1]

			}

			featureSelectionOut[featureSelectionOut$PROTEIN==levels(l$PROTEIN)[i]&featureSelectionOut$Lambda==res.L,"Decision"][-which(feature.rank<=accept.rank)]<-"remove"


		}

		print(paste(i, 'out of ', nlevels(l$PROTEIN), 'proteins are done for relarization for every lambda.', sep=' '))
		
	}    #End of the for loop of protein


	featureSelectionOut$PRSSprot<-featureSelectionOut$PRSS
		
	write.csv(featureSelectionOut,file=paste(address,"FeatureSelectionOutTable.csv",sep=""))
	save(featureSelectionOut,file=paste(address,"FeatureSelectionOutTable.RData",sep=""))

	# (1) select the best lambda #(1) allProtein - lambda
	#bestLambda<-as.numeric(names(sort(tapply(featureSelectionOut$PRSSprot,featureSelectionOut[,c("Lambda")],function(x) mean(x))))[1])

	temp<-tapply(featureSelectionOut$PRSSprot,featureSelectionOut[,c("Lambda")],function(x) median(x,na.rm=TRUE))
	temp1<-abs(diff(temp))
	d <- median(abs(diff(lambda)))   #Mike consider the slope since ratio is affected by outliers
	ratio <- temp1/d 
	#ratio<-temp1/sum(temp1,na.rm=TRUE)

	# criteria: has the smallest descreasing and the absolute cutoff of 0.1
	#eta<-0.1
	#eta<-0.2
	set1<-which(ratio<eta)
	set2<-which(rank(temp,ties.method="min")==1)[1]
	sign<-sum(set1%in%ifelse(set2==length(temp),(set2-1),set2),na.rm=TRUE)>0

	if(!sign) bestLambda<-unique(featureSelectionOut$Lambda)[set2]	


	########Suggestions by Mike########
	##Not necessary to always choose the second last lambda if there are fluctuations in PRSS##
	if(sign) {
	# make sure the TRUE is continues	
		if(all(diff(set1)==1)){
			bestLambda<-unique(featureSelectionOut$Lambda)[set1[1]]
		}else{
			ly <- which(diff(set1)!=1)
			ly2 <- ly[length(ly)]     #Find the last fluctuation
			bestLambda<-unique(featureSelectionOut$Lambda)[set1[ly2+1]]
		}	

	}


	# (2) select the per-protein lambda

	featureSelectionOut <- featureSelectionOut[order(featureSelectionOut$Lambda, featureSelectionOut$FEATURE, featureSelectionOut$PROTEIN),]
	temp<-unique(featureSelectionOut[,c("PROTEIN","Lambda","PRSSprot")])
	temp1<-matrix(temp$PRSSprot,ncol=length(unique(featureSelectionOut$Lambda)))
	colnames(temp1)<-unique(featureSelectionOut$Lambda)
	rownames(temp1)<-levels(featureSelectionOut$PROTEIN)

	temp2<-t(apply(temp1,1,function(x) abs(diff(x))))

	d <- median(abs(diff(lambda)))   #Mike consider the slope since ratio is affected by outliers

	ratio<-t(apply(temp2,1,function(x) x/d))

	#ratio<-t(apply(temp2,1,function(x) x/ifelse(sum(x,na.rm=TRUE)==0,1,sum(x,na.rm=TRUE))))

	# criteria: has the smallest descreasing and the absolute cutoff of 0.1
	#eta<-0.2
	set1<-apply(ratio,1,function(x) which(x<eta))
	set2<-apply(temp1,1,function(x) which(rank(x,ties.method="min")==1)[1])

	if(is.matrix(set1)){     #In case set1 is not a list, we force set1 be a list
		set1 <- as.list(as.data.frame(set1))
	}

	temp4<-rep(0,dim(ratio)[1])

	for (i in 1:dim(ratio)[1]){
		sign<-sum(set1[[i]]%in%ifelse(set2[i]==dim(temp1)[2],(set2[i]-1),set2[i]),na.rm=TRUE)>0	
		if(!sign) temp4[i]<-set2[i]	
		if(sign) {
			# make sure the TRUE is continues	
			if(all(diff(set1[[i]])==1)){
				temp4[i]<-set1[[i]][1]
			}else{
				ly <- which(diff(set1[[i]])!=1)
				ly2 <- ly[length(ly)]     #Find the last fluctuation
				temp4[i]<-set1[[i]][ly2+1]
			}	
		}

	}
	
	perProteinLambda<-data.frame(Protein=levels(featureSelectionOut$PROTEIN),perProtLambda=unique(featureSelectionOut$Lambda)[temp4],perProtPosition=temp4)


	# generate PRSS plot

	pdf(file=paste(address,"PRSSPlot.pdf",sep=""))
	
	colLambda<-ifelse(unique(featureSelectionOut$Lambda)==bestLambda,7,0)
	borderLambda<-colLambda
	borderLambda[borderLambda==0]<-1
	borderLambda[borderLambda!=1]<-2
	sub<-unique(featureSelectionOut[,c("PRSSprot","Lambda")])
	
	boxplot(PRSSprot~Lambda,data=sub,main="All proteins",xlab="Lambda",ylab="PRSS",col=colLambda,border=borderLambda)

	upPRSS<-quantile(sub[,"PRSSprot"],0.90)
	#upPRSS<-max(sub[,"PRSSprot"], na.rm=TRUE)
	#upPRSS<-0.2

	for (i in 1:nlevels(featureSelectionOut$PROTEIN)){
	
		sub<-featureSelectionOut[featureSelectionOut$PROTEIN==levels(featureSelectionOut$PROTEIN)[i],]
		numKeep<-tapply(sub$Decision,sub[,"Lambda"],function(x) sum(x=="keep",na.rm=TRUE))
		temp<-unique(sub[,c("Lambda","PRSSprot")])
		
		plot(temp,type="b",xaxt="n",ylim=c(0,upPRSS),ylab="PRSS")
		coltemp<-rep(1,dim(temp)[1])
		coltemp[perProteinLambda$perProtPosition[i]]<-2
		points(temp,pch=19,col=coltemp)
		axis(1,at=lambda,labels=lambda)
		axis(3,at=lambda,labels=numKeep)
		mtext("# of selected features",at=lambda[1],line=2,cex=.7)
		mtext(levels(featureSelectionOut$PROTEIN)[i],line=2.5,cex=1.2)
	}

	dev.off()


	# define final feature

	#(1) allProtein - lambda

	featureSelectionFinal<-featureSelectionOut[featureSelectionOut$Lambda==bestLambda,]
	featureSelectionFinal<-featureSelectionFinal[,c("PROTEIN","PEPTIDE","FEATURE","weights","PRSSprot","Decision","Lambda")]
	colnames(featureSelectionFinal)<-c("PROTEIN","PEPTIDE","FEATURE","weights","PRSS","Decision","Lambda")

	write.csv(featureSelectionFinal,file=paste(address,"featureSelectionFinalTable.csv",sep=""))
	save(featureSelectionFinal,file=paste(address,"featureSelectionFinalTable.RData",sep=""))

	#(2) per-Protein - lambda

	featureSelectionPerProtFinal<-NULL
	
	for (i in 1:nlevels(featureSelectionOut$PROTEIN)){
		temp<-featureSelectionOut[featureSelectionOut$PROTEIN==levels(featureSelectionOut$PROTEIN)[i]&featureSelectionOut$Lambda==perProteinLambda$perProtLambda[i],]	
		featureSelectionPerProtFinal<-rbind(featureSelectionPerProtFinal,temp)	
	
	}

	featureSelectionPerProtFinal<-featureSelectionPerProtFinal[,c("PROTEIN","PEPTIDE","FEATURE","weights","PRSSprot","Decision","Lambda")]
	colnames(featureSelectionPerProtFinal)<-c("PROTEIN","PEPTIDE","FEATURE","weights","PRSS","Decision","Lambda")

	write.csv(featureSelectionPerProtFinal,file=paste(address,"featureSelectionPerProtFinalTable.csv",sep=""))
	save(featureSelectionPerProtFinal,file=paste(address,"featureSelectionPerProtFinalTable.RData",sep=""))


	# perProtLambda
	
	workPerProt<-work[work$FEATURE%in%featureSelectionPerProtFinal[featureSelectionPerProtFinal$Decision=="keep","FEATURE"],]
	workPerProt$PROTEIN<-factor(workPerProt$PROTEIN)
	workPerProt$PEPTIDE<-factor(workPerProt$PEPTIDE)
	workPerProt$FEATURE<-factor(workPerProt$FEATURE)

	temp<-featureSelectionPerProtFinal[featureSelectionPerProtFinal$Decision=="keep",c("FEATURE","weights")]
	workPerProt<-merge(workPerProt,temp,by="FEATURE",sort=FALSE)	
	at<-which(colnames(workPerProt)=="TRANSITION")
names<-c(colnames(workPerProt)[2:at],"FEATURE",colnames(workPerProt)[(at+1):dim(workPerProt)[2]])
	workPerProt<-workPerProt[,names]
	colnames(workPerProt)[dim(workPerProt)[2]]<-"WEIGHTS"
	workPerProt<-workPerProt[with(workPerProt,order(PROTEIN,PEPTIDE,FEATURE,LABEL,RUN)),]
	
	# allProtLambda

	work<-work[work$FEATURE%in%featureSelectionFinal[featureSelectionFinal$Decision=="keep","FEATURE"],]
	work$PROTEIN<-factor(work$PROTEIN)
	work$PEPTIDE<-factor(work$PEPTIDE)
	work$FEATURE<-factor(work$FEATURE)

	temp<-featureSelectionFinal[featureSelectionFinal$Decision=="keep",c("FEATURE","weights")]
	work<-merge(work,temp,by="FEATURE",sort=FALSE)	
	at<-which(colnames(work)=="TRANSITION")
	names<-c(colnames(work)[2:at],"FEATURE",colnames(work)[(at+1):dim(work)[2]])
	work<-work[,names]
	colnames(work)[dim(work)[2]]<-"WEIGHTS"
	work<-work[with(work,order(PROTEIN,PEPTIDE,FEATURE,LABEL,RUN)),]

	return(work)
}	

# end-inHouse (The end of the feature selection algorithm)



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



##########################################################################################

.make.contrast.run.quantification<-function(fit,contrast.matrix,sub1,labeled){

	if(class(fit)=="lm"){
		coef.name<-names(coef(fit))
	}else{
		coef.name<-names(fixef(fit))
	}

#####
# intercept
#####
	temp<-coef.name[grep("Intercept",coef.name)]
	
	if(length(temp)>0){
		intercept.c<-rep(1,length(temp))
		names(intercept.c)<-temp
	}else{
		intercept.c<-NULL
	}

#####
# feature
#####
	temp<-coef.name[grep("FEATURE",coef.name)[!grep("FEATURE",coef.name)%in%grep(":",coef.name)]]

	if(length(temp)>0){
		tempSub<-unique(sub1[,c("FEATURE","RUN")])
		tempSub1<-xtabs(~RUN+FEATURE,data=tempSub)
		tempSub2<-tempSub1[contrast.matrix==1,]

		feature.c<-as.numeric(tempSub2[-1])/sum(tempSub2)
		names(feature.c)<-temp
		
	}else{
		feature.c<-NULL
	}
	

#####
# run : different with other quantification - first try
#####

	if(!labeled){ ## label-free
		temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	
		if(length(temp)>0){
			run.c<-contrast.matrix[-1]
			names(run.c)<-temp
		}else{
			run.c<-NULL
		}
	}else{ ## label-based

		temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	
		if(length(temp)>0){
			run.c<-rep(1/nlevels(sub1$RUN),length(temp))
			names(run.c)<-temp
		}else{
			run.c<-NULL
		}
	}
	


#####
# ref
#####
	temp<-coef.name[grep("ref",coef.name)]
	
	if(length(temp)>0){
		if (nlevels(sub1$LABEL)==2){
			levels<-levels(sub1$ref)
			ref.c<-contrast.matrix[1:(length(levels)-1)]
		}

		names(ref.c)<-temp
	}else{
		ref.c<-NULL
	}

####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	tempSub<-dim(unique(sub1[,c("RUN","FEATURE")]))[1]
	if(length(temp)>0){
		rf.c<-rep(1/tempSub,length(temp))
		names(rf.c)<-temp
	}else{
		rf.c<-NULL
	}
	
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
	}else{
		subjectNested.c<-NULL
	}

	
	contrast<-c(intercept.c,feature.c,run.c,ref.c,rf.c, subjectNested.c)

	if(class(fit)=="lm"){
		contrast1<-contrast[!is.na(coef(fit))]
	}else{
		contrast1<-contrast[!is.na(fixef(fit))]
	}
	
	return(contrast1)
}



##########################################################################################
.make.contrast.run.quantification.reference<-function(fit,contrast.matrix,sub1){

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

	if(length(temp)>0){
		tempSub<-unique(sub1[,c("FEATURE","RUN")])
		tempSub1<-xtabs(~RUN+FEATURE,data=tempSub)
		tempSub2<-tempSub1[contrast.matrix==1,]

		feature.c<-as.numeric(tempSub2[-1])/sum(tempSub2)
		names(feature.c)<-temp
		
	}else{
		feature.c<-NULL
	}


#####
# run
#####
	temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	
	if(length(temp)>0){
		run.c<-rep(1/nlevels(sub1$RUN),length(temp))
		names(run.c)<-temp
	}else{
		run.c<-NULL
	}


#####
# ref
#####
	temp<-coef.name[grep("ref",coef.name)]
	if(length(temp)>0){
		ref.c<-rep(0,length(temp))	
		names(ref.c)<-temp
	}else{
		ref.c<-NULL
	}


	contrast<-c(intercept.c,feature.c,run.c,ref.c)

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
#	if(length(temp)>0){
#		if(class(fit)=="lm"){
#			feature.c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
#		}else{
#			### need to fix fit$omdel
#			tempfeature<-eval(getCall(fit)$data)
#			tempfeature$FEATURE<-factor(tempfeature$FEATURE)
#			feature.c<-rep(1/nlevels(tempfeature$FEATURE),length(temp))
#		}
#		names(feature.c)<-temp
#	}else{
#		feature.c<-NULL
#	}

	if(length(temp)>0){
		tempSub<-unique(sub1[,c("FEATURE","GROUP")])
		tempSub1<-xtabs(~GROUP+FEATURE,data=tempSub)

		if (nlevels(sub1$LABEL)==2){
			tempSub1<-tempSub1[-1,]
			tempSub2<-tempSub1[contrast.matrix==1,]
		}
		# free
		if (nlevels(sub1$LABEL)==1){
			tempSub2<-tempSub1[contrast.matrix==1,]
		}

		feature.c<-as.numeric(tempSub2[-1])/sum(tempSub2)
		names(feature.c)<-temp
		
	}else{
		feature.c<-NULL
	}
		
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
	}else{
		subjectNested.c<-NULL
	}

#####
# subject_original
#####
	temp<-coef.name[grep("SUBJECT_ORIGINAL",coef.name)[!grep("SUBJECT_ORIGINAL",coef.name)%in%grep(":",coef.name)]]

	if(length(temp)>0){
#		if(class(fit)=="lm"){
#			subjectOrig.c<-rep(1/nlevels(fit$model$SUBJECT_ORIGINAL),length(temp))
#
#		}else{
#			### need to fix fit$omdel
#			subjectOrig.c<-rep(1/nlevels(eval(getCall(fit)$data)$SUBJECT_ORIGINAL),length(temp))
#		}
		subjectOrig.c<-rep(1/nlevels(sub1$SUBJECT_ORIGINAL),length(temp))
		names(subjectOrig.c)<-temp
	}else{
		subjectOrig.c<-NULL
	}

#####
# subject_original : group
#####
	temp<-coef.name[intersect(grep("SUBJECT_ORIGINAL",coef.name),grep("GROUP",coef.name))]

	if(length(temp)>0){
		tempSub<-unique(sub1[,c("SUBJECT_ORIGINAL","GROUP")])
		tempSub1<-xtabs(~GROUP,data=tempSub)
		tempSub2<-tempSub1[-1]
	
		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\:")),nrow=2))
		temp2<-as.vector(xtabs(~temp1[,2]))	

		if (nlevels(sub1$LABEL)==2){
			subjectOrigGroup.c<-rep(as.numeric(contrast.matrix)/(tempSub2),temp2)
		}
		# free
		if (nlevels(sub1$LABEL)==1){
			subjectOrigGroup.c<-rep(contrast.matrix[-1]/(tempSub2),temp2)
		}
		names(subjectOrigGroup.c)<-temp
	}else{
		subjectOrigGroup.c<-NULL
	}
	
#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]

	if(length(temp)>0){
		if (nlevels(sub1$LABEL)==2){
			group.c<-contrast.matrix
		}
		# free
		if (nlevels(sub1$LABEL)==1){
			group.c<-contrast.matrix[-1]
		}
		names(group.c)<-temp
		
		## if some group's coef is NA, need to remove
		#tempname<-rownames(summary(fit)$coefficients)
		#tempname1<-tempname[grep("GROUP",tempname)[!grep("GROUP",tempname)%in%grep(":",tempname)]]
		#group.c<-group.c[names(group.c)==tempname1]
		
	}else{
		group.c<-NULL
	}
	
#####
# run
#####
	temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	
	if(length(temp)>0){
#		if(class(fit)=="lm"){
#			run.c<-rep(1/nlevels(fit$model$RUN),length(temp))
#		}else{
#			### need to fix fit$omdel
#			run.c<-rep(1/nlevels(eval(getCall(fit)$data)$RUN),length(temp))
#		}
			run.c<-rep(1/nlevels(sub1$RUN),length(temp))

		names(run.c)<-temp
	}else{
		run.c<-NULL
	}

	
####
# feature by group : different from labeled
#####
	temp<-coef.name[intersect(grep("GROUP",coef.name),grep("FEATURE",coef.name))]

	if(length(temp)>0){
		tempSub<-unique(sub1[,c("GROUP","FEATURE")])
		tempSub1<-xtabs(~GROUP,data=tempSub)
		tempSub2<-tempSub1[-1]
	
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
	}else{
		gf.c<-NULL
	}
	
####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	tempSub<-dim(unique(sub1[,c("RUN","FEATURE")]))[1]
	if(length(temp)>0){
		rf.c<-rep(1/tempSub,length(temp))
		names(rf.c)<-temp
	}else{
		rf.c<-NULL
	}

	contrast<-c(intercept.c,feature.c,subjectNested.c,subjectOrig.c, subjectOrigGroup.c, group.c,run.c,gf.c,rf.c)

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
#	if(length(temp)>0){
#		if(class(fit)=="lm"){
#			feature.c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
#		}else{
#			### need to fix fit$omdel
#			tempfeature<-eval(getCall(fit)$data)
#			tempfeature$FEATURE<-factor(tempfeature$FEATURE)
#			feature.c<-rep(1/nlevels(tempfeature$FEATURE),length(temp))
#		}
#		names(feature.c)<-temp
#	}else{
#		feature.c<-NULL
#	}

	if(length(temp)>0){
		tempSub<-unique(sub1[,c("FEATURE","GROUP")])
		tempSub1<-xtabs(~GROUP+FEATURE,data=tempSub)
		
		if (nlevels(sub1$LABEL)==2){
			tempSub1<-tempSub1[-1,]
			tempSub2<-tempSub1[contrast.matrix==1,]
		}
		# free
		if (nlevels(sub1$LABEL)==1){
			tempSub2<-tempSub1[contrast.matrix==1,]
		}

		feature.c<-as.numeric(tempSub2[-1])/sum(tempSub2)
		names(feature.c)<-temp
		
	}else{
		feature.c<-NULL
	}
		


#####
# subject_nested
#####
	temp<-coef.name[grep("SUBJECT_NESTED",coef.name)]
	
	if(length(temp)>0){
		subjectNested.c<-rep(0,length(temp))	
		names(subjectNested.c)<-temp
	}else{
		subjectNested.c<-NULL
	}
	
#####
# subject_original
#####
	temp<-coef.name[grep("SUBJECT_ORIGINAL",coef.name)[!grep("SUBJECT_ORIGINAL",coef.name)%in%grep(":",coef.name)]]
	
	if(length(temp)>0){
		subjectOrig.c<-rep(0,length(temp))	
		names(subjectOrig.c)<-temp
	}else{
		subjectOrig.c<-NULL
	}
	
#####
# subject_original : group
#####
	temp<-coef.name[intersect(grep("SUBJECT_ORIGINAL",coef.name),grep("GROUP",coef.name))]

	if(length(temp)>0){
		subjectOrigGroup.c<-rep(0,length(temp))	
		names(subjectOrigGroup.c)<-temp
	}else{
		subjectOrigGroup.c<-NULL
	}
	
#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]

	if(length(temp)>0){
		group.c<-rep(0,length(temp))	
		names(group.c)<-temp
	}else{
		group.c<-NULL
	}

#####
# run
#####
	temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]

	if(length(temp)>0){
#		if(class(fit)=="lm"){
#			run.c<-rep(1/nlevels(fit$model$RUN),length(temp))
#		}else{
#			### need to fix fit$omdel
#			temprun<-eval(getCall(fit)$data)
#			tempfeature$FEATURE<-factor(tempfeature$FEATURE)
#			feature.c<-rep(1/nlevels(tempfeature$FEATURE),length(temp))
#			run.c<-rep(1/nlevels(eval(getCall(fit)$d)$RUN),length(temp))
#		}
		run.c<-rep(1/nlevels(sub1$RUN),length(temp))
		names(run.c)<-temp
	}else{
		run.c<-NULL
	}
	
####
# feature by group : different from labeled
#####
	temp<-coef.name[intersect(grep("GROUP",coef.name),grep("FEATURE",coef.name))]
	if(length(temp)>0){
		gf.c<-rep(0,length(temp))	
		names(gf.c)<-temp
	}else{
		gf.c<-NULL
	}

####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	tempSub<-dim(unique(sub1[,c("RUN","FEATURE")]))[1]
	if(length(temp)>0){
		rf.c<-rep(1/tempSub,length(temp))
		names(rf.c)<-temp
	}else{
		rf.c<-NULL
	}

	contrast<-c(intercept.c,feature.c,subjectNested.c,subjectOrig.c,subjectOrigGroup.c,group.c,run.c,gf.c,rf.c)

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
#	if(length(temp)>0){
#		if(class(fit)=="lm"){
#			feature.c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
#		}else{
#			### need to fix fit$omdel
#			tempfeature<-eval(getCall(fit)$data)
#			tempfeature$FEATURE<-factor(tempfeature$FEATURE)
#			feature.c<-rep(1/nlevels(tempfeature$FEATURE),length(temp))
#		}
#		names(feature.c)<-temp
#	}else{
#		feature.c<-NULL
#	}

	if(length(temp)>0){
		tempSub<-unique(sub1[,c("FEATURE","SUBJECT_NESTED")])
		tempSub1<-xtabs(~SUBJECT_NESTED+FEATURE,data=tempSub)

		if (nlevels(sub1$LABEL)==2){
			tempSub1<-tempSub1[-1,]
			tempSub2<-tempSub1[contrast.matrix==1,]
		}
		# free
		if (nlevels(sub1$LABEL)==1){
			tempSub2<-tempSub1[contrast.matrix==1,]
		}

		feature.c<-as.numeric(tempSub2[-1])/sum(tempSub2)
		names(feature.c)<-temp
		
	}else{
		feature.c<-NULL
	}

		


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
	}else{
		subjectNested.c<-NULL
	}

#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]
	if(length(temp)>0){
		if (nlevels(sub1$LABEL)==2){
			group.c<-contrastGroup
		}
		# free
		if (nlevels(sub1$LABEL)==1){
			group.c<-contrastGroup[-1]
		}
		names(group.c)<-temp
	}else{
		group.c<-NULL
	}

#####
# run
#####
	temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]

	if(length(temp)>0){
#run.c<-rep(1/nlevels(fit$model$RUN),length(temp))
## when no technical replicate : subject_nested = run
		run.c<-contrast.matrix[-1]

## however, with technical replicate : subject_nested != run, need others

		names(run.c)<-temp
	}else{
		run.c<-NULL
	}

####
# feature by group : different from labeled
#####
	temp<-coef.name[intersect(grep("GROUP",coef.name),grep("FEATURE",coef.name))]

	if(length(temp)>0){
		tempSub<-unique(sub1[,c("GROUP","FEATURE")])
		tempSub1<-xtabs(~GROUP,data=tempSub)
		tempSub2<-tempSub1[-1]
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
	}else{
		gf.c<-NULL
	}

####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	tempSub<-dim(unique(sub1[,c("RUN","FEATURE")]))[1]
	if(length(temp)>0){
		rf.c<-rep(1/tempSub,length(temp))
		names(rf.c)<-temp
	}else{
		rf.c<-NULL
	}
	
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
		contrastList<-unique(sub1[,c("GROUP","SUBJECT_ORIGINAL","SUBJECT_NESTED")])
		contrastList<-contrastList[contrastList$GROUP!="0",] ## remove GROUP==0
		contrastList$GROUP<-factor(contrastList$GROUP) ## remove '0' group
		contrastList$SUBJECT_ORIGINAL<-factor(contrastList$SUBJECT_ORIGINAL) ## remove '0' group

		contrastGroup<-rep(0,(nlevels(sub1$GROUP)-1))
		contrastGroup[as.numeric(contrastList[contrast.matrix==1,"GROUP"])]<-1
		
		contrastSubjectOriginal<-rep(0,(nlevels(sub1$SUBJECT_ORIGINAL)))
		contrastSubjectOriginal[as.numeric(contrastList[contrast.matrix==1,"SUBJECT_ORIGINAL"])]<-1

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
	
	if(length(temp)>0){
		intercept.c<-rep(1,length(temp))
		names(intercept.c)<-temp
	}else{
		intercept.c<-NULL
	}

#####
# feature
#####
	temp<-coef.name[grep("FEATURE",coef.name)[!grep("FEATURE",coef.name)%in%grep(":",coef.name)]]
#	if(length(temp)>0){
#		if(class(fit)=="lm"){
#			feature.c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
#		}else{
#			### need to fix fit$omdel
#			tempfeature<-eval(getCall(fit)$data)
#			tempfeature$FEATURE<-factor(tempfeature$FEATURE)
#			feature.c<-rep(1/nlevels(tempfeature$FEATURE),length(temp))
#		}
#		names(feature.c)<-temp
#	}else{
#		feature.c<-NULL
#	}

	if(length(temp)>0){
		tempSub<-unique(sub1[,c("FEATURE","SUBJECT_NESTED")])
		tempSub1<-xtabs(~SUBJECT_NESTED+FEATURE,data=tempSub)

		if (nlevels(sub1$LABEL)==2){
			tempSub1<-tempSub1[-1,]
			tempSub2<-tempSub1[contrast.matrix==1,]
		}
		# free
		if (nlevels(sub1$LABEL)==1){
			tempSub2<-tempSub1[contrast.matrix==1,]
		}

		feature.c<-as.numeric(tempSub2[-1])/sum(tempSub2)
		names(feature.c)<-temp
		
	}else{
		feature.c<-NULL
	}
		

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
	}else{
		subjectNested.c<-NULL
	}
	
#####
# subject_original
#####
	temp<-coef.name[grep("SUBJECT_ORIGINAL",coef.name)[!grep("SUBJECT_ORIGINAL",coef.name)%in%grep(":",coef.name)]]
	
	if(length(temp)>0){
		subjectOrig.c<-contrastSubjectOriginal[-1]
		names(subjectOrig.c)<-temp
	}else{
		subjectOrig.c<-NULL
	}


#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]

	if(length(temp)>0){
		if (nlevels(sub1$LABEL)==2){
			group.c<-contrastGroup
		}
		# free
		if (nlevels(sub1$LABEL)==1){
			group.c<-contrastGroup[-1]
		}
		names(group.c)<-temp
	}else{
		group.c<-NULL
	}

#####
# run
#####
	temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	
	if(length(temp)>0){
#		if(class(fit)=="lm"){
#			run.c<-rep(1/nlevels(fit$model$RUN),length(temp))
#		}else{
#			### need to fix fit$omdel
#			run.c<-rep(1/nlevels(eval(getCall(fit)$data)$RUN),length(temp))
#		}
		run.c<-rep(1/nlevels(sub1$RUN),length(temp))
		names(run.c)<-temp
	}else{
		run.c<-NULL
	}

####
# feature by group : different from labeled
#####
	temp<-coef.name[intersect(grep("GROUP",coef.name),grep("FEATURE",coef.name))]

	if(length(temp)>0){
		tempSub<-unique(sub1[,c("GROUP","FEATURE")])
		tempSub1<-xtabs(~GROUP,data=tempSub)
		tempSub2<-tempSub1[-1]
	
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
	}else{
		gf.c<-NULL
	}
	
####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	tempSub<-dim(unique(sub1[,c("RUN","FEATURE")]))[1]
	if(length(temp)>0){
		rf.c<-rep(1/tempSub,length(temp))
		names(rf.c)<-temp
	}else{
		rf.c<-NULL
	}

	contrast<-c(intercept.c,feature.c,subjectNested.c,subjectOrig.c,group.c,run.c,gf.c,rf.c)

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
#	if(length(temp)>0){
#		if(class(fit)=="lm"){
#			feature.c<-rep(1/nlevels(fit$model$FEATURE),length(temp))
#		}else{
#			### need to fix fit$omdel
#			tempfeature<-eval(getCall(fit)$data)
#			tempfeature$FEATURE<-factor(tempfeature$FEATURE)
#			feature.c<-rep(1/nlevels(tempfeature$FEATURE),length(temp))
#		}
#		names(feature.c)<-temp
#	}else{
#		feature.c<-NULL
#	}

	if(length(temp)>0){
		tempSub<-unique(sub1[,c("FEATURE","SUBJECT_NESTED")])
		tempSub1<-xtabs(~SUBJECT_NESTED+FEATURE,data=tempSub)
		
		if (nlevels(sub1$LABEL)==2){
			tempSub1<-tempSub1[-1,]
			tempSub2<-tempSub1[contrast.matrix==1,]
		}
		# free
		if (nlevels(sub1$LABEL)==1){
			tempSub2<-tempSub1[contrast.matrix==1,]
		}

		feature.c<-as.numeric(tempSub2[-1])/sum(tempSub2)
		names(feature.c)<-temp
		
	}else{
		feature.c<-NULL
	}


#####
# subject_nested
#####
	temp<-coef.name[grep("SUBJECT_NESTED",coef.name)]
	if(length(temp)>0){
		subjectNested.c<-rep(0,length(temp))	
		names(subjectNested.c)<-temp
	}else{
		subjectNested.c<-NULL
	}

#####
# subject_original
#####
	temp<-coef.name[grep("SUBJECT_ORIGINAL",coef.name)[!grep("SUBJECT_ORIGINAL",coef.name)%in%grep(":",coef.name)]]
	
	if(length(temp)>0){
		subjectOrig.c<-rep(0,length(temp))
		names(subjectOrig.c)<-temp
	}else{
		subjectOrig.c<-NULL
	}
	
#####
# group : different from labeled
#####
	temp<-coef.name[grep("GROUP",coef.name)[!grep("GROUP",coef.name)%in%grep(":",coef.name)]]

	if(length(temp)>0){
		group.c<-rep(0,length(temp))	
		names(group.c)<-temp
	}else{
		group.c<-NULL
	}

#####
# run
#####
	temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	
	if(length(temp)>0){
#		if(class(fit)=="lm"){
#			run.c<-rep(1/nlevels(fit$model$RUN),length(temp))
#		}else{
#			### need to fix fit$omdel
#			run.c<-rep(1/nlevels(eval(getCall(fit)$data)$RUN),length(temp))
#		}
		run.c<-rep(1/nlevels(sub1$RUN),length(temp))
		names(run.c)<-temp
	}else{
		run.c<-NULL
	}
	
####
# feature by group : different from labeled
#####
	temp<-coef.name[intersect(grep("GROUP",coef.name),grep("FEATURE",coef.name))]
	if(length(temp)>0){
		gf.c<-rep(0,length(temp))	
		names(gf.c)<-temp
	}else{
		gf.c<-NULL
	}

####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	tempSub<-dim(unique(sub1[,c("RUN","FEATURE")]))[1]
	if(length(temp)>0){
		rf.c<-rep(1/tempSub,length(temp))
		names(rf.c)<-temp
	}else{
		rf.c<-NULL
	}
	
	contrast<-c(intercept.c,feature.c,subjectNested.c,subjectOrig.c,group.c,run.c,gf.c,rf.c)

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
	levels<-levels(tempdata$GROUP)
	labels<-paste("GROUP", levels, sep="")
	patients<-NULL
	for(i in 1:length(levels)){
		sub<-tempdata[tempdata$GROUP==levels[i],]
		sub.patients<-cbind(GROUP=paste("GROUP", levels[i], sep=""), SUBJECT=paste("SUBJECT", as.character(levels(sub$SUBJECT)), sep=""), Value=contrast.matrix[as.numeric(as.character(levels[i]))])
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
	tempSub<-as.numeric(as.character(levels(sub1[,c("GROUP")])))
	tempcontrast<-contrast.matrix[tempSub]

	group.c<-tempcontrast[-1] ## for label-free, need to remove first
	names(group.c)<-temp
	if(length(temp)==0) group.c<-NULL


#####
# subject_nested 
#####
	temp<-coef.name[grep("SUBJECT_NESTED",coef.name)]
	if(length(temp)>0){
		temp1<-t(matrix(unlist(strsplit(as.character(temp),"\\.")),nrow=2))
		temp2<-as.vector(xtabs(~temp1[,1]))	
		
		tempdata<-fit$model
		levels<-levels(tempdata$GROUP)
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
# subject by group : only for time-course - SUBJECT and GROUP (order) even GROUP:SUBJECT in model
#####
	temp<-coef.name[intersect(grep("SUBJECT",coef.name),grep("GROUP",coef.name))]
	
	if(length(temp)>0){
#		subject.c<-rep(0,length(temp))
#		names(subject.c)<-temp
	
	tempdata<-fit$model
	levels<-levels(tempdata$GROUP)
	labels<-paste("GROUP", levels, sep="")
	patients<-NULL
	for(i in 1:length(levels)){
		sub<-tempdata[tempdata$GROUP==levels[i],]
		sub.patients<-cbind(GROUP=paste("GROUP", levels[i], sep=""), SUBJECT=paste("SUBJECT", as.character(levels(sub$SUBJECT)), sep=""), Value=contrast.matrix[as.numeric(as.character(levels[i]))])
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
	contrast<-c(intercept.c,group.c,subjectNested.c, subject.c, gs.c)

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
	### todo : for unbalanced case, variance is weighted by degree of freedom	
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

.estimableFixedQuantification<-function(cf,cm){

	cm <- matrix(cm, nrow=1)
	ct <- cm %*% cf[, 1]
	result<-cbind(est=ct)
	colnames(result)<-c("log-intensities")
	
	return(result)
}

##########################################################################################

.estimableRandomQuantification<-function(cf,cm){

	cm <- matrix(cm, nrow=1)
	ct <- cm %*% cf
	result<-cbind(est=ct)
	colnames(result)<-c("log-intensities")
	
	return(result)
}

##########################################################################################

.estimableFixedQuantificationSurvival<-function(cf,cm){

	cm <- matrix(cm, nrow=1)
	ct <- cm %*% cf
	result<-cbind(est=ct)
	colnames(result)<-c("log-intensities")
	
	return(result)
}


##########################################################################################

.make.contrast.run.quantification.Survival<-function(fit,contrast.matrix,sub1,labeled){

		coef.name<-names(fit$coefficients)

#####
# intercept
#####
	temp<-coef.name[grep("Intercept",coef.name)]
	
	if(length(temp)>0){
		intercept.c<-rep(1,length(temp))
		names(intercept.c)<-temp
	}else{
		intercept.c<-NULL
	}

#####
# feature
#####
	temp<-coef.name[grep("FEATURE",coef.name)[!grep("FEATURE",coef.name)%in%grep(":",coef.name)]]

	if(length(temp)>0){
		tempSub<-unique(sub1[,c("FEATURE","RUN")])
		tempSub1<-xtabs(~RUN+FEATURE,data=tempSub)
		tempSub2<-tempSub1[contrast.matrix==1,]

		feature.c<-as.numeric(tempSub2[-1])/sum(tempSub2)
		names(feature.c)<-temp
		
	}else{
		feature.c<-NULL
	}
	

#####
# run : different with other quantification - first try
#####

	if(!labeled){ ## label-free
		temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	
		if(length(temp)>0){
			run.c<-contrast.matrix[-1]
			names(run.c)<-temp
		}else{
			run.c<-NULL
		}
	}else{ ## label-based

		temp<-coef.name[grep("RUN",coef.name)[!grep("RUN",coef.name)%in%grep(":",coef.name)]]
	
		if(length(temp)>0){
			run.c<-rep(1/nlevels(sub1$RUN),length(temp))
			names(run.c)<-temp
		}else{
			run.c<-NULL
		}
	}
	


#####
# ref
#####
	temp<-coef.name[grep("ref",coef.name)]
	
	if(length(temp)>0){

		levels<-levels(sub1$ref)
		ref.c<-contrast.matrix[1:(length(levels)-1)]
	
		names(ref.c)<-temp
		
	}else{
		ref.c<-NULL
	}

####
# run by feature
#####
	temp<-coef.name[intersect(grep("RUN",coef.name),grep("FEATURE",coef.name))]
	tempSub<-dim(unique(sub1[,c("RUN","FEATURE")]))[1]
	if(length(temp)>0){
		rf.c<-rep(1/tempSub,length(temp))
		names(rf.c)<-temp
	}else{
		rf.c<-NULL
	}
	
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
	}else{
		subjectNested.c<-NULL
	}

	
	contrast<-c(intercept.c,feature.c,run.c,ref.c,rf.c, subjectNested.c)

		contrast1<-contrast[!is.na(fit$coefficients)]

	
	return(contrast1)
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

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
.multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

