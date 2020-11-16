# ##================================
# ## label-free, single
# ##================================
# .make.contrast.free.single <- function(fit, contrast.matrix, sub1) {
#     
#     if (class(fit) == "lm") {
#         coef_name <- names(coef(fit))
#     } else {
#         coef_name <- names(fixef(fit))
#     }
#     
#     ## change contrast.matrix without Group1
#     # cons=1
#     # if(contrast.matrix[1]==0) contrast.free<-contrast.matrix[-1]
#     # if(contrast.matrix[1]<0){ 
#     # 	contrast.free<-contrast.matrix[-1]/abs(contrast.matrix[1])
#     # 	cons=abs(contrast.matrix[1])
#     # 	}
#     # if(contrast.matrix[1]>0){ 
#     # 	contrast.free<-contrast.matrix[-1]*(-1)/contrast.matrix[1]
#     # 	cons=contrast.matrix[1]*(-1)
#     # 	}
#     # if(class(fit)=="lm"){
#     # 	coef_name<-names(coef(fit))
#     # }else{
#     # 	coef_name<-names(fixef(fit))
#     # }
#     
#     ## intercept
#     temp <- coef_name[grep("Intercept", coef_name)]
#     intercept_c <- rep(0, length(temp))
#     names(intercept_c) <- temp
#     if (length(temp) == 0) {
#         intercept_c <- NULL
#     }
#     
#     ## subject
#     temp <- coef_name[setdiff(grep("SUBJECT", coef_name), 
#                               grep(":|NESTED", coef_name))]
#     if (length(temp) > 0) {
#         
#         # subject_c<-rep(0,length(temp))
#         # names(subject_c)<-temp
#         
#         tempdata <- fit$model
#         group_levels <- levels(tempdata$GROUP)
#         labels <- paste("GROUP", group_levels, sep="")
#         patients <- NULL
#         for (i in 1:length(group_levels)) {
#             sub <- tempdata[tempdata$GROUP == group_levels[i], ]
#             sub_patients <- cbind(
#                 GROUP=paste("GROUP", group_levels[i], sep=""), 
#                 SUBJECT=paste("SUBJECT", as.character(group_levels(sub$SUBJECT)), sep=""), 
#                 Value=contrast.matrix[as.numeric(as.character(group_levels[i]))])
#             patients <- data.frame(rbind(patients, sub_patients))
#         }
#         patient_count <- tapply(patients$SUBJECT, patients$GROUP, 
#                                 function(x) length(unique(x)))
#         patient_seq <- rep(0, length(temp))
#         for (i in 1:length(as.character(patients$SUBJECT))) {
#             match <- any(temp == as.character(patients$SUBJECT)[i])
#             if (match & as.numeric(as.character(patients$Value[i])) != 0) {
#                 res <- temp == as.character(patients$SUBJECT)[i]
#                 index <- which(res == TRUE)
#                 group <- as.character(patients[i, ]$GROUP)
#                 count <- as.numeric(patient_count[names(patient_count) == group])
#                 value <- as.numeric(as.character(patients[i, ]$Value))
#                 patient_value <- c(rep(0, index-1), 
#                                    value / count, 
#                                    rep(0, length(temp) - index))
#             } else {
#                 patient_value <- rep(0, length(temp))
#             }
#             patient_seq <- patient_value + patient_seq
#         }
#         subject_c <- patient_seq
#         names(subject_c) <- temp
#     } else if (length(temp) == 0) {
#         subject_c <- NULL
#     }
#     
#     ## group: different from labeled
#     temp <- coef_name[setdiff(grep("GROUP", coef_name), grep(":", coef_name))]
#     ## when there are some groups which are all missing
#     tempSub <- as.numeric(as.character(levels(sub1[, c("GROUP")])))
#     tempcontrast <- contrast.matrix[tempSub]
#     ## for label-free, need to remove first
#     group_c <- tempcontrast[-1] 
#     names(group_c) <- temp
#     if (length(temp) == 0) {
#         group_c<-NULL
#     }
#     
#     ## subject_nested
#     temp <- coef_name[grep("SUBJECT_NESTED", coef_name)]
#     if (length(temp) > 0) {
#         temp1 <- t(matrix(unlist(strsplit(as.character(temp), "\\.")), nrow=2))
#         temp2 <- as.vector(xtabs(~ temp1[, 1]))
#         tempdata <- fit$model
#         group_levels <- levels(tempdata$GROUP)
#         sub_contrast <- contrast.matrix[as.numeric(as.character(group_levels))]
#         
#         ## the base is alway be the first SUBJECT_NESTED
#         ## (SUBJECT_NESTED1.1)
#         temp3 <- temp2
#         if (length(temp2) == length(sub_contrast)) {
#             ## this line first, otherwise length of temp3 >1
#             temp3[1] <- temp2[1] + 1
#         } else {
#             temp3 <- c(1, temp3)
#         }
#         
#         ## in case of unequal sample per group, wrong
#         # subjectNested_c<-rep(contrast.matrix/(temp3),temp2) 
#         
#         subjectNested_c <- rep(sub_contrast / temp3, temp3)[-1]
#         names(subjectNested_c) <- temp
#     } else if (length(temp) == 0) {
#         subjectNested_c <- NULL
#     }
#     
#     ## subject by group: only for time-course - SUBJECT and GROUP (order) even GROUP:SUBJECT in model
#     temp <- coef_name[intersect(grep("SUBJECT", coef_name), 
#                                 grep("GROUP", coef_name))]
#     if (length(temp) > 0) {
#         
#         # subject_c<-rep(0,length(temp))
#         # names(subject_c)<-temp
#         
#         tempdata <- fit$model
#         group_levels <- levels(tempdata$GROUP)
#         labels <- paste("GROUP", group_levels, sep="")
#         patients <- NULL
#         for (i in 1:length(group_levels)) {
#             sub <- tempdata[tempdata$GROUP == group_levels[i], ]
#             sub_patients <- cbind(
#                 GROUP=paste("GROUP", group_levels[i], sep=""), 
#                 SUBJECT=paste("SUBJECT", as.character(levels(sub$SUBJECT)), sep=""), 
#                 Value=contrast.matrix[as.numeric(as.character(group_levels[i]))])
#             patients <- data.frame(rbind(patients, sub_patients))
#         }
#         patient_count <- tapply(patients$SUBJECT, patients$GROUP, 
#                                 function(x) length(unique(x)))
#         interaction_seq <- rep(0, length(temp))
#         interaction_labels <- paste(as.character(patients$GROUP), 
#                                     as.character(patients$SUBJECT), sep=":")
#         for (i in 1:length(as.character(patients$SUBJECT))) {
#             match <- any(temp == interaction_labels[i])
#             if (match & as.numeric(as.character(patients$Value[i])) != 0) {
#                 res <- temp == interaction_labels[i]
#                 index <- which(res == TRUE)
#                 group <- as.character(patients[i, ]$GROUP)
#                 count <- as.numeric(patient_count[names(patient_count) == group])
#                 value <- as.numeric(as.character(patients[i, ]$Value))
#                 interaction_value <- c(rep(0, index - 1), 
#                                        value / count, 
#                                        rep(0, length(temp) - index))
#             } else {
#                 interaction_value <- rep(0, length(temp))
#             }
#             interaction_seq <- interaction_value + interaction_seq
#         }
#         gs_c <- interaction_seq
#         names(gs_c) <- temp
#     } else if (length(temp) == 0) {
#         gs_c <- NULL
#     }
#     
#     ## combine all
#     contrast <- c(intercept_c, group_c, subjectNested_c, subject_c, gs_c)
#     if (class(fit) == "lm") {
#         contrast1 <- contrast[!is.na(coef(fit))]
#     } else {
#         contrast1 <- contrast[!is.na(fixef(fit))]
#     }
#     
#     return(contrast1)
# }
# 
# 
# ########### estimate   ###########
# .getParameterFixed <- function(obj) {
#     temp1 <- summary.lm(obj)
#     cf <- temp1$coefficients
#     vcv <- temp1$cov.unscaled * temp1$sigma ^ 2
#     ## TODO (): for unbalanced case, variance is weighted by degree of freedom	
#     df <- obj$df.residual
#     parameter <- list(cf=cf, vcv=vcv, df=df)
#     
#     return(parameter)
# }
# 
# .getParameterRandom <- function(obj, df.full) {
#     cf <- as.matrix(fixef(obj))
#     vcv <- as.matrix(vcov(obj))
#     df <- df.full
#     parameter <- list(cf=cf, vcv=vcv, df=df)
#     
#     return(parameter)
# }
# 
# 
# .estimableFixedRandom <- function(parameter, cm) {
#     cm <- matrix(cm, nrow=1)
#     ct <- cm %*% parameter$cf[, 1]
#     vc <- sqrt(diag(cm %*% parameter$vcv %*% t(cm)))
#     prob <- 2 * (1 - pt(abs(ct / vc), parameter$df))
#     result <- cbind(est=ct, stderr=vc, t=ct / vc, df=parameter$df, prob=prob)
#     colnames(result) <- c("logFC", "SE", "Tvalue", "DF", "pvalue")
#     
#     return(result)
# }