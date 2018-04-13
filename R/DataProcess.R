#' dataProcess() is a large function which does what it says on the tin.
#'
#' Ideally, this function should perform a series of normalizations on a msnset
#' and figure out how much stuff is missing among other operations, then perform
#' comparisons across the experimental conditions.
#'
#' @param raw  A Raw msnset to analyze, in my hands, this is provided by
#'   SWATH2stats.
#' @param logTrans  Transform the data on a log2 scale?
#' @param normalization  Normalization method to perform.  (Wow this feels like
#'   rnaseq now)
#' @param nameStandards  Global standards to help normalization.
#' @param betweenRunInterferenceSCore  That variable name is too long, I forgot
#'   what it does while typing it out.
#' @param address  the only place this exists is in a paste0() call, I have no
#'   clue.
#' @param fillIncompleteRows  If experimental design is missing, fill it in.
#' @param featureSubset  Subset the data by feature (topn or 3 likely, or all)
#' @param remove_proteins_with_interference  Do we?
#' @param n_top_feature  How many features to subset? This should be used by
#'   featuresubset
#' @param summaryMethod  Which summary method to use?
#' @param equalFeatureVar  Calculate variance for equal features (I guess?)
#' @param censoredInt  What value for censored rows?
#' @param cutoffCensored  How to decide a censored row
#' @param MBimpute  Impute Impute damned spot!
#' @param original_scale  Return the results on the same scale as the input data (I assume)?
#' @param remove50missing  wtf?
#' @param maxQuantileforCensored  yeah, whatever
#' @param clusters The size of a SNOW cluster to speak up the calculations.
#' @export
dataProcess <- function(raw, logTrans=2, normalization="equalizeMedians", nameStandards=NULL,
                        betweenRunInterferenceScore=FALSE, address="", fillIncompleteRows=TRUE,
                        featureSubset="all", remove_proteins_with_interference=FALSE,
                        n_top_feature=3, summaryMethod="TMP", equalFeatureVar=TRUE,
                        censoredInt="NA", cutoffCensored="minFeature", MBimpute=TRUE,
                        original_scale=FALSE, remove50missing=FALSE,
                        maxQuantileforCensored=0.999, clusters=1) {
  ## Standardize the column names.
  colnames(raw) <- toupper(colnames(raw))
  requiredInput <- c("PROTEINNAME", "PEPTIDESEQUENCE", "PRECURSORCHARGE",
                     "FRAGMENTION", "PRODUCTCHARGE", "ISOTOPLABELCHARGE",
                     "CONDITION", "BIOREPLICATE", "RUN", "INTENSITY")

  ## [THT: disambiguation for PeptideSequence & PeptideModifiedSequence - begin]
  ## PeptideModifiedSequence is also allowed.
  providedInput <- colnames(raw)
  if (all(requiredInput %in% providedInput)) {
    logging::loginfo("The required inputs were provided.")
  } else if (all(setdiff(requiredInput, "PEPTIDESEQUENCE") %in%
                 providedInput) && "PEPTIDEMODIFIEDSEQUENCE" %in%
             providedInput) {
    logging::loginfo("The required inputs were provided.")
    ## if PeptideModifiedSequence is provided instead of PeptideSequence,
    ## change the column name as PeptideSequence
  } else {
    missedInput <- which(!(requiredInput %in% providedInput))
    logstring <- paste0(
      "The required inputs: ", toString(requiredInput[missedInput]),
      " were not provided. the required inputs are: \n",
      "Proteinname, PeptideSequence/PeptideModifiedSequence, PrecursorCharge, FragmentIon,
ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity")
    logging::logwarn(logstring)
  }

  ## check logTrans is 2,10 or not
  if (logTrans != 2 & logTrans != 10) {
    logging::logwarn(paste0("The transformation is: ", logTrans,
                            ". Setting transformation to log2."))
    logTrans <- 2
  }
  ## check no row for some feature : balanced structure or not
  if (!(fillIncompleteRows == TRUE |
        fillIncompleteRows == FALSE) |
      !is.logical(fillIncompleteRows)) {
    logging::logwarn(paste0("The option to fill incomplete rows is: ",
                            fillIncompleteRows, " setting it to TRUE."))
    fillIncompleteRows <- TRUE
  }
  ## check input for summaryMethod
  if (sum(summaryMethod == c("linear", "TMP")) == 0) {
    logging::logwarn(paste0("The summary method is: ",
                            summaryMethod, " setting it to 'TMP'."))
    summaryMethod <- "TMP"
  } else {
    logging::loginfo(paste0("The summary method is: ", summaryMethod))
  }
  ## check input for cutoffCensored
  if (sum(cutoffCensored == c("minFeature", "minRun", "minFeatureNRun")) == 0) {
    logging::logwarn(paste0("The cutoff censoring method is: ",
                            cutoffCensored, " setting it to minFeature."))
    cutoffCensored <- "minFeature"
  } else {
    logging::loginfo(paste0("The cutoff censor method is: ", cutoffCensored, "."))
  }
  ## check input for censoredInt
  if (sum(censoredInt == c("0", "NA")) == 0 & !is.null(censoredInt)) {
    logging::logwarn(paste0("The censored int is: ", censoredInt,
                            " setting it to 0."))
    censoredInt <- 0
  } else {
    logging::loginfo(paste0("The censored int is: ", as.character(censoredInt)))
  }
  ## check input for censoredInt and MBimpute
  if (summaryMethod == "TMP" & MBimpute & is.null(censoredInt)) {
    logging::logwarn("The provided combination of censoredInt and MBimpute is invalid.
Imputation will not be performed.")
    MBimpute <- FALSE
  }

  normalization <- toupper(normalization)
  ## [THT: if (!all(normalization %in% c("NONE", "FALSE",
  ## "EQUALIZEMEDIANS", "QUANTILE", "GLOBALSTANDARDS")))]
  ## [THT: send a warning message if the user mixes "NONE" with any of the last three choices]
  if (!(normalization == "NONE" |
        normalization == "FALSE" |
        normalization == "EQUALIZEMEDIANS" |
        normalization == "QUANTILE" |
        normalization == "GLOBALSTANDARDS")) {
    logging::logwarn(paste0("An invalid normalization was provided, skipping normalization."))
    normalization <- "NONE"
  }
  ## need the names of global standards
  if (!is.element("NONE", normalization) &
      !is.element("FALSE", normalization) &
      is.element("GLOBALSTANDARDS", normalization) &
      is.null(nameStandards)) {
    logging::logwarn(paste0("Global standard normalization was chosen, but no standards were provided.
Disabling normalization, if you want to use them, add the nameStandards option.
Sleeping for 5 seconds to give you a chance to stop and reconsider your life choices."))
    Sys.sleep(5)
    normalization <- "NONE"
  }
  ## check whether class of intensity is factor or chaterer, if yes, neec to chage as numeric
  if (is.factor(raw[["INTENSITY"]]) | is.character(raw[["INTENSITY"]])) {
    suppressWarnings(raw[["INTENSITY"]] <- as.numeric(as.character(raw[["INTENSITY"]])))
  }

  ## check whether the intensity has 0 value or negative value
  ##    if (length(which(raw$Intensity<=0))>0 & !skylineReport) {
  ##    if (is.null(censoredInt)) {
  ##        stop("Intensity has 0 or negative values. Please check these intensities
  ##and change them. \n")
  ##    } else if (censoredInt=="NA") {
  ##        stop("Intensity has 0 or negative values. Please check these intensities and
  ##change them. \n")
  ##        }
  ##}

  ## here, need to get standard protein name
  ## column name : standardtype..
  ## what value it has, normzalition, unique(proteinname)
  ## if normalition== "standard" & no normalizaion selection, error message
  ## For Skyline
  ## required cols : ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge,
  ##                 IsotopeLabelType, and Condition, BioReplicate, Run, Intensity

  ## annotation information :
  if (any(is.na(raw[["RUN"]]))) {
    missing_runs <- is.na(raw[["RUN"]])
    logging::logwarn("There is missing information in the 'Run' column, filling it with 'undefined'.")
    raw[missing_runs, "RUN"] <- "undefined"
  }
  if (any(is.na(raw[["BIOREPLICATE"]]))) {
    missing_replicates <- is.na(raw[["BIOREPLICATE"]])
    logging::logwarn("There is missing information in the 'BioReplicate' column, filling it with 'undefined'.")
    raw[missing_replicates, "BIOREPLICATE"] <- "undefined"
  }
  if (any(is.na(raw[["CONDITION"]]))) {
    missing_conditions <- is.na(raw[["CONDITION"]])
    logging::logwarn("There is missing information in the 'Condition' column, filling it with 'undefined'.")
    raw[missing_replicates, "CONDITION"] <- "undefined"
  }

  ## make letters case-insensitive

  if (any(is.element(colnames(raw), "FRACTION"))) {
    fraction <- "FRACTION"
  } else {
    fraction <- NULL
  }
  if (any(is.element(colnames(raw), "TECHREPLICATE"))) {
    tech.rep <- "TECHREPLICATE"
  } else {
    tech.rep <- NULL
  }

  require.col <- c("PROTEINNAME", "PEPTIDESEQUENCE", "PRECURSORCHARGE",
                   "FRAGMENTION", "PRODUCTCHARGE", "ISOTOPELABELTYPE",
                   "CONDITION", "BIOREPLICATE", "RUN", "INTENSITY", fraction, tech.rep)
  raw_subset <- raw[, require.col]
  ## before remove, get PeptideSequence and combination of PeptideSequence and precursorcharge
  ## for global standard normalization
  tempPeptide <- unique(raw[, c("PEPTIDESEQUENCE", "PRECURSORCHARGE")])
  tempPeptide$PEPTIDE <- paste(tempPeptide$PEPTIDESEQUENCE, tempPeptide$PRECURSORCHARGE, sep="_")
  rm(raw)
  ## assign peptide, transition
  raw_subset <- data.frame(
    raw_subset,
    "PEPTIDE" = paste(raw_subset[["PEPTIDESEQUENCE"]], raw_subset[["PRECURSORCHARGE"]], sep="_"),
    "TRANSITION" = paste(raw_subset[["FRAGMENTION"]], raw_subset[["PRODUCTCHARGE"]], sep="_"))

  if (length(unique(raw_subset$ISOTOPELABELTYPE)) > 2) {
    logging::logerror(
               "There are more than two levels of labeling.  Only label-free and reference-labeled are supported.")
    stop()
  }
  ## change light, heavy -> L,H
  ## [THT: should check if users really provide light/heavy, L/H, l/h, or something else ]
  ## [THT: should also check if users provide only H (instead of L)]
  raw_subset[["ISOTOPELABELTYPE"]] <- factor(raw_subset[["ISOTOPELABELTYPE"]])
  if (nlevels(raw_subset[["ISOTOPELABELTYPE"]]) == 2) {
    levels(raw_subset[["ISOTOPELABELTYPE"]]) <- c("H", "L")
  }
  if (nlevels(raw_subset[["ISOTOPELABELTYPE"]]) == 1) {
    levels(raw_subset[["ISOTOPELABELTYPE"]]) <- c("L")
  }
  if (any(is.element(colnames(raw_subset), "FRACTION"))) {
    fraction <- "FRACTION"
  } else {
    fraction <- NULL
  }
  if (any(is.element(colnames(raw_subset), "TECHREPLICATE"))) {
    tech.rep <- "TECHREPLICATE"
  } else {
    tech.rep <- NULL
  }
  require.col <- c("PROTEINNAME", "PEPTIDE", "TRANSITION", "ISOTOPELABELTYPE",
                   "CONDITION", "BIOREPLICATE","RUN", "INTENSITY", fraction, tech.rep)

  raw_subset <- raw_subset[, require.col]
  if (ncol(raw_subset) == 10 &
      any(is.element(colnames(raw_subset), "FRACTION")) &
      any(is.element(colnames(raw_subset), "TECHREPLICATE"))) {
    ## If you choose to uc everything, then stick to it!
    colnames(raw_subset) <- c("PROTEIN", "PEPTIDE", "TRANSITION", "LABEL",
                              "CONDITION", "SAMPLE", "RUN", "INTENSITY",
                              "FRACTION", "TECHREPLICATE")
  } else if (ncol(raw_subset) == 9 &
             any(is.element(colnames(raw_subset), "FRACTION"))) {
    colnames(raw_subset) <- c("PROTEIN", "PEPTIDE", "TRANSITION", "LABEL",
                              "CONDITION", "SAMPLE", "RUN", "INTENSITY", "FRACTION")
  } else {
    colnames(raw_subset) <- c("PROTEIN", "PEPTIDE", "TRANSITION", "LABEL",
                              "CONDITION", "SAMPLE", "RUN", "INTENSITY")
  }

  ## create work data for quant analysis
  ## -----------------------------------
  raw_subset <- raw_subset[!is.na(raw_subset[["PROTEIN"]]), ]
  raw_subset <- raw_subset[raw_subset[["PROTEIN"]] != "", ]
  ## This data structure receives a lot of abuse!
  ## I am thinking there should just be a function to validate this and then invoke it
  ## periodically through this process, but whatever.
  ## Furthermore, this could be filled in automagically without these if() and whatnot.
  work <- data.frame(
    "PROTEIN" = raw_subset[["PROTEIN"]],
    "PEPTIDE" = raw_subset[["PEPTIDE"]],
    "TRANSITION" = raw_subset[["TRANSITION"]],
    "FEATURE" = paste(raw_subset[["PEPTIDE"]], raw_subset[["TRANSITION"]], sep="_"),
    "LABEL" = raw_subset[["LABEL"]],
    "GROUP_ORIGINAL" = raw_subset[["CONDITION"]],
    "SUBJECT_ORIGINAL" = raw_subset[["SAMPLE"]],
    "RUN" = raw_subset[["RUN"]],
    "GROUP" = 0,
    "SUBJECT" = 0,
    "INTENSITY" = raw_subset[["INTENSITY"]])
  work$GROUP_ORIGINAL <- factor(work$GROUP_ORIGINAL)
  work$SUBJECT_ORIGINAL <- factor(work$SUBJECT_ORIGINAL, levels=unique(work$SUBJECT_ORIGINAL))
  work$LABEL <- factor(work$LABEL, levels=levels(work$LABEL))
  work[work$LABEL=="L", "GROUP"] <- work[work$LABEL=="L", "GROUP_ORIGINAL"]
  work[work$LABEL=="L", "SUBJECT"] <- work[work$LABEL=="L", "SUBJECT_ORIGINAL"]
  work <- data.frame(work, SUBJECT_NESTED=paste(work$GROUP, work$SUBJECT, sep="."))

  if (any(is.element(colnames(raw_subset), "FRACTION"))) {
    work <- data.frame(work, "FRACTION" = raw_subset[["FRACTION"]])
  }
  if (any(is.element(colnames(raw_subset), "TECHREPLICATE"))) {
    work <- data.frame(work, "TECHREPLICATE" = raw_subset[["TECHREPLICATE"]])
  }
  logging::loginfo("Data successfully reformatted for further analyses.")

  ## 2016. 08.29 : replace <1 with zero for log2(intensity)
  work[["ABUNDANCE"]] <- work[["INTENSITY"]]
  if (length(which(!is.na(work[["ABUNDANCE"]]) & work[["ABUNDANCE"]] < 1)) > 0) {
    logging::loginfo("Replacing abundance values < 1 with 1, as they will be put on the log scale.")
    ## Why not just set the NAs to 0, too?  I bet there is a good reason but I do not know it.
    work[!is.na(work[["ABUNDANCE"]]) & work[["ABUNDANCE"]] < 1, "ABUNDANCE"] <- 1
  }

  ## now, INTENSITY keeps original values.
  ## NA means no observation. assume that spectral tools are not report if no observation.
  ## zero means detected but zero.
  ## considered intenseity <1 -> intensity = 1

  ## based on logTrans option, assign log transformation
  ## remove log2 or log10 intensity
  ## [THT: add one more conidtion to have the program complain if a user
  ## provide unexpected value for logTrans]
  if (logTrans == 2) {
    work[["ABUNDANCE"]] <- log2(work[["ABUNDANCE"]])
  } else if (logTrans == 10) {
    work[["ABUNDANCE"]] <- log10(work[["ABUNDANCE"]])
  } else {
    ## Above there was a check for only log 2 and 10, but we can do e if we want.
    ## I might go back up there and remove that check. Long live e! 2.718282 rules!
    work[["ABUNDANCE"]] <- log(work[["ABUNDANCE"]]) / log(logTrans)
  }
  logging::loginfo(paste0("Log ", logTrans, " transformation complete."))

  ## Check multi-method or not : multiple run for a replicate
  work[["RUN"]] <- as.factor(work[["RUN"]])

  ## dot functions offend me, just don't export it!
  checkMultirun <- countMultiRun(work)
  if (checkMultirun$is.risky) {
    logging::logwarn("MSstats suspects that there are fractionations and potentially
technical replicates too. Please add Fraction column in the input.")
  }
  ## Check the various fraction information and set it appropriately.
  work <- check_fractions(work, checkMultirun)

  ## check missingness for multirun
  ## check no value for some feature: balanced structure or not
  ## need to separate label-free or label-based
  logging::loginfo(paste0("The fillincomplete rows option is: ", fillIncompleteRows))
  ## [THT: better to write a function for single method, and call that function
  ## here and for the case with multuple methods]
  ## only 1 method
  ## I agree!

  if (!checkMultirun$out | length(unique(work$FRACTION)) == 1) {
    ## label-free experiments
    if (nlevels(work$LABEL) == 1) {
      work <- extract_singlefrac_labelfree_data(work, fillIncompleteRows=fillIncompleteRows)
    } else {
      work <- extract_singlefrac_labeled_data(work, fillIncompleteRows=fillIncompleteRows)
    } # end 1 method
  } else { # multiple fractionations
    work <- extract_multifrac_data(work)
  } # end multiple fractionations

  ## factorize GROUP, SUBJECT, GROUP_ORIGINAL, SUBJECT_ORIGINAL, SUBJECT_ORIGINAL_NESTED, FEATURE, RUN
  ## -------------------------------------------------------------------------------------------------
  work[["PROTEIN"]] <- as.factor(work[["PROTEIN"]])
  work[["PEPTIDE"]] <- as.factor(work[["PEPTIDE"]])
  work[["TRANSITION"]] <- as.factor(work[["TRANSITION"]])
  work <- work[with(work, order(LABEL, GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN,
                                PROTEIN, PEPTIDE, TRANSITION)), ]
  work[["GROUP"]] <- as.factor(work[["GROUP"]])
  work[["SUBJECT"]] <- as.factor(work[["SUBJECT"]])
  ## SUBJECT_ORIGINAL_NESTED will sorted as GROUP_ORIGINAL, SUBJECT_ORIGINAL
  work[["SUBJECT_NESTED"]] <- factor(work[["SUBJECT_NESTED"]], levels=unique(work[["SUBJECT_NESTED"]]))
  ## FEATURE will sorted as PROTEIN, PEPTIDE, TRANSITION
  work[["FEATURE"]] <- factor(work[["FEATURE"]], levels=unique(work[["FEATURE"]]))
  ## RUN will sorted as GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN
  work[["originalRUN"]] <- work[["RUN"]]
  work[["RUN"]] <- factor(work[["RUN"]],
                          levels=unique(work[["RUN"]]),
                          labels=seq(1, length(unique(work[["RUN"]]))))
  logging::loginfo("Recast the following columns as factors:
group, subject, group_original, subject_original, subject_original_nested, feature, and run.")

  ## Normalization ##
  ## ------------- ##
  ## Normalization : option 0. none
  if (is.element("NONE", normalization) | is.element("FALSE", normalization)) {
    ## after "toupper", FALSE becomes character.
    ## To my mind, that peculiar magical recasting of FALSE by toupper() is either a semi-convincing
    ## reason to do it, or an abomination that should be killed with fire -- I'm not sure which.
    logging::loginfo("Not performing normalization.")
  }

  ## Normalization : option 1. constant normalization , equalize medians ##
  ## -------------------------------------------------------------------
  if (!is.element("NONE", normalization) & !is.element("FALSE", normalization) &
      is.element("EQUALIZEMEDIANS", normalization)) {
    if (nlevels(work[["LABEL"]]) == 1) {
      work <- normalize_singlefrac_by_medians(work)
    }

    if (nlevels(work[["LABEL"]]) == 2 ) {
      work <- normalize_singlefrac_by_heavy_standard(work)
    } # for labe-based
    logging::loginfo("Performing equalize medians normalization.")
  } ## end equaliemedian normalization

  ## Normalization : option 2. quantile normalization ##
  ## ------------------------------------------------ ##
  if (!is.element("NONE", normalization) & !is.element("FALSE", normalization)
      & is.element("QUANTILE", normalization)) {
    if (nlevels(work[["LABEL"]]) == 1) {
      work <- normalize_singlefrac_by_quantiles(work)
    }
    if (nlevels(work[["LABEL"]]) == 2) {
      work <- normalize_singlefrac_by_quantile_heavy(work)
    }
    logging::loginfo("Performing quantile normalization.")
  }

  ## Normalization : option 3. global standards - for endogenous ##
  ## ----------------------------------------------------------- ##
  if (!is.element("NONE", normalization) & !is.element("FALSE", normalization) &
      is.element("GLOBALSTANDARDS", normalization)) {
    work <- normalize_by_globalstandards(work)
  }

  ## ----------------------------------------------------------- ##
  ## if there are multiple method, need to merge after normalization + before feature selection ##
  if (length(unique(work[["FRACTION"]])) > 1) {
    work <- multifraction_normalize_confusing(work)
  }

  ## BetweenRunInterferenceScore
  ## need to make new function
  if (betweenRunInterferenceScore) {
    ## only output light
    l <- subset(work, LABEL=="L")
    ## add ProtFeature and ProtPeptide, because the shared peptides appear in multiple proteins
    l[["ProtFeature"]] <- paste(l[["PROTEIN"]], l[["FEATURE"]], sep="/")
    l[["ProtPeptide"]] <- paste(l[["PROTEIN"]], l[["PEPTIDE"]], sep="/")
    temp <- tapply(
      l[["ABUNDANCE"]],
      l[, c("RUN","ProtPeptide")],
      function(x) {
        mean(x,na.rm=TRUE)
      })

    temp1 <- data.frame(
      "ProtPeptide" = rep(colnames(temp), each=dim(temp)[1]),
      "RUN" = rep(rownames(temp), dim(temp)[2]),
      "meanPEPTIDE" = as.numeric(unlist(temp)))
    temp2 <- merge(l[, c("PROTEIN", "PEPTIDE", "FEATURE", "ProtPeptide",
                         "ProtFeature", "RUN", "ABUNDANCE")],
                   temp1, by=c("ProtPeptide","RUN"))
    temp3 <- temp2[!is.na(temp2[["ABUNDANCE"]]), ]
    temp4 <- tapply(
      rownames(temp3),
      temp3[, c("ProtFeature")],
      function(x) {
        cor(temp3[x, "ABUNDANCE"], temp3[x, "meanPEPTIDE"])
      })
    names <- unique(temp2[, c("PROTEIN", "PEPTIDE", "FEATURE", "ProtFeature")])
    names <- names[with(names, order(ProtFeature)), ]
    BetweenRunInterferenceFile <- data.frame(names[, c("PROTEIN","PEPTIDE","FEATURE")],
                                             "BetweenRunInterferenceScore" = temp4)
    BetweenRunInterferenceFile <- BetweenRunInterferenceFile[with(BetweenRunInterferenceFile,
                                                                  order(PROTEIN,PEPTIDE,FEATURE)), ]
    logging::loginfo("Between run interference score was calculated and saved to a csv file.")
  } else {
    logging::loginfo("Between run interference score was not calculated.")
  }

  ## Below two lines were merely for in-house testing and comparisons when needed
  ## work.NoImpute <- work
  ## AbundanceAfterImpute <- Imputation(work, cutoffCensored, censoredInt,
  ##         remove50missing, MBimpute, original_scale)

  ## ------------- ##
  ## how to decide censored or not
  ## ------------- ##

  ## after normalization, zero intensity could be negative
  work[!is.na(work[["ABUNDANCE"]]) & work[["ABUNDANCE"]] < 0, "ABUNDANCE"] <- 0
  work[!is.na(work[["INTENSITY"]]) & work[["INTENSITY"]] == 1, "ABUNDANCE"] <- 0
  ## If imputation=TRUE and there is any value for maxQuantileforCensored, apply cutoff for censored missing
  if (summaryMethod == "TMP" & MBimpute) {
    work <- summarize_tmp_mbimpute(work, maxQuantileforCensored=maxQuantileforCensored,
                                   censoredInt=censoredInt)
  }

  ## ------------- ##
  ## featureSubset ##
  ## ------------- ##
  ## !! need to decide how to present : keep original all data and make new column to mark, or just present selected subset

  ## Take a copy of the data before removing anything, this will be used later.
  full_dataset <- work

  if (featureSubset == "all") {
    logging::loginfo("Feature Subset: using all features in the data set.")
  }
  if (featureSubset == "highQuality") {
    logging::loginfo("Feature Subset: selecting high quality features defaults to top3 for the moment.")
    featureSubset <- "top3"

    ##message("* Use feature selection algorithm in order to remove features with interference.")
    ## 2016.04.25. MC
    ## there is the possibility to remain features which have completely missing in the certain condition after imputation
    ## Therefore, remove the features which are completely missing in the certain condition before imputation
    ## !! need to check for zero for skylineReport=TRUE
    ##work$remove <- FALSE

    ##Impute the missing valuess before feature selection
    ## It should be handle within feature_selection function using runQuantification function
    ##AbundanceAfterImpute <- Imputation(work, cutoffCensored, censoredInt, MBimpute, remove50missing)
    ##work <- AbundanceAfterImpute
    ##
    ##selection_list <- feature_selection(work,
    ##                   cutoffCensored,
    ##                   censoredInt,
    ##                   remove_proteins_with_interference)
    ##SelectionAfterImpute <- work
    ## 20160425-MC : after selecting feature, original ABUNCANCE should be used.
    ##work$ABUNDANCE <- work$ABUNDANCE.O
    ##work <- work[, -which(colnames(work) %in% c("ABUNDANCE.O", "feature.label", "run.label", "cen", "pred", "ref", "Protein_Peptide"))]
    ## if (cutoff_peptidelist == "distance.perProtein") {
    ## work[which(work$PEPTIDE %in% c(selection_list$PeptideList.remove.perProteinCutoff, selection_list$PeptideList.issue)), "remove"] <- TRUE
    ##}
    ##if (cutoff_peptidelist == "distance.amongAllProteins") {
    ## work[which(work$PEPTIDE %in% c(selection_list$PeptideList.remove.allProteinCutoff, selection_list$PeptideList.issue)), "remove"] <- TRUE
    ##}
    ##if (cutoff_peptidelist == "error.amongAllProteins") {
    ## tmp <- selection_list$Model.Based.Error
    ## tmp <- tmp[tmp$Model.Based.Error > 0, ]
    ## cutoff_error <- quantile(tmp$Model.Based.Error, prob=cutoff_percentile)
    ## removePeptide <- tmp[tmp$Model.Based.Error > cutoff_error, "Peptide"]
    ## work[which(work$PEPTIDE %in% c(removePeptide)), "remove"] <- TRUE
    ##}
  }

  if (featureSubset == "top3") {
    work <- subset_top3_result(work)
  }

  if (featureSubset == "topN") {
    work <- subset_topn_result(work)
  }

  ## check missingness
  ## transitions are completely missing in one condition : missingness ##
  final.decision <- decide_missingness(work)

  ## output : summary ##
  ## ---------------- ##
  ## output for label
  logging::loginfo(paste0(
             "Summary: ", length(unique(work[["LABEL"]])), " levels of isotope labeling were observed."))
  temp <- data.frame("Summary of Features :")
  colnames(temp) <- " "
  rownames(temp) <- " "
  summary.f <- matrix(NA, nrow=3)
  summary.f[1] <- nlevels(work[["PROTEIN"]])
  temp <- unique(work[, c("PROTEIN", "PEPTIDE")])
  temp1 <- xtabs(~ PROTEIN, data=temp)
  temp2 <- summary(as.numeric(temp1))
  summary.f[2] <- paste(temp2["Min."], temp2["Max."], sep="-")
  temp <- unique(work[, c("PEPTIDE", "FEATURE")])
  temp1 <- xtabs(~ PEPTIDE, data=temp)
  temp2 <- summary(as.numeric(temp1))
  summary.f[3] <- paste(temp2["Min."], temp2["Max."], sep="-")
  colnames(summary.f) <- "count"
  rownames(summary.f) <- c("Number of Proteins", "Number of Peptides/Protein", "Number of Transitions/Peptide")

  ## output for process
  logging::loginfo("Summary of Features: ")
  logging::loginfo(paste0(rownames(summary.f)[1]," : ", summary.f[1]))
  logging::loginfo(paste0(rownames(summary.f)[2]," : ", summary.f[2]))
  logging::loginfo(paste0(rownames(summary.f)[3]," : ", summary.f[3]))

  ## protein list with 1 feature
  temp <- unique(work[, c("PROTEIN", "FEATURE")])
  temp1 <- xtabs(~ PROTEIN, data=temp)
  temp2 <- as.data.frame(temp1[temp1 == 1])
  if (nrow(temp2) > 0) {
    logging::loginfo(paste0("Proteins: ", toString(rownames(temp2)),
                            " have only single transitions, consider removing them from the dataset."))
  }
  temp <- data.frame("Summary of Samples :")
  colnames(temp) <- " "
  rownames(temp) <- " "
  print(temp)
  summary.s <- matrix(NA, ncol=nlevels(work[["GROUP_ORIGINAL"]]), nrow=3)
  ## # of MS runs
  temp <- unique(work[, c("GROUP_ORIGINAL", "RUN")])
  temp1 <- xtabs(~ GROUP_ORIGINAL, data=temp)
  summary.s[1, ] <- temp1
  ## # of biological replicates
  temp <- unique(work[, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL")])
  temp1 <- xtabs(~ GROUP_ORIGINAL, data=temp)
  summary.s[2, ] <- temp1
  ## # of technical replicates
  c.tech <- round(summary.s[1, ] / (summary.s[2, ] * length(unique(work[["FRACTION"]]))))
  ##summary.s[3, ] <- ifelse(c.tech==1, 0, c.tech)
  summary.s[3, ] <- c.tech
  colnames(summary.s) <- unique(work[["GROUP_ORIGINAL"]])
  rownames(summary.s) <- c("Number of MS runs",
                           "Number of Biological Replicates",
                           "Number of Technical Replicates")
  logging::loginfo(knitr::kable(summary.s))
  logging::loginfo(paste0("Missingness summary: ",
                          sum(final.decision != 0),
                          " are missing completely in one condition."))
  if (sum(final.decision != 0) != 0) {
    tmp.final <- final.decision[final.decision != 0]
    if (length(tmp.final) > 5) {
      logging::loginfo(paste0("  -> ", toString(names(tmp.final[1:5])), " ..."))
    } else {
      logging::loginfo(paste0("  -> ", toString(names(tmp.final)), " ..."))
    }
    rm(tmp.final)
  }

  without <- xtabs(~ RUN, work)
  withall <- xtabs(~ RUN, full_dataset)
  run.missing <- without / withall
  logging::loginfo(paste0("Missingness summary: ", sum(run.missing < 0.25), " are missing 75% observations."))
  if (sum(run.missing < 0.25) != 0) {
    logging::loginfo(paste0("  -> ", paste("RUN", names(without[run.missing < 0.25]), sep=" ")))
  }

  ## I am not sure how this is different from what is above?
  if (sum(final.decision != 0) != 0) {
    tmp.final <- final.decision[final.decision != 0]
    if (length(tmp.final) > 5 ) {
      logging::loginfo(paste0("  -> ", toString(names(tmp.final[1:5])), " ..."))
    } else {
      logging::loginfo(paste0("  -> ", toString(names(tmp.final)), " ..."))
    }
    rm(tmp.final)
  }

  if (sum(run.missing < 0.25) != 0) {
    logging::loginfo(paste0("  -> ", paste("RUN", names(without[run.missing < 0.25]), sep=" ")))
  }

  ## check any protein has only light for labeled-experiment
  if (nlevels(work[["LABEL"]]) == 2) {
    temp <- unique(work[, c("PROTEIN", "LABEL")])
    temp1 <- xtabs(~ PROTEIN, data=temp)
    if (any(temp1 != 2)) {
      ## check that is L or H
      namepro <- names(temp1[temp1 != 2])
      for (j in 1:length(namepro)) {
        if (unique(work[work[["PROTEIN"]] == namepro[j], "LABEL"]) == "L") {
          logging::loginfo(paste0(namepro[j], "has only endogeneous intensities in label-based experiment.
Please check this protein or remove it."))
        }
        if (unique(work[work[["PROTEIN"]] == namepro[j], "LABEL"]) == "H") {
          logging::loginfo(paste0(namepro[j], " has only reference intensities in label-based experiment.
Please check this protein or remove it."))
        }
      }
    }
  } ## End checking the nlevels for work$LABEL

  logging::loginfo("Processing data for analysis is complete.")
  ## get the summarization per subplot (per RUN)
  ## -------------------------------------------
  rqresult <- try(runQuantification(work, summaryMethod=summaryMethod, equalFeatureVar=equalFeatureVar,
                                    cutoffCensored=cutoffCensored, censoredInt=censoredInt,
                                    remove50missing=remove50missing, MBimpute=MBimpute,
                                    original_scale=original_scale, logsum=FALSE,
                                    featureSubset=featureSubset, message.show=FALSE,
                                    clusters=clusters))
  if (class(rqresult) == "try-error") {
    logging::logwarn(paste0("Cannot summarize per subplot with method: ", summaryMethod, "."))
    rqall <- NULL
    rqmodelqc <- NULL
    workpred <- NULL
  } else {
    label <- nlevels(work[["LABEL"]]) == 2
    if (sum(is.element(colnames(rqresult[["rqdata"]]), "RUN")) == 0) {
      ## logsum is summarization per subject
      lab <- unique(work[, c("GROUP", "GROUP_ORIGINAL", "SUBJECT_ORIGINAL",
                             "SUBJECT_NESTED", "SUBJECT")])
      if (label) {
        lab <- lab[lab[["GROUP"]] != 0, ]
      }
      rqall <- merge(rqresult[["rqdata"]], lab, by="SUBJECT_ORIGINAL")
    } else {
      lab <- unique(work[, c("RUN", "originalRUN", "GROUP", "GROUP_ORIGINAL", "SUBJECT_ORIGINAL",
                             "SUBJECT_NESTED", "SUBJECT")])
      if (label) {
        lab <- lab[lab[["GROUP"]] != 0, ]
      }
      rqall <- merge(rqresult[["rqdata"]], lab, by="RUN")
    }
    rqall[["GROUP"]] <- factor(rqall[["GROUP"]])
    rqall[["Protein"]] <- factor(rqall[["Protein"]])
    rqmodelqc <- rqresult[["ModelQC"]]
    ##MC : can"t use this predicted value.
    ##workpred <- rqresult$PredictedBySurvival
    workpred <- NULL
    logging::loginfo(paste0("The summarization per subplot by method: ", summaryMethod, " is finished."))
  }

  ## return work data.frame and run quantification
  ##Align the run quantification data
  if (any(is.element(colnames(rqall), "RUN"))) {
    rqall <- rqall[order(rqall$Protein, as.numeric(as.character(rqall[["RUN"]]))), ]
    rownames(rqall) <- NULL
  }

  ##Mike: Below is for in-house verification occasionally
  ##processedquant <- list(ProcessedData=work.NoImpute, RunlevelData=rqall,
  ##  SummaryMethod=summaryMethod, ModelQC=rqmodelqc, PredictBySurvival=workpred,
  ## ImputedData=AbundanceAfterImpute)
  processedquant <- list(
    "ProcessedData" = work,
    "RunlevelData" = rqall,
    "SummaryMethod" = summaryMethod,
    "ModelQC" = rqmodelqc,
    "PredictBySurvival" = workpred)
  return(processedquant)
}

#' Manual function allowing foreach to return a list of multiple variables
#' Neat!
#'
#' @param x  A foreach result to recast
#' @param ... The rest of the stuff to recast, I am guessing.
#' @return  A list of variables rather than foreach's mess
resultsAsLists <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))

  ## I think the above may be written as the following, which is the limit
  ## of my intelligence to understand.  Smart people making things too short for stupid people (me)
  ## to understand makes me sad.
  ## lapply(seq_along(x),
  ##        function(i) {
  ##           c(x[[i]], lapply(list(...), function(y) {
  ##             y[[i]]
  ##           })
  ##           ) ## Completes the c(), how interesting. I would never think of this.
  ##         })
}

#' Perform the quantification for the sanitized data
#'
#' I haven't read this function yet, so I have no clue what it does, that should be amended.
#' So for the moment I am going to put stupid strings into the parameter list
#' docstrings until I understand what is going on.  Holy crap, I thought
#' dataProcess() was long and obnoxious, this function is 1,500 lines of dense
#' stuff.  Sometimes I hate smart people, a stupider person would have split
#' this into smaller, understandable pieces.
#'
#' @param data  Yay!
#' @param summaryMethod  Which method to use for summarizing
#' @param equalFeatureVar  Set vars for equal features?
#' @param cutoffCensored  Remove censored values?
#' @param censoredInt  Int for censored!
#' @param remove50missing  Yeah I hate those 50 missing guys, too.
#' @param MBimpute  Impute Impute!
#' @param original_scale  ok, so this bugs me, camelCase, snake_case,
#'   dotted.names, pick one to use most of the time!
#' @param logsum oh yeah, or just don't do anything to separate words because
#'   consistency is the hobgoblin of little minds.
#' @param featureSubset My eyes don't focus well enough to do camelCase
#'   efficiently.  So I am sad that you came back this way.
#' @param message.show  and I think using a '.' as a space is nasty.  ick.
#' @param clusters  for my money, I would go with cluster_size, or cpus or
#'   something.
#' @return some new data?
runQuantification <- function(data, summaryMethod="TMP", equalFeatureVar=TRUE,
                              cutoffCensored="minFeature", censoredInt="NA", remove50missing=FALSE,
                              MBimpute=TRUE, logsum=FALSE, featureSubset="all",
                              message.show=FALSE, clusters=1, original_scale=FALSE,
                              pct_adjust=0.99) {

  ##Since the imputation has been done before feature selection, delete the columns of censoring indicator to avoid imputing the same intensity again
  ##if(featureSubset == "highQuality") {
  ##    data$cen <- NULL; data$pred <- NULL; data$INTENSITY <- 2^data$ABUNDANCE
  ##}
  ##If we want to impute again after the feature selection
  ##if(featureSubset == "highQuality" & ImputeAgain==TRUE) {
  ##    data$ABUNDANCE <- data$ABUNDANCE.O
  ##}
  ## if there is "remove" column, remove TRUE
  ## 2016. 08.29 should change column name for this remove variable. from feature selection??
  if (any(is.element(colnames(data), "remove"))) {
    data <- data[!data[["remove"]], ]
  }
  data[["LABEL"]] <- factor(data[["LABEL"]])
  label <- nlevels(data[["LABEL"]]) == 2

  ## set ref which is distinguish reference and endogenous. any reference=0. endogenous is the same as RUN
  if (label) {
    data[["ref"]] <- 0
    data$ref[data$LABEL != "H"] <- data$RUN[data$LABEL != "H"]
    data$ref <- factor(data$ref)
    ##      unique(data[,c("RUN","LABEL","GROUP","ref")])
  }

  ##  finalresult <- data.frame(Protein=rep(levels(data$PROTEIN),each=nlevels(data$RUN)),
  ##      RUN=rep(c(levels(data$RUN)),nlevels(data$PROTEIN)),Condition=NA,
  ##      BioReplicate=NA,LogIntensities=NA,NumFeature=NA,NumPeaks=NA)
  ## for saving predicting value for impute option
  predAbundance <- NULL
  ##
  ## method 1 : model based summarization
  ##
  if (summaryMethod == "linear" & is.null(censoredInt)) {
    data <- data[!is.na(data$ABUNDANCE), ]
    data$PROTEIN <- factor(data$PROTEIN)
    data$RUN <- factor(data$RUN)
    result <- NULL
    dataafterfit <- NULL
    for (i in 1: nlevels(data$PROTEIN)) {
      sub <- data[data$PROTEIN==levels(data$PROTEIN)[i],]
      sub$SUBJECT_NESTED <- factor(sub$SUBJECT_NESTED)
      sub$FEATURE <- factor(sub$FEATURE)
      sub$RUN <- factor(sub$RUN)
      if (!label) {
        temp <- data.frame(xtabs(~ RUN, data=sub))
        sub.result <- data.frame(
          "Protein" = rep(unique(sub$PROTEIN), each=nlevels(sub$RUN)),
          "RUN" = rep(c(levels(sub$RUN)),1),
          "LogIntensities" = NA,
          "NumFeature" = length(unique(sub$FEATURE)),
          "NumPeaks" = temp$Freq)
      } else {
        sub$ref <- factor(sub$ref)
        temp <- data.frame(xtabs(~ ref, data=sub))
        sub.result <- data.frame(
          "Protein" = rep(levels(data$PROTEIN)[i], each=nlevels(sub$ref)),
          "RUN" = rep(c(levels(sub$ref)[-1], "Ref"), 1),
          "LogIntensities" = NA,
          "NumFeature" = length(unique(sub$FEATURE)),
          "NumPeaks" = c(temp[-1, "Freq"], temp[1, "Freq"]))
      }
      singleFeature <- checkSingleFeature(sub)
      singleSubject <- checkSingleSubject(sub)
      TechReplicate <- checkTechReplicate(sub) ## use for label-free model
      ## fit the model
      ##if (message.show) {
      message(paste0("Getting the summarization per subplot for protein ",
                     unique(sub$PROTEIN), "(",i," of ",length(unique(data$PROTEIN)),")"))
      ##}

      fit <- try(.fit.quantification.run(sub, singleFeature, singleSubject, TechReplicate,
                                         labeled=label, equalFeatureVar), silent=TRUE)
      if (class(fit) == "try-error") {
        message(paste0("*** error : can't fit the model for ", levels(data$PROTEIN)[i]))
        result <- rbind(result, sub.result)
        if (nrow(sub) != 0) {
          sub$residuals <- NA
          sub$fitted <- NA
        }
      } else {
        if (class(fit) == "lm") {
          cf <- summary(fit)$coefficients
        } else {
          cf <- fixef(fit)
        }

        ## calculate sample quantification for all levels of sample
        a <- 1
        for (j in 1:nlevels(sub$RUN)) {
          contrast.matrix <- rep(0, nlevels(sub$RUN))
          contrast.matrix[j] <- 1
          contrast <- make.contrast.run.quantification(fit, contrast.matrix, sub, labeled=label)
          if (class(fit) == "lm") {
            sub.result[a, 3] <- estimableFixedQuantification(cf, contrast)
          } else {
            sub.result[a, 3] <- estimableRandomQuantification(cf, contrast)
          }
          a <- a + 1
        }

        ## for label-based case, need reference quantification
        if (label) {
          contrast <- make.contrast.run.quantification.reference(fit, contrast.matrix, sub)
          if (class(fit) == "lm") {
            sub.result[a, 3] <- estimableFixedQuantification(cf, contrast)
          } else {
            sub.result[a, 3] <- estimableRandomQuantification(cf, contrast)
          }
        }

        result <- rbind(result, sub.result)
        if (class(fit) == "lm") { ### lm model
          sub$residuals <- fit$residuals
          sub$fitted <- fit$fitted.values
        } else {  ### lmer model
          sub$residuals <- resid(fit)
          sub$fitted <- fitted(fit)
        }

        dataafterfit <- rbind(dataafterfit, sub)
      }
    } ## end-loop for each protein
  } ## for linear model summary

  ##
  ## Method 2 : Tukey Median Polish
  ##
  if (summaryMethod == "TMP") {
    ##data <- data[!is.na(data$ABUNDANCE),]
    data[["PROTEIN"]] <- factor(data[["PROTEIN"]])
    data[["RUN"]] <- factor(data[["RUN"]])

    result_list <- perform_tukey_polish_median(
      data,
      label=label,
      pct_adjust=0.99,
      clusters=clusters,
      MBimpute=MBimpute,
      censoredInt=censoredInt,
      message.show=message.show,
      cutoffCensored=cutoffCensored,
      original_scale=original_scale,
      remove50missing=remove50missing)

    result <- data.frame()
    for (i in 1:length(result_list)) {
      result_class <- class(result_list[[i]])
      if (result_class == "data.frame") {
        result <- rbind(result, result_list[[i]])
      }
    }
    ##predAbundance <- predAbundance[-which(duplicated(predAbundance))] # remove duplicates
    dataafterfit <- NULL
  } else if (summaryMethod == "linear") {
    ##
    ## method 4 : survival model for censored missing values
    ##
    result <- perform_survival_censored(data, summaryMethod=summaryMethod,
                                        equalFeatureVar=equalFeatureVar,
                                        cutoffCensored=cutoffCensored,
                                        censoredInt=censoredInt,
                                        remove50missing=remove50missing, MBimpute=MBimpute,
                                        logsum=logsum, featureSubset=featureSubset,
                                        message.show=message.show, clusters=clusters,
                                        original_scale=original_scale, pct_adjust=0.99)
  }
  ## final result
  finalout <- list(
    "rqdata" = result,
    "ModelQC" = dataafterfit,
    "PredictedBySurvival" = predAbundance)
  return(finalout)
}

#' Use lm to perform the actual quantification
#'
#' @param sub yeah I dunnot yet, I need to read the context of this function.
#' @param singleFeature  A single feature, I am guessing.
#' @param singleSubject and a subject
#' @param TechReplicate  look, if you are going to do camelCase, then why
#'   CamelCase?  I bet you have some smart reason, but it seems like fuckery to me.
#' @param labeled whatever
#' @param equalFeatureVar set var for the equal features?
#' @return a fit from lm
fit.quantification.run <- function(sub, singleFeature, singleSubject,
                                   TechReplicate, labeled, equalFeatureVar) {
  if (!labeled) { ### label-free case
    ## for single Feature, original value is the run quantification
    if (singleFeature) {
      fit.full <- lm(ABUNDANCE ~ RUN, data=sub)
    } else {
      fit.full <- lm(ABUNDANCE ~ FEATURE + RUN, data=sub)
    }
  } else { ### labeled-based case
    ## update v3
    if (singleFeature) {
      fit.full <- lm(ABUNDANCE ~ RUN + ref, data=sub)
    }else{ ## several subjects
      fit.full <- lm(ABUNDANCE ~ FEATURE + RUN + ref, data=sub)
    }
  }
  ## make equal variance for feature : need to be updated
  if (!equalFeatureVar) {
    fit.full <- iter.wls.fit.model(data=sub, fit=fit.full, nrepeats=1)
  }
  return(fit.full)
}

#' Examine the work data structure for viable fractions.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @param multi_run  The result from checkmultirun()
#' @return Send the work data back to dataProcess()
check_fractions <- function(work, multi_run) {
  if (multi_run[["out"]]) {
    if (any(is.element(colnames(work), "FRACTION"))) {
      logging::loginfo(paste0("Multiple fractions were found; ",
                              length(unique(work[["FRACTION"]])), " fractions per MS replicate."))
    } else {
      ## need to make new column "Fraction"
      ## each sample has no technical replicate, all runs are fractionated MS runs.
      work[["FRACTION"]] <- NA
      info <- unique(work[, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "RUN")])
      info$condition <- paste(info$GROUP_ORIGINAL, info$SUBJECT_ORIGINAL, sep="_")
      tmp <- work[!is.na(work$ABUNDANCE), ]

      ## get on one sample first
      info.sample1 <- info[info$condition == unique(info$condition)[1], ]
      ## assign fraction info first
      info.sample1$FRACTION <- seq(1, nrow(info.sample1))
      for (k in 1:length(unique(info.sample1$RUN))) {
        ## then fine the same fraction for next sample
        unique.feature <- unique(tmp[tmp$RUN %in% info.sample1$RUN[k], "FEATURE"])
        ## I like silly tmp variables as much as the next guy, but really!?
        tmptmp <- tmp[which(tmp$FEATURE %in% unique.feature), ]
        tmptmp$condition <- paste(tmptmp$GROUP_ORIGINAL, tmptmp$SUBJECT_ORIGINAL, sep="_")
        count.feature <- reshape2::dcast(RUN ~ GROUP_ORIGINAL + SUBJECT_ORIGINAL,
                                         data=tmptmp, fun.aggregate=length, value.var="ABUNDANCE")

        ## !! get one run which has maximum overlapped feature by each sample
        same.frac <- apply(count.feature[,-which(colnames(count.feature) %in% c("RUN"))], 2,
                           function(x) count.feature[which.max(x), "RUN"])
        work[which(work$RUN %in% same.frac), "FRACTION"] <-
          info.sample1[which(info.sample1$RUN %in% info.sample1$RUN[k]), "FRACTION"]
      }
      rm(tmp)

      ## final check up
      checkup <- sum(is.na(unique(work$FRACTION))) > 0
      if (!checkup) {
        logging::loginfo(paste0("Multiple fractions detected; ",
                                length(unique(work$FRACTION)), " fractions per MS replicate."))
      } else {
        logging::logwarn(paste0(
                   "It is hard to find the same fractionation across samples.
Please add a fraction column to the input. Sleeping for 5 seconds to let you reconsider, then
filling fraction with 1 and presumably failing later."))
        Sys.sleep(5)
        work[["FRACTION"]] <- 1
      }
    }

    ##
    ## need an additional step to remove overlapped features across several fractions.
    ##
    if (length(unique(work$FRACTION)) > 1) {
      tmp <- work[!is.na(work$ABUNDANCE) & work$ABUNDANCE > 1, ]
      count.feature <- reshape2::dcast(FEATURE ~ FRACTION, data=tmp,
                                       fun.aggregate=length, value.var="ABUNDANCE")
      ## 1. first, keep features which are measured in one fraction
      count.fraction <- apply(
        count.feature[, -which(colnames(count.feature) %in% c("FEATURE"))],
        1,
        function(x) {
          sum(x > 0)
        })
      keep.feature <- count.feature[count.fraction == 1, "FEATURE"]
      ## 2. second, if features are measured in multiple fractionations,
      ## use the fractionation with maximum number of measurements.
      ## if there are multiple maximum number of measurements, remove features completely.
      first_feature_count <- count.feature[count.fraction > 1, ]

      if (nrow(first_feature_count) > 0) {
        count.fraction <- apply(
          first_feature_count[, -which(colnames(first_feature_count) %in% c("FEATURE"))],
          1,
          function(x) {
            sum(x == max(x))
          })
        count.feature2 <- first_feature_count[count.fraction == 1, ] ## except maximum,

        ##  otherone replace with na
        if (nrow(count.feature2) > 0) {
          remove.fraction <- apply(
            count.feature2,
            1,
            function(x) {
              paste(x[1], names(x[-1])[x[-1] != max(x[-1]) & x[-1] != 0], sep="_")
            })
          remove.fraction <- unlist(remove.fraction)
          work[work$tmp %in% remove.fraction, "INTENSITY"] <- NA
          work[work$tmp %in% remove.fraction, "ABUNDANCE"] <- NA
        }

        ## 3. third, Then check whether there are multiple maximum number
        ## of measurements across fractionation
        count.feature3 <- first_feature_count[count.fraction > 1, ]
        if (nrow(count.feature3) > 0) {
          ## 3.1 : maximum number of measurement == 1, remove that feature
          max.feature <- apply(
            count.feature3[, -which(colnames(count.feature3) %in% c("FEATURE"))],
            1,
            function(x) {
              max(x)
            })
          max.feature.1 <- count.feature3[max.feature == 1, "FEATURE"]
          work <- work[-which(work$FEATURE %in% max.feature.1), ]
          count.feature3 <- count.feature3[-which(count.feature3$FEATURE %in% max.feature.1), ]
          if (nrow(count.feature3) > 0) {
            ## 3.2 : remove fractionations which have not maximum number of measurements
            remove.fraction <- apply(
              count.feature3,
              1,
              function(x) {
                paste(x[1], names(x[-1])[x[-1] != max(x[-1]) & x[-1] != 0], sep="_" )
              })
            remove.fraction <- unlist(remove.fraction)
            work[work$tmp %in% remove.fraction, "INTENSITY"] <- NA
            work[work$tmp %in% remove.fraction, "ABUNDANCE"] <- NA

            ## 3.3 : among fractionations, keep one fractionation which has maximum average
            tmptmp <- work[which(work$FEATURE %in% count.feature3$FEATURE), ]
            tmptmp <- tmptmp[!is.na(tmptmp$ABUNDANCE), ]
            mean.frac.feature <- aggregate(ABUNDANCE ~ FEATURE + tmp, data=tmptmp, FUN=mean)
            mean.frac.feature$FEATURE <- factor(mean.frac.feature$FEATURE)
            temp2 <- split(mean.frac.feature, mean.frac.feature$FEATURE)
            temp3 <- lapply(
              temp2,
              function(x) {
                x <- x[order(x$ABUNDANCE, decreasing=TRUE),]
                x <- x$tmp[2:nrow(x)] # get the fractionation which should be removed.
              })

            remove.fraction <- unlist(temp3, use.names=FALSE)
            work[work$tmp %in% remove.fraction, "INTENSITY"] <- NA
            work[work$tmp %in% remove.fraction, "ABUNDANCE"] <- NA
            rm(tmp)
            rm(tmptmp)
            rm(temp2)
            rm(temp3)
          }
        }
      }
      work <- work[, -which(colnames(work) %in% c("tmp"))]
    }
  } else { ## no fractionation
    work$FRACTION <- 1
  }
  return(work)
}

#'  This counts how many runs are in a data set.
#'
#'  If there are multiple runs, then we have some normalization options available.
#'
#'  @param data  The data to count.
#'  @return  A list containing information whether this is risky or not.
countMultiRun <- function(data) {
  ## if some feature are missing for this spedific run, it could be error. that is why we need balanced design.
  ## with balanced design (fill in NA in each row), it should have different unique number of measurments per fractionation
  ## change it
  ## 2017.05 24 : however, after going through converter functions,
  ## with balanced design, impossible to detect fractionation
  is.risky <- FALSE
  ## if there is fraction info and multiple value for fraction column, we don"t need to count here.
  if(any(is.element(colnames(data), "FRACTION"))) {
    ## already fraction info are available. there are multiple runs.
    out <- TRUE
  } else {
    ## there is no fraction information. First chec, whether there are tech replicates or not.
    info <- unique(data[, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "RUN")])
    info$condition <- paste(info$GROUP_ORIGINAL, info$SUBJECT_ORIGINAL, sep="_")
    count.tech <- xtabs(~ condition, info)
    is.multiplerun <- any(count.tech > 1)
    if (!is.multiplerun) {
      ## only one run for condition*bioreplicate -> no tech replicate at all, no multiple run.
      out <- FALSE
    } else {
      ## need to distinguish whether technical replicate or multiple runs.
      ## For one specific sample, most of features are measured across runs -> tech replicate
      ## if not, multiple runs.
      tmp <- data[!is.na(data$ABUNDANCE), ]
      ## get on sample information
      info.sample1 <- info[info$condition == unique(info$condition)[1], ]
      tmp.sample1 <- tmp[tmp$GROUP_ORIGINAL == unique(info.sample1$GROUP_ORIGINAL) &
                         tmp$SUBJECT_ORIGINAL == unique(info.sample1$SUBJECT_ORIGINAL), ]
      standardFeature <- unique(tmp.sample1[tmp.sample1$RUN == unique(tmp.sample1$RUN[1]),
                                            "FEATURE"])
      tmp.sample1$RUN <- factor(tmp.sample1$RUN)
      ## get overlapped feature ID
      countdiff <- tapply(tmp.sample1$FEATURE,
                          tmp.sample1$RUN,
                          function(x) length(intersect(unique(x), standardFeature)))
      per.overlap.feature <- (countdiff)[-1] / max(unique(countdiff))
      ## first, technical replicate, no fraction :
      ## all runs should have more than 50% features should be the same.
      if (all(per.overlap.feature > 0.5)) { ## then there are technical replicates
        out <- FALSE
      } else if (all(per.overlap.feature < 0.5)) {
        out <- TRUE
      } else {
        ## hard to distinguish fractionation automatically. need information
        ## if both fractionation and technical replicates are there, can"t distinguish.
        ## need fractionation info. + even technical replicate information is needed.
        out <- FALSE
        is.risky <- TRUE
      }
    }
  }
  result <- list(
    "out" = out,
    "is.risky" = is.risky)
  return(result)
}

#' Examine the work data structure single fraction labeled data.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @param fillIncompleteRows  Fill those rows if they are not finished.
#' @return Send the work data back to dataProcess()
extract_singlefrac_labeled_data <- function(work, fillIncompleteRows=TRUE) {
  ## label-based experiment
  ## count the reference and endobenous separately
  work.l <- work[work$LABEL == "L", ]
  work.h <- work[work$LABEL == "H", ]

  ## get feature by Run count of data
  structure.l <- tapply(work.l$ABUNDANCE,
                        list(work.l$FEATURE, work.l$RUN),
                        function(x) length(x))
  structure.h <- tapply(work.h$ABUNDANCE,
                        list(work.h$FEATURE, work.h$RUN),
                        function(x) length(x))

  ## first, check some features which completely missing across run
  missingcomplete.l <- NULL
  missingcomplete.h <- NULL
  ## 1. reference peptides
  featurestructure.h <- apply(structure.h, 1, function(x) sum(is.na(x)))
  ## get feature ID of reference which are completely missing across run
  featureID.h <- names(featurestructure.h[featurestructure.h == ncol(structure.h)])
  if (length(featureID.h) > 0) {
    ## print message
    print_string <- paste0(
      "CAUTION : some REFERENCE features have missing intensities in all
 the runs. The completely missing REFERENCE features are ",
 paste(featureID.h, collapse=", "),
 ". Please check whether features in the list are correctly
 generated from spectral processing tool. \n")
    logging::logwarn(print_string)

    ## add missing rows if option is TRUE
    if (fillIncompleteRows) {
      ## get unique Run information
      nameID <- unique(work.h[, c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT",
                                  "SUBJECT_NESTED", "RUN", "FRACTION")])
      ## get PROTEIN and FEATURE information
      ## here use whole work dataset
      tempTogetfeature <- work[which(work$FEATURE %in% featureID.h), ]
      tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])

      ## then generate data.frame for missingness,
      ##for(j in 1:nrow(nameID)) {
      ##    ## merge feature info and run info as "work" format
      ##    tempmissingwork <- data.frame(tempfeatureID, LABEL="H",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL[j], SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL[j], RUN=nameID$RUN[j], GROUP=nameID$GROUP[j], SUBJECT=nameID$SUBJECT[j], SUBJECT_NESTED=nameID$SUBJECT_NESTED[j], INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD[j])
      ##    ## merge with tempary space, missingwork
      ##    missingcomplete.h <- rbind(missingcomplete.h, tempmissingwork)
      ##}

      ## MC : 2016.04.21 : use merge for simplicity
      tmp <- merge(nameID, tempfeatureID, by=NULL)
      missingcomplete.h <- data.frame(
        "PROTEIN" = tmp$PROTEIN,
        "PEPTIDE" = tmp$PEPTIDE,
        "TRANSITION" = tmp$TRANSITION,
        "FEATURE" = tmp$FEATURE,
        "LABEL" = "H",
        "GROUP_ORIGINAL" = tmp$GROUP_ORIGINAL,
        "SUBJECT_ORIGINAL" = tmp$SUBJECT_ORIGINAL,
        "RUN" = tmp$RUN,
        "GROUP" = tmp$GROUP,
        "SUBJECT" = tmp$SUBJECT,
        "SUBJECT_NESTED" = tmp$SUBJECT_NESTED,
        "INTENSITY" = NA,
        "ABUNDANCE" = NA,
        "FRACTION" = tmp$FRACTION)
      rm(tmp)
    }   # end fillIncompleteRows option
  } # end for reference peptides

  ## 2. endogenous peptides
  featurestructure.l <- apply(structure.l, 1, function(x) sum(is.na(x)))
  ## get feature ID of reference which are completely missing across run
  featureID.l <- names(featurestructure.l[featurestructure.l == ncol(structure.l)])
  if (length(featureID.l) > 0) {
    ## print message
    print_string <- paste0(
      strwrap(prefix=" ", initial="",
              "CAUTION : some ENDOGENOUS features have missing intensities in
all the runs. The completely missing ENDOGENOUS features are "),
paste(featureID.l, collapse=", "),
strwrap(prefix=" ", initial="",
        ". Please check whether features in the list are correctly generated
 from spectral processing tool."))
    logging::logwarn(print_string)

    ## add missing rows if option is TRUE
    if (fillIncompleteRows) {
      ## get unique Run information
      nameID <- unique(work.l[, c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT",
                                  "SUBJECT_NESTED", "RUN", "FRACTION")])
      ## get PROTEIN and FEATURE information
      ## here use whole work dataset
      tempTogetfeature <- work[which(work$FEATURE %in% featureID.l), ]
      tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
      ## then generate data.frame for missingness,
      ##for (j in 1:nrow(nameID)) {
      ##    ## merge feature info and run info as "work" format
      ##    tempmissingwork <- data.frame(tempfeatureID, LABEL="L",GROUP_ORIGINAL=nameID$GROUP_ORIGINAL[j], SUBJECT_ORIGINAL=nameID$SUBJECT_ORIGINAL[j], RUN=nameID$RUN[j], GROUP=nameID$GROUP[j], SUBJECT=nameID$SUBJECT[j], SUBJECT_NESTED=nameID$SUBJECT_NESTED[j], INTENSITY=NA, ABUNDANCE=NA, METHOD=nameID$METHOD[j])
      ##    ## merge with tempary space, missingwork
      ##    missingcomplete.l <- rbind(missingcomplete.l, tempmissingwork)
      ##}

      ## MC : 2016.04.21 : use merge for simplicity
      tmp <- merge(nameID, tempfeatureID, by=NULL)
      missingcomplete.l <- data.frame(
        "PROTEIN" = tmp$PROTEIN,
        "PEPTIDE" = tmp$PEPTIDE,
        "TRANSITION" = tmp$TRANSITION,
        "FEATURE" = tmp$FEATURE,
        "LABEL" = "L",
        "GROUP_ORIGINAL" = tmp$GROUP_ORIGINAL,
        "SUBJECT_ORIGINAL" = tmp$SUBJECT_ORIGINAL,
        "RUN" = tmp$RUN,
        "GROUP" = tmp$GROUP,
        "SUBJECT" = tmp$SUBJECT,
        "SUBJECT_NESTED" = tmp$SUBJECT_NESTED,
        "INTENSITY" = NA,
        "ABUNDANCE" = NA,
        "FRACTION" = tmp$FRACTION)
      rm(tmp)
    } # end fillIncompleteRows option
  } # end endogenous peptides

  ## second, check other some missingness
  ## for missign row, need to assign before looping. need to assign at the beginning because it need either cases, with missingness or not
  missingwork.l <- NULL
  missingwork.h <- NULL

  ## structure value should be 1 for reference and endogenous separately, if not there are missingness. if more there are duplicates.
  ## if count of NA is not zero and not number of run (excluding complete missingness across runs)
  missing.l <- names(featurestructure.l[featurestructure.l != ncol(structure.l) & featurestructure.l != 0])
  missing.h <- names(featurestructure.h[featurestructure.h != ncol(structure.h) & featurestructure.h != 0])
  flagmissing.l <- length(missing.l) > 0
  flagmissing.h <- length(missing.h) > 0

  ## structure value is greater than 1, there are duplicates
  flagduplicate.l <- sum(structure.l[!is.na(structure.l)] > 1) > 0
  flagduplicate.h <- sum(structure.h[!is.na(structure.h)] > 1) > 0
  ## if there is missing rows for endogenous
  if (flagmissing.l | flagmissing.h) {
    print_string <- strwrap(prefix=" ", initial="",
                            "CAUTION: the input dataset has incomplete rows. If missing
 peaks occur they should be included in the dataset as separate rows, and the missing intensity
 values should be indicated with 'NA'. The incomplete rows are listed below.")
    logging:logwarn(print_string)
    ## endogenous intensities
    if (flagmissing.l) {
      if (length(missing.l) > 1){
        runstructure <- apply(structure.l[which(rownames(structure.l) %in% missing.l), ],
                              2, function(x) sum(is.na(x))) > 0
      } else if (length(missing.l) == 1) {
        runstructure <- is.na(structure.l[which(rownames(structure.l) %in% missing.l), ]) > 0
      }

      ## get the name of Run
      runID <- names(runstructure[runstructure==TRUE])
      ## then for each run, which features are missing,
      for (j in 1:length(runID)) {
        ## get subject, group information for this run
        nameID <- unique(work.l[work.l$RUN==runID[j],
                                c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT",
                                  "SUBJECT_NESTED", "RUN", "FRACTION")])
        ## MC : 2016/04/21. if there is one row, can"t catch up data.frame
        ## get feature ID
        if (length(missing.l) > 1) {
          featureID <- structure.l[which(rownames(structure.l) %in% missing.l),
                                   colnames(structure.l) == runID[j]]

          ## get feature ID which has no measuremnt.
          finalfeatureID <- names(featureID[is.na(featureID)])
          ## print features ID
          print_string <- paste0("*** Subject : ",
                                 as.character(nameID[,"SUBJECT_ORIGINAL"]) ,
                                 ", Condition : ",
                                 as.character(nameID[,"GROUP_ORIGINAL"]),
                                 " has incomplete rows for some ENDOGENOUS features (",
                                 paste(finalfeatureID, collapse=", "),")")
          logging::loginfo(print_string)
        } else if (length(missing.l) == 1) {
          finalfeatureID <- missing.l
          ## print features ID
          print_string <- paste0("*** Subject : ",
                                 as.character(nameID[,"SUBJECT_ORIGINAL"]),
                                 ", Condition : ",
                                 as.character(nameID[,"GROUP_ORIGINAL"]),
                                 " has incomplete rows for some ENDOGENOUS features (",
                                 finalfeatureID,")")
          logging::loginfo(print_string)
        }

        ## add missing rows if option is TRUE
        if (fillIncompleteRows) {
          tempTogetfeature <- work.l[which(work.l$FEATURE %in% finalfeatureID), ]
          ## get PROTEIN and FEATURE infomation
          tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
          ## merge feature info and run info as "work" format
          tempmissingwork <- data.frame(
            tempfeatureID,
            "LABEL" = "L",
            "GROUP_ORIGINAL" = nameID$GROUP_ORIGINAL,
            "SUBJECT_ORIGINAL" = nameID$SUBJECT_ORIGINAL,
            "RUN" = nameID$RUN,
            "GROUP" = nameID$GROUP,
            "SUBJECT" = nameID$SUBJECT,
            "SUBJECT_NESTED" = nameID$SUBJECT_NESTED,
            "INTENSITY" = NA,
            "ABUNDANCE" = NA,
            "FRACTION" = nameID$FRACTION)
          ## merge with tempary space, missingwork
          missingwork.l <- rbind(missingwork.l,tempmissingwork)
        } # end fillIncompleteRows options
      } # end loop for run ID
    } # end for endogenous

    ## reference intensities
    if (flagmissing.h) {
      ## first, which run has missing
      if (length(missing.h) > 1){
        runstructure <- apply(structure.h[which(rownames(structure.h) %in% missing.h), ], 2,
                              function(x) sum(is.na(x))) > 0
      } else if (length(missing.h) == 1) {
        runstructure <- is.na(structure.h[which(rownames(structure.h) %in% missing.h), ]) > 0
      }

      ## get the name of Run
      runID <- names(runstructure[runstructure == TRUE])
      ## then for each run, which features are missing,
      for (j in 1:length(runID)) {
        ## get subject, group information for this run
        nameID <- unique(work.h[work.h$RUN == runID[j],
                                c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT",
                                  "SUBJECT_NESTED", "RUN", "FRACTION")])
        ## MC : 2016/04/21. if there is one row, can"t catch up data.frame
        ## get feature ID
        if (length(missing.h) > 1) {
          featureID <- structure.h[which(rownames(structure.h) %in% missing.h), colnames(structure.h) == runID[j] ]
          ## get feature ID which has no measuremnt.
          finalfeatureID <- names(featureID[is.na(featureID)])
          ## print features ID
          print_string <- paste0("*** Subject : ",
                                 as.character(nameID[,"SUBJECT_ORIGINAL"]),
                                 ", Condition : ",
                                 as.character(nameID[,"GROUP_ORIGINAL"]),
                                 " has incomplete rows for some REFERENCE features (",
                                 paste(finalfeatureID, collapse=", "))
          logging::loginfo(print_string)
        } else if (length(missing.h) == 1) {
          finalfeatureID <- missing.h
          ## print features ID
          print_string <- paste0("*** Subject : ",
                                 as.character(nameID[, "SUBJECT_ORIGINAL"]),
                                 ", Condition : ",
                                 as.character(nameID[,"GROUP_ORIGINAL"]),
                                 " has incomplete rows for some REFERENCE features (",
                                 finalfeatureID,")")
          logging::loginfo(print_string)
        }

        ## add missing rows if option is TRUE
        if (fillIncompleteRows) {
          tempTogetfeature <- work.h[which(work.h$FEATURE %in% finalfeatureID), ]
          ## get PROTEIN and FEATURE infomation
          tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE")])
          ## merge feature info and run info as "work" format
          tempmissingwork <- data.frame(
            tempfeatureID,
            "LABEL" = "H",
            "GROUP_ORIGINAL" = nameID$GROUP_ORIGINAL,
            "SUBJECT_ORIGINAL" = nameID$SUBJECT_ORIGINAL,
            "RUN" = nameID$RUN,
            "GROUP" = nameID$GROUP,
            "SUBJECT" = nameID$SUBJECT,
            "SUBJECT_NESTED" = nameID$SUBJECT_NESTED,
            "INTENSITY" = NA,
            "ABUNDANCE" = NA,
            "FRACTION" = nameID$FRACTION)
          ## merge with tempary space, missingwork
          missingwork.h <- rbind(missingwork.h, tempmissingwork)
        } # end fillIncompleteRows options
      } # end loop for run ID
    } # end for endogenous
  } # end for flag missing

  ## merge missing rows if fillIncompleteRows=TRUE or message.
  if (fillIncompleteRows) {
    ## merge with work
    ## in future, use rbindlist?? rbindlist(list(work, missingwork))
    work <- rbind(work, missingcomplete.l, missingcomplete.h, missingwork.l, missingwork.h)
    print_string <- "DONE : Incomplete rows for missing peaks are added with intensity values=NA."
    logging::loginfo(print_string)
  } else if (!is.null(missingcomplete.l) | !is.null(missingcomplete.h) |
             !is.null(missingwork.l) | !is.null(missingwork.l) ) {
    print_string <- strwrap(prefix=" ", initial="",
                            "Please check whether features in the list are generated from
 spectral processing tool or not. Or the option, fillIncompleteRows=TRUE, will add incomplete
 rows for missing peaks with intensity=NA.")
    logging::loginfo(print_string)
  }

  ## if there are duplicates measurements
  if (flagduplicate.h) {
    ## first, which run has duplicates
    runstructure <- apply(structure.h, 2, function(x) sum(x[!is.na(x)] > 1) > 0)
    runID <- names(runstructure[runstructure == TRUE])
    ## then for each run, which features have duplicates,
    for (j in 1:length(runID)) {
      nameID <- unique(work[work$RUN==runID[j], c("SUBJECT_ORIGINAL",
                                                  "GROUP_ORIGINAL",
                                                  "GROUP",
                                                  "SUBJECT",
                                                  "SUBJECT_NESTED",
                                                  "RUN",
                                                  "FRACTION")])

      featureID <- structure.h[, colnames(structure.h) == runID[j]]
      finalfeatureID <- featureID[!is.na(featureID) & featureID > 1]
      print_string <- paste0("*** Subject : ",
                             as.character(nameID[,"SUBJECT_ORIGINAL"]),
                             ", Condition : ",
                             as.character(nameID[,"GROUP_ORIGINAL"]),
                             " has multiple rows (duplicate rows) for some REFERENCE features (",
                             head(paste(names(finalfeatureID), collapse=", ")), ")")
      logging::logwarn(print_string)
    }
    ## save in process file.
    print_string <- "Please remove duplicate rows in the list above. "
    logging::logwarn(print_string)
  } # end flag duplicate for reference

  if (flagduplicate.l) {
    ## first, which run has duplicates
    runstructure <- apply(structure.l, 2, function(x) sum(x[!is.na(x)] > 1) > 0)
    runID <- names(runstructure[runstructure == TRUE])
    ## then for each run, which features have duplicates,
    for (j in 1:length(runID)) {
      nameID <- unique(work[work$RUN==runID[j], c("SUBJECT_ORIGINAL",
                                                  "GROUP_ORIGINAL",
                                                  "GROUP",
                                                  "SUBJECT",
                                                  "SUBJECT_NESTED",
                                                  "RUN",
                                                  "FRACTION")])

      featureID <- structure.l[, colnames(structure.l) == runID[j]]
      finalfeatureID <- featureID[!is.na(featureID) & featureID > 1]
      print_string <- paste0("*** Subject : ",
                             head(as.character(nameID[,"SUBJECT_ORIGINAL"])),
                             ", Condition : ",
                             head(as.character(nameID[,"GROUP_ORIGINAL"])),
                             " has multiple rows (duplicate rows) for some ENDOGENOUS features (",
                             head(paste(names(finalfeatureID), collapse=", ")),")")
      message(print_string)
      logging::logwarn(print_string)
    }
    print_string <- "ERROR : Please remove duplicate rows in the list above. "
    logging::logwarn(print_string)
  } # end flag duplicate for endogenous

  ## no missing and no duplicates
  if (!flagmissing.h & !flagmissing.l & !flagduplicate.h & !flagduplicate.l) {
    print_string <- "Balanced data format with NA for missing feature intensities - okay"
    logging::loginfo(print_string)
  }

  return(work)
}

#' Examine the work data structure single fraction unlabeled data.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @param fillIncompleteRows  Fill those rows!!
#' @return Send the work data back to dataProcess()
extract_singlefrac_labelfree_data <- function(work, fillIncompleteRows=TRUE) {
  ## get feature by Run count of data
  structure = tapply(work$ABUNDANCE, list(work$FEATURE, work$RUN),
                     function(x) {
                       length (x)
                     })
  ## structure value should be 1 for label-free, if not there are missingness. if more there are duplicates.
  flagmissing = sum(is.na(structure)) > 0
  flagduplicate = sum(structure[!is.na(structure)] > 1) > 0
  ## if there is missing rows
  if (flagmissing) {
    print_string <- "CAUTION: the input dataset has incomplete rows. If missing peaks
 occur they should be included in the dataset as separate rows, and the missing intensity
 values should be indicated with 'NA'. The incomplete rows are listed below."
    logging::logwarn(print_string)
    ## first, which run has missing
    runstructure <- apply(structure, 2, function(x) sum(is.na (x))) > 0
    ## get the name of Run
    runID <- names(runstructure[runstructure == TRUE])
    ## for missign row, need to assign before looping
    missingwork <- NULL
    ## then for each run, which features are missing,
    for (j in 1:length(runID)) {
      ## get subject, group information for this run
      nameID <- unique(work[work$RUN == runID[j],
                            c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT",
                              "SUBJECT_NESTED", "RUN", "FRACTION")])
      ## get feature ID
      featureID <- structure[, colnames(structure) == runID[j]]
      ## get feature ID which has no measuremnt.
      finalfeatureID <- featureID[is.na(featureID)]
      ## print features ID
      print_string <- strwrap(prefix=" ", initial="",
                              paste0("*** Subject : ",
                                     head(as.character(nameID[["SUBJECT_ORIGINAL"]])),
                                     ", Condition : ",
                                     head(as.character(nameID[["GROUP_ORIGINAL"]])),
                                     " has incomplete rows for some features (",
                                     toString(head(names(finalfeatureID)))))
      logging::loginfo(print_string)
      ## add missing rows if option is TRUE
      if (fillIncompleteRows) {
        tempTogetfeature <- work[which(work$FEATURE %in% names(finalfeatureID)), ]
        ## get PROTEIN and FEATURE infomation
        tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE",
                                                     "TRANSITION", "FEATURE")])
        ## merge feature info and run info as "work" format
        ## One of my personal style choices is a little strange, I like to quote the
        ## data.frame() and list() ... material which is used to create the columns/elements.
        ## But that can wait for another day, there are far far bigger problems with this code
        ## and one can easily say that my style choice in this is inconsistent.
        tempmissingwork <- data.frame(
          tempfeatureID,
          "LABEL" = "L",
          "GROUP_ORIGINAL" = nameID$GROUP_ORIGINAL,
          "SUBJECT_ORIGINAL" = nameID$SUBJECT_ORIGINAL,
          "RUN" = nameID$RUN,
          "GROUP" = nameID$GROUP,
          "SUBJECT" = nameID$SUBJECT,
          "SUBJECT_NESTED" = nameID$SUBJECT_NESTED,
          "INTENSITY" = NA,
          "ABUNDANCE" = NA,
          "FRACTION" = nameID$FRACTION)
        ## merge with tempary space, missingwork
        missingwork <- rbind(missingwork, tempmissingwork)
      } # end fillIncompleteRows options
    } # end loop for run ID

    ## [THT: this part can probably be merged into the above.
    ## Also, it might be better to check fillIncompleteRows earlier
    ## and terminate the process when it"s FALSE]
    if (fillIncompleteRows) {
      ## merge with work
      ## in future, use rbindlist?? rbindlist(list(work, missingwork))
      work <- rbind(work, missingwork)
      ## print message
      print_string <- "Incomplete rows for missing peaks are added with intensity values=NA."
      message(print_string)
      logging::loginfo(print_string)
    } else {
      print_string <- strwrap(prefix=" ", initial="",
                              "Please check whether features in the list are generated from
 spectral processing tool. Or the option, fillIncompleteRows=TRUE, will add incomplete rows
 for missing peaks with intensity=NA.")
      logging::logwarn(print_string)
    }
  } # end for flag missing

  ## if there are duplicates measurements
  if (flagduplicate) {
    ## first, which run has duplicates
    runstructure <- apply(structure, 2, function(x) sum(x[!is.na(x)] > 1) > 0)
    runID <- names(runstructure[runstructure == TRUE])
    ## then for each run, which features have duplicates,
    for (j in 1:length(runID)) {
      nameID <- unique(work[work$RUN == runID[j],
                            c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT",
                              "SUBJECT_NESTED", "RUN", "FRACTION")])
      featureID <- structure[, colnames(structure) == runID[j]]
      finalfeatureID <- featureID[!is.na(featureID) & featureID > 1]
      print_string <- paste0("*** Subject : ",
                             head(as.character(nameID[["SUBJECT_ORIGINAL"]])),
                             ", Condition : ",
                             head(as.character(nameID[["GROUP_ORIGINAL"]])),
                             " has multiple rows (duplicate rows) for some features (",
                             toString(head(names(finalfeatureID))), ")")
      message(print_string)
      logging::logwarn(print_string)
    }
    print_string <- "Please remove duplicate rows in the list above."
    logging::logwarn(print_string)
  } # end flag duplicate

  ## no missing and no duplicates
  if (!flagmissing & !flagduplicate) {
    print_string <- "Balanced data format with NA for missing feature intensities."
    logging::loginfo(print_string)
  }
  ## end label-free
  return(work)
}

#' Examine the work data structure for multi fraction data.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @return Send the work data back to dataProcess()
extract_multifrac_data <- function(work) {
  ## This large block handles multiple fractionations.
  ## It appears to duplicate the labelfree/labeled code which I moved into
  ## the functions extract_singlfrac_labelfree_data()
  ## and extract_singlefrac_labeled_data()
  allflagmissing <- NULL
  allflagduplicate <- NULL
  ## check each method
  for (k in 1:length(unique(work$FRACTION))) {
    worktemp <- work[work$FRACTION == k, ]
    worktemp$RUN <- factor(worktemp$RUN)
    worktemp$FEATURE <- factor(worktemp$FEATURE)
    structure <- tapply(worktemp$ABUNDANCE, list(worktemp$FEATURE, worktemp$RUN), function(x) length(x))
    ## structure value should be 2 for labeled, 1 for label-free, if not there are missingness
    if (nlevels(worktemp$LABEL) == 2) { ## label-based
      flag <- sum(is.na(structure)) > 0 | sum(structure[!is.na(structure)] < 2) > 0
    } else { ## label-free
      flag <- sum(is.na(structure)) > 0
    }
    allflagmissing <- c(allflagmissing,flag)

    ## for duplicate
    if (nlevels(worktemp$LABEL) == 2) { # label-based
      worktemp.h <- worktemp[worktemp$LABEL == "H", ]
      worktemp.l <- worktemp[worktemp$LABEL == "L", ]
      structure.h <- tapply(worktemp.h$ABUNDANCE,
                            list(worktemp.h$FEATURE, worktemp.h$RUN),
                            function(x) length(x))
      structure.l <- tapply(worktemp.l$ABUNDANCE,
                            list(worktemp.l$FEATURE, worktemp.l$RUN),
                            function(x) length(x))
      flagduplicate <- sum(structure.h[!is.na(structure.h)] > 1) > 0 |
        sum(structure.l[!is.na(structure.l)] > 1) > 0
    } else { # label-free
      flagduplicate <- sum(structure[!is.na(structure)] > 1) > 0
    }
    allflagduplicate <- c(allflagduplicate, flag)
  } # end to check any flag among methods

  if (sum(allflagmissing) != 0) {
    print_string <- strwrap(prefix=" ", initial="",
                            "CAUTION: the input dataset has incomplete rows. Missing feature
 intensities should be present in the dataset, and their intensities should be indicated with
 'NA'. The incomplete rows are listed below.")
    logging::logwarn(print_string)

    ## for missign row, need to assign before looping
    missingwork <- NULL
    missingcomplete.h <- NULL
    missingcomplete.l <- NULL
    missingwork.h <- NULL
    missingwork.l <- NULL
    for (k in 1:length(unique(work$FRACTION))) {
      ## see which method has missing rows
      if (allflagmissing[k]) {
        worktemp <- work[work$FRACTION==k, ]
        worktemp$RUN <- factor(worktemp$RUN)
        worktemp$FEATURE <- factor(worktemp$FEATURE)
        if (nlevels(worktemp$LABEL) == 1) { ## label-free
          structure = tapply(worktemp$ABUNDANCE,
                             list(worktemp$FEATURE, worktemp$RUN),
                             function(x) length(x))
          ## first, which run has missing
          runstructure <- apply(structure, 2, function(x) sum(is.na(x))) > 0
          ## get the name of Run
          runID <- names(runstructure[runstructure == TRUE])

          ## then for each run, which features are missing,
          for (j in 1:length(runID)) {
            nameID <- unique(worktemp[worktemp$RUN == runID[j],
                                      c("SUBJECT_ORIGINAL", "GROUP_ORIGINAL", "GROUP", "SUBJECT",
                                        "SUBJECT_NESTED", "RUN", "FRACTION")])
            ## get feature ID
            featureID <- structure[, colnames(structure) == runID[j]]
            ## get feature ID which has no measuremnt.
            finalfeatureID <- featureID[is.na(featureID)]
            ## print features ID
            print_string <- paste0("*** Subject : ",
                                   as.character(nameID[,"SUBJECT_ORIGINAL"]),
                                   ", Condition : ",
                                   as.character(nameID[,"GROUP_ORIGINAL"]),
                                   " has incomplete rows for some features (",
                                   paste(names(finalfeatureID), collapse=", "),")")
            logging::loginfo(print_string)

            ## add missing rows if option is TRUE
            if (fillIncompleteRows) {
              tempTogetfeature <- work[which(work$FEATURE %in% names(finalfeatureID)), ]
              ## get PROTEIN and FEATURE infomation
              tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE",
                                                           "TRANSITION", "FEATURE")])
              ## merge feature info and run info as "work" format
              tempmissingwork <- data.frame(
                tempfeatureID,
                "LABEL" = "L",
                "GROUP_ORIGINAL" = nameID$GROUP_ORIGINAL,
                "SUBJECT_ORIGINAL" = nameID$SUBJECT_ORIGINAL,
                "RUN" = nameID$RUN,
                "GROUP" = nameID$GROUP,
                "SUBJECT" = nameID$SUBJECT,
                "SUBJECT_NESTED" = nameID$SUBJECT_NESTED,
                "INTENSITY" = NA,
                "ABUNDANCE" = NA,
                "FRACTION" = nameID$FRACTION)
              ## merge with tempary space, missingwork
              missingwork <- rbind(missingwork, tempmissingwork)
            } # end fillIncompleteRows options
          } # end loop for run
        } else { # end label-free
          ## label-based
          ## count the reference and endobenous separately
          work.l <- worktemp[worktemp$LABEL == "L", ]
          work.h <- worktemp[worktemp$LABEL == "H", ]
          ## get feature by Run count of data
          structure.l <- tapply(work.l$ABUNDANCE,
                                list(work.l$FEATURE, work.l$RUN),
                                function(x) length(x))
          structure.h <- tapply(work.h$ABUNDANCE,
                                list(work.h$FEATURE, work.h$RUN),
                                function(x) length(x))
          ## 1. reference peptides
          featurestructure.h <- apply(structure.h, 1, function (x) sum(is.na(x)))
          ## get feature ID of reference which are completely missing across run
          featureID.h <- names(featurestructure.h[featurestructure.h == ncol(structure.h)])

          if (length(featureID.h) > 0) {
            ## print message
            print_string <- paste0(strwrap(prefix=" ", initial="",
                                           "CAUTION : some REFERENCE features have missing intensities
 in all the runs. The completely missing REFERENCE features are "),
 paste(featureID.h, collapse=", "),
 strwrap(prefix=" ", initial="",
         ". Please check whether features in the list are correctly
 generated from spectral processing tool."))
            logging::logwarn(print_string)
            ## add missing rows if option is TRUE
            if (fillIncompleteRows) {
              if (nrow(work.h) == 0) {
                work.h <- work[work$LABEL=="H", ]
                ## get unique Run information
                nameID <- unique(work.h[, c("SUBJECT_ORIGINAL",
                                            "GROUP_ORIGINAL",
                                            "GROUP",
                                            "SUBJECT",
                                            "SUBJECT_NESTED",
                                            "RUN",
                                            "FRACTION")])
                nameID$FRACTION <- k
              } else {
                ## get unique Run information
                nameID <- unique(work.h[, c("SUBJECT_ORIGINAL",
                                            "GROUP_ORIGINAL",
                                            "GROUP",
                                            "SUBJECT",
                                            "SUBJECT_NESTED",
                                            "RUN",
                                            "FRACTION")])
              }

              ## get PROTEIN and FEATURE information
              ## here use whole worktemp dataset
              tempTogetfeature <- worktemp[which(worktemp$FEATURE %in% featureID.h), ]
              tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE",
                                                           "TRANSITION", "FEATURE")])
              ## then generate data.frame for missingness,
              for (j in 1:nrow(nameID)) {
                ## merge feature info and run info as "work" format
                tempmissingwork <- data.frame(
                  tempfeatureID,
                  "LABEL" = "H",
                  "GROUP_ORIGINAL" = nameID$GROUP_ORIGINAL[j],
                  "SUBJECT_ORIGINAL" = nameID$SUBJECT_ORIGINAL[j],
                  "RUN" = nameID$RUN[j],
                  "GROUP" = nameID$GROUP[j],
                  "SUBJECT" = nameID$SUBJECT[j],
                  "SUBJECT_NESTED" = nameID$SUBJECT_NESTED[j],
                  "INTENSITY" = NA,
                  "ABUNDANCE" = NA,
                  "FRACTION" = nameID$FRACTION[j])
                ## merge with tempary space, missingwork
                missingcomplete.h <- rbind(missingcomplete.h, tempmissingwork)
              }
            } # end fillIncompleteRows option
          } # end for reference peptides

          ## 2. endogenous peptides
          featurestructure.l <- apply(structure.l, 1, function (x) sum(is.na(x)))
          ## get feature ID of reference which are completely missing across run
          featureID.l <- names(featurestructure.l[featurestructure.l == ncol(structure.l)])
          if (length(featureID.l) > 0) {
            ## print message
            print_string <- paste0(strwrap(prefix=" ", initial="",
                                           "CAUTION : some ENDOGENOUS features have missing intensities
 in all the runs. The completely missing ENDOGENOUS features are "),
 toString(featureID.l),
 strwrap(prefix = " ", initial = "",
         ". Please check whether features in the list are correctly
 generated from spectral processing tool."))
            logging::logwarn(print_string)
            ## add missing rows if option is TRUE
            if (fillIncompleteRows) {
              ## get unique Run information
              nameID <- unique(work.l[, c("SUBJECT_ORIGINAL",
                                          "GROUP_ORIGINAL",
                                          "GROUP",
                                          "SUBJECT",
                                          "SUBJECT_NESTED",
                                          "RUN",
                                          "FRACTION")])

              ## get PROTEIN and FEATURE information
              ## here use whole worktemp dataset
              tempTogetfeature <- worktemp[which(worktemp$FEATURE %in% featureID.l), ]
              tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE",
                                                           "TRANSITION", "FEATURE")])
              ## then generate data.frame for missingness,
              for (j in 1:nrow(nameID)) {
                ## merge feature info and run info as "work" format
                tempmissingwork <- data.frame(
                  tempfeatureID,
                  "LABEL" = "L",
                  "GROUP_ORIGINAL" = nameID$GROUP_ORIGINAL[j],
                  "SUBJECT_ORIGINAL" = nameID$SUBJECT_ORIGINAL[j],
                  "RUN" = nameID$RUN[j],
                  "GROUP" = nameID$GROUP[j],
                  "SUBJECT" = nameID$SUBJECT[j],
                  "SUBJECT_NESTED" = nameID$SUBJECT_NESTED[j],
                  "INTENSITY" = NA,
                  "ABUNDANCE" = NA,
                  "FRACTION" = nameID$FRACTION[j])
                ## merge with tempary space, missingwork
                missingcomplete.l <- rbind(missingcomplete.l, tempmissingwork)
              }
            } # end fillIncompleteRows option
          } # end endogenous peptides

          ## second, check other some missingness
          ## structure value should be 1 for reference and endogenous separately, if not
          ## there are missingness. if more there are duplicates.
          ## if count of NA is not zero and not number of run (excluding complete
          ## missingness across runs)
          missing.l <- names(featurestructure.l[featurestructure.l!=ncol(structure.l) &
                                                featurestructure.l != 0])
          missing.h <- names(featurestructure.h[featurestructure.h!=ncol(structure.h) &
                                                featurestructure.h != 0])
          flagmissing.l <- length(missing.l) > 0
          flagmissing.h <- length(missing.h) > 0
          ## structure value is greater than 1, there are duplicates
          flagduplicate.l <- sum(structure.l[!is.na(structure.l)] > 1) > 0
          flagduplicate.h <- sum(structure.h[!is.na(structure.h)] > 1) > 0

          ## if there is missing rows for endogenous
          if (flagmissing.l | flagmissing.h) {
            print_string <- strwrap(prefix = " ", initial = "",
                                    "CAUTION: the input dataset has incomplete rows. If
 missing peaks occur they should be included in the dataset as separate rows, and the missing
 intensity values should be indicated with 'NA'. The incomplete rows are listed below.")
            logging::logwarn(print_string)
            ## endogenous intensities
            if (flagmissing.l) {
              ## first, which run has missing
              runstructure <- apply(structure.l[-which(rownames(structure.l) %in% featureID.l),],
                                    2, function(x) sum(is.na(x))) > 0
              ## get the name of Run
              runID <- names(runstructure[runstructure==TRUE])
              ## then for each run, which features are missing,
              for (j in 1:length(runID)) {
                ## get subject, group information for this run
                nameID <- unique(work.l[work.l$RUN==runID[j], c("SUBJECT_ORIGINAL",
                                                                "GROUP_ORIGINAL",
                                                                "GROUP",
                                                                "SUBJECT",
                                                                "SUBJECT_NESTED",
                                                                "RUN",
                                                                "FRACTION")])
                ## get feature ID
                featureID <- structure.l[-which(rownames(structure.l) %in% featureID.l),
                                         colnames(structure.l) == runID[j]]
                ## get feature ID which has no measuremnt.
                finalfeatureID <- featureID[is.na(featureID)]
                ## print features ID
                print_string <- paste0("*** Subject : ",
                                       as.character(nameID[,"SUBJECT_ORIGINAL"]) ,
                                       ", Condition : ", as.character(nameID[, "GROUP_ORIGINAL"]),
                                       " has incomplete rows for some ENDOGENOUS features (",
                                       paste(names(finalfeatureID), collapse=", "),")")
                logging::logwarn(print_string)
                ## add missing rows if option is TRUE
                if (fillIncompleteRows) {
                  tempTogetfeature <- work.l[which(work.l$FEATURE %in% names(finalfeatureID)), ]
                  ## get PROTEIN and FEATURE infomation
                  tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE",
                                                               "TRANSITION", "FEATURE")])
                  ## merge feature info and run info as "work" format
                  tempmissingwork <- data.frame(
                    tempfeatureID,
                    "LABEL" = "L",
                    "GROUP_ORIGINAL" = nameID$GROUP_ORIGINAL,
                    "SUBJECT_ORIGINAL" = nameID$SUBJECT_ORIGINAL,
                    "RUN" = nameID$RUN,
                    "GROUP" = nameID$GROUP,
                    "SUBJECT" = nameID$SUBJECT,
                    "SUBJECT_NESTED" = nameID$SUBJECT_NESTED,
                    "INTENSITY" = NA,
                    "ABUNDANCE" = NA,
                    "FRACTION" = nameID$FRACTION)
                  ## merge with tempary space, missingwork
                  missingwork.l <- rbind(missingwork.l, tempmissingwork)
                } # end fillIncompleteRows options
              } # end loop for run ID
            } # end for endogenous

            ## reference intensities
            if (flagmissing.h) {
              ## first, which run has missing
              runstructure <- apply(structure.h[-which(rownames(structure.h) %in% featureID.h),],
                                    2, function(x) sum(is.na(x))) > 0
              ## get the name of Run
              runID <- names(runstructure[runstructure==TRUE])
              ## then for each run, which features are missing,
              for (j in 1:length(runID)) {
                ## get subject, group information for this run
                nameID <- unique(work.h[work.h$RUN==runID[j], c("SUBJECT_ORIGINAL",
                                                                "GROUP_ORIGINAL",
                                                                "GROUP",
                                                                "SUBJECT",
                                                                "SUBJECT_NESTED",
                                                                "RUN",
                                                                "FRACTION")])

                ## get feature ID
                featureID <- structure.h[-which(rownames(structure.h) %in% featureID.h),
                                         colnames(structure.h)==runID[j]]
                ## get feature ID which has no measuremnt.
                finalfeatureID <- featureID[is.na(featureID)]
                ## print features ID
                print_string <- paste0("*** Subject : ",
                                       as.character(nameID[,"SUBJECT_ORIGINAL"]) ,
                                       ", Condition : ",
                                       as.character(nameID[,"GROUP_ORIGINAL"]),
                                       " has incomplete rows for some REFERENCE features (",
                                       paste(names(finalfeatureID), collapse=", "),")")
                logging::logwarn(print_string)
                ## add missing rows if option is TRUE
                if (fillIncompleteRows) {
                  tempTogetfeature <- work.h[which(work.h$FEATURE %in% names(finalfeatureID)), ]
                  ## get PROTEIN and FEATURE infomation
                  tempfeatureID <- unique(tempTogetfeature[, c("PROTEIN", "PEPTIDE",
                                                               "TRANSITION", "FEATURE")])
                  ## merge feature info and run info as "work" format
                  tempmissingwork <- data.frame(
                    tempfeatureID,
                    "LABEL" = "H",
                    "GROUP_ORIGINAL" = nameID$GROUP_ORIGINAL,
                    "SUBJECT_ORIGINAL" = nameID$SUBJECT_ORIGINAL,
                    "RUN" = nameID$RUN,
                    "GROUP" = nameID$GROUP,
                    "SUBJECT" = nameID$SUBJECT,
                    "SUBJECT_NESTED" = nameID$SUBJECT_NESTED,
                    "INTENSITY" = NA,
                    "ABUNDANCE" = NA,
                    "FRACTION" = nameID$FRACTION)
                  ## merge with tempary space, missingwork
                  missingwork.h <- rbind(missingwork.h, tempmissingwork)
                } # end fillIncompleteRows options
              } # end loop for run ID
            } # end for endogenous
          } # end any missingness
        } # end label-based
      } # if only any flag for method
    } # end loop for methods

    if (fillIncompleteRows) {
      ## merge with work
      ## in future, use rbindlist?? rbindlist(list(work, missingwork))
      if (nlevels(worktemp$LABEL) == 1) {
        work <- rbind(work, missingwork)
      } else {
        work <- rbind(work, missingcomplete.l, missingcomplete.h, missingwork.l, missingwork.h)
      }

      ## print message
      print_string <- "Incomplete rows for missing peaks are added with intensity values=NA."
      logging::loginfo(print_string)
    } else if (!is.null(missingcomplete.l) |
               !is.null(missingcomplete.h) |
               !is.null(missingwork.l) |
               !is.null(missingwork.l) |
               !is.null(missingwork)) {
      print_string <- strwrap(prefix=" ", initial="",
                              "Please check whether features in the list are generated from
 spectral processing tool. Or the option, fillIncompleteRows=TRUE, will add incomplete rows
 for missing peaks with intensity=NA.")
      logging::logwarn(print_string)
    }
  } else {
    print_string <- "Balanced data format with NA for missing feature intensities."
    logging::loginfo(print_string)
  }
  ## for duplicate, in future
  return(work)
}

#' Examine the work data structure and quantile normalize single fraction data.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @return Send the work data back to dataProcess()
normalize_singlefrac_by_quantiles <- function(work) {
  ## for label-free, just use endogenous
  nmethod <- unique(work$FRACTION)
  quantileall <- NULL
  ## ABUNDANCE=0 replace with 1, in order to distinguish later.
  work[!is.na(work$ABUNDANCE) & work$ABUNDANCE == 0, "ABUNDANCE"] <- 1
  for (j in 1:length(nmethod)) {
    namerun <- unique(work[work$FRACTION == nmethod[j], "RUN"])
    worktemp <- work[which(work$RUN %in% namerun & !is.na(work$INTENSITY)), ]
    worktemp$RUN <- factor(worktemp$RUN)
    worktemp$FEATURE <- factor(worktemp$FEATURE)
    quantiletemp <- as.matrix(xtabs(ABUNDANCE ~ FEATURE + RUN, data=worktemp))
    ## need to put NA for missing value in endogenous
    quantiletemp[quantiletemp == 0] <- NA
    ## using preprocessCore library
    quantiledone <- preprocessCore::normalize.quantiles(quantiletemp)
    rownames(quantiledone) <- rownames(quantiletemp)
    colnames(quantiledone) <- colnames(quantiletemp)
    ## get quantiled to long format for apply difference endogenous
    quantilelong <- reshape2::melt(quantiledone, id=rownames(quantiledone))
    colnames(quantilelong) <- c("FEATURE", "RUN", "ABUNDANCE_quantile")
    rm(quantiledone)
    ## quantileall <- rbindlist(list(quantileall,quantilelong))
    quantileall <- rbind(quantileall, quantilelong)
    rm(quantilelong)
  }
  work <- merge(work, quantileall, by=c("FEATURE", "RUN"))
  rm(quantileall)
  ## reorder
  work <- data.frame(
    "PROTEIN" = work$PROTEIN,
    "PEPTIDE" = work$PEPTIDE,
    "TRANSITION" = work$TRANSITION,
    "FEATURE" = work$FEATURE,
    "LABEL" = work$LABEL,
    "GROUP_ORIGINAL" = work$GROUP_ORIGINAL,
    "SUBJECT_ORIGINAL" = work$SUBJECT_ORIGINAL,
    "RUN" = work$RUN,
    "GROUP" = work$GROUP,
    "SUBJECT" = work$SUBJECT,
    "SUBJECT_NESTED" = work$SUBJECT_NESTED,
    "INTENSITY" = work$INTENSITY,
    "ABUNDANCE" = work$ABUNDANCE_quantile,
    "FRACTION" = work$FRACTION,
    "originalRUN" = work$originalRUN)

  work <- work[with(work,
                    order(LABEL, GROUP_ORIGINAL, SUBJECT_ORIGINAL, RUN, PROTEIN, PEPTIDE,
                          TRANSITION)), ]
  ## for skyline case, separate 1 and zero
  work[!is.na(work$INTENSITY) & work$INTENSITY == 1, "ABUNDANCE"] <- 0

  return(work)
}

#' Examine the work data structure and normalize by global standards.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @return Send the work data back to dataProcess()
normalize_by_globalstandards <- function(work) {
  work$RUN <- factor(work$RUN)
  combine <- data.frame(RUN=levels(work$RUN))
  allPeptide <- unique(work$PEPTIDE)
  allProtein <- unique(work$PROTEIN)
  for (i in 1:length(nameStandards)) {
    ## if Peptides
    ## namePeptide <- allPeptide[grep(nameStandards[i],allPeptide)]
    ## cannot grep for modified peptide sequence, [,],+ sign
    namePeptide <- tempPeptide[tempPeptide$PEPTIDESEQUENCE == nameStandards[i], "PEPTIDE"]
    if (length(namePeptide) != 0) {
      tempStandard <- work[work$PEPTIDE == namePeptide,]
    } else {
      ## if Proteins
      nameProtein <- allProtein[allProtein == nameStandards[i]]
      ## if we use "grep", can" find the proteins name with some symbol, such as "sp|P30153|2AAA_HUMAN"
      if (length(nameProtein) != 0) {
        tempStandard <- work[work$PROTEIN==nameProtein, ]
      } else {
        print_string <- paste0("global standard peptides or proteins, ",
                               nameStandards[i] ,
                               ", is not in dataset. Check whether 'nameStandards' is correct.")
        logging::logwarn(print_string)
      }
    }

    ##\ here, by RUN, but need to check !!!
    tempStandard <- tempStandard[tempStandard$GROUP != "0" ,]
    tempStandard$RUN <- factor(tempStandard$RUN)
    tempStandard <- tempStandard[!is.na(tempStandard$ABUNDANCE), ]
    meanStandard <- tapply(tempStandard$ABUNDANCE, tempStandard$RUN,
                           function(x) mean(x, na.rm=TRUE))
    meanStandard <- data.frame(RUN=names(meanStandard), meanStandard)
    combine <- merge(combine, meanStandard, by="RUN", all=TRUE)
    colnames(combine)[i + 1] <- paste0("meanStandard", i)
  }
  rownames(combine) <- combine$RUN
  combine <- subset(combine, select=-c(RUN))
  ## get mean among global standards
  allmean <- apply(combine, 1,
                   function(x) mean(x, na.rm=TRUE))
  ## allmean[is.na(allmean)] <- 0
  allmeantemp <- data.frame(RUN=names(allmean), allmean)
  allrun <- unique(work[, c("RUN", "FRACTION")])
  allmeantemp <- merge(allmeantemp, allrun, by="RUN")
  median.all <- tapply(allmeantemp$allmean, allmeantemp$FRACTION,
                       function(x) median(x,na.rm=TRUE))
  ## adjust
  nmethod <- unique(work$FRACTION)
  for (j in 1:length(nmethod)) {
    namerun <- unique(work[work$FRACTION==nmethod[j], "RUN"])
    for (i in 1:length(namerun)) {
      ## ABUNDANCE is normalized
      if (!is.na(allmean[names(allmean)==namerun[i]])) {
        work[work$RUN == namerun[i] & work$LABEL=="L","ABUNDANCE"] <-
          work[work$RUN == namerun[i] & work$LABEL == "L","ABUNDANCE"] -
          allmean[names(allmean) == namerun[i]]+median.all[j]
      }
    }
  } # end loop method
  print_string <- "Normalization : normalization with global standards protein - okay"
  logging::loginfo(print_string)
  return(work)
}

#' Examine the work data structure and normalize multifraction data
#' This code is confusering.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @return Send the work data back to dataProcess()
multifraction_normalize_confusing <- function(work) {
  ## check any features measured across all runs.
  tmp <- work[!is.na(work$ABUNDANCE), ]
  check.multiple.run <- xtabs(~ FEATURE + FRACTION, tmp)
  check.multiple.run.TF <- check.multiple.run != 0
  check.multiple.run.feature <- apply(check.multiple.run.TF, 1, sum)
  ## each feature should be measured only in one method
  overlap.feature <- names(check.multiple.run.feature[check.multiple.run.feature > 1 ])
  if (length(overlap.feature) > 0) {
    print_string <- paste0("** Please check the listed featurues (",
                           paste(overlap.feature, collapse=", "),
                           ") \n Those features are measured across all fractionations.")
    logging::logwarn(print_string)
    stop_string <- "Keep only one intensity of listed features among fractinations from one sample."
    logging::logerror(stop_string)
    stop()
  }
  ## merge ##
  ## get which Run id should be merged
  ## decide which two runs should be merged
  if (any(is.element(colnames(work), "TECHREPLICATE"))) {
    runid.multiple <- unique(work[, c("GROUP_ORIGINAL",
                                      "SUBJECT_ORIGINAL",
                                      "RUN",
                                      "originalRUN",
                                      "FRACTION",
                                      "TECHREPLICATE")])
    ## if there are technical replicates from the same group and subject, can"t match.
    run.match <- try(reshape2::dcast(GROUP_ORIGINAL + SUBJECT_ORIGINAL + TECHREPLICATE ~ FRACTION,
                                     data=runid.multiple, value.var = "originalRUN"), silent=TRUE)
    ## While it is true, I like try(), why would you immediately
    ## follow up an unsuccessful try() with a stop!?
    if (class(run.match) == "try-error") {
      stop_string <- "Cannot figure out which runs come from the same sample."
      logging::logerror(stop_string)
      stop()
    } else {
      work$newRun <- NA
      run.match$GROUP_ORIGINAL <- as.character(run.match$GROUP_ORIGINAL)
      run.match$SUBJECT_ORIGINAL <- as.character(run.match$SUBJECT_ORIGINAL)
      for (k in 1:nrow(run.match)) {
        work[which(work$originalRUN %in% run.match[k, 4:ncol(run.match)]), "newRun"] <- paste(
          paste(run.match[k, 1:4], collapse="_"), "merged", sep="_")
      }
      ## remove extra run NAs
      tmp <- work[is.na(work$ABUNDANCE), ]
      na.count <- reshape2::dcast(FEATURE ~ FRACTION, data=tmp,
                                  fun.aggregate=length, value.var="ABUNDANCE")
      na.count.long <- melt(na.count, id.vars=c("FEATURE"))
      na.count.long <- na.count.long[na.count.long$value == length(unique(work$newRun)), ]
      na.count.long$tmp <- paste(na.count.long$FEATURE, na.count.long$variable, sep="_")
      work$tmp <- paste(work$FEATURE, work$FRACTION, sep="_")
      work <- work[-which(work$tmp %in% na.count.long$tmp), ]
      ##
      work$originalRUN <- work$newRun
      ## update RUN based on new originalRUN
      work$RUN <- work$originalRUN
      work$RUN <- factor(work$RUN, levels=unique(work$RUN),
                         labels=seq(1, length(unique(work$RUN))))
      work <- work[, -which(colnames(work) %in% c("tmp","newRun"))]
    }
  } else { ## Fraction, but no tech replicate
    runid.multiple <- unique(work[, c("GROUP_ORIGINAL",
                                      "SUBJECT_ORIGINAL",
                                      "RUN",
                                      "originalRUN",
                                      "FRACTION")])
    ## if there are technical replicates from the same group and subject, can"t match.
    run.match <- try(reshape2::dcast(GROUP_ORIGINAL + SUBJECT_ORIGINAL ~ FRACTION,
                                     data=runid.multiple, value.var="originalRUN"), silent=TRUE)
    if (class(run.match) == "try-error") {
      stop_string <- "Cannot figure out which runs come from the same sample."
      logging::logerror(stop_string)
      stop()
    } else {
      work$newRun <- NA
      run.match$GROUP_ORIGINAL <- as.character(run.match$GROUP_ORIGINAL)
      run.match$SUBJECT_ORIGINAL <- as.character(run.match$SUBJECT_ORIGINAL)
      for (k in 1:nrow(run.match)) {
        work[which(work$originalRUN %in% run.match[k, 3:ncol(run.match)]), "newRun"] <- paste(
          paste(run.match[k, 1:3], collapse="_"), "merged", sep="_")
      }
      ## remove extra run NAs
      tmp <- work[is.na(work$ABUNDANCE), ]
      na.count <- reshape2::dcast(FEATURE ~ FRACTION, data=tmp,
                                  fun.aggregate=length, value.var="ABUNDANCE")
      na.count.long <- melt(na.count, id.vars=c("FEATURE"))
      na.count.long <- na.count.long[na.count.long$value == length(unique(work$newRun)), ]
      na.count.long$tmp <- paste(na.count.long$FEATURE, na.count.long$variable, sep="_")
      work$tmp <- paste(work$FEATURE, work$FRACTION, sep="_")
      work <- work[-which(work$tmp %in% na.count.long$tmp), ]
      ##
      work$originalRUN <- work$newRun
      ## update RUN based on new originalRUN
      work$RUN <- work$originalRUN
      work$RUN <- factor(work$RUN, levels=unique(work$RUN), labels=seq(1, length(unique(work$RUN))))
      work <- work[, -which(colnames(work) %in% c("tmp","newRun"))]
    }
  }
  return(work)
}

#' Examine the work data structure and quantile normalize by the heavy fraction.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @return Send the work data back to dataProcess()
normalize_singlefrac_by_quantile_heavy <- function(work) {
  nmethod <- unique(work$FRACTION)
  quantileall <- NULL
  for (j in 1:length(nmethod)) {
    namerun <- unique(work[work$FRACTION == nmethod[j], "RUN"])
    ## for label-based, make quantile normalization for reference
    ##worktemp <- work[which(work$RUN %in% namerun & work$LABEL=="H" &
    ##          !is.na(work$INTENSITY)),] ## because for sparse of reference
    worktemp <- work[which(work$RUN %in% namerun & work$LABEL == "H"),]
    worktemp$RUN <- factor(worktemp$RUN)
    worktemp$FEATURE <- factor(worktemp$FEATURE)
    quantiletemp <- as.matrix(xtabs(ABUNDANCE ~ FEATURE + RUN, data=worktemp))
    rm(worktemp)
    ## need to put NA for missing value in endogenous
    quantiletemp[quantiletemp == 0] <- NA
    ## using preprocessCore library
    quantiledone <- preprocessCore::normalize.quantiles(quantiletemp)
    rownames(quantiledone) <- rownames(quantiletemp)
    colnames(quantiledone) <- colnames(quantiletemp)
    ## get quantiled to long format for apply difference endogenous
    quantilelong.h <- melt(quantiledone, id=rownames(quantiledone))
    colnames(quantilelong.h) <- c("FEATURE","RUN","ABUNDANCE_quantile")
    quantilelong.h <- data.frame(quantilelong.h, LABEL="H")
    ## endogenous, in order to applying
    ##worktemp.l <- work[which(work$RUN %in% namerun & work$LABEL=="L" &
    ##       !is.na(work$INTENSITY)),] ## because for sparse of reference
    worktemp.l <- work[which(work$RUN %in% namerun & work$LABEL=="L"),]
    worktemp.l$RUN <- factor(worktemp.l$RUN)
    worktemp.l$FEATURE <- factor(worktemp.l$FEATURE)
    quantiletemp.l <- as.matrix(xtabs(ABUNDANCE ~ FEATURE + RUN, data=worktemp.l))
    rm(worktemp.l)
    ## need to put NA for missing value in endogenous
    quantiletemp.l[quantiletemp.l == 0] <- NA
    ## apply the difference from reference
    quantiledone.l <- quantiletemp.l-(quantiletemp-quantiledone)
    ## get quantiled to long format for apply difference endogenous
    quantilelong.l <- melt(quantiledone.l, id=rownames(quantiledone.l))
    colnames(quantilelong.l) <- c("FEATURE", "RUN", "ABUNDANCE_quantile")
    quantilelong.l <- data.frame(quantilelong.l, LABEL="L")
    rm(quantiletemp)
    rm(quantiledone)
    rm(quantiletemp.l)
    rm(quantiledone.l)
    ## quantileall <- rbindlist(list(quantileall,quantilelong.h, quantilelong.l))
    quantileall <- rbind(quantileall,quantilelong.h, quantilelong.l)
  }

  ## merge with original data
  work <- merge(work, quantileall, by=c("FEATURE","RUN","LABEL"))

  ## reorder
  work <- data.frame(
    "PROTEIN" = work$PROTEIN,
    "PEPTIDE" = work$PEPTIDE,
    "TRANSITION" = work$TRANSITION,
    "FEATURE" = work$FEATURE,
    "LABEL" = work$LABEL,
    "GROUP_ORIGINAL" = work$GROUP_ORIGINAL,
    "SUBJECT_ORIGINAL" = work$SUBJECT_ORIGINAL,
    "RUN" = work$RUN,
    "GROUP" = work$GROUP,
    "SUBJECT" = work$SUBJECT,
    "SUBJECT_NESTED" = work$SUBJECT_NESTED,
    "INTENSITY" = work$INTENSITY,
    "ABUNDANCE" = work$ABUNDANCE_quantile,
    "FRACTION" = work$FRACTION,
    "originalRUN" = work$originalRUN)

  work <- work[with(work,
                    order(LABEL, GROUP_ORIGINAL, SUBJECT_ORIGINAL,
                          RUN, PROTEIN, PEPTIDE, TRANSITION)),]
  return(work)
}

#' Examine the work data structure and normalize single fraction data by medians.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @return Send the work data back to dataProcess()
normalize_singlefrac_by_medians <- function(work) {
  ## Constant normalization by endogenous per method
  ## [MC : use median of medians]
  median.run.method <- aggregate(ABUNDANCE ~ RUN + FRACTION,
                                 data=work,
                                 median,
                                 na.rm=TRUE)
  median.method <- tapply(median.run.method$ABUNDANCE, median.run.method$FRACTION,
                          median, na.rm = TRUE)
  nmethod <- unique(work$FRACTION)
  for(j in 1:length(nmethod)) {
    namerun <- unique(work[work$FRACTION == nmethod[j], "RUN"])
    for (i in 1:length(namerun)) {
      ## ABUNDANCE is normalized
      namerun.idx <- which(work$RUN == namerun[i])
      work[namerun.idx, "ABUNDANCE"] <- work[namerun.idx, "ABUNDANCE"] -
        median.run.method[median.run.method$RUN == namerun[i], "ABUNDANCE"] +
        median.method[j]
    }
  }
  return(work)
}

#' Examine the work data structure and normalize single fraction data by heavy standards.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @return Send the work data back to dataProcess()
normalize_singlefrac_by_heavy_standard <- function(work) {
  ## Constant normalization by heavy standard per method
  h <- work[work$LABEL == "H", ]
  ## [MC : use median of medians]
  median.run.method <- aggregate(ABUNDANCE ~ RUN + FRACTION,
                                 data=h,
                                 median,
                                 na.rm=TRUE)
  median.method <- tapply(median.run.method$ABUNDANCE, median.run.method$FRACTION,
                          median, na.rm = TRUE)
  nmethod <- unique(work$FRACTION)
  for(j in 1:length(nmethod)) {
    namerun <- unique(work[work$FRACTION == nmethod[j],"RUN"])
    for (i in 1:length(namerun)) {
      ## ABUNDANCE is normalized
      namerun.idx <- which(work$RUN == namerun[i])
      work[namerun.idx, "ABUNDANCE"] <- work[namerun.idx, "ABUNDANCE"] -
        median.run.method[median.run.method$RUN == namerun[i], "ABUNDANCE"] +
        median.method[j]
    }
  } # end loop method
  return(work)
}

#' Examine the work data structure and summarize it by mbimpute.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @param maxQuantileforCensored  What is the censoring cutoff by quantile?
#' @param censoredInt "NA!" not NA?
#' @return Send the work data back to dataProcess()
summarize_tmp_mbimpute <- function(work, maxQuantileforCensored=0.999, censoredInt="NA") {
  work$LABEL <- factor(work$LABEL)
  label <- nlevels(work$LABEL) == 2
  work$censored <- FALSE
  ## if intensity = 1, but abundance > cutoff after normalization, it also should be censored.
  if (!is.null(maxQuantileforCensored)) {
    ## label-free
    if (!label) {
      ## calculate outlier cutoff
      ## only consider intensity > 1
      tmp <- work[!is.na(work$INTENSITY) & work$INTENSITY > 1, "ABUNDANCE"]
      ## or
      ##tmp <- work[!is.na(work$INTENSITY), "ABUNDANCE"]
      log2int.prime.quant <- quantile(tmp, prob=c(0.01, 0.25, 0.5, 0.75,
                                                  maxQuantileforCensored), na.rm = TRUE)
      iqr <- log2int.prime.quant[4] - log2int.prime.quant[2]

      ## need to decide the multiplier from high intensities
      multiplier <- (log2int.prime.quant[5] - log2int.prime.quant[4]) / iqr
      cutoff.lower <- (log2int.prime.quant[2] - multiplier * iqr)
      work[!is.na(work$INTENSITY) &
           work$ABUNDANCE < cutoff.lower, "censored"] <- TRUE
      print_string <- paste0("** Log2 intensities under cutoff =",
                             format(cutoff.lower, digits=5),
                             " were considered as censored missing values.")
      logging::loginfo(print_string)

      ## if censoredInt == "0, and cutoff is negative, still zero should becensored
      if (cutoff.lower <= 0 & !is.null(censoredInt) & censoredInt == "0") {
        work[!is.na(work$INTENSITY) & work$INTENSITY == 1, "censored"] <- TRUE
        work[!is.na(work$ABUNDANCE) & work$ABUNDANCE <= 0, "censored"] <- TRUE
        print_string <- "** Log2 intensities = 0 were considered as censored missing values."
        logging::loginfo(print_string)
      }

      ## if censoredInt == NA, original NA also shoule be "censored"
      if (!is.null(censoredInt) & censoredInt == "NA") {
        work[is.na(work$INTENSITY), "censored"] <- TRUE
        print_string <- "** Log2 intensities = NA were considered as censored missing values."
        logging::loginfo(print_string)
      }
    }

    ## labeled : only consider light. Assume that missing in heavy is random.
    if(label) {
      work.tmp <- work[which(work$LABEL %in% "L"), ]
      ## calculate outlier cutoff
      ## only consider intensity > 1
      tmp <- work.tmp[!is.na(work.tmp$INTENSITY) & work.tmp$INTENSITY > 1, "ABUNDANCE"]
      log2int.prime.quant <- quantile(tmp,
                                      prob=c(0.01, 0.25, 0.5, 0.75,
                                             maxQuantileforCensored), na.rm = TRUE)
      iqr <- log2int.prime.quant[4] - log2int.prime.quant[2]
      ## need to decide the multiplier from high intensities
      multiplier <- (log2int.prime.quant[5] - log2int.prime.quant[4]) / iqr
      cutoff.lower <- (log2int.prime.quant[2] - multiplier * iqr)
      ##work$censored <- FALSE
      work[work$LABEL == "L" & !is.na(work$INTENSITY) &
           work$ABUNDANCE < cutoff.lower, "censored"] <- TRUE
      print_string <- paste0("** Log2 endogenous intensities under cutoff =",
                             format(cutoff.lower, digits=5),
                             " were considered as censored missing values.")
      logging::loginfo(print_string)

      ## if censoredInt == "0, and cutoff is negative, still zero should becensored
      if (cutoff.lower <= 0 & !is.null(censoredInt) & censoredInt == "0") {
        work[work$LABEL == "L" &
             !is.na(work$INTENSITY) & work$INTENSITY == 1, "censored"] <- TRUE
        work[work$LABEL == "L" &
             !is.na(work$ABUNDANCE) & work$ABUNDANCE <= 0, "censored"] <- TRUE
        print_string <- "** Log2 endogenous intensities = 0 were considered as censored missing values."
        logging::loginfo(print_string)
      }

      ## if censoredInt == NA, original NA also shoule be "censored"
      if (!is.null(censoredInt) & censoredInt == "NA") {
        work[work$LABEL == "L" & is.na(work$INTENSITY), "censored"] <- TRUE
        print_string <- "** Log2 endogenous intensities = NA were considered as censored missing values."
        logging::loginfo(print_string)
      }
    }

  } else { ## will MBimpute, but not apply algorithm for cutoff
    if (censoredInt == "0") {
      work[work$LABEL == "L" & !is.na(work$INTENSITY) & work$INTENSITY == 1, "censored"] <- TRUE
      work[work$LABEL == "L" & !is.na(work$ABUNDANCE) & work$ABUNDANCE <= 0, "censored"] <- TRUE
    }
    if (censoredInt == "NA") {
      work[work$LABEL == "L" & is.na(work$ABUNDANCE), "censored"] <- TRUE
    }
  } ## End of MBimpute testing
  return(work)
}

#' Extract the top3 for each result.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @return Send the work data back to dataProcess()
subset_top3_result <- function(work) {
  logging::loginfo("** Use top3 features that have highest average of log2(intensity) across runs.")
  ## INTENSITY vs ABUNDANCE? [THT: make more sense to use ABUNDANCE]
  ## how to decide top3 for DIA?
  work$remove <- FALSE
  temp1 <- aggregate(INTENSITY ~ PROTEIN + FEATURE,
                     data=work,
                     function(x) mean(x, na.rm=TRUE))
  temp2 <- split(temp1, temp1$PROTEIN)
  temp3 <- lapply(temp2, function(x) {
    x <- x[order(x$INTENSITY, decreasing=TRUE),]
    x <- x$FEATURE[1:3]
  })
  selectfeature <- unlist(temp3, use.names=FALSE)
  selectfeature <- selectfeature[!is.na(selectfeature)]
  ## get subset
  work[-which(work$FEATURE %in% selectfeature), "remove"] <- TRUE
  return(work)
}

#' take the top n for each result.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @return Send the work data back to dataProcess()
subset_topn_result <- function(work) {
  ## check whether there is the input for "N"
  logging::loginfo(paste0("** Use top", n_top_feature,
                          " features that have highest average of log2(intensity) across runs."))
  ## INTENSITY vs ABUNDANCE? [THT: make more sense to use ABUNDANCE]
  ## how to decide top3 for DIA?
  work$remove <- FALSE
  worktemp <- work[!is.na(work$ABUNDANCE) & work$ABUNDANCE != 0, ]
  temp1 <- aggregate(INTENSITY ~ PROTEIN + FEATURE,
                     data=worktemp,
                     function(x) mean(x, na.rm=TRUE))
  temp2 <- split(temp1, temp1$PROTEIN)
  temp3 <- lapply(temp2, function(x) {
    x <- x[order(x$INTENSITY, decreasing=TRUE), ]
    x <- x$FEATURE[1:n_top_feature]
  })
  selectfeature <- unlist(temp3, use.names=FALSE)
  selectfeature <- selectfeature[!is.na(selectfeature)]
  ## get subset
  work[-which(work$FEATURE %in% selectfeature), "remove"] <- TRUE
  return(work)
}

#' Examine the work data structure and figure out how much stuff is missing.
#'
#' This was just kill/yanked from dataProcess(), hopefully it does not depend on
#' anything there.
#'
#' @param work  The work data structure from dataProcess()
#' @return Send the work data back to dataProcess()
decide_missingness <- function(work) {
  if (nlevels(work$LABEL) == 1) {
    ##Use the data frame before imputation to summarize the missingness
    all.work <- work
    test <- tapply(is.na(work[, "ABUNDANCE"]),
                   work[, c("GROUP_ORIGINAL", "FEATURE")],
                   function(x) sum(x, na.rm=TRUE))
    numObs <- tapply(work[, "ABUNDANCE"],
                     work[, c("GROUP_ORIGINAL", "FEATURE")],
                     function(x) length(x))
    test1 <- test == numObs
    test2 <- apply(test1, 2, function(x) sum(x, na.rm=TRUE))
    filterList <- names(test2)[test2 > 0]
    final.decision <- ifelse(test2 > 0, 1, 0)
  }

  if (nlevels(work$LABEL) == 2) {
    ##Use the data frame before imputation to summarize the missingness
    ## first, remove NA
    all.work <- work  # with all NA observations
    work.miss <- na.omit(work)
    ## draw table
    light <- subset(work.miss, LABEL=="L")
    heavy <- subset(work.miss, LABEL=="H")
    ## use FEATURE because the name of transition can be used in other peptide
    count.light <- xtabs(~ FEATURE + GROUP_ORIGINAL, light)
    count.heavy <- xtabs(~ FEATURE + GROUP_ORIGINAL, heavy)
    count.light <- count.light == 0
    count.heavy <- count.heavy == 0
    count.light <- as.data.frame(count.light)
    count.heavy <- as.data.frame(count.heavy)
    ## summary of missingness
    decision <- count.light
    decision[] <- 0
    for (i in 1:ncol(decision)) {
      for (j in 1:nrow(decision)) {
        ## either light or heavy has no obs -> subject to filter
        if (count.light[j, i]==TRUE || count.heavy[j, i]==TRUE) {
          decision[j, i] <- 1
        }
      }
    }

    final.decision <- apply(decision, 1, sum)
    ## assign "subject to filter" column
    work <- data.frame(work, "SuggestToFilter"=0)
    for (i in 1:length(final.decision)) {
      ## assign subject_to_filter=1 for entire transition
      if (final.decision[i] != 0) {
        work[work$FEATURE == names(final.decision[i]), "SuggestToFilter"] <- 1
      }
    }
  }
  return(final.decision)
}

#' I need to read up on this
#'
#' I am familiar with some statistics tests from Tukey, but not the polish median.
#' I know for a fact that not all of these options are required, but I am not
#' ready to remove them yet.
#'
#' @param modified_data  The data provided to msstats after some sanitization.
#' @param label  Is this data labeled?
#' @param pct_adjust  This is a peculiar feature of this code, this adjustment
#'   is applied liberally throughout, and I am not sure yet which it is not done
#'   just once at the beginning.  But until I figure that out I will not mess
#'   with it.
#' @param clusters  Currently not used until I clean it up enough that I feel
#'   confident re-invoking parallel
#' @param MBimpute  Perform imputation?
#' @param censoredInt  "NA", not NA!
#' @param message.show  Only print the peptides which fail when false.
#' @param cutoffCensored  Where is the cutoff for censoring?
#' @param original_scale  Keep the data on its original scale?
#' @param remove50missing  If coverage is < 50%, drop the protein?
#' @param ...  I might add more stuff, so drop it into arglist.
#' @return a list of tukeys!
perform_tukey_polish_median <- function(modified_data, label=1,
                                        pct_adjust=0.99, clusters=1,
                                        MBimpute=TRUE, censoredInt="NA",
                                        message.show=FALSE,
                                        cutoffCensored="minFeature",
                                        original_scale=FALSE,
                                        remove50missing=FALSE, ...) {
  arglist <- list(...)
  ## create cluster for paralleled workflow
  ## Bizarrely, the default cluster size is 1, so there is exactly
  ## 0 utility in using SNOW.  Thus for debugging I will just do a for loop, as I suspect
  ## I will probably smash the lower level function shortly.

  ## FIXME: I am breaking the parallelization of this function!
  ## message(paste0("Cluster Size: ", clusters, "\n"))
  ## doSNOW::registerDoSNOW(parallel::makeCluster(clusters, type="SOCK"))
  ##
  ##pb <- txtProgressBar(max=nlevels(data$PROTEIN), style=3)
  ##progress <- function(n) setTxtProgressBar(pb, n)
  ##opts <- list(progress=progress)
  ## Emacs does not handle indenting well for lines with foreach(stuff) %dopar% {
  ##MS_results <- foreach(i=1: nlevels(data[["PROTEIN"]]), .combine="resultsAsLists", .options.snow=opts, .multicombine=TRUE, .init=list(list(), list())) %dopar% {

  tukey_result <- list()
  end <- nlevels(modified_data[["PROTEIN"]])
  bar <- utils::txtProgressBar(max=end, style=3)
  for (i in 1:end) {
    pct_done <- i / end
    setTxtProgressBar(bar, pct_done)
    sub <- modified_data[modified_data[["PROTEIN"]] == levels(modified_data[["PROTEIN"]])[i], ]
    if (message.show) {
      message(paste0("Getting the summary by Tukey's median polish per subplot for protein ",
                     unique(sub$PROTEIN), "(",
                     i, " of ", length(unique(modified_data$PROTEIN)), ")"))
    }
    sub[["FEATURE"]] <- factor(sub[["FEATURE"]])
    sub[["feature.label"]] <- paste(sub[["FEATURE"]], sub[["LABEL"]], sep="_")
    sub[["run.label"]] <- paste(sub[["RUN"]], sub[["LABEL"]], sep="_")

    ## how to decide censored or not
    if (isTRUE(MBimpute)) {
      if (!is.null(censoredInt)) {
        ## 1. censored
        if (censoredInt == "0") {
          sub[sub[["censored"]] == TRUE, "ABUNDANCE"] <- 0
          sub[["cen"]] <- ifelse(sub[["censored"]], 0, 1)
        }
        ## 2. all censored missing
        if (censoredInt == "NA") {
          sub[sub[["censored"]] == TRUE, "ABUNDANCE"] <- NA
          sub[["cen"]] <- ifelse(sub[["censored"]], 0, 1)
        }
      }
    }

    ## if all measurements are NA,
    if (nrow(sub) == (sum(is.na(sub[["ABUNDANCE"]])) +
                      sum(!is.na(sub[["ABUNDANCE"]]) &
                          sub[["ABUNDANCE"]] == 0))) {
      message(paste0("Can't summarize for ",
                     unique(sub[["PROTEIN"]]), "(",i,
                     " of ", length(unique(modified_data[["PROTEIN"]])),
                     ") because all measurements are NAs."))
      next()
      ## return(NULL)
    }

    ## remove features which are completely NAs
    if (isTRUE(MBimpute)) {
      if (!is.null(censoredInt)) {
        ## 1. censored
        if (censoredInt == "0") {
          subtemp <- sub[sub[["LABEL"]] == "L" &
                         !is.na(sub[["ABUNDANCE"]]) &
                         sub[["ABUNDANCE"]] != 0, ]
        }

        ## 2. all censored missing
        if (censoredInt=="NA") {
          subtemp <- sub[sub[["LABEL"]] == "L" & !is.na(sub[["ABUNDANCE"]]), ]
        }
      }
    } else {
      subtemp <- sub[sub[["LABEL"]] == "L" &
                     !is.na(sub[["ABUNDANCE"]]) &
                     sub[["ABUNDANCE"]] != 0, ]
    }

    countfeature <- xtabs(~ FEATURE, subtemp)
    namefeature <- names(countfeature)[countfeature == 0]
    if (length(namefeature) != 0) {
      sub <- sub[-which(sub[["FEATURE"]] %in% namefeature), ]
      if (nrow(sub) == 0) {
        message(paste0("Can't summarize for ",
                       unique(sub[["PROTEIN"]]), "(", i,
                       " of ", length(unique(modified_data[["PROTEIN"]])),
                       ") because all measurements are NAs."))
        next()
        ## return(NULL)
      } else {
        sub[["FEATURE"]] <- factor(sub[["FEATURE"]])
      }
    }

    ## remove features which have only 1 measurement.
    namefeature1 <- names(countfeature)[countfeature == 1]
    if (length(namefeature1) != 0) {
      sub <- sub[-which(sub[["FEATURE"]] %in% namefeature1), ]
      if (nrow(sub) == 0) {
        message(paste0("Can't summarize for ",
                       unique(sub[["PROTEIN"]]), "(", i, " of ",
                       length(unique(modified_data[["PROTEIN"]])),
                       ") because features have only one measurement across MS runs."))
        next()
        ## return(NULL)
      } else {
        sub[["FEATURE"]] <- factor(sub[["FEATURE"]])
      }
    }

    ## check one more time
    ## if all measurements are NA,
    if (nrow(sub) == (sum(is.na(sub[["ABUNDANCE"]])) +
                      sum(!is.na(sub[["ABUNDANCE"]]) & sub[["ABUNDANCE"]] == 0))) {
      message(paste0("After removing features which has only 1 measurement, Can't summarize for ",
                     unique(sub[["PROTEIN"]]), "(", i,
                     " of ", length(unique(modified_data[["PROTEIN"]])),
                     ") because all measurements are NAs."))
      next()
    }

    ## remove run which has no measurement at all
    ## remove features which are completely NAs
    if (isTRUE(MBimpute)) {
      if (!is.null(censoredInt)) {
        ## 1. censored
        if (censoredInt == "0") {
          subtemp <- sub[sub[["LABEL"]] == "L" &
                         !is.na(sub[["ABUNDANCE"]]) &
                         sub[["ABUNDANCE"]] != 0, ]
        }
        ## 2. all censored missing
        if (censoredInt == "NA") {
          subtemp <- sub[sub[["LABEL"]] == "L" &
                         !is.na(sub[["ABUNDANCE"]]), ]
        }
      }
    } else {
      subtemp <- sub[sub[["LABEL"]] == "L" &
                     !is.na(sub[["ABUNDANCE"]]) &
                     sub[["ABUNDANCE"]] != 0, ]
    }
    count <- aggregate(ABUNDANCE ~ RUN, data=subtemp, length)
    norun <- setdiff(unique(modified_data[["RUN"]]), count[["RUN"]])
    if (length(norun) != 0 &
        length(intersect(norun, as.character(unique(sub[["RUN"]]))))) {
      ## removed NA rows already, if there is no overlapped run, error
      sub <- sub[-which(sub[["RUN"]] %in% norun), ]
      sub[["RUN"]] <- factor(sub[["RUN"]])
    }

    if (isTRUE(remove50missing)) {
      ## count # feature per run
      if (!is.null(censoredInt)) {
        if (censoredInt == "NA") {
          subtemp <- sub[sub[["LABEL"]] == "L" &
                         !is.na(sub[["INTENSITY"]]), ]
        }
        if (censoredInt == "0") {
          subtemp <- sub[sub[["LABEL"]] == "L" &
                         !is.na(sub[["INTENSITY"]]) &
                         sub[["INTENSITY"]] != 0, ]
        }
      }

      numFea <- xtabs(~ RUN, subtemp) ## RUN or run.label?
      numFea <- numFea / length(unique(subtemp[["FEATURE"]]))
      numFea <- numFea <= 0.5
      removerunid <- names(numFea)[numFea]
      ## if all measurements are NA,
      if (length(removerunid) == length(numFea)) {
        message(paste0("Can't summarize for ",
                       unique(sub[["PROTEIN"]]), "(",
                       i," of ",
                       length(unique(modified_data[["PROTEIN"]])),
                       ") because all runs have more than 50% NAs and are removed with the option, remove50missing=TRUE."))
        next()
        ## return(NULL)
      }
    }

    ## This section is disturbingly long and I think should be
    ## split off into one or more functions.
    ## check whether we need to impute or not.
    if (sum(sub[["cen"]] == 0) > 0) {
      ## cutoffCensored
      ## 1. put 0 to censored
      ##if (cutoffCensored=="0") {
      ##  if (censoredInt=="NA") {
      ##      sub[is.na(sub$INTENSITY),"ABUNDANCE"] <- 0
      ##  }
      ##  if (censoredInt=="0") {
      ##      sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"] <- 0
      ##  }
      ##}

      ## 2. put minimum in feature level to NA
      if (cutoffCensored == "minFeature") {
        if (censoredInt == "NA") {
          cut <- aggregate(ABUNDANCE ~ feature.label,
                           data=sub,
                           function(x) {
                             min(x, na.rm=TRUE)
                           })
          ## cutoff for each feature is little less than minimum abundance in a run.
          ## I am thinking this 0.99 should be parameterized.
          cut[["ABUNDANCE"]] <- pct_adjust * cut[["ABUNDANCE"]]
          ## remove runs which has more than 50% missing values
          if (isTRUE(remove50missing)) {
            if (length(removerunid) != 0) {
              sub <- sub[-which(sub[["RUN"]] %in% removerunid), ]
              sub[["RUN"]] <- factor(sub[["RUN"]])
            }
          }
          for (j in 1:length(unique(cut[["feature.label"]]))) {
            ## The following line is entirely more difficult than it needs to be
            ## I would suggest splitting into two easier to read: 1 where you take the index
            ## and another where you set it.  friggin' smart people.
            sub[is.na(sub[["ABUNDANCE"]]) &
                sub[["feature.label"]] == cut[["feature.label"]][j], "ABUNDANCE"] <- cut[["ABUNDANCE"]][j]
          }
        } else if (censoredInt == "0") {
          subtemptemp <- sub[!is.na(sub[["ABUNDANCE"]]) & sub[["ABUNDANCE"]] != 0, ]
          cut <- aggregate(ABUNDANCE ~ feature.label, data=subtemptemp, FUN=min)
          ## cutoff for each feature is little less than minimum abundance in a run.
          cut[["ABUNDANCE"]] <- pct_adjust * cut[["ABUNDANCE"]]

          ## remove runs which has more than 50% missing values
          if (isTRUE(remove50missing)) {
            if (length(removerunid) != 0) {
              sub <- sub[-which(sub[["RUN"]] %in% removerunid), ]
              sub[["RUN"]] <- factor(sub[["RUN"]])
            }
          }
          for (j in 1:length(unique(cut[["feature.label"]]))) {
            sub[!is.na(sub[["ABUNDANCE"]]) &
                sub[["ABUNDANCE"]] == 0 &
                sub[["feature.label"]] == cut[["feature.label"]][j],
                "ABUNDANCE"] <- cut[["ABUNDANCE"]][j]
          }
        }  ## End if censoredInt is 0
      } ## End of minfeature cutoffcensored

      ## 3. put minimum in RUN to NA
      if (cutoffCensored == "minRun") {
        ## remove runs which has more than 50% missing values
        if (remove50missing) {
          if (length(removerunid) != 0) {
            sub <- sub[-which(sub[["RUN"]] %in% removerunid), ]
            sub[["RUN"]] <- factor(sub[["RUN"]])
          }
        }
        if (censoredInt == "NA") {
          cut <- aggregate(ABUNDANCE ~ run.label,
                           data=sub,
                           function(x) {
                             min(x, na.rm=TRUE)
                           })
          ## cutoff for each Run is little less than minimum abundance in a run.
          ## Maybe instead of doing this multiplication in like 8 places, just do it once
          ## at the beginning of this code block?
          cut[["ABUNDANCE"]] <- pct_adjust * cut[["ABUNDANCE"]]
          for (j in 1:length(unique(cut[["run.label"]]))) {
            sub[is.na(sub[["ABUNDANCE"]]) &
                sub[["run.label"]] == cut[["run.label"]][j], "ABUNDANCE"] <- cut[["ABUNDANCE"]][j]
          }
        }

        if (censoredInt == "0") {
          subtemptemp <- sub[!is.na(sub[["ABUNDANCE"]]) & sub[["ABUNDANCE"]] != 0, ]
          cut <- aggregate(ABUNDANCE ~ run.label,
                           data=subtemptemp,
                           FUN=min)
          cut[["ABUNDANCE"]] <- pct_adjust * cut[["ABUNDANCE"]]
          for (j in 1:length(unique(cut[["run.label"]]))) {
            sub[!is.na(sub[["ABUNDANCE"]]) &
                sub[["ABUNDANCE"]] == 0 &
                sub[["run.label"]] == cut[["run.label"]][j], "ABUNDANCE"] <- cut[["ABUNDANCE"]][j]
          }
        }
      } ## End of minRun block

      ## 20150829 : 4. put minimum RUN and FEATURE
      if (cutoffCensored == "minFeatureNRun") {
        if (censoredInt == "NA") {
          ## cutoff for each feature is little less than minimum abundance in a run.
          cut.fea <- aggregate(ABUNDANCE ~ feature.label,
                               data=sub,
                               function(x) {
                                 min(x, na.rm=TRUE)
                               })
          cut.fea[["ABUNDANCE"]] <- pct_adjust * cut.fea[["ABUNDANCE"]]
          ## remove runs which has more than 50% missing values
          ## before removing, need to contribute min feature calculation
          if (isTRUE(remove50missing)) {
            if (length(removerunid) != 0) {
              sub <- sub[-which(sub[["RUN"]] %in% removerunid), ]
              sub[["RUN"]] <- factor(sub[["RUN"]])
            }
          }
          ## cutoff for each Run is little less than minimum abundance in a run.
          cut.run <- aggregate(ABUNDANCE ~ run.label,
                               data=sub,
                               function(x) {
                                 min(x, na.rm=TRUE)
                               })
          cut.run[["ABUNDANCE"]] <- pct_adjust * cut.run[["ABUNDANCE"]]
          if (length(unique(cut.fea[["feature.label"]])) > 1) {
            for (j in 1:length(unique(cut.fea[["feature.label"]]))) {
              for (k in 1:length(unique(cut.run[["run.label"]]))) {
                ## get smaller value for min Run and min Feature
                finalcut <- min(cut.fea[["ABUNDANCE"]][j], cut.run[["ABUNDANCE"]][k])
                sub[is.na(sub[["ABUNDANCE"]]) &
                    sub[["feature.label"]] == cut.fea[["feature.label"]][j] &
                    sub[["run.label"]] == cut.run[["run.label"]][k], "ABUNDANCE"] <- finalcut
              }
            }
          }
          ## if single feature, not impute
        } ## End of minFeatureNRun and censoredInt is NA

        if (censoredInt == "0") {
          subtemptemp <- sub[!is.na(sub[["ABUNDANCE"]]) & sub[["ABUNDANCE"]] != 0,]
          cut.fea <- aggregate(ABUNDANCE ~ feature.label, data=subtemptemp, FUN=min)
          cut.fea[["ABUNDANCE"]] <- pct_adjust * cut.fea[["ABUNDANCE"]]
          ## remove runs which has more than 50% missing values
          ## before removing, need to contribute min feature calculation
          if (isTRUE(remove50missing)) {
            if (length(removerunid) != 0) {
              sub <- sub[-which(sub[["RUN"]] %in% removerunid),]
              sub[["RUN"]] <- factor(sub[["RUN"]])
            }
          }

          cut.run <- aggregate(ABUNDANCE ~ run.label, data=subtemptemp, FUN=min)
          cut.run[["ABUNDANCE"]] <- pct_adjust * cut.run[["ABUNDANCE"]]
          if (length(unique(cut.fea[["feature.label"]])) > 1) {
            for (j in 1:length(unique(cut.fea[["feature.label"]]))) {
              for (k in 1:length(unique(cut.run[["run.label"]]))) {
                ## get smaller value for min Run and min Feature
                finalcut <- min(cut.fea[["ABUNDANCE"]][j], cut.run[["ABUNDANCE"]][k])

                sub[!is.na(sub[["ABUNDANCE"]]) &
                    sub[["ABUNDANCE"]] == 0 &
                    sub[["feature.label"]] == cut.fea[["feature.label"]][j] &
                    sub[["run.label"]] == cut.run[["run.label"]][k], "ABUNDANCE"] <- finalcut
              }
            }
          } else { # single feature
            sub[!is.na(sub[["ABUNDANCE"]]) &
                sub[["ABUNDANCE"]] == 0, "ABUNDANCE"] <- cut.fea[["ABUNDANCE"]]
          }
        }
      } ## End of cutoff == minFeatureNRun

      if (isTRUE(MBimpute)) {
        if (!label) { ## label-free
          if (nrow(sub[sub[["cen"]] == 0, ]) > 0) {
            ## impute by survival model
            subtemp <- sub[!is.na(sub[["ABUNDANCE"]]), ]
            countdf <- nrow(subtemp) <
              (length(unique(subtemp[["FEATURE"]])) + length(unique(subtemp[["RUN"]])) - 1)

            ## fit the model
            if (length(unique(sub[["FEATURE"]])) == 1) {
              fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type="left") ~ RUN,
                                           data=sub,
                                           dist="gaussian")
            } else {
              if (countdf) {
                fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type="left") ~ RUN,
                                             data=sub,
                                             dist="gaussian")
              } else {
                fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type="left") ~ FEATURE + RUN,
                                             data=sub,
                                             dist="gaussian")
              }
            }

            ## get predicted value from survival
            sub <- data.frame(sub,
                              "pred" = predict(fittest, new_data=sub, type="response"))
            ## the replace censored value with predicted value
            sub[sub[["cen"]] == 0, "ABUNDANCE"] <- sub[sub[["cen"]] == 0, "pred"]

            ## save predicted value
            ## predAbundance <- c(predAbundance,predict(fittest, newmodified_data=sub, type="response"))
            ##predAbundance <- c(predict(fittest, newmodified_data=sub, type="response"))
          }
        } else { ## label-based
          ## only endogenous will be imputed
          sub.h <- sub[sub[["LABEL"]] == "H", ]
          sub.l <- sub[sub[["LABEL"]] == "L", ]

          if (nrow(sub.l[sub.l[["cen"]] == 0, ]) > 0) {
            ## impute by survival model
            subtemp <- sub.l[!is.na(sub.l[["ABUNDANCE"]]), ]
            countdf <- nrow(subtemp) < (length(unique(subtemp[["FEATURE"]])) +
                                        length(unique(subtemp[["RUN"]])) - 1)
            ##set.seed(100)
            ## fit the model
            if (length(unique(sub.l[["FEATURE"]])) == 1) {
              fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type="left") ~ RUN,
                                           data=sub.l,
                                           dist="gaussian")
            } else {
              if (countdf) {
                fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type="left") ~ RUN,
                                             data=sub.l,
                                             dist="gaussian")
              } else {
                fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type="left") ~ FEATURE + RUN,
                                             data=sub.l,
                                             dist="gaussian")
              }
            }

            ## get predicted value from survival
            sub.l <- data.frame(sub.l,
                                "pred" = predict(fittest, new_data=sub.l, type="response"))

            ## predAbundance <- c(predAbundance,predict(fittest, newmodified_data=sub, type="response"))
            ## predAbundance <- c(predict(fittest, newmodified_data=sub.l, type="response"))
            ## the replace censored value with predicted value
            sub.l[sub.l[["cen"]] == 0, "ABUNDANCE"] <- sub.l[sub.l[["cen"]] == 0, "pred"]
            sub.h[["pred"]] <- NA
            ## for label-based, need to merge again
            sub <- rbind(sub.h, sub.l)
          }
        } ## No more mbimpute of label-based
      }
    }  ## end if sub$cen == 0

    ## then, finally remove NA in abundance
    sub <- sub[!is.na(sub[["ABUNDANCE"]]), ]
    if (nlevels(sub[["FEATURE"]]) > 1) { ## for more than 1 features
      if (!label) { ## label-free
        modified_data_w <- reshape2::dcast(RUN ~ FEATURE,
                                           data=sub,
                                           value.var="ABUNDANCE",
                                           keep=TRUE,
                                           ## I am concerned that the following line is incorrect!!
                                           ## FIXME FIXME!
                                           fun.aggregate=sum)
        ## RIGHT ABOVE, IS THAT CORRECT!?!?! FIXME FIXME FIXME FIXME
        rownames(modified_data_w) <- modified_data_w[["RUN"]]
        modified_data_w <- modified_data_w[, -1]
        modified_data_w[modified_data_w == 1] <- NA
        if (isTRUE(original_scale)) {
          modified_data_w <- 2 ^ modified_data_w
          medmodified_data <- medpolish(modified_data_w, na.rm=TRUE, trace.iter=FALSE)
          tmpresult <- medmodified_data[["overall"]] + medmodified_data[["row"]]
          tmpresult <- log2(tmpresult)
        } else {
          medmodified_data <- medpolish(modified_data_w, na.rm=TRUE, trace.iter=FALSE)
          tmpresult <- medmodified_data[["overall"]] + medmodified_data[["row"]]
          ## if fractionated sample, need to get per sample run
          ## ?? if there are technical replicates, how to match sample and MS run for different fractionation??
          ##if( length(unique(sub$METHOD)) > 1 ) {
          ## runinfo <- unique(sub[, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "RUN", "METHOD")])
          ## runinfo$uniquesub <- paste(runinfo$GROUP_ORIGINAL, runinfo$SUBJECT_ORIGINAL, sep="_")
          ##}
        }
        ## count # feature per run
        if (!is.null(censoredInt)) {
          if (censoredInt=="NA") {
            subtemp <- sub[!is.na(sub[["INTENSITY"]]), ]
            subtempimpute <- sub[is.na(sub[["INTENSITY"]]), ]
            subtempimpute <- subtempimpute[!is.na(subtempimpute[["ABUNDANCE"]]), ]
          }
          if (censoredInt == "0") {
            subtemp <- sub[!is.na(sub[["INTENSITY"]]) & sub[["INTENSITY"]] != 0, ]
            subtempimpute <- sub[!is.na(sub[["INTENSITY"]]) & sub[["INTENSITY"]] == 0, ]
            subtempimpute <- subtempimpute[!is.na(subtempimpute[["ABUNDANCE"]]) &
                                           subtempimpute[["ABUNDANCE"]] != 0, ]
          }

          subtemp[["RUN"]] <- factor(subtemp[["RUN"]], levels=rownames(modified_data_w))
          numFea <- xtabs(~ RUN, subtemp)
          numFeaPercentage <- 1 - numFea / length(unique(subtemp[["FEATURE"]]))
          numFeaTF <- numFeaPercentage >= 0.5
          subtempimpute[["RUN"]] <- factor(subtempimpute[["RUN"]], levels=rownames(modified_data_w))
          numimpute <- xtabs(~ RUN, subtempimpute)
          sub.result <- data.frame(
            "Protein" = unique(sub[["PROTEIN"]]),
            "LogIntensities" = tmpresult,
            "RUN" = names(tmpresult),
            "NumMeasuredFeature" = as.vector(numFea),
            "MissingPercentage" = as.vector(numFeaPercentage),
            "more50missing" = numFeaTF,
            "NumImputedFeature" = as.vector(numimpute))
        } else {
          subtemp <- sub[!is.na(sub[["INTENSITY"]]), ]
          subtemp[["RUN"]] <- factor(subtemp[["RUN"]], levels=rownames(modified_data_w))
          numFea <- xtabs(~ RUN, subtemp)
          numFeaPercentage <- 1 - numFea / length(unique(subtemp[["FEATURE"]]))
          numFeaTF <- numFeaPercentage >= 0.5
          sub.result <- data.frame(
            "Protein" = unique(sub[["PROTEIN"]]),
            "LogIntensities" = tmpresult,
            "RUN" = names(tmpresult),
            "NumMeasuredFeature" = as.vector(numFea),
            "MissingPercentage" = as.vector(numFeaPercentage),
            "more50missing" = numFeaTF)
        }
        ## result <- rbind(result, sub.result)
      } else { ## labeled
        modified_data_w = reshape2::dcast(run.label ~ FEATURE,
                                          data=sub,
                                          value.var="ABUNDANCE",
                                          keep=TRUE,
                                          fun.aggregate=sum)
        ### SAME QUESTION HERE FIXME FIXME FIXME
        rownames(modified_data_w) <- modified_data_w[["run.label"]]
        modified_data_w <- modified_data_w[, -1]
        ##modified_data_w[modified_data_w==1] <- NA
        medmodified_data <- medpolish(modified_data_w, na.rm=TRUE, trace.iter=FALSE)
        tmpresult <- medmodified_data[["overall"]] + medmodified_data[["row"]]

        reformresult <- data.frame(tmpresult)
        end <- nchar(rownames(reformresult))
        reformresult[["LABEL"]] <- substr(rownames(reformresult), end, end)
        reformresult[["RUN"]] <- substr(rownames(reformresult), 1, end - 2)
        colnames(reformresult)[1] <- "ABUNDANCE"
        ## now single feature, adjust reference feature difference
        h <- reformresult[reformresult[["LABEL"]] == "H", ]
        allmed <- median(h[["ABUNDANCE"]], na.rm=TRUE)
        for (k in 1:length(unique(h[["RUN"]]))) {
          ## ABUNDANCE is normalized
          reformresult.logical <- reformresult[["RUN"]] == unique(h[["RUN"]])[k]
          reformresult.idx <- which(reformresult.logical)
          reformresult[reformresult.idx, "ABUNDANCE"] <- reformresult[reformresult.idx, "ABUNDANCE"] -
            reformresult[reformresult.logical & reformresult[["LABEL"]] == "H", "ABUNDANCE"] + allmed
        }

        reformresult <- reformresult[reformresult[["LABEL"]] == "L", ]
        subtemp <- reformresult[!is.na(reformresult[["ABUNDANCE"]]), ]
        ## count # feature per run
        if (!is.null(censoredInt)) {
          if (censoredInt == "NA") {
            subtemp <- sub[sub[["LABEL"]] == "L" & !is.na(sub[["INTENSITY"]]), ]
            subtempimpute <- sub[sub[["LABEL"]] == "L" & is.na(sub[["INTENSITY"]]), ]
            subtempimpute <- subtempimpute[!is.na(subtempimpute[["ABUNDANCE"]]), ]
          }

          if (censoredInt == "0") {
            subtemp <- sub[sub[["LABEL"]] == "L" &
                           !is.na(sub[["INTENSITY"]]) &
                           sub[["INTENSITY"]] != 0, ]
            subtempimpute <- sub[sub[["LABEL"]] == "L" &
                                 !is.na(sub[["INTENSITY"]]) &
                                 sub[["INTENSITY"]] == 0, ]
            subtempimpute <- subtempimpute[!is.na(subtempimpute[["ABUNDANCE"]]) &
                                           subtempimpute[["ABUNDANCE"]] != 0, ]
          }

          numFea <- xtabs(~ RUN, subtemp)
          numFeaPercentage <- 1 - numFea / length(unique(subtemp[["FEATURE"]]))
          numFeaTF <- numFeaPercentage >= 0.5
          numimpute <- xtabs(~ RUN, subtempimpute)
          sub.result <- data.frame(
            "Protein" = unique(sub[["PROTEIN"]]),
            "LogIntensities" = reformresult[["ABUNDANCE"]],
            "RUN" = reformresult[["RUN"]],
            "NumMeasuredFeature" = as.vector(numFea),
            "MissingPercentage" = as.vector(numFeaPercentage),
            "more50missing" = numFeaTF,
            "NumImputedFeature" = as.vector(numimpute))
        } else {
          subtemp <- sub[sub[["LABEL"]] == "L" & !is.na(sub[["INTENSITY"]]), ]
          numFea <- xtabs(~ RUN, subtemp)
          numFeaPercentage <- 1 - numFea / length(unique(subtemp[["FEATURE"]]))
          numFeaTF <- numFeaPercentage >= 0.5
          sub.result <- data.frame(
            "Protein" = unique(sub[["PROTEIN"]]),
            "LogIntensities" = reformresult[["ABUNDANCE"]],
            "RUN" = reformresult[["RUN"]],
            "NumMeasuredFeature" = as.vector(numFea),
            "MissingPercentage" = as.vector(numFeaPercentage),
            "more50missing" = numFeaTF)
        }
        ## result <- rbind(result, sub.result)
      }
    } else { ## single feature
      if (label) { ## label-based
        ## single feature, adjust reference feature difference
        h <- sub[sub[["LABEL"]] == "H", ]
        allmed <- median(h[["ABUNDANCE"]], na.rm=TRUE)
        for (k in 1:length(unique(h[["RUN"]]))) {
          ## ABUNDANCE is normalized
          subrun.logical <- sub[["RUN"]] == unique(h[["RUN"]])[k]
          subrun.idx <- which(subrun.logical)
          sub[subrun.idx, "ABUNDANCE"] <- sub[subrun.idx, "ABUNDANCE"] -
            sub[subrun.logical & sub[["LABEL"]] == "H", "ABUNDANCE"] + allmed
        }
        sub <- sub[sub[["LABEL"]] == "L", ]
      }
      ## single feature, use original values
      subtemp <- sub[!is.na(sub[["ABUNDANCE"]]), ]
      if (!is.null(censoredInt)) {
        if (censoredInt == "NA") {
          subtempcount <- sub[!is.na(sub[["INTENSITY"]]), ]
          subtempimpute <- sub[is.na(sub[["INTENSITY"]]), ]
          subtempimpute <- subtempimpute[!is.na(subtempimpute[["ABUNDANCE"]]), ]
        }

        if (censoredInt == "0") {
          subtempcount <- sub[!is.na(sub[["INTENSITY"]]) & sub[["INTENSITY"]] != 0, ]
          subtempimpute <- sub[!is.na(sub[["INTENSITY"]]) & sub[["INTENSITY"]] == 0, ]
          subtempimpute <- subtempimpute[!is.na(subtempimpute[["ABUNDANCE"]]) &
                                         subtempimpute[["ABUNDANCE"]] != 0, ]
        }

        numFea <- xtabs(~ RUN, subtempcount)
        numFeaPercentage <- 1 - numFea / length(unique(subtemp[["FEATURE"]]))
        numFeaTF <- numFeaPercentage >= 0.5
        numimpute <- xtabs(~ RUN, subtempimpute)
        sub.result <- data.frame(
          "Protein" = subtemp[["PROTEIN"]],
          "LogIntensities" = subtemp[["ABUNDANCE"]],
          "RUN" = subtemp[["RUN"]],
          "NumMeasuredFeature" = as.vector(numFea),
          "MissingPercentage" = as.vector(numFeaPercentage),
          "more50missing" = numFeaTF,
          "NumImputedFeature" = as.vector(numimpute))
      } else {
        subtempcount <- subtemp
        numFea <- xtabs(~ RUN, subtempcount)
        numFeaPercentage <- 1 - numFea / length(unique(subtemp[["FEATURE"]]))
        numFeaTF <- numFeaPercentage >= 0.5
        sub.result <- data.frame(
          "Protein" = subtemp[["PROTEIN"]],
          "LogIntensities" = subtemp[["ABUNDANCE"]],
          "RUN" = subtemp[["RUN"]],
          "NumMeasuredFeature" = as.vector(numFea),
          "MissingPercentage" = as.vector(numFeaPercentage),
          "more50missing" = numFeaTF)
      }
      ## result <- rbind(result, sub.result)
    }
    tukey_result[[i]] <- sub.result
  } ## loop for proteins
  close(bar)

  ## FIXME: when I put the clustering code back together, don't forget this.
  ##stopCluster(cl) # foreach autocloses
  ## Clean up the parallelized results
  ##results.list <- list()
  ##predAbundance.list <- list()
  ##for (j in 1:length(MS_results[[1]])) {
  ##  ## deal with the "results" first
  ##  results.list[[j]] <- MS_results[[1]][[j]]
  ##  predAbundance.list[[j]] <- MS_results[[2]][[j]]
  ##}

  return(tukey_result)
}

#' Use Rs Survival package to quantify the censored data.
#'
#' I have not yet carefully read this code yet, so I dunno what it does properly.
#'
#' @param modified_data the MS dataset after filtering and such.
#' @param equalFeatureVar  I think this variable is not used, it is top of my
#'   list to remove.
#' @param cutoffCensored  What criterion for censoring?
#' @param censoredInt "NA" not NA!
#' @param remove50missing  Drop proteins for which we have < 50% coverage (I
#'   think?)
#' @param MBimpute  Impute Impute damn spot!
#' @param logsum  I am not sure.
#' @param featureSubset  Keep them all!  I think this parameter is not used.
#' @param message.show  Print only the duds which this is false.
#' @param clusters  I think this is not used for this function.
#' @param original_scale  Put the data back on the original scale?
#' @param pct_adjust  Adjust the data by this factor?
#' @param ...  Extra arguments dropped into arglist for future reference.
#' @return a dataframe of quantified proteins/peptides.
perform_survival_censored <- function(modified_data,
                                      equalFeatureVar=TRUE, cutoffCensored="minFeature",
                                      censoredInt="NA", remove50missing=FALSE,
                                      MBimpute=TRUE, logsum=FALSE,
                                      featureSubset="all", message.show=FALSE, clusters=1,
                                      original_scale=FALSE, pct_adjust=0.99, ...) {
  ##data <- data[!is.na(data$ABUNDANCE),]
  data[["PROTEIN"]] <- factor(data[["PROTEIN"]])
  data[["RUN"]] <- factor(data[["RUN"]])
  if (label) {
    result <- NULL
    for(i in 1:length(unique(data[["PROTEIN"]]))) {
      sub <- data[data$PROTEIN==unique(data[["PROTEIN"]])[i], ]
      if (message.show) {
        message(paste0("Getting the summarization for censored missing values per subplot for protein ",
                       unique(sub[["PROTEIN"]]), "(",i," of ",length(unique(data[["PROTEIN"]])),")"))
      }
      sub[["FEATURE"]] <- factor(sub[["FEATURE"]])
      sub[["feature.label"]] <- paste(sub[["FEATURE"]], sub[["LABEL"]], sep="_")
      sub[["run.label"]] <- paste(sub[["RUN"]], sub[["LABEL"]], sep="_")
      ## if all measurements are NA,
      if (nrow(sub)==sum(is.na(sub[["ABUNDANCE"]]))) {
        message(paste0("Can't summarize for ",
                       unique(sub[["PROTEIN"]]), "(",
                       i," of ", length(unique(datafeature[["PROTEIN"]])),
                       ") because all measurements are NAs."))
        next()
      }
      ## remove run which has no measurement at all
      subtemp <- sub[sub$LABEL == "L" & !is.na(sub$INTENSITY), ]
      count <- aggregate(ABUNDANCE ~ RUN, data=subtemp, length)
      norun <- setdiff(unique(data$RUN), count$RUN)

      if (length(norun) != 0 & length(intersect(norun, as.character(unique(sub$RUN))))) {
        ## removed NA rows already, if there is no overlapped run, error
        sub <- sub[-which(sub$RUN %in% norun), ]
        sub$RUN <- factor(sub$RUN)
      }
      if (length(unique(sub$RUN)) == 1) {
        message(paste0("* Only 1 MS run in ",
                       levels(data$PROTEIN)[i],
                       " has measurement. Can't summarize with censored intensities."))
        next
      }

      ## remove features which are completely NAs or zero
      subtemp <- sub[sub$LABEL == "L" & !is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
      countfeature <- xtabs(~ FEATURE, subtemp)
      namefeature <- names(countfeature)[countfeature == 0]
      if (length(namefeature) != 0) {
        sub <- sub[-which(sub$FEATURE %in% namefeature), ]
        sub$FEATURE <- factor(sub$FEATURE)
      }

      ## how to decide censored or not
      ## 1. censored
      if (censoredInt == "0") {
        sub$cen <- ifelse(!is.na(sub$INTENSITY) & sub$INTENSITY == 0, 0, 1)
      }
      ## 2. all censored missing
      if (censoredInt == "NA") {
        sub$cen <- ifelse(is.na(sub$INTENSITY), 0, 1)
      }

      ## cutoffCensored
      ## 1. put minimum in protein level to NA
      ##if (cutoffCensored=="minEachProtein") {
      ##  if (censoredInt=="NA") {
      ##      cut <- min(sub$ABUNDANCE, na.rm=TRUE)
      ##      sub[is.na(sub$INTENSITY),"ABUNDANCE"] <- cut
      ##  }
      ##  if (censoredInt=="0") {
      ##      cut <- min(sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,"ABUNDANCE"])
      ##      sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"] <- cut
      ##  }
      ##}

      ## 2. put minimum in feature level to NA
      if (cutoffCensored == "minFeature") {
        if (censoredInt == "NA") {
          cut <- aggregate(ABUNDANCE ~ feature.label,
                           data=sub, function(x) min(x, na.rm=TRUE))
          ## cutoff for each Run is little less than minimum abundance in a run.
          cut$ABUNDANCE <- pct_adjust * cut$ABUNDANCE
          for (j in 1:length(unique(cut$feature.label))) {
            sub[is.na(sub$INTENSITY) &
                sub$feature.label == cut$feature.label[j], "ABUNDANCE"] <- cut$ABUNDANCE[j]
          }
        }
        if (censoredInt == "0") {
          subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
          cut <- aggregate(ABUNDANCE ~ feature.label, data=subtemptemp, FUN=min)
          ## cutoff for each Run is little less than minimum abundance in a run.
          cut$ABUNDANCE <- pct_adjust * cut$ABUNDANCE
          for (j in 1:length(unique(cut$feature.label))) {
            sub[!is.na(sub$INTENSITY) & sub$INTENSITY == 0 &
                sub$feature.label == cut$feature.label[j], "ABUNDANCE"] <- cut$ABUNDANCE[j]
          }
        }
      }
      ## 3. put minimum in RUN to NA
      if (cutoffCensored == "minRun") {
        if (censoredInt == "NA") {
          cut <- aggregate(ABUNDANCE ~ run.label, data=sub, function(x) min(x, na.rm=TRUE))
          ## cutoff for each Run is little less than minimum abundance in a run.
          cut$ABUNDANCE <- pct_adjust * cut$ABUNDANCE
          for (j in 1:length(unique(cut$run.label))) {
            sub[is.na(sub$INTENSITY) &
                sub$run.label == cut$run.label[j], "ABUNDANCE"] <- cut$ABUNDANCE[j]
          }
        }
        if (censoredInt == "0") {
          subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
          cut <- aggregate(ABUNDANCE ~ run.label, data=subtemptemp, FUN=min)
          ## cutoff for each Run is little less than minimum abundance in a run.
          cut$ABUNDANCE <- pct_adjust * cut$ABUNDANCE
          for (j in 1:length(unique(cut$run.label))) {
            sub[!is.na(sub$INTENSITY) &
                sub$INTENSITY == 0 &
                sub$run.label == cut$run.label[j],
                "ABUNDANCE"] <- cut$ABUNDANCE[j]
          }
        }
      }

      ## 20150829 : 4. put minimum RUN and FEATURE
      if (cutoffCensored == "minFeatureNRun") {
        if (censoredInt == "NA") {
          ## cutoff for each feature is little less than minimum abundance in a run.
          cut.fea <- aggregate(ABUNDANCE ~ feature.label,
                               data=sub, function(x) min(x, na.rm=TRUE))
          cut.fea$ABUNDANCE <- pct_adjust * cut.fea$ABUNDANCE
          ## cutoff for each Run is little less than minimum abundance in a run.
          cut.run <- aggregate(ABUNDANCE ~ run.label,
                               data=sub, function(x) min(x, na.rm=TRUE))
          cut.run$ABUNDANCE <- pct_adjust * cut.run$ABUNDANCE

          if (length(unique(sub$feature.label)) > 1) {
            for (j in 1:length(unique(sub$feature.label))) {
              for (k in 1:length(unique(sub$run.label))) {
                ## get smaller value for min Run and min Feature
                finalcut <- min(cut.fea$ABUNDANCE[j], cut.run$ABUNDANCE[k])
                sub[is.na(sub$INTENSITY) &
                    sub$feature.label == cut.fea$feature.label[j] &
                    sub$run.label == cut.run$run.label[k],
                    "ABUNDANCE"] <- finalcut
              }
            }
          }
          ## if single feature, not impute
        }

        if (censoredInt == "0") {
          subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
          cut.fea <- aggregate(ABUNDANCE ~ feature.label,
                               data=subtemptemp, FUN=min)
          cut.fea$ABUNDANCE <- pct_adjust * cut.fea$ABUNDANCE
          ## remove runs which has more than 50% missing values
          ## before removing, need to contribute min feature calculation
          if (remove50missing) {
            if (length(removerunid) != 0) {
              sub <- sub[-which(sub$RUN %in% removerunid), ]
              sub$RUN <- factor(sub$RUN)
            }
          }

          cut.run <- aggregate(ABUNDANCE ~ run.label,
                               data=subtemptemp, FUN=min)
          cut.run$ABUNDANCE <- pct_adjust * cut.run$ABUNDANCE
          if (length(unique(sub$feature.label)) > 1) {
            for (j in 1:length(unique(sub$feature.label))) {
              for (k in 1:length(unique(sub$run.label))) {
                ## get smaller value for min Run and min Feature
                finalcut <- min(cut.fea$ABUNDANCE[j], cut.run$ABUNDANCE[k])
                sub[!is.na(sub$INTENSITY) &
                    sub$INTENSITY == 0 &
                    sub$feature.label == cut.fea$feature.label[j] &
                    sub$run.label == cut.run$run.label[k],
                    "ABUNDANCE"] <- finalcut
              }
            }
          } else{ # single feature
            sub[!is.na(sub$INTENSITY) & sub$INTENSITY == 0, "ABUNDANCE"] <- cut.fea$ABUNDANCE
          }
        }
      }

      ## when number of measurement is less than df, error for fitting
      subtemp <- sub[!is.na(sub$ABUNDANCE), ]
      countdf <- nrow(subtemp) < (length(unique(subtemp$FEATURE)) +
                                  length(unique(subtemp$RUN)) - 1)
      ## fit the model
      if (length(unique(sub$FEATURE)) == 1) {
        ## with single feature, not converge, wrong intercept
        ## need to check
        fittest <- survival::survreg(
                               survival::Surv(
                                           ABUNDANCE,
                                           cen,
                                           type="left") ~ RUN + ref,
                               data=sub, dist="gaussian")
      } else {
        if (countdf) {
          fittest <- survival::survreg(
                                 survival::Surv(
                                             ABUNDANCE,
                                             cen,
                                             type="left") ~ RUN + ref,
                                 data=sub, dist="gaussian")
        } else {
          fittest <- survival::survreg(
                                 survival::Surv(
                                             ABUNDANCE,
                                             cen,
                                             type="left") ~ FEATURE + RUN + ref,
                                 data=sub, dist="gaussian")
        }
      }

      sub.result <- data.frame(
        "Protein" = unique(sub$PROTEIN),
        "RUN" = rep(c(levels(sub$RUN)), 1),
        "LogIntensities" = NA)
      ## get the parameters
      cf <- summary(fittest)$coefficients
      ## calculate sample quantification for all levels of sample
      a <- 1

      for (j in 1:nlevels(sub$RUN)) {
        contrast.matrix <- rep(0,nlevels(sub$RUN))
        contrast.matrix[j] <- 1
        contrast <- make.contrast.run.quantification.Survival(
          fittest, contrast.matrix, sub, labeled=TRUE)

        sub.result[a, 3] <- estimableFixedQuantificationSurvival(cf, contrast)
        a <- a + 1
      }
      result <- rbind(result, sub.result)
    }

    datamat <- reshape2::dcast(Protein ~ RUN, data=result, value.var="LogIntensities", keep=TRUE)
    datamat <- melt(datamat, id.vars=c("Protein"))
    colnames(datamat) <- c("Protein", "RUN", "LogIntensities")
    result <- datamat
  } else {
    result <- NULL
    for (i in 1:length(unique(data$PROTEIN))) {
      sub <- data[data$PROTEIN == unique(data$PROTEIN)[i], ]
      if (message.show) {
        message(paste0("Getting the summarization for censored missing values per subplot for protein ",
                       unique(sub$PROTEIN), "(",
                       i," of ", length(unique(data$PROTEIN)),")"))
      }

      sub$FEATURE <- factor(sub$FEATURE)
      ## if all measurements are NA,
      if (nrow(sub) == sum(is.na(sub$ABUNDANCE))) {
        message(paste0("Can't summarize for ",
                       unique(sub$PROTEIN), "(",
                       i," of ",length(unique(data$PROTEIN)),
                       ") because all measurements are NAs."))
        next
      }
      ## remove run which has no measurement at all
      subtemp <- sub[!is.na(sub$INTENSITY), ]
      count <- aggregate(ABUNDANCE ~ RUN, data=subtemp, length)
      norun <- setdiff(unique(data$RUN), count$RUN)
      if (length(norun) != 0 &
          length(intersect(norun, as.character(unique(sub$RUN)))) != 0) {
        ## removed NA rows already, if there is no overlapped run, error
        sub <- sub[-which(sub$RUN %in% norun), ]
        sub$RUN <- factor(sub$RUN)
      }

      if (length(unique(sub$RUN)) == 1) {
        message(paste0("* Only 1 MS run in ",
                       levels(data$PROTEIN)[i],
                       " has measurement. Can't summarize with censored intensities."))
        next
      }

      ## remove features which are (completely NAs or zero)
      subtemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
      countfeature <- xtabs(~ FEATURE, subtemp)
      namefeature <- names(countfeature)[countfeature == 0]
      if (length(namefeature) != 0) {
        sub <- sub[-which(sub$FEATURE %in% namefeature), ]
        sub$FEATURE <- factor(sub$FEATURE)
      }
      if (nrow(sub) == 0) {
        message(paste0("* All measurements are NAs or only one measurement per feature in ",
                       levels(data$PROTEIN)[i],
                       ". Can't summarize with censored intensities."))
        next
      }

      ## how to decide censored or not
      ## 1. censored
      if (censoredInt == "0") {
        sub$cen <- ifelse(!is.na(sub$INTENSITY) & sub$INTENSITY == 0, 0, 1)
      }
      ## 2. all censored missing
      if (censoredInt == "NA") {
        sub$cen <- ifelse(is.na(sub$INTENSITY), 0, 1)
      }
      ## 3. decide random above some point
      ## get runs which has all features
      ##subtemp <- sub[!is.na(sub$INTENSITY),]
      ##count <- aggregate(ABUNDANCE~RUN,data=subtemp, length)
      ##completerun <- count[count$ABUNDANCE==length(unique(sub$FEATURE)),"RUN"]
      ##if (length(completerun)!=0) {
      ##  subtemp <- sub[which(sub$RUN %in% completerun),]
      ##}else{
      ##  subtemp <- sub[which(sub$RUN %in% count[which.max(count$ABUNDANCE),"RUN"]),]
      ##}
      ## get feature mean and make order of feature
      ## mean or median?
      ##featureorder <- aggregate(ABUNDANCE~FEATURE,data=subtemp, mean)
      ##featureorder <- featureorder[with(featureorder, order(ABUNDANCE, decreasing=T)),]

      ## runs which has any missing
      ##if (length(completerun)!=0) {
      ##  incompleterun <- count[count$ABUNDANCE!=length(unique(sub$FEATURE)),"RUN"]
      ##}else{
      ##  incompleterun <- count[-which.max(count$ABUNDANCE),"RUN"]
      ##}

      ##if (length(incompleterun)!=0) {
      ##  for(j in 1:length(incompleterun)) {
      ##      temp <- sub[sub$RUN==incompleterun[j],]
      ##      temptemp <- temp[!is.na(temp$INTENSITY),]
      ##      minfeature <- temptemp[which.min(temptemp$ABUNDANCE),"FEATURE"]
      ##      abovefeature <- featureorder[1:which(featureorder$FEATURE==minfeature),"FEATURE"]
      ##      sub[which(sub$RUN==incompleterun[j] & sub$FEATURE %in% abovefeature & is.na(sub$INTENSITY)),"ABUNDANCE"] <- NA
      ##      sub[which(sub$RUN==incompleterun[j] & sub$FEATURE %in% abovefeature & is.na(sub$INTENSITY)),"cen"] <- 1
      ##  }
      ##}

      ## cutoffCensored
      ## 1. put minimum in protein level to NA
      ##if (cutoffCensored=="minEachProtein") {
      ##  if (censoredInt=="NA") {
      ##      cut <- min(sub$ABUNDANCE, na.rm=TRUE)
      ##      sub[is.na(sub$INTENSITY),"ABUNDANCE"] <- cut
      ##  }

      ##  if (censoredInt=="0") {
      ##      cut <- min(sub[!is.na(sub$INTENSITY) & sub$INTENSITY!=0,"ABUNDANCE"])
      ##      sub[!is.na(sub$INTENSITY) & sub$INTENSITY==0,"ABUNDANCE"] <- cut
      ##  }
      ##}

      ## 2. put minimum in feature level to NA
      if (cutoffCensored == "minFeature") {
        if (censoredInt == "NA") {
          cut <- aggregate(ABUNDANCE ~ FEATURE,
                           data=sub, function(x) min(x, na.rm=TRUE))
          ## cutoff for each Run is little less than minimum abundance in a run.
          cut$ABUNDANCE <- pct_adjust * cut$ABUNDANCE
          for (j in 1:length(unique(cut$FEATURE))) {
            sub[is.na(sub$INTENSITY) &
                sub$FEATURE == cut$FEATURE[j], "ABUNDANCE"] <- cut$ABUNDANCE[j]
          }
        }
        if (censoredInt == "0") {
          subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0 ,]
          cut <- aggregate(ABUNDANCE ~ FEATURE,
                           data=subtemptemp, FUN=min)
          ## cutoff for each Run is little less than minimum abundance in a run.
          cut$ABUNDANCE <- pct_adjust * cut$ABUNDANCE
          for (j in 1:length(unique(cut$FEATURE))) {
            sub[!is.na(sub$INTENSITY) &
                sub$INTENSITY == 0 &
                sub$FEATURE == cut$FEATURE[j],
                "ABUNDANCE"] <- cut$ABUNDANCE[j]
          }
        }
      }
      ## 3. put minimum in RUN to NA
      if (cutoffCensored == "minRun") {
        if (censoredInt == "NA") {
          cut <- aggregate(ABUNDANCE ~ RUN,
                           data=sub, function(x) min(x, na.rm=TRUE))
          ## cutoff for each Run is little less than minimum abundance in a run.
          cut$ABUNDANCE <- pct_adjust * cut$ABUNDANCE
          for (j in 1:length(unique(cut$RUN))) {
            sub[is.na(sub$INTENSITY) &
                sub$RUN == cut$RUN[j],
                "ABUNDANCE"] <- cut$ABUNDANCE[j]
          }
        }
        if (censoredInt == "0") {
          subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
          cut <- aggregate(ABUNDANCE ~ RUN, data=subtemptemp, FUN=min)
          ## cutoff for each Run is little less than minimum abundance in a run.
          cut$ABUNDANCE <- pct_adjust * cut$ABUNDANCE
          for (j in 1:length(unique(cut$RUN))) {
            sub[!is.na(sub$INTENSITY) &
                sub$INTENSITY == 0 &
                sub$RUN == cut$RUN[j],
                "ABUNDANCE"] <- cut$ABUNDANCE[j]
          }
        }
      }
      ## 20150829 : 4. put minimum RUN and FEATURE
      if (cutoffCensored == "minFeatureNRun") {
        if (censoredInt == "NA") {
          ## cutoff for each feature is little less than minimum abundance in a run.
          cut.fea <- aggregate(ABUNDANCE ~ FEATURE,
                               data=sub, function(x) min(x, na.rm=TRUE))
          cut.fea$ABUNDANCE <- pct_adjust * cut.fea$ABUNDANCE
          ## cutoff for each Run is little less than minimum abundance in a run.
          cut.run <- aggregate(ABUNDANCE ~ RUN,
                               data=sub, function(x) min(x, na.rm=TRUE))
          cut.run$ABUNDANCE <- pct_adjust * cut.run$ABUNDANCE
          if (length(unique(sub$FEATURE)) > 1) {
            for (j in 1:length(unique(sub$FEATURE))) {
              for (k in 1:length(unique(sub$RUN))) {
                ## get smaller value for min Run and min Feature
                finalcut <- min(cut.fea$ABUNDANCE[j], cut.run$ABUNDANCE[k])
                sub[is.na(sub$INTENSITY) &
                    sub$FEATURE == cut.fea$FEATURE[j] &
                    sub$RUN == cut.run$RUN[k], "ABUNDANCE"] <- finalcut
              }
            }
          }
          ## if single feature, not impute
        }
        if (censoredInt == "0") {
          subtemptemp <- sub[!is.na(sub$INTENSITY) & sub$INTENSITY != 0, ]
          cut.fea <- aggregate(ABUNDANCE ~ FEATURE, data=subtemptemp, FUN=min)
          cut.fea$ABUNDANCE <- pct_adjust * cut.fea$ABUNDANCE
          cut.run <- aggregate(ABUNDANCE ~ RUN, data=subtemptemp, FUN=min)
          cut.run$ABUNDANCE <- pct_adjust * cut.run$ABUNDANCE
          if (length(unique(sub$FEATURE)) > 1) {
            for (j in 1:length(unique(sub$FEATURE))) {
              for (k in 1:length(unique(sub$RUN))) {
                ## get smaller value for min Run and min Feature
                finalcut <- min(cut.fea$ABUNDANCE[j], cut.run$ABUNDANCE[k])
                sub[!is.na(sub$INTENSITY) &
                    sub$INTENSITY == 0 &
                    sub$FEATURE == cut.fea$FEATURE[j] &
                    sub$RUN == cut.run$RUN[k],
                    "ABUNDANCE"] <- finalcut
              }
            }
          } else { # single feature
            sub[!is.na(sub$INTENSITY) &
                sub$INTENSITY == 0, "ABUNDANCE"] <- cut.fea$ABUNDANCE
          }
        }
      }

      ## when number of measurement is less than df, error for fitting
      subtemp <- sub[!is.na(sub$ABUNDANCE), ]
      countdf <- nrow(subtemp) < (length(unique(subtemp$FEATURE)) +
                                  length(unique(subtemp$RUN)) - 1)
      ## fit the model
      if (length(unique(sub$FEATURE)) == 1) {
        fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type="left") ~ RUN,
                                     data=sub, dist="gaussian")
      } else {
        if (countdf) {
          fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type="left") ~ RUN,
                                       data=sub, dist="gaussian")
        } else {
          fittest <- survival::survreg(survival::Surv(ABUNDANCE, cen, type="left") ~ FEATURE + RUN,
                                       data=sub, dist="gaussian")
        }
      }

      sub.result <- data.frame(
        "Protein" = unique(sub$PROTEIN),
        "RUN" = rep(c(levels(sub$RUN)), 1),
        "LogIntensities" = NA)
      ## get the parameters
      cf <- summary(fittest)$coefficients
      ## calculate sample quantification for all levels of sample
      a <- 1
      for (j in 1:nlevels(sub$RUN)) {
        contrast.matrix <- rep(0, nlevels(sub$RUN))
        contrast.matrix[j] <- 1
        contrast <- make.contrast.run.quantification.Survival(
          fittest, contrast.matrix, sub, labeled=FALSE)
        sub.result[a, 3] <- estimableFixedQuantificationSurvival(cf, contrast)
        a <- a + 1
      }
      result <- rbind(result, sub.result)
    }

    datamat <- reshape2::dcast(Protein ~ RUN,
                               data=result, value.var="LogIntensities", keep=TRUE)
    datamat <- melt(datamat, id.vars=c("Protein"))
    colnames(datamat) <- c("Protein", "RUN", "LogIntensities")
    result <- datamat
  }
  return(result)
}
