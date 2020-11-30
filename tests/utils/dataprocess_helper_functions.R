invoke_dataprocess<-function(input_file, summary_method="linear", mb_impute=FALSE,
                             censored_int="0"){
  # parameterized data process function
  
  #against stable version
  output_stable = MSstats::dataProcess(as.data.frame(unclass(input_file)), 
                                       summaryMethod=summary_method, 
                                       MBimpute=mb_impute, 
                                       censoredInt=censored_int)
  #against dev version
  output_dev = MSstatsdev::dataProcess(input_file, 
                                       summaryMethod=summary_method, 
                                       MBimpute=mb_impute, 
                                       censoredInt=censored_int)
  return(list(stable=output_stable, dev=output_dev))
}

compare_values_processed_data <- function(input_df){
  #couldn't make this logic work. Saving for later
  #this logic is overkilled for now. Need to handle columns more flexibly
  input_df$match.orgint <- input_df$INTENSITY.x == input_df$INTENSITY.y
  input_df$match.abn <- input_df$ABUNDANCE.x == input_df$ABUNDANCE.y
  input_df$match.censored <- input_df$censored.x == input_df$censored.y
  input_df$match.ftrslct <- input_df$feature_quality.x == input_df$feature_quality.y
  input_df$match.otr <- input_df$is_outlier.x == input_df$is_outlier.y
  return(input_df)
}

compare_values_runlevel_data <- function(input_df){
  #this logic is overkilled for now. Need to handle columns more flexibly
  input_df$match.int <- input_df$LogIntensities.x == input_df$LogIntensities.y
  input_df$match.perc <- input_df$MissingPercentage.x == input_df$MissingPercentage.y
  input_df$match.numimpf <- input_df$NumImputedFeature.x == input_df$NumImputedFeature.y
  input_df$match.nummsrf <- input_df$NumMeasuredFeature.x == input_df$NumMeasuredFeature.y
  return(input_df)
}

handle_na_values_processed_data <- function(input_df ){
  ## if both NA, they also match.
  input_df[is.na(
    input_df$INTENSITY.x) & is.na(input_df$INTENSITY.y), 'match.orgint'] <- TRUE
  input_df[is.na(
    input_df$ABUNDANCE.x) & is.na(input_df$ABUNDANCE.y), 'match.abn'] <- TRUE
  input_df[is.na(
    input_df$censored.x) & is.na(input_df$censored.y), 'match.censored'] <- TRUE
  input_df[is.na(
    input_df$feature_quality.x) & is.na(input_df$feature_quality.y), 'match.ftrslct'] <- TRUE
  input_df[is.na(
    input_df$is_outlier.x) & is.na(input_df$is_outlier.y), 'match.otr'] <- TRUE
  
  ## others, should be FALSE : if one of them NA, match column NA -> change to FALSE
  input_df[is.na(input_df$match.orgint), 'match.orgint'] <- FALSE
  input_df[is.na(input_df$match.abn), 'match.abn'] <- FALSE
  input_df[is.na(input_df$match.censored), 'match.censored'] <- FALSE
  input_df[is.na(input_df$match.ftrslct), 'match.ftrslct'] <- FALSE
  input_df[is.na(input_df$match.otr), 'match.otr'] <- FALSE
  return(input_df)
}

handle_na_values_runlevel_data <- function(input_df){
  input_df[is.na(input_df$LogIntensities.x) & is.na(
    input_df$LogIntensities.y), 'match.int'] <- TRUE
  input_df[is.na(input_df$MissingPercentage.x) & is.na(
    input_df$MissingPercentage.y), 'match.perc'] <- TRUE
  input_df[is.na(input_df$NumImputedFeature.x) & is.na(
    input_df$NumImputedFeature.y), 'match.numimpf'] <- TRUE
  input_df[is.na(input_df$NumMeasuredFeature.x) & is.na(
    input_df$NumMeasuredFeature.y), 'match.nummsrf'] <- TRUE
  
  ## others, should be FALSE : if one of them NA, match column NA -> change to FALSE
  input_df[is.na(input_df$match.int), 'match.int'] <- FALSE
  input_df[is.na(input_df$match.perc), 'match.perc'] <- FALSE
  input_df[is.na(input_df$match.numimpf), 'match.numimpf'] <- FALSE
  input_df[is.na(input_df$match.nummsrf), 'match.nummsrf'] <- FALSE
  return(input_df)
}