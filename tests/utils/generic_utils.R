library(data.table)
library(xlsx)

rename_column <- function(df, old_names, new_names){
  setnames(df, old = old_names, new = new_names)
  return(df)
}

merge_dataframes <- function(df1, df2, col_names){
  return(merge(df1, df2, by = setdiff(colnames(df1), col_names), 
               all.x = T, all.y = T))
}

generate_xlsx <- function(dataframes, file_name){
  write.xlsx(dataframes[[1]], file=file_name, sheetName=dataframes[[1]]$notes, row.names=FALSE)
  for (each_dataframe in dataframes[-1] ){
    write.xlsx(each_dataframe, file=file_name, sheetName=each_dataframe$notes, append=TRUE, row.names=FALSE)
  }
}