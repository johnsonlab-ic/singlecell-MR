
#' Check_colnames
#'
#' This function checks the column names. Used in cellCOLOC.
#'
#' @export

check_colnames=function(df,colnames){
    
    df_colnames=colnames(df)
    input_colnames=colnames
    
    #check if all input colnames exist in df colnames
    missing=setdiff(input_colnames,df_colnames)
    
    if(length(missing)==0){
        message("All column names present.")
    } else {
        stop(paste0("Columns missing from GWAS: ",missing,". Please check formatting."))
    }
    
}




