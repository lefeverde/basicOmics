#' Split names by delim character
#'
#' @param df a data.frame
#' @param delim character to split names by
#' @param split_item index of split string to keep
#' @param remove_space should whitespace be replaced with "_" default: TRUE
#'
#' @return vector of names
#' @export
#'R
#' @examples
split_fix_names <- function(df, delim=':', split_item=1, remove_space=TRUE){
  nms <- limma::strsplit2(names(df), delim)[,split_item]
  if(remove_space){
    nms <- gsub(' ', '_',nms)
  }
  return(nms)
}

#' Turn summary object from lm into a data.frame
#'
#' @param lmlist lm summary object
#' @param coef the coefficient to extract
#'
#' @return data.frame with coefficient
#' @export
#'
#' @examples
lmlist_to_df <- function(lmlist, coef=2){
  rnms <- limma::strsplit2(names(lmlist), ' ')[,2]
  cnms <- c("Estimate","Std_Error","z_or_t_value", "p_value" )

  coef_list <- lapply(lmlist, function(x){
    x <- x$coefficients[coef,]
  })
  out_df <- data.frame(do.call(rbind, coef_list))
  row.names(out_df) <- rnms
  colnames(out_df) <- cnms
  return(out_df)

}

