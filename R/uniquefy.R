#' Return a unique data.frame by max value
#'
#' @param mat data.frame with 1st column as string and second as numeric
#' @param group_col string denoting group column
#' @param val_col string denoting value column
#'
#' @return unique data.frame
#' @export
#'
#' @examples
#' set.seed(42)
#' mat <- data.frame(group=c('a', 'a', 'a', 'b', 'c', 'c'), vals=runif(6))
#' mat$vals[6] <- -4
uniquefy_by_abs_max <- function(mat, group_col=NULL, val_col=NULL){
  #TODO make this workhorse functions which
  # and split the variance and other uniquefiers
  # into compatible functions
  if(is.null(group_col)){
    group_col <- names(mat)[1]
  }
  if(is.null(val_col)){
    val_col <- names(mat)[2]
  }
  mat <- data.frame(mat)
  mat <- mat[order(mat[group_col], -abs(mat[val_col])),]
  mat <- mat[!duplicated(mat[group_col]),]
}

#### TODO Functions ####

#' Uniquefy by variance
#'
#' Idea for using max var to de-duplicate probes came
#' form here: https://www.biostars.org/p/51756/#51875
#' Probe set needs to be a data frame with one column
#' the rownames are the probe ID and the single column is
#' refseq or equivalent gene annotation
#'
#' @param expr_mat
#' @param probe_set
#'
#' @return
#' @keywords internal
#'
#' @examples
uniquefy_by_variance <- function(expr_mat, probe_set){
  # TODO Refactor to follow uniquefy_by_abs_max syntax
  names(probe_set) <- 'annots'

  # filtering so annots and expr have same
  # probes in both
  if(class(expr_mat)[[1]] == 'ExpressionSet'){
    expr_mat <- expr_mat@assayData$exprs

  } else {
    expr_mat <- as.matrix(expr_mat)
  }

  expr_mat <- expr_mat[row.names(expr_mat) %in% row.names(probe_set),]
  probe_set <- probe_set[row.names(probe_set) %in% row.names(expr_mat),,drop=FALSE]

  ordered_probes <- probe_set[order(probe_set$annots),,drop=FALSE]
  dup_annots <- ordered_probes$annots[duplicated(ordered_probes$annots)]
  dup_probes <- ordered_probes[ordered_probes$annots %in% dup_annots,,drop=FALSE]

  unique_expr <- expr_mat[! row.names(expr_mat) %in% row.names(dup_probes),]
  dup_expr <- expr_mat[row.names(expr_mat) %in% row.names(dup_probes),]
  gene_annot <- dup_probes[match(row.names(dup_expr), row.names(dup_probes)),, drop=FALSE]
  rv <- data.frame(prob_id = row.names(gene_annot), annots=gene_annot, prob_variance=matrixStats::rowVars(as.matrix(dup_expr)))
  max_var <- rv[order(rv$annots, -rv$prob_variance),]
  max_var <- max_var[!duplicated(max_var$annots),]
  probes_to_keep <- c(row.names(max_var), row.names(unique_expr))

  fixed_expr_mat <- expr_mat[row.names(expr_mat) %in% probes_to_keep,]
  return(fixed_expr_mat)

}
