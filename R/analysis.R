#' get overlap ratio of 2 gene sets
#'
#' Returns the overlap divded by the length of the smaller set
#'
#' @param x gene vector
#' @param y gene vector
#'
#' @return ratio
#' @noRd
#'
#'
overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}


#' Creates matrix of overlap ratios
#'
#' This function creates a matrix of the ratio
#' of overlapping genesets. It takes a list of genes
#' or alternatively 2 gene lists
#'
#' @param l1 gene list
#' @param l2 optional 2nd gene list
#'
#' @return top diagonal matrix of overlap ratios
#' @export
#'
#' @examples
#' \dontrun{
#' # returns a self (l1 vs l1) overlap matrix
#' genes <- paste0('gene', seq(1, 100))
#' set.seed(42)
#' n_genes_to_sample <- rnbinom(10, 10, .5)
#' l1 <- lapply(seq_along(n_genes_to_sample), function(i){
#'  n <- n_genes_to_sample[i]
#'  sample(genes, n)
#' }) %>%
#'  setNames(., paste0('L1_', LETTERS[1:10]))
#' omat <- get_overlap_matix(l1)
#' }
#'
get_overlap_matix <- function(l1, l2=NULL){
  if(is.null(l2)){
    l2 <- l1
  }

  l1_nms <- names(l1)
  l2_nms <- names(l2)


  w <- matrix(NA, nrow=length(l1_nms), ncol=length(l2_nms))
  rownames(w) <- l1_nms
  colnames(w) <- l2_nms

  for (i in 1:nrow(w)) {
    for (j in 1:ncol(w)) {
      w[i,j] = overlap_ratio(l1[[i]], l2[[j]])
    }
  }
  return(w)
}


#' Returns summaries of gene expression
#'
#'
#' This just summarizes the number of significant up & down
#' regulated genes per group
#'
#' @param long_res_table res df in long format
#' @param log2_thresh threshold to genes by lfc
#' @param logFC_col str of logFC_col
#' @param padj_col str of padj_col
#'
#' @return table summary
#' @export
#'
#' @examples
expression_summaries <- function(long_res_table, log2_thresh=0, logFC_col = NULL, padj_col=NULL){
  #TODO fix hardcoding of coefficient
  if(is.null(logFC_col)){
    potential_cols <- c('log2FoldChange','logFC')
    logFC_col <- na.omit(match(potential_cols, colnames(long_res_table)) )
    colnames(long_res_table)[logFC_col] <- 'log2FoldChange'
  }
  if(is.null(padj_col)){
    potential_cols <- c('padj','adj.P.Val')
    logFC_col <- na.omit(match(potential_cols, colnames(long_res_table)) )
    colnames(long_res_table)[logFC_col] <- 'padj'
  }
  group_sums <-
    long_res_table %>%
    dplyr::filter(padj < .05) %>%
    dplyr::filter(.,abs(log2FoldChange) >= log2_thresh) %>%
    group_by(., coefficient) %>%
    summarise(.,
              up=sum(log2FoldChange > 0),
              down=sum(log2FoldChange < 0))
  return(group_sums)
}
