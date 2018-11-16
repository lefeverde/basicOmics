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

  rnms <- try(limma::strsplit2(names(lmlist), ' ')[,2], silent = T) # gets rid of extranous name
  if(class(rnms) == 'try-error'){
    rnms <- names(lmlist)
  }
  cnms <- c("Estimate","Std_Error","z_or_t_value", "p_value" )

  coef_list <- lapply(lmlist, function(x){
    x <- x$coefficients[coef,]
  })
  out_df <- data.frame(do.call(rbind, coef_list))
  row.names(out_df) <- rnms
  colnames(out_df) <- cnms
  return(out_df)

}

#' Return a unique data.frame by max value
#'
#' @param mat data.frame with 1st column as string and second as numeric
#'
#' @return unique data.frame
#' @export
#'
#' @examples
#' set.seed(42)
#' mat <- data.frame(group=c('a', 'a', 'a', 'b', 'c', 'c'), vals=runif(6))
#' mat$vals[6] <- -4
uniquefy_by_abs_max <- function(mat){
  names(mat) <- c('group', 'vals')
  mat <- mat[order(mat$group, -abs(mat$vals)),]
  mat <- mat[!duplicated(mat$group),]
  return(mat)
}

#' Title
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
#' @export
#'
#' @examples
uniquefy_by_variance <- function(expr_mat, probe_set){

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

get_deseq2_results <- function(dds_object, gene_annots){
  results_list <- list()
  res_names <- resultsNames(dds_object)
  for(i in  res_names[2:length(res_names)]){
    fixed_name <- limma::strsplit2(i, '_')[,2]
    #results_list[[i]] <-results(dds_object, name = i)
    temp_res <- data.frame(coefficient=i,results(dds_sub, name=i, tidy = T))
    results_list[[i]] <- temp_res
    #temp_res <- lfcShrink(dds_object, )
  }
  names(results_list) <- NULL
  res_df <- do.call(rbind,results_list)
  rn <- res_df$row
  res_df$row <- NULL
  sorted_annots <- gene_annots[match(rn, row.names(gene_annots)),]
  out_df <- data.frame(ENSEMBL_ID=rn, sorted_annots, res_df)
  row.names(out_df) <- NULL
  return(out_df)
}

uniquefy_results_by_lfc <- function(res_df){

  #res_df <- res_df[row.names(res_df) %in% row.names(res_df),]
  res_df <- res_df[order(res_df$external_gene_name, -res_df$log2FoldChange),]
  res_df <- res_df[!duplicated(res_df$external_gene_name),]
  return(res_df)
}


# get_deseq2_results <- function(dds_object, gene_annots){
#   results_list <- list()
#   res_names <- DESeq2::resultsNames(dds_object)
#   for(i in  res_names[2:length(res_names)]){
#     fixed_name <- limma::strsplit2(i, '_')[,2]
#     #results_list[[i]] <-results(dds_object, name = i)
#     temp_res <- data.frame(coefficient=i,DESeq2::results(dds_object, name=i, tidy = T))
#     results_list[[i]] <- temp_res
#     #temp_res <- lfcShrink(dds_object, )
#   }
#   names(results_list) <- NULL
#
#   res_df <- do.call(rbind,results_list)
#   rn <- res_df$row
#   res_df$row <- NULL
#   sorted_annots <- gene_annots[match(rn, row.names(gene_annots)),]
#   out_df <- data.frame(ENSEMBL_ID=rn, sorted_annots, res_df)
#   return(out_df)
#
# }


#' Title
#'
#' Takes a limma fit object and returns a data
#' frame with each each coefficient summarized
#' @param limma_fit
#' @param skip_intercept
#'
#' @return
#' @export
#'
#' @examples
get_limma_results <- function(limma_fit, coefs=NULL){
  all_coefs <- topTable(limma_fit, number = Inf)
  gene_order <- row.names(all_coefs)
  if(is.null(coefs)){
    coefs <- colnames(limma_fit$design)
    coefs <- coefs[2:length(coefs)]
  }


  fc_and_fdr <- list()
  fc_and_fdr[['gene_id']] <- gene_order
  for(i in coefs){
    cur_top <- topTable(limma_fit,
                        coef = i,
                        confint = TRUE,
                        number = Inf,
                        adjust='fdr')
    cur_top <- cur_top[gene_order,]
    fc_and_fdr[[(paste0(i, '_pval'))]] <- cur_top$P.Value
    fc_and_fdr[[(paste0(i, '_fdr_pval'))]] <- cur_top$adj.P.Val
    fc_and_fdr[[(paste0(i, '_log2_fold_change'))]] <- cur_top$logFC
  }
  out_df <- as.data.frame(fc_and_fdr,
                          check.names=F,
                          row.names = fc_and_fdr[['gene_id']])
  out_df$gene_id <- NULL
  gene_annot_cols <- !sapply(all_coefs, is.numeric)
  if(sum(gene_annot_cols) >= 1){
    cur_annots <- all_coefs[,gene_annot_cols, drop=F]
    cur_annots <- cur_annots[row.names(out_df),]
    out_df <- cbind(cur_annots, out_df)
  }
  return(out_df)
}


#' Gene expression filter by quantile
#' This function first removes rows in which
#' the max value equals the min value and
#' then filters the resulting genes by row
#' medians. Rows with a median greater than
#' the resulting value of the quantile are
#' retained.
#'
#' @param expr_mat an expression matrix of class matrix
#' @param quant_val quantile threshold
#'
#' @return filtered by quantile  expresion matrix
#' @export
#'
#' @examples
gene_quantile_filter <- function(expr_mat, quant_val=.05){
  rm <- Biobase::rowMax(expr_mat)
  qthresh <- quantile(rm, quant_val)
  filt_mat <- expr_mat[rm > qthresh,]
  filt_mat <- filt_mat[Biobase::rowMedians(filt_mat) > qthresh,]
  return(filt_mat)
}
