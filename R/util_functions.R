#' Extract genes sets from a pathway results data.frame
#'
#' pathway results data.frame with a column
#' of gene sets seperated by a delimiter
#'
#' @param path_df data.frame with gene sets
#' @param delim delimeter seperating genes (default: '/')
#'
#' @return list with names as pathway ids
#' @export
#'
#' @examples
get_geneSets <- function(path_df, delim="/"){
  nms <- path_df$ID
  gene_sets <-  strsplit(path_df$geneID, delim, fixed=TRUE)
  names(gene_sets) <- nms
  return(gene_sets)
}


#' Creates a much nicer formatted design matrix than the default model.matrix
#'
#' I made this function because the way the model.matrix function
#' creates a desing matrix, can cause issues with differential
#' gene expression tools, such as edgeR, DESeq2,
#' and limma. If the intercept column is included,
#' it first removes the parens and parens from the
#' intercept column and then removing the variable
#' name from the colnames and replaces them with
#' factor level vs reference level syntax. If the
#' design is specified without an intercept, variable
#' concatentation is simply removed.
#'
#' @param ... formula with outcome as first variable
#' @param data data.frame
#' @param show_warnings logical whether to show wanrings
#'
#' @return model.matrix
#' @export
#'
#' @examples
better_model_matrix <- function(..., data, show_warnings=TRUE){
  m <- model.matrix(..., data=data)

  # Getting variables (columns) from formula
  mod_str <- deparse(...)

  var_str <-
    limma::strsplit2(mod_str, '\\+') %>%
    drop(.) %>%
    trimws(.) %>%
    gsub('~', '', .) #%>%


  # cleaning model matrix
  if(var_str[1] == "0"){
    # checks if intercept is 1st term
    outcome_var <- var_str[2]
  } else {
    outcome_var <- var_str[1]
  }
  if(show_warnings){
    w <- paste0('\n', outcome_var, ' used as outcome variable.\nIf this is incorrect, reorder the formula so that the outcome variable comes first.'   )
    warning(w)
  }
  df_str <- deparse(substitute(data))
  colnames(m) <-
    gsub(outcome_var, '', colnames(m)) %>%
    # gsub(df_str, '', .) %>%
    gsub('\\(', '', .) %>%
    gsub('\\)', '', .)

  model_has_intercept <-
    colnames(m) %>%
    grepl('Intercept', ., ignore.case = TRUE) %>%
    any(.)

  if(model_has_intercept){
    colnames_to_prettify <- levels(data[,outcome_var])
    ref_level <- colnames_to_prettify[1]
    colnames_to_prettify <- colnames_to_prettify[-1]
    stopifnot(all(colnames_to_prettify %in% colnames(m)))
    pretty_colnames <- sapply(colnames_to_prettify, function(x){
      paste0(x, '_vs_', ref_level)
    })
    col_idxs <- match(colnames_to_prettify, colnames(m))
    colnames(m)[col_idxs] <- pretty_colnames

  }
  return(m)
}




#' get results from limma fit object
#'
#' Takes a limma fit object created using either \link[limma]{eBayes} or
#' \link[limma]{treat} and returns a long data.frame of the results.
#' If the `fit` object contains `lods`, then results are obtained using
#' \link[limma]{topTable}. Else, if the fit object contains a `treat.lfc`
#' topTreat is used. If no coefficients are passed, it tries to guess
#' which are the relevant ones checking which are made up of 1's and 0's.
#'
#' @param limma_fit fit object produced by limma
#' @param coefs coefficients to return
#' @param print_summary logical, whether print summary of the results
#'
#' @return data.frame of results in long(ish) format
#' @export
#'
#' @importFrom limma topTable
#'
#'
#' @examples
get_limma_results <- function(limma_fit, coefs=NULL, print_summary=TRUE){
  #TODO add tests
  #TODO think about handling levels with whitespace
  #TODO maybe automagically figure type of gene annots
  d <- limma_fit$design
  if(is.null(coefs)){
    d <- limma_fit$design
    one_hot_cols <- sapply(seq(1,ncol(d)), function(i){
      all(d[,i] %in% c(1,0))
    })
    coefs <- colnames(d)[one_hot_cols]
    coefs <- coefs[!grepl('intercept', coefs, ignore.case = T)]
  }

  res_list <- lapply(coefs, function(x){
    if('lods' %in% names(limma_fit)){
      cur_res <-
        topTable(limma_fit, coef = x, number = Inf, confint=TRUE)
    }
    if('treat.lfc' %in% names(limma_fit)){
      cur_res <-
        topTreat(limma_fit, coef = x, number = Inf, confint=TRUE)
    }
      cur_res <- cur_res %>%
        tibble::rownames_to_column(., 'gene_id') %>%
        cbind(coefficient=x, .)
  }) %>% do.call(rbind, .)

  if(print_summary){
    expression_summaries(res_list) %>%
      print
  }
  return(res_list)
}





#' Returns DESeq2 results in long format
#'
#' @param dds_object DDS object after running DESeq
#' @param use_shrinkage logical whether to perform lfcshrinkage
#' @param gene_annots Gene annotations which if not given are assumed
#'
#' @return data.frame of results in long(ish) format
#' @export
#'
#' @examples
get_deseq2_results <- function(dds_object, gene_annots=NULL, use_shrinkage=FALSE){
  if(is.null(gene_annots)){
    gene_annots <-
      SummarizedExperiment::mcols(dds_object) %>%
      data.frame(.) %>%
      dplyr::select_if(., is.character)

    warning_msg <- paste0('gene annots not given.\nUsing: ',
                          paste(names(gene_annots), collapse = ', '),
                          ' for gene annotations')
    warning(warning_msg)

  }

  # Get results by coefs which
  # have `_vs_` in them
  results_list <- list()
  res_names <- DESeq2::resultsNames(dds_object)
  res_names <- res_names[grepl('_vs_', res_names)]

  for(i in  res_names){
    # Checks whether to get results via shrinkage
    if(use_shrinkage){
      temp_res <- DESeq2::lfcShrink(dds_object,
                                    coef = i,
                                    type = 'normal',
                                    parallel=TRUE)
    } else if(!use_shrinkage){
      temp_res <- DESeq2::results(dds_object, name = i)
    }
    temp_res <-
      temp_res %>%
        data.frame(.) %>%
        tibble::rownames_to_column(., var = "ensembl_gene_id") %>%
        data.frame(coefficient=i, .)

    results_list[[i]] <- temp_res
  }
  names(results_list) <- NULL
  res_df <- do.call(rbind,results_list)
  rn <- res_df$ensembl_gene_id
  res_df$row <- NULL
  sorted_annots <- gene_annots[match(rn, row.names(gene_annots)),,drop=FALSE]
  out_df <- data.frame(ensembl_gene_id=rn, sorted_annots, res_df)
  row.names(out_df) <- NULL
  return(out_df)
}



#' Selects the top n genes by variance
#'
#' @param expr_mat expression mat
#' @param n_genes number of genes to extract
#'
#' @return gene matrix
#' @export
#'
#' @examples
select_genes_by_variance <- function(expr_mat, n_genes=5000){
  rv_idx <- matrixStats::rowVars(expr_mat) %>%
    order(., decreasing = TRUE)
  rv_idx <- rv_idx[1:n_genes]
  filt_expr_mat <- expr_mat[rv_idx,]
  return(filt_expr_mat)
}


#' Gene expression filter by quantile
#'
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
#' @return filtered by quantile expresion matrix
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

#### TODO Functions ####

#' Split names by delim character
#'
#' @param df a data.frame
#' @param delim character to split names by
#' @param split_item index of split string to keep
#' @param remove_space should whitespace be replaced with "_" default: TRUE
#'
#' @return vector of names
#' @keywords internal
#'
#' @examples
split_fix_names <- function(df, delim=':', split_item=1, remove_space=TRUE){
  nms <- limma::strsplit2(names(df), delim)[,split_item]
  if(remove_space){
    nms <- gsub(' ', '_',nms)
  }
  return(nms)
}



#' Title
#'
#' @param n_genes
#' @param m_samples
#' @param seed
#' @param groups
#'
#' @return
#' @keywords internal
#'
#' @examples
create_test_dds <- function(n_genes=500, m_samples=60, seed=42, groups=NULL){
  set.seed(seed)
  dds <- DESeq2::makeExampleDESeqDataSet(n=n_genes, m = m_samples)
  mcols(dds) <- NULL
  mcols(dds) <- data.frame(gene_id=names(dds),
                           row.names=names(dds))
  if(is.null(groups)){
    groups <-
      c(rep('A', m_samples/3),
        rep('B', m_samples/3),
        rep('C', m_samples/3)) %>%
      factor(.)
  }
  colData(dds) <- DataFrame(condition=groups,
                            row.names = colnames(dds))

  return(dds)

}

#' Turn summary object from lm into a data.frame
#'
#' @param lmlist lm summary object
#' @param coef the coefficient to extract
#'
#' @return data.frame with coefficient
#' @keywords internal
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

#' Gets lfc from a pathway res df
#'
#' This returns the lfc from a genelist as a
#' named vector of genes by lfc
#'
#' @param gene_list list of gene vectors
#' @param res_df long res df table
#'
#' @return
#' @keywords internal
#'
#' @examples
get_lfc <- function(gene_list, res_df){
  ol <- lapply(gene_list, function(x){
    lfc <- res_df$logFC[match(x, res_df$external_gene_name)] %>%
      setNames(., x)
  })
  names(ol) <- names(gene_list)
  return(ol)
}


#' Split a named vector of DEGs by fold change
#'
#' @param lfc_list DEG vector of FCs
#'
#' @return
#' @keywords internal
#'
#' @examples
split_by_reg <- function(lfc_list){
  split_degs <- list(
    sort(lfc_list[lfc_list > 0, drop=F], decreasing = TRUE),
    sort(lfc_list[lfc_list < 0, drop=F]))
  names(split_degs) <- c('upreg', 'downreg')
  return(split_degs)
}


