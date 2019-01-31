#' Split a named vector of DEGs by fold change
#'
#' @param lfc_list DEG vector of FCs
#'
#' @return
#' @export
#'
#' @examples
split_by_reg <- function(lfc_list){
  split_degs <- list(
    sort(lfc_list[lfc_list > 0, drop=F], decreasing = TRUE),
    sort(lfc_list[lfc_list < 0, drop=F]))
  names(split_degs) <- c('upreg', 'downreg')
  #cat(paste0(names(split_degs[[1]]), ' '))
  #cat('\n',paste0(names(split_degs[[1]]), ' '))
  return(split_degs)
}

#' rbinds a named list of data.frames
#'
#' rbinds together data.frames in which the list name
#' is turned into a column and bound to the data.frame
#' becoming the first column of the data.frame
#'
#' @param df_list named list of data.frames
#' @param sum_col_name colname of new column bound to data.frame
#'
#' @return a single data.frame row bound data.frame
#' @export
#'
#' @examples
rbind_named_df_list <- function(df_list, sum_col_name='col_from_list'){
  nms <- names(df_list)
  if(is.null(nms)){
    stop('df_list needs to be a named list')
  }
  out_df <- list()
  empty_dfs <- NULL
  for(x in nms){
    temp_pth <- data.frame(df_list[[x]])
    if(nrow(temp_pth) > 0){
      temp_pth <- cbind(x, temp_pth)
      row.names(temp_pth) <- NULL
      names(temp_pth)[1] <- sum_col_name
      out_df[[x]] <- temp_pth
    } else {
      empty_dfs <- c(empty_dfs, x)
      next
    }
    if(!is.null(empty_dfs)){
      warning(paste('Some items contained empty or otherwise disagreeable data.frames and were not included in the bound data.frame.
                    These include:\n'),
              paste0(head(empty_dfs), sep='\n'))
    }

  }
  out_df <- do.call(rbind, out_df)
  row.names(out_df) <- NULL
  return(out_df)

  }


#' Extract genes sets from a pathway results data.frame
#'
#' pathway results data.frame with a column
#' of gene sets seperated by a delimiter
#'
#' @param path_df data.frame with gene sets
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
#' I made this function because the way the
#' model.matrix creates issues with differential
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




#' Split names by delim character
#'
#' @param df a data.frame
#' @param delim character to split names by
#' @param split_item index of split string to keep
#' @param remove_space should whitespace be replaced with "_" default: TRUE
#'
#' @return vector of names
#' @export
#'
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
uniquefy_by_abs_max <- function(mat, group_col=NULL, val_col=NULL){
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
#' @export
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

#' Returns DESeq2 results in long format
#'
#' @param dds_object DDS object after running DESeq
#' @param use_shrinkage logical whether to perform lfcshrinkage
#' @param gene_annots Gene annotations which if not given are assumed
#'
#' @return results in long format
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
    #temp_res <- data.frame(coefficient=i, results(dds_object, name=i, tidy = T))
    #DESeq2::lfcShrink(dds_object, coef = i, type = 'normal', parallel=TRUE) %>%
    #temp_res <- lfcShrink(dds_object, )
  }
  names(results_list) <- NULL
  res_df <- do.call(rbind,results_list)
  rn <- res_df$row
  res_df$row <- NULL
  sorted_annots <- gene_annots[match(rn, row.names(gene_annots)),,drop=FALSE]
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

#' get results from limma fit object
#'
#' Takes a limma fit object and returns a long data.frame
#' of the results. If no coefficients are passed, it tries
#' to guess which are the relevant ones checking which
#' are made up of 1's and 0's.
#'
#' @param limma_fit
#' @param skip_intercept
#'
#' @return
#' @export
#'
#' @examples
get_limma_results <- function(limma_fit, coefs=NULL){

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
    cur_res <- limma::topTable(limma_fit, coef = x, number = Inf) %>%
      tibble::as_tibble(.) %>%
      # tibble::rownames_to_column(., 'ensembl_gene_id') %>%
      # data.frame(coefficient=x, .)
      data.frame(coefficient=x, .) %>%
      tibble::rownames_to_column(., 'ensembl_gene_id')

  }) %>% do.call(rbind, .)
  return(res_list)
}
#' Selects the top n genes by variance
#'
#' @param expr_mat expression mat
#' @param n_genes number of genes to extract
#'
#' @return
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

#' Title
#'
#' @param n_genes
#' @param m_samples
#' @param seed
#' @param groups
#'
#' @return
#' @export
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
