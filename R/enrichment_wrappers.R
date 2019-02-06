#' Over Representation analysis (ORA)
#'
#' ORA wrapper which splits the long res list into
#' combined, up-regulated, and down-regulated list of
#' genes. It then performs ORA using clusterProfiler
#'
#' @param res_object res object in long(ish) format
#' @param lfc_thresh lfc threshold  absolute value
#'
#' @return data.frame of enriched pathways
#' @export
#'
#' @examples
enrich_wrapper <- function(res_object, lfc_thresh=0){
  res_object2 <- res_object[abs(res_object$logFC) >= lfc_thresh,] %>%
    dplyr::filter(., adj.P.Val < .05)
  temp_reslist <- list(
    combined=res_object2,
    upreg=res_object2[res_object2$logFC > 0,],
    downreg=res_object2[res_object2$logFC < 0,]
  )
  for(i in temp_reslist){
    print( paste(lfc_thresh,paste0(dim(i), collapse = ' ')))
  }
  ol <- lapply(temp_reslist, function(temp_res){
    temp_list <-  list(
      try(DOSE::enrichDGN(temp_res$entrezgene,  readable = TRUE, pvalueCutoff = .1, qvalueCutoff = .1)),
      try(DOSE::enrichDO(temp_res$entrezgene, readable = TRUE, pvalueCutoff = .1, qvalueCutoff = .1)),
      try(ReactomePA::enrichPathway(temp_res$entrezgene,readable = TRUE, pvalueCutoff = .1, qvalueCutoff = .1 )),
      try(clusterProfiler::enrichKEGG(temp_res$entrezgene, pvalueCutoff = .1, qvalueCutoff = .1))) %>% setNames(., c('DGN', 'DO','Reactome', 'KEGG' )) %>%
      lapply(., data.frame)
  }) %>% lapply(., function(x){
    rbind_named_df_list(x, 'database')
  })# %>%
    #rbind_named_df_list(., 'regulation')
  return(ol)
}

#### TODO Functions ####

#' Converts ensembl to entrez
#'
#' @param gene_vector vector of human ensembl genes
#'
#' @return data.frame of with ensembl gene and entrez gene ids
#' @keywords internal
#'
#' @examples
ensembl_to_entrez_genes <- function(gene_vector){
  # if(!exists(mart)){
  #   mart <- biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  # }
  atts <- c('ensembl_gene_id', 'entrezgene')
  filts <- c('ensembl_gene_id')
  entrez_genes <- biomaRt::getBM(atts, filts, gene_vector, mart)
  return(entrez_genes)
}


#' Performs enrichment analysis accross all GO domains
#'
#' @param ent_genes vector of entrez genes
#'
#' @return pathway list by database
#' @keywords internal
#'
#'
#' @examples
enrich_go_wrapper <- function(ent_genes){
  l <- list(enrichGO(ent_genes, ont = 'CC', readable=TRUE, OrgDb=org.Hs.eg.db), enrichGO(ent_genes, ont = 'MF', readable=TRUE, OrgDb=org.Hs.eg.db),enrichGO(ent_genes, ont = 'BP', readable=TRUE, OrgDb=org.Hs.eg.db))

  # ent_genes <- bitr(gene_vector, 'SYMBOL', toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
  dbs <- c( 'GO_CC', 'GO_MF', 'GO_BP')
  names(l) <- dbs
  l <- lapply(l, as.data.frame)
  l <- do.call(rbind, l)
  l$database <- limma::strsplit2(row.names(l), '\\.')[,1]
  row.names(l) <- NULL
  return(l)
}

#' Performs enrichment analysis using KEGG, REACTOME, DO, and DisGenNet
#'
#' @param ent_genes vector of entrez genes
#'
#' @return pathway list by database
#' @keywords internal
#'
#' @examples
enrich_pathway_wrapper <- function(ent_genes){

  # ent_genes <- bitr(gene_vector, 'SYMBOL', toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
  dbs <- c('KEGG', 'REACTOME', 'DO', 'DisGenNet')
  # l <- list(enrichKEGG(ent_genes, pvalueCutoff = 1), enrichPathway(ent_genes, pvalueCutoff = 1), enrichDO(ent_genes), enrichDGN(ent_genes))
  l <- list(clusterProfiler::enrichKEGG(ent_genes),  ReactomePA::enrichPathway(ent_genes), DOSE::enrichDO(ent_genes), DOSE::enrichDGN(ent_genes))
  names(l) <- dbs
  l <- lapply(l, as.data.frame)
  l <- do.call(rbind, l)
  l$database <- limma::strsplit2(row.names(l), '\\.')[,1]
  row.names(l) <- NULL

  return(l)


}

