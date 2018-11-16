
#' Title
#'
#' @param gene_vector vector of human ensembl genes
#'
#' @return data.frame of with ensembl gene and entrez gene ids
#' @export
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
#' @param gene_vector
#'
#' @return
#' @export
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
#' @param gene_vector
#'
#' @return
#' @export
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

