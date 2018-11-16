
main <- function(){
  load('/Users/daniellefever/Desktop/Lans_Taylor_Lab/NAFLD_MPS/data/published_data/microarray/microarray_meta/deg_metanalysis/fixed_pdata_fibrosis_list.RData')
  load('~/Desktop/Lans_Taylor_Lab/NAFLD_MPS/data/published_data/microarray/microarray_meta/ica/fastica_iter10_GSE48452_res.RData')
  ica_list <- list()
  ica_list[['GSE48452']] <- res
  load('~/Desktop/Lans_Taylor_Lab/NAFLD_MPS/data/published_data/microarray/microarray_meta/ica/fastica_iter10_GSE49541_res.RData')
  ica_list[['GSE49541']] <- res
  load('~/Desktop/Lans_Taylor_Lab/NAFLD_MPS/data/published_data/microarray/microarray_meta/ica/fastica_iter10_GSE89632_res.RData')
  ica_list[['GSE89632']] <- res

  p <- fixed_plist[["GSE48452"]]
  cur_a <- res$A
  cont_assoc <- pdata_continuous_lm_wrapper(p, cur_a)
  p_vals <- lapply(cont_assoc, function(x){
    x <- x['p_value']
  })
  nms <- names(p_vals)
  p_vals <- do.call(cbind, p_vals)
  colnames(p_vals) <- nms
}



