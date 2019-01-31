# library(MineICA)
# library(doParallel)
# load('./filtered_expression_matrixes.RData')
#
# make_ica_data <- function(geo_id){
#   registerDoParallel(cores=8)
#   e <- t(filt_expr[[geo_id]])
#   #e <- data.frame(t(apply(e,1,scale,scale=FALSE)))
#   n_samples <- ncol(e)
#   res <- clusterFastICARuns(e, n_samples, alg.type="parallel", funClus="hclust", method="ward")
#   #res <- runICA(method= 'fastICA', X = e, nbComp = n_samples, alg.type = 'parallel')
#   out_handle <- paste0('fastica_iter10_', geo_id, '_res.RData')
#   save(res, file = out_handle)
# }
# make_ica_data('GSE48452')
# make_ica_data('GSE49541')
# make_ica_data('GSE89632')
