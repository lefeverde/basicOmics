
main <- function(){
  load('/Users/daniellefever/Desktop/Lans_Taylor_Lab/NAFLD_MPS/data/published_data/microarray/microarray_meta/deg_metanalysis/fixed_pdata_fibrosis_list.RData')
  load('~/Desktop/Lans_Taylor_Lab/NAFLD_MPS/data/published_data/microarray/microarray_meta/ica/fastica_iter10_GSE89632_res.RData')
  p <- fixed_plist[["GSE89632"]]
  cur_a <- res$A
  cont_assoc <- pdata_continuous_lm_wrapper(p, cur_a)
  p_vals <- lapply(cont_assoc, function(x){
    x <- x['p_value']
  })
  nms <- names(p_vals)
  p_vals <- do.call(cbind, p_vals)
  colnames(p_vals) <- nms
}



#' Compare the sample contributions according to their annotation level across the components.
#'
#' Wilcoxon or Kruskal-Wallis tests are performed depending on the number of levels in the considered annotation.
#' @title Comparison of distributions of sample groups
#' @param A A matrix of dimensions 'samples x components' containing the sample contributions
#' @param annot A matrix of dimensions 'samples x variables' containing the sample annotations
#' @param colAnnot The name of the column of \code{annot} to be considered
#' @return A vector of p-values
#' @author Anne Biton
#' @seealso \code{wilcox.test}, \code{kruskal.test}
#' @keywords internal
wilcoxOrKruskalOnA <- function (A, colAnnot, annot) {

    comp <- NULL
    A <- A[rownames(annot),]

    res_tests <- foreach(comp=as.list(A),.combine = c, .errorhandling = "stop") %dopar% {
      annotComp <- data.frame(comp=comp)
      annotComp[[colAnnot]] <- as.factor(annot[[colAnnot]])
      if (length(levels(annotComp[[colAnnot]])) == 2)
        res.test <- wilcox.test(comp~eval(as.name(colAnnot)), data = annotComp, na.action = "na.omit")
      else
        res.test <- kruskal.test(comp~eval(as.name(colAnnot)), data = annotComp, na.action = "na.omit")
      return(res.test$p.value)
    }
    return(unlist(res_tests))

  }
