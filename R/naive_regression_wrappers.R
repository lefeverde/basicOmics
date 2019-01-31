#' Function performs lm over accross a number of genes
#'
#' @param expr_mat matrix of expression values
#' @param p_data matrix with sample data
#' @param model_formula a string of the form: y ~ x + o
#'
#' @return
#' @export
#'
#' @examples
basic_linear_regression <- function(expr_mat, p_data, model_formula){
  if(class(model_formula) != 'character'){
    stop('model_formula needs to be a string of the form y ~ + x1 + x2...')
  }
  expr_mat <- data.frame(expr_mat[row.names(p_data),], check.names = F)
  l <- lapply(expr_mat, function(x){
    tmp <- summary(lm(as.formula(model_formula), data=p_data))
  })
  out_df <- lmlist_to_df(l)
  # nms <- names(l)
  # x_list <- lapply(l, function(x){
  #   x <- x$coefficients[2,]
  # })
  # out_df <- do.call(rbind, x_list)
  # row.names(out_df) <- nms
  return(out_df)
}

#' Performs logistic regression with numerous Y vals
#'
#' @param expr_mat
#' @param p_data
#' @param model_formula
#'
#' @return
#' @export
#'
#' @examples
basic_logistic_regression <- function(expr_mat, p_data, model_formula){
  ### Args:
  ### expr_mat: matrix of expression values
  ### p_data: matrix with sample data
  ### model_formula: needs to be y ~ x + o
  ### y is outcome var, x is gene, o is other covariates

  if(class(model_formula) != 'character'){
    stop('model_formula needs to be a string of the form y ~ + x1 + x2...')
  }
  expr_mat <- data.frame(expr_mat, check.names = F)
  l <- lapply(expr_mat, function(x){
    tmp <- summary(glm(as.formula(model_formula), family = 'binomial', data=p_data))
    tmp <- coefficients(tmp)
    if(nrow(tmp) == 1){
      temp_row <- matrix(rep(NA, 4), nrow = 1)
      colnames(temp_row) <- colnames(tmp)
      x <- temp_row
    } else {
      x <- tmp[2,]
    }
  })
  res_logreg <- data.frame(do.call(rbind, l), check.names = F)
  colnames(res_logreg) <- c("Estimate","Std_Error","z_or_t_value", "p_value" )
  row.names(res_logreg) <- names(l)
  return(res_logreg)
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
#' @keywords external
#' @import foreach
#' @import doParallel
wilcoxOrKruskalOnA <- function (A, colAnnot, annot) {

  comp <- NULL
  A <- A[rownames(annot),]

  res_tests <- foreach::foreach(comp=as.list(A),.combine = c, .errorhandling = "stop") %dopar% {
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

