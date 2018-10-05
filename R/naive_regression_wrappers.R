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
  forms <- paste()
  l <- lapply(expr_mat, function(x){
    tmp <- summary(lm(as.formula(model_formula), data=p_data))
  })
  nms <- names(l)
  x_list <- lapply(l, function(x){
    x <- x$coefficients[2,]
  })
  out_df <- do.call(rbind, x_list)
  row.names(out_df) <- nms
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
basic_generalized_regression <- function(expr_mat, p_data, model_formula){
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

}

#' Performs linear regression over all continuous variables
#' in phenodata data.frame. I made this to test associations
#' between continuous variables and the independent components
#' from ICA.
#'
#' @param p_data phenodata data frame
#' @param ic_mat matrix of independent components
#'
#' @return list of linear regression sums
#' @export
#'
#' @examples
pdata_continuous_lm_wrapper <- function(p_data, ic_mat){
  big_lmlist <- list()
  for(i in names(p_data)){
    cur_var <- p_data[,i]
    if(is.numeric(cur_var)){
      lm_sum <- summary(lm(ic_mat ~ cur_var))
      big_lmlist[[i]] <- lmlist_to_df(lm_sum)
    }

  }
  return(big_lmlist)
}

#' Basically wrapper to apply Kruskal-Wallis test to
#' the categorical sample variables. I made this to
#' test associations between categorical variables
#' and the independent components from ICA. Why
#' Kruskal-Wallis? See MineICA
#'
#'
#' @param p_data
#' @param ic_mat
#'
#' @return
#' @export
#'
#' @examples
pdata_categorical_wrapper <- function(p_data, ic_mat){
  big_catlist <- list()
  for(i in names(p_data)){
    cur_var <- p_data[,i]
    if(! is.numeric(cur_var)){
      #lm_sum <- summary(lm(ic_mat ~ cur_var))

      big_catlist[[i]] <- lmlist_to_df(lm_sum)
    }

  }
  return(big_lmlist)
}


