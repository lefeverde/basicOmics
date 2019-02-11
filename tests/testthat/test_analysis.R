context("analysis.R tests ")
library(basicOmics)


test_that('overlap_ratio calculates overlap correctly',{
  x <- LETTERS[1:5]
  y <- LETTERS[2:10]
  expect_equal(overlap_ratio(x,y), .4)
})

test_that('get_overlap_matix produces correct overlap matrix',{
  set.seed(42)
  n_genes_to_sample <- rnbinom(3, 5, .2)
  genes <- paste0('gene', seq(1, max(n_genes_to_sample)))
  l1 <- lapply(seq_along(n_genes_to_sample), function(i){
    n <- n_genes_to_sample[i]
    sample(genes, n)
  }) %>%
    setNames(., paste0('L1_', LETTERS[1:length(n_genes_to_sample)]))
  e_mat <-  c( 1.0000000, 0.4642857, 0.5000000, 0.4642857, 1.0000000, 0.2857143, 0.5000000, 0.2857143, 1.0000000) %>%
    matrix(data=., nrow=3, ncol=3, dimnames = list(names(l1), names(l1)))
  omat <- get_overlap_matix(l1)

  expect_equal(e_mat, omat)


})
