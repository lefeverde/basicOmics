context("analysis.R tests ")
library(basicOmics)


get_example_data <- function(){
  rds_to_load <-
    c('example_dge_object.rds',
      'example_voom_object.rds',
      'example_fit_object.rds',
      'expected_limma_res.rds') %>%
    lapply(., function(x){
      system.file('extdata',
                  x,
                  package = 'basicOmics',
                  mustWork = TRUE)
    }) %>% setNames(., c('dge', 'v', 'fit', 'e_res'))
}


test_that('overlap_ratio calculates overlap correctly',{
  x <- LETTERS[1:5]
  y <- LETTERS[2:10]
  expect_equal(overlap_ratio(x,y), .8)
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
  e_mat <-  c(1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 0.4615385, 1.0000000, 0.4615385, 1.0) %>%
    matrix(data=., nrow=3, ncol=3, dimnames = list(names(l1), names(l1)))
  omat <- get_overlap_matix(l1)
  expect_equivalent(e_mat, omat, tolerance=1e-6)
})



