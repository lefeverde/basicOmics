make_example_pdata <- function(){
  example_pdata <- data.frame(A=rep(5, 5), B=as.factor(c('a','a','b','c','e')), group=as.factor(paste0('G', seq(1,5))))
  example_pdata$group <- relevel(example_pdata$group, 'G3')
  example_pdata2 <-
    rbind(example_pdata, example_pdata) %>%
    dplyr::arrange(., group)
  example_pdata2$sex <- factor(c(rep('male', 4), rep('female', 6)))
  return(example_pdata2)
}

