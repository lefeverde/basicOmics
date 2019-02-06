library(edgeR)
library(dplyr)
library(tximport)
ex_p <- p_sub[,c('file_name', 'sex','age', 'bmi_surg', 'diagnosis', 'group')]

set.seed(42)
ex_p <-
  ex_p[ex_p$group != 'Lob',] %>%
  rownames_to_column(., var='sample_name') %>%
  group_by(., group) %>%
  sample_n(., 10) %>%
  as.data.frame(.)
row.names(ex_p) <- ex_p$sample_name


files  <- paste0('/Users/daniellefever/Desktop/R/thesisAnalysis/data-raw/DiStefano_raw_data/kallisto/', ex_p$file_name, '.gz')
names(files) <- row.names(ex_p)
ex_lengthscaledtpm_txi <- tximport(files, type = 'kallisto', tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = 'lengthScaledTPM')



