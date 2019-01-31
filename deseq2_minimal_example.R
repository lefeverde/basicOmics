library(DESeq2)
library(dplyr)



#colData(dds) <- DataFrame(condition=factor(c(rep('A', 12),rep('B', 18),rep('C', 30))),row.names = colnames(dds))


#counts(dds)[1,] <- as.integer(c(rnbinom(m_samples/2, 1, mu=10), rnbinom(m_samples/2, 10, mu=50)))

set.seed(42)
n_genes <- 500
m_samples <- 60
dds <- makeExampleDESeqDataSet(n=n_genes, m = m_samples)
colData(dds) <- DataFrame(condition=factor(c(rep('A', m_samples/4 - 3),rep('B', m_samples/4 + 3),rep('C', m_samples/2))),row.names = colnames(dds))
counts(dds)[1,] <- as.integer(c(rnbinom(m_samples/2, 5, mu=15) + rnbinom(m_samples/2, 5, mu=15), rnbinom(m_samples/2, 1, mu=100) + rnbinom(m_samples/2, 1, mu=100)))
counts(dds)[1, 1] <- as.integer(3*max(counts(dds)[1,]))
counts(dds) <- matrix(sapply(counts(dds), as.integer), dim(counts(dds)))
mcols(dds) <- NULL
dds <- DESeq(dds)
mcols(dds)[["maxCooks"]] <- apply(assays(dds)[["cooks"]],1,max)
plotCounts(dds, gene = 'gene1')
dds <- DESeq(dds)

res <- data.frame(results(dds, cooksCutoff = 1.5))



mout <- data.frame(mcols(dds_out))
gene1 <-
  do.call(rbind, lapply(assays(dds_out), function(x){data.frame(x)['gene1',]})) %>%
  t(.) %>%
  data.frame(.)
gene1 <- cbind(data.frame(colData(dds))[row.names(gene1),], gene1)



dds_out <-
  estimateSizeFactors(dds) %>%
  estimateDispersions(.) %>%
  nbinomWaldTest(.) %>%
  replaceOutliersWithTrimmedMean(., cooksCutoff = 3) %>%
  nbinomWaldTest(.)





counts(dds)[1,] <- c(as.integer(5000), rbinom(m_samples/2 -1, 10, .1), rbinom(m_samples/2, 50, .25))

res <-
  DESeq(dds) %>%
  results(., name='condition_B_vs_A') %>%
  data.frame(.)

count_mat <-
  rbinom(n_genes*m_samples, 1000, .5) %>%
  matrix(., ncol=m_samples) %>%
  data.frame(.,row.names=paste0('gene', seq(n_genes)))
names(count_mat) <- paste0('sample', seq(m_samples))


