#' Plots PC1 by PC2 using ggplot
#'
#' @param transformed_data
#' @param sample_map
#' @param leg_row_num
#' @param gene_num
#' @param return_data
#'
#' @return
#' @export
#'
#' @examples
pca_plotter <- function(transformed_data, sample_map,leg_row_num=3, gene_num=Inf, return_data=FALSE){
  library(ggplot2)
  library(matrixStats)
  # Filter out any samples not listed in sample_map
  cur_subset_mat <- data.frame(transformed_data)
  cur_subset_mat <- transformed_data[,colnames(transformed_data) %in% rownames(sample_map)]
  #return(cur_subset_mat)
  # Taken from DESeq2 plotPCA function
  # Calculates the row wise variance
  rv <- rowVars(as.matrix(cur_subset_mat))


  # select the gene_num genes by variance
  # the seq_len thing looks weird, but it was in DESeq2 function
  # so leaving it.
  select <- order(rv, decreasing=TRUE)[seq_len(min(gene_num, length(rv)))]
  #return(select)
  # perform a PCA on the data in assay(x) for the selected genes
  fnmt_pcomp <- prcomp(t((cur_subset_mat)[select,]))

  var_exp <- (fnmt_pcomp$sdev^2)/sum(fnmt_pcomp$sdev^2)

  # Puts first 3 pcs into df
  plot_data <- data.frame(pc1=fnmt_pcomp$x[,1],
                          pc2=fnmt_pcomp$x[,2],
                          pc3=fnmt_pcomp$x[,3])
  # sorting because paranoia

  plot_data <- plot_data[sort(rownames(plot_data)),]
  # Getting PCs
  plot_data <- data.frame(pc1=fnmt_pcomp$x[,1],
                          pc2=fnmt_pcomp$x[,2],
                          pc3=fnmt_pcomp$x[,3])
  # merges sample metadata into df by rowname
  plot_data <- merge(plot_data, sample_map, by=0)
  # gets rid of extraneous column
  rownames(plot_data) <- plot_data$Row.names
  plot_data$Row.names <- NULL
  # changes metadata to column name to group
  colnames(plot_data)[4] <- 'group'
  plot_data$group <- as.factor(plot_data$group)
  #eturn(plot_data)
  # This just makes the labels for axises
  axlab.1 <- paste("PC1 (", signif(var_exp[1]*100, digits=4),"%)", sep="")
  axlab.2 <- paste("PC2 (", signif(var_exp[2]*100, digits=4), "%)", sep="")
  axlab.3 <- paste("PC3 (", signif(var_exp[3]*100, digits=4), "%)", sep="")

  if(return_data){
    return(plot_data)
  }
  # And here comes the plot!
  plt1 <- ggplot(data=plot_data,
                 aes(x=pc1,
                     y=pc2,
                     fill=group,
                     colour=group,
                     shape=group,
                     label=row.names(plot_data))) +
    # The geom point aes specificies colouring by group
    # and changes point shape by group as well
    geom_point(size = rel(1.95), aes(shape=factor(group), colour=factor(group))) +
    # geom_point(size = rel(1.5)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    # gets a pretty colour set
    # stat_ellipse(alpha=.15, geom = "polygon") +
    scale_colour_brewer(palette="Set1") +
    # scale_fill_brewer(palette="Set1") +
    labs(x=axlab.1, y=axlab.2) + theme_bw()

  # This all just setting the themes the way I like it

  plt2 <- plt1 + theme(plot.margin = unit(c(1,1,1,1), "cm"),panel.background = element_blank(), axis.title.y=element_text(size=rel(1.75), face="bold", margin=margin(0,7.5,0,0)), axis.title.x=element_text(size=rel(1.75), face="bold",margin=margin(7.5,0,0,0)),axis.text.y=element_text(size=rel(1.5),colour="black"),axis.text.x=element_text(size=rel(1.5), colour="black"), legend.title=element_blank(),legend.key = element_blank(),legend.text=element_text(size=rel(1.25)),legend.position = 'bottom',panel.border=element_rect(fill=NA,colour="black", size=rel(1.9)))
  #title=element_text(size=22,


  # this just splits the legend into two rows
  # when there is more than 3 groups because of
  # ugly formatting
  # if(length(levels(factor(plot_data$group))) > 3){
  #   plt2 <- plt2 + guides(col=guide_legend(nrow = 2))
  # }
  group_num <- length(levels(factor(plot_data$group)))
  if (group_num > 6){
    plt2 <- plt2 + scale_shape_manual(values = seq(1,group_num))
  }
  plt2 <- plt2 + guides(col=guide_legend(nrow = leg_row_num))
  plt2 <- plt2 + scale_x_continuous(breaks = pretty(plot_data$pc1, n=7)) + scale_y_continuous(breaks = pretty(plot_data$pc2, n=7))


  return(plt2)
}



#' Title
#'
#' @param transformed_data
#' @param sample_map
#' @param gene_num
#'
#' @return
#' @export
#'
#' @examples
make_pca_plot_data <- function(transformed_data, sample_map, gene_num=Inf){
  library(ggplot2)
  library(matrixStats)
  #TODO make plotting functions work with continous variables
  if(all(row.names(sample_map) != colnames(transformed_data))){
    stop('row.names of sample_map need to equal colnames of transformed_data')
  }
  # Filter out any samples not listed in sample_map
  cur_subset_mat <- data.frame(transformed_data)
  cur_subset_mat <- transformed_data[,colnames(transformed_data) %in% rownames(sample_map)]
  #return(cur_subset_mat)
  # Taken from DESeq2 plotPCA function
  # Calculates the row wise variance
  cur_subset_mat <- as.matrix(cur_subset_mat)
  rv <- rowVars(cur_subset_mat)
  names(rv) <- row.names(cur_subset_mat)
  cur_subset_mat <- cur_subset_mat[rv != 0,]
  rv <- rv[row.names(cur_subset_mat)]
  # select the gene_num genes by variance
  # the seq_len thing looks weird, but it was in DESeq2 function
  # so leaving it.
  select <- order(rv, decreasing=TRUE)[seq_len(min(gene_num, length(rv)))]
  #return(select)
  # perform a PCA on the data in assay(x) for the selected genes
  fnmt_pcomp <- prcomp(t((cur_subset_mat)[select,]), center = TRUE, scale = FALSE)

  var_exp <- (fnmt_pcomp$sdev^2)/sum(fnmt_pcomp$sdev^2)

  # Puts first 3 pcs into df
  plot_data <- data.frame(pc1=fnmt_pcomp$x[,1],
                          pc2=fnmt_pcomp$x[,2],
                          pc3=fnmt_pcomp$x[,3])
  # sorting because paranoia

  plot_data <- plot_data[sort(rownames(plot_data)),]
  # Getting PCs
  plot_data <- data.frame(pc1=fnmt_pcomp$x[,1],
                          pc2=fnmt_pcomp$x[,2],
                          pc3=fnmt_pcomp$x[,3])
  plot_data <- data.frame(sample_map[row.names(plot_data),,drop=FALSE], plot_data)

  # This just makes the labels for axises

  ax_labels <- c(paste("PC1 (", signif(var_exp[1]*100, digits=4),"%)", sep=""), paste("PC2 (", signif(var_exp[2]*100, digits=4), "%)", sep=""), paste("PC3 (", signif(var_exp[3]*100, digits=4), "%)", sep=""))

  # returns a list with both plot data
  # and the axes lables
  ol <- list()
  ol[['plot_data']] <- plot_data
  ol[['ax_labels']] <- ax_labels
  return(ol)
}




pca_data_wrapper <- function(tmat, sample_map, gene_num=Inf){
  # This is a wrapper for the make_pca_plot_data()
  # function.
  pca_dat <- make_pca_plot_data(tmat, sample_map[,1,drop=F], gene_num)

  base_pca <- list()
  for(i in colnames(sample_map)){
    temp_data <- pca_dat[[1]]
    temp_data$group <- sample_map[,i, drop=T]
    base_pca[[i]] <- temp_data
  }
  ol <- list(base_pca, pca_dat[[2]])

  return(ol)
}
