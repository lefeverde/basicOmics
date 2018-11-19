#TODO make plotting functions work with continous variables

make_pca_plot_data <- function(transformed_data, sample_map, gene_num=Inf){
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
  # changes metadata to column name to 'group'
  colnames(plot_data)[4] <- 'group'
  #plot_data$group <- as.factor(plot_data$group)
  #eturn(plot_data)
  # This just makes the labels for axises 
  ax_labels <- c(paste("PC1 (", signif(var_exp[1]*100, digits=4),"%)", sep=""), paste("PC2 (", signif(var_exp[2]*100, digits=4), "%)", sep=""), paste("PC3 (", signif(var_exp[3]*100, digits=4), "%)", sep=""))
  
  # returns a list with both plot data 
  # and the axes lables 
  ol <- list()
  ol[['plot_data']] <- plot_data
  ol[['ax_labels']] <- ax_labels
  return(ol)
}


pca_categorical_triplot <- function(plot_data, ax_labels, plot_ellipse=F, plot_title=''){
  ### THis uses ggplot facet rapping to 
  ### plot PC1 vs PC2, PC2 vs PC3, and
  ### PC1 vs PC3
  require(gridExtra)
  require(ggplot2)
  require(RColorBrewer)
  require(grid)
  require(lemon)
  # easiest way to do this
  pc1 <- plot_data$pc1 
  pc2 <- plot_data$pc2
  pc3 <- plot_data$pc3
  
  data_list <- list() # couldn't think of better way
  #ax_df <- rbind(c(ax_labels[1], ax_labels[2]), c(ax_labels[2], ax_labels[3]), c(ax_labels[1], ax_labels[3]))
  ax_df <- rbind(c(ax_labels[1], ax_labels[2]), c(ax_labels[1], ax_labels[3]), c(ax_labels[2], ax_labels[3]))
  
  # Kludge to add number per group by 
  # counting the number and appending
  # the values in the group col to include this
  
  # This feels like the worst R hack I've done
  # in atleast a 6 months. All this is 
  # basically trying to emulate a python dict
  # where the 'key' is the row.name and
  # resulting value in the 1st col is the 
  # is the output 'val'.
  # TODO think about turning this into a proper object 
  # where it can act like a python dict
  freq_df <- as.data.frame(table(plot_data$group))
  freq_df <- data.frame(freq=freq_df$Freq,row.names = freq_df$Var1)
  #plot_data$group <- as.character(plot_data$group[[1]])
  plot_data$group <- as.character(plot_data$group)
  #plot_groups <- as.character(plot_data$group)
  for(g in row.names(freq_df)){
    cur_freq <- as.character(freq_df[g,])
    append_group_name <- paste0(g, '\n(n=', cur_freq, ')\n')
    #append_group_name <- paste0(g, ' (n= ', cur_freq, ' ) ')
    #plot_data$group[plot_data$group == g] <- as.factor(append_group_name)
    plot_data$group[plot_data$group == g] <- append_group_name
    #plot_groups[plot_groups == g] <- as.character(append_group_name)
    
  }
  
  data_list[['pc1_vs_pc2']] <- data.frame(plotvars='PC1 vs. PC2',group=plot_data$group,x=pc1, y=pc2)
  #data_list[['pc2_vs_pc3']] <- data.frame(plotvars='PC2 vs. PC3',group=plot_data$group,x=pc2, y=pc3)
  data_list[['pc1_vs_pc3']] <- data.frame(plotvars='PC1 vs. PC3',group=plot_data$group,x=pc1, y=pc3)
  data_list[['pc2_vs_pc3']] <- data.frame(plotvars='PC2 vs. PC3',group=plot_data$group,x=pc2, y=pc3)
  rel_size <- 1.05 # Controls the rel size of text in theme
  
  
  #### This is just my theme ####
  
  # Makes the colors 
  # If there are more than 12,
  # color ramppallette is used to make more
  
  my_theme <- theme(plot.margin = unit(c(1,1,1,1), "cm"),panel.background = element_blank(), axis.title.y=element_text(size=rel(rel_size), face="bold", margin=margin(0,7.5,0,0)), axis.title.x=element_text(size=rel(rel_size), face="bold",margin=margin(7.5,0,0,0)),axis.text.y=element_text(size=rel(rel_size - .05),colour="black"),axis.text.x=element_text(size=rel(rel_size - .05), colour="black"), legend.title=element_text(size=rel(rel_size-.05)),legend.title.align=0.5,legend.key = element_blank(),legend.text=element_text(size=rel(rel_size-.25)),legend.position = 'bottom',panel.border=element_rect(fill=NA,colour="black", size=rel(1)))
  
  #### Basically puts the plots into a list ####
  plot_list <- list() # again, best solution I could think of
  idx <- 0 # this is so the list can be increment
  for(i in data_list){
    idx <- idx + 1 # increments 
    x_lab <- ax_df[idx,][1]
    y_lab <- ax_df[idx,][2]
    temp_plot <- ggplot(data = i, aes(x=x, y=y, colour=group)) + geom_point() + labs(x=x_lab, y=y_lab) + my_theme + geom_hline(yintercept=0) + geom_vline(xintercept=0) + scale_colour_brewer(palette="Set1") + theme(legend.position = 'right') 
    temp_plot$labels$colour <- plot_title # changes legend tit
    leg <- g_legend(temp_plot)
    if(plot_ellipse){
      temp_plot <- temp_plot + stat_ellipse(alpha=.15, geom = "polygon") 
    }
    plot_list[[idx]] <- temp_plot + theme(legend.position = 'none')
    
    
    #plot_list[[idx]] <- grob(plot_list[[idx]]) 
  }
  # lmat <- t(matrix(c(1,2,3, 4,4,4), nrow = 3))
  lmat <- t(matrix(c(1,2)))
  pca_plots <- arrangeGrob(grobs = plot_list, nrow = 1)
  # pca_and_legend_grid <- grid.arrange(pca_plots, legend=leg,  ncol=2, layout_matrix=lmat, widths = c(10,1), heights = c(10,1))
  # pca_and_legend_grid <- grid.arrange(pca_plots, leg,  ncol=2,  layout_matrix=lmat, widths = c(10,1), heights = c(10,.1))
  pca_and_legend_grid <- arrangeGrob(pca_plots, leg,  ncol=2,  layout_matrix=lmat, widths = c(10,1), heights = c(10,.1))
  # pca_and_legend_grid <- grid.arrange(pca_plots, legend=leg, ncol=2,  width = c(10,5))
  return(pca_and_legend_grid)
  # return(plot_list)
  
  # Super kludge to get it to plot all on one row
  #return(grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], plot_list[[3]], nrow = 1, position = 'top'))
  #return(plot_list)
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