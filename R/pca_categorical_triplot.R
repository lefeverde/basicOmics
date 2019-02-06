#### TODO Functions ####


#' PCA categorical triplot
#'
#' Currently this is too unweildy to use out in the wild
#' but it's not a bad overall, just needs some work.
#'
#' @param plot_data
#' @param ax_labels
#' @param plot_ellipse
#' @param plot_title
#' @param plot_legend
#'
#' @return
#' @keywords internal
#'
#' @examples
pca_categorical_triplot <- function(plot_data, ax_labels, plot_ellipse=F, plot_title='', plot_legend=TRUE){
  ### THis uses ggplot facet rapping to
  ### plot PC1 vs PC2, PC2 vs PC3, and
  ### PC1 vs PC3

  # easiest way to do this
  pc1 <- plot_data$pc1
  pc2 <- plot_data$pc2
  pc3 <- plot_data$pc3

  data_list <- list() # couldn't think of better way
  ax_df <- rbind(c(ax_labels[1], ax_labels[2]), c(ax_labels[1], ax_labels[3]), c(ax_labels[2], ax_labels[3]))


  # This counts the occurances the groups
  # to be added to legend
  group_cnts <- plyr::count(plot_data$group)
  g <- as.character(group_cnts[,1])
  freqs <- group_cnts[,2]
  append_group_name <- paste0(g, '\n(n=', freqs, ')\n')



  data_list[['pc1_vs_pc2']] <- data.frame(plotvars='PC1 vs. PC2',group=plot_data$group,x=pc1, y=pc2)

  data_list[['pc1_vs_pc3']] <- data.frame(plotvars='PC1 vs. PC3',group=plot_data$group,x=pc1, y=pc3)
  data_list[['pc2_vs_pc3']] <- data.frame(plotvars='PC2 vs. PC3',group=plot_data$group,x=pc2, y=pc3)
  rel_size <- 1.05 # Controls the rel size of text in theme


  #### This is just my theme ####

  # Makes the colors
  # If there are more than 12,
  # color ramppallette is used to make more


  my_theme <- theme(panel.background = element_blank(), axis.title.y=element_text(size=rel(rel_size-.1), face="bold"), axis.title.x=element_text(size=rel(rel_size-.1), face="bold"),axis.text.y=element_text(size=rel(rel_size - .25),colour="black"),axis.text.x=element_text(size=rel(rel_size - .25), colour="black"),legend.direction = "horizontal", legend.title=element_blank(),legend.title.align=0.5,legend.key = element_blank(),legend.text=element_text(size=rel(rel_size-.25)),legend.position = 'bottom',panel.border=element_rect(fill=NA,colour="black", size=rel(1)))

  #### Basically puts the plots into a list ####
  plot_list <- list() # again, best solution I could think of
  idx <- 0 # this is so the list can be increment
  for(i in data_list){
    idx <- idx + 1 # increments
    x_lab <- ax_df[idx,][1]
    y_lab <- ax_df[idx,][2]
    temp_plot <- ggplot(data = i, aes(x=x, y=y, colour=group)) + geom_point(size=rel(.5), alpha=.75) + labs(x=x_lab, y=y_lab) + my_theme + geom_hline(yintercept=0) + geom_vline(xintercept=0)
    temp_plot <- temp_plot + scale_x_continuous(breaks = pretty(i$x, n=5)) + scale_y_continuous(breaks = pretty(i$y, n=5))
    #temp_plot$labels$colour <- plot_title # changes legend tit
    #leg <- g_legend(temp_plot)
    if(plot_ellipse){
      temp_plot <- temp_plot + stat_ellipse(alpha=.15, geom = "polygon")
    }
    # Re labels legend with frequency counts
    if(length(g) < 7){
      temp_plot <- temp_plot + scale_colour_manual(values=ggsci::pal_d3()(length(g)),breaks=g, labels=append_group_name)
    } else {
      temp_plot <- temp_plot + scale_colour_discrete(breaks=g, labels=append_group_name)
    }
    temp_plot <- temp_plot + guides(colour = guide_legend(override.aes = list(size=rel(1.375), alpha=1)))
    plot_list[[idx]] <- temp_plot

  }
  #return(plot_list)
  if(plot_legend== TRUE){
    out_plot <- ggpubr::ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]], ncol=3, nrow=1,common.legend = TRUE,legend="bottom")

  }
  if(plot_legend == FALSE){
    out_plot <- ggpubr::ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]], ncol=3, nrow=1,common.legend = TRUE,legend="none")
  }

  return(out_plot)


}
