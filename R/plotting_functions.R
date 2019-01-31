#' Creates a volcano plot using data in long format
#'
#' This function creates a volc plot using the returned results
#' in long format. Meant to be used with my results wappers
#'
#' @param plot_data plot data in long format
#' @param plot_labels whether labels should be drawn for some genes
#' @param label_num number of genes to draw labels
#'
#' @return
#' @export
#'
#' @example
volcano_plotter <- function(plot_data, plot_labels=TRUE, label_num=10){

  require(ggplot2) #TODO consider fixing explicit `library` load
  library("ggrepel") #Avoid overlapping labels
  split_names <- limma::strsplit2(plot_data$coefficient, '_')
  name_len <- length(split_names[1,])
  plot_data$group <- do.call(paste, data.frame(split_names[,2:name_len]))


  plot_data$fill_fact <- as.factor(ifelse((abs(plot_data$log2FoldChange)>= 1.3 & plot_data$padj < .05), 2, 1))
  plt1 <- ggplot(data=plot_data,
                 aes(x=log2FoldChange,
                     y=-log10(padj),
                     label=external_gene_name,
                     group=group,
                     colour=fill_fact)) +
    scale_colour_manual(values=c('black', '#e31a1c')) +
    theme_bw() +
    geom_point(size=rel(.375), alpha=.5, position = position_jitterdodge()) +
    geom_hline(yintercept = 1.30103, colour='#1f78b4', size=rel(.5), linetype='dashed') +
    geom_vline(xintercept = 1.3, colour='#1f78b4', size=rel(.5), linetype='dashed') +
    geom_vline(xintercept = -1.3, colour='#1f78b4', size=rel(.5), linetype='dashed') +
    scale_x_continuous(breaks = pretty(plot_data$log2FoldChange, n=7)) +
    scale_y_continuous(breaks = pretty(-log10(plot_data$padj), n=7))

  plt2 <- plt1  + theme(plot.margin = unit(c(1,1,1,1), "cm"),
                        panel.background = element_blank(),
                        axis.title.y=element_text(size=rel(1.25),face="bold",
                                                  margin=margin(0,7.5,0,0)),
                        axis.title.x=element_text(size=rel(1.25),
                                                  #vjust=-.5,
                                                  #hjust=.25,
                                                  face="bold",
                                                  margin=margin(7.5,0,0,0)),
                        axis.text.y=element_text(size=rel(1.25),
                                                 colour="black"),
                        axis.text.x=element_text(size=rel(1.25),
                                                 colour="black"),
                        legend.title=element_blank(),
                        legend.key = element_blank(),
                        legend.text=element_text(size=rel(1)),
                        legend.position = 'none',
                        panel.border=element_rect(fill=NA,
                                                  colour="black",
                                                  size=rel(1)),
                        plot.title=element_text(size=rel(1.5),
                                                hjust=.5,
                                                colour="black",
                                                face="bold"))
  plt3 <- plt2 + labs(x='Log2 Fold Change',
                      y='-Log10 FDR p-value',
                      title=plot_data$group[1])
  if(plot_labels){
    sig_genes <- plot_data[plot_data$padj < .05,]

    # Breaks genes into negative and positive log2fc
    # This is so both pos/neg sides of volc get annots
    #browser()
    neg_genes <- sig_genes[sig_genes$log2FoldChange < 0,]
    neg_genes <- neg_genes[order(neg_genes$log2FoldChange,neg_genes$padj),][1:label_num,]
    # same as above except for the pos genes
    pos_genes <- sig_genes[sig_genes$log2FoldChange > 0,]
    pos_genes <- pos_genes[order(-pos_genes$log2FoldChange, pos_genes$padj),][1:label_num,]
    top_genes <- rbind(neg_genes, pos_genes)
    plt3 <- plt3 + geom_text_repel(data=top_genes,
                                   aes(x=log2FoldChange,
                                       y=-log10(padj),
                                       label=external_gene_name),
                                   size=rel(1.05),
                                   segment.size=rel(.05))
  }
  # TODO fix this hardcoded kludge
  cur_df <- na.omit(temp_split[[i]])
  cur_df <- cur_df[abs(cur_df$log2FoldChange) >= 1.3 & cur_df$padj < .05,]
  upreg_num <- nrow(cur_df[cur_df$log2FoldChange > 0,])
  downreg_num <- nrow(cur_df[cur_df$log2FoldChange < 0,])
  out_str <- paste0('\nUpregulated genes (log2FC >= 1.3): ', upreg_num, '\nDownregulated genes (log2FC <= 1.3): ', downreg_num)
  plt3 <- plt3  + annotate(geom = 'text', label=out_str, x = -Inf, y = Inf, hjust = -.5, vjust = 1)
  #plt3 <- plt3 + facet_wrap(.~group, ncol = 1, scales='free')


  # plt4 <- plt3 + facet_wrap(~ treatment, nrow = 2) + theme(legend.position = 'none',
  #                                                          strip.text=element_text(size=20),
  #                                                          panel.spacing = unit(1.75, "lines"))
  return(plt3)
}

