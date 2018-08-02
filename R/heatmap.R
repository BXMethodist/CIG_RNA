library('ggplot2')
library('pheatmap')
library('grid')
library('gridExtra')

setwd('/Users/boxia/PycharmProjects/CIG_RNA/results/')

for (cell_type in c('HUVEC', 'CD34', 'GM12878', 'HSMM', 'HMEC', 'NHLF', 'MSC', 'H1-hESC',  'MRG', 'neural')){
# for (cell_type in c('HUVEC')){
  for (marker in c('h3k4me3', 'h3k27me3')) {
    for (target in c('k27m3', 'k4m3')) {
      plot_list = list()
      df = read.csv(paste(paste(cell_type, marker, 'lncRNA', target, sep="_"), 'sig.xls', sep=""), sep='\t', row.names=1)
      df = df+  0.01
      bk = unique(c(seq(0, 10, length=100),seq(10,40,length=100),seq(40,as.integer(max(df))+2, length=100)))
      colors1=colorRampPalette(c("blue","black"))(100)
      colors2=colorRampPalette(c("black","yellow"))(100)
      colors3=colorRampPalette(c("yellow","yellow"))(as.integer(max(df))+2-40)
      # colors = colorRampPalette(c("blue", "black","yellow"))(length(bk))
      colors = c(colors1, colors2, colors3)
      x = pheatmap(df,scale= "none",cluster_rows=FALSE,cluster_cols=FALSE,dendrogram=('none'),color=colors,breaks=bk,border_color=NA, show_colnames=FALSE, show_rownames=FALSE, legend=FALSE)[[4]]
      plot_list[[paste(cell_type, 'pos')]] = x
      ggsave(paste(paste("/Users/boxia/PycharmProjects/CIG_RNA/results/", cell_type, marker, target, sep="_"), "lncRNA_heatmap.png", sep="_"), arrangeGrob(grobs=plot_list, ncol=1))
    }
  }
}


## Coding genes
for (cell_type in c('HUVEC', 'CD34', 'GM12878', 'HSMM', 'HMEC', 'NHLF', 'MSC', 'H1-hESC',  'MRG', 'neural')){
# for (cell_type in c('HUVEC')){
  for (marker in c('h3k4me3', 'h3k27me3')) {
    for (target in c('k4m3', 'k27m3')) {
    # for (target in c('k4m3')) {
      plot_list = list()
      df = read.csv(paste(paste(cell_type, marker, 'gene', target, sep="_"), 'sig.xls', sep=""), sep='\t')
      df$gene_id = NULL
      df = df+  0.01
      bk = unique(c(seq(0, 10, length=100),seq(10,40,length=100),seq(40,as.integer(max(df))+2, length=100)))
      colors1=colorRampPalette(c("blue","black"))(100)
      colors2=colorRampPalette(c("black","yellow"))(100)
      colors3=colorRampPalette(c("yellow","yellow"))(as.integer(max(df))+2-40)
      # colors = colorRampPalette(c("blue", "black","yellow"))(length(bk))
      colors = c(colors1, colors2, colors3)
      x = pheatmap(df,scale= "none",cluster_rows=FALSE,cluster_cols=FALSE,dendrogram=('none'),color=colors,breaks=bk,border_color=NA, show_colnames=FALSE, show_rownames=FALSE, legend=FALSE)[[4]]
      plot_list[[paste(cell_type, 'pos')]] = x
      ggsave(paste(paste("/Users/boxia/PycharmProjects/CIG_RNA/results/", cell_type, marker, target, sep="_"), "gene_heatmap.png", sep="_"), arrangeGrob(grobs=plot_list, ncol=1))
    }
  }
}



