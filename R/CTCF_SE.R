# Title     : TODO
# Objective : TODO
# Created by: boxia
# Created on: 7/12/18

setwd('/Users/boxia/PycharmProjects/CIG_RNA/CTCF_SE/')

celltypes = c('CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte', 'cardiac-muscle',
                     'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
                     'neural-cell', 'neutrophil', 'osteoblast', 'smooth-muscle-cell','H1-hESC')
#celltypes = c('CD4', 'HUVEC', 'NHLF', 'HMEC')
pdf("CTCF_SE_vs_non_CTCF_SE_SE_width_cmp.pdf")
par(mfrow = c(2, 4))
for (celltype in celltypes){
    df = read.csv(paste(celltype, '_CTCF_H3K27ac_SE_width_cmp.xls', sep=''), sep='\t')
    boxplot(df[, 'CTCF'], df[,'none.CTCF'], outline=FALSE,names=c('CSE','OSE'),col=c("red","cyan"), main=celltype, frame=FALSE)
    box(bty="l")
}
dev.off()


celltypes = c('CD4', 'HSMM', 'HMEC', 'HUVEC', 'B-cell', 'NHLF', 'astrocyte',
              'monocyte', 'hepatocyte', 'fibroblast-dermis', 'myotube', 'keratinocyte',
              'neural-cell', 'osteoblast', 'smooth-muscle-cell')
celltypes = c('CD4', 'HUVEC', 'NHLF', 'HMEC')
pdf("gene_CTCF_SE_vs_non_CTCF_SE_SE_width_cmp.pdf")
par(mfrow = c(2, 4))
for (celltype in celltypes){
  df = read.csv(paste(celltype, '_CSE_vs_OSE_genes_SEwidth.xls', sep=''), sep='\t')
  boxplot(df[, 'CSE'], df[,'OSE'], outline=FALSE,names=c('CSE','OSE'),col=c("red","cyan"), main=celltype, frame=FALSE)
  box(bty="l")
}
dev.off()

setwd('/Users/boxia/PycharmProjects/CIG_RNA/RNAexp/')
celltypes = c('HSMM', 'HMEC', 'HUVEC',  'NHLF')
pdf("CSE_OSE.pdf")
par(mfrow = c(2, 4))
for (celltype in celltypes){
  df = read.csv(paste(celltype, '_CSE_vs_OSE.xls', sep=''), sep='\t', row.names =1)
  boxplot(df[, 'CSE'], df[,'OSE'], outline=FALSE,names=c('CSE','OSE'),col=c("red","cyan"), main=celltype, frame=FALSE)
  box(bty="l")
}
dev.off()

celltypes = c('HMEC', 'HUVEC')
pdf("stem_CSE.pdf")
par(mfrow = c(2, 4))
for (celltype in celltypes){
  df = read.csv(paste(celltype, '_stem_CSE.xls', sep=''), sep='\t', row.names =1)
  boxplot(df[, 'H1.hESC_FPKM'], df[,paste(celltype, '_FPKM', sep='')], outline=FALSE,names=c('H1.hESC_FPKM',paste(celltype, '_FPKM', sep='')),col=c("purple","yellow"), main=celltype, frame=FALSE)
  box(bty="l")
}
dev.off()

celltypes = c('HMEC', 'HUVEC')
pdf("stem_OSE.pdf")
par(mfrow = c(2, 4))
for (celltype in celltypes){
  df = read.csv(paste(celltype, '_stem_OSE.xls', sep=''), sep='\t', row.names =1)
  boxplot(df[, 'H1.hESC_FPKM'], df[,paste(celltype, '_FPKM', sep='')], outline=FALSE,names=c('H1.hESC_FPKM',paste(celltype, '_FPKM', sep='')),col=c("purple","yellow"), main=celltype, frame=FALSE)
  box(bty="l")
}
dev.off()



