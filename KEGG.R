
rm(list=ls())
gc()
options(stringsAsFactors=FALSE)

# Load libraries
{
  library(ggplot2)
  library(GSEABase)
  library(GSVA)
  library(gplots)
  library(pheatmap)
  library(gdata)
  library(geneplotter)
}

# LOad the method of phylogenetic tree
{
  cor.dist <- function (x) 
  {
    as.dist(1-cor(t(x), use="pairwise.complete.obs",method="kendall"))
  }
  cor.dist.spearman <- function (x) 
  {
    as.dist(1-cor(t(x), use="pairwise.complete.obs",method="spearman"))
  }
  dist.manhattan <- function(x){dist(x, method="manhattan")}
  hclust.ward.d <- function(x){hclust(x, method="ward.D")}
  hclust.ward.d2 <- function(x){hclust(x, method="ward.D2")}
}


# Read Hallmark Gene Sets, please obtain it from original authors' website: MSigDB
hm.sets <- getGmt("h.all.v7.0.symbols.gmt")

# Load primary dataset
if(!file.exists("m4T1.primary.Rdata")){
  load("m4T1.D4.sceset.qc.new.Rdata")
  exprs <- log2(calculateCPM(counts(sce.all)) +1)
  exprs_norm <- t(scale(t(exprs)))
  save(ann, exprs_norm, file = "m4T1.primary.Rdata")
}
load("m4T1.primary.Rdata")

# Isolate D4 breast cancer cells
ann.cellname_D4 <- ann[ann$CellType=="BoneMet_D4", ]
exprs_D4 <- exprs_norm[ , colnames(exprs_norm)==row.names(ann.cellname_D4)]
rownames(exprs_D4) <- toupper(rownames(exprs_D4))

# Run GSVA on the Breast Cancer data using Hallmark Gene Sets
datasets.hm <- gsva(as.matrix(exprs_D4), hm.sets, method="gsva")
rownames(datasets.hm) <- gsub("HALLMARK_","",rownames(datasets.hm))

# Generate GSVA analysis
pheatmap(datasets.hm, treeheight_col = 30, treeheight_row = 0, show_colnames = F, clustering_method="ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(50), fontsize_row = 7)
# heatmap.2(datasets.hm, dendrogram = "col", lmat=rbind( c(0, 3, 4), c(2,1,0 ) ), lwid=c(1.5, 4, 2 ),
          # col=c(rep("blue",20), dChip.colors(100), rep("red",20)), trace="none", keysize = 1,
          # scale=("row"), dist=cor.dist, hclust=hclust.ward.d)$colInd


#============================#
#       Musician: Resonance  #
#           Date: 2019/10/05 #
# Revised author: Resonance  #
#           Time: 2019/10/07 #
#============================#