rm(list = ls())
library(dplyr)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

####挑选出ATMi/IRAK1i加入后突变数目增加的gene####
####Select genes where the number of mutations increases after the addition of ATMi/IRAK1i inhibitors####
load("2.mutant_gene.Rdata")
total.gene$log2FC_EtoIvsEto <- log2((total.gene$EtoI+1)/(total.gene$Eto+1))
total.gene$log2FC_EtoAvsEto <- log2((total.gene$EtoA+1)/(total.gene$Eto+1))
EtoI.up.gene <- total.gene[total.gene$log2FC_EtoIvsEto>=log2(4),]$SYMBOL
EtoA.up.gene <- total.gene[total.gene$log2FC_EtoAvsEto>=log2(4),]$SYMBOL
others.gene <- total.gene[total.gene$log2FC_EtoIvsEto<0 & 
                            total.gene$log2FC_EtoAvsEto<0,]$SYMBOL

####导入其它ChIPseq在基因TSS±5kb范围的counts矩阵####
####Import ChIP-seq counts matrices within ±5 kb around gene TSS####
counts <- read.table("TSS_ALL_seq_counts.txt",header=F)
colnames(counts) <- c("chr","start","end","gene_id","none","strand_direction",
                      "DDX39AKO_EdUHU_RNA","DRIP_R-loop","FANCD2_ChIP","ATACseq","EUseq",
                      "H3K4me1_ChIP","H3K4me3_ChIP","H3K9ac_ChIP","H3K27ac_ChIP",
                      "H3K27me3_ChIP","H3K36me3_ChIP","H3K9me3_ChIP",
                      "Canonical_MiDAS","siRAD51_atypical_MiDAS","siBRCA2_atypical_MIDAS",
                      "PolII_ChIP","REPLI_earlyS","REPLI_midS","REPLI_lateS",
                      "EZH2_ChIP","HDAC1_ChIP","HDAC2_ChIP","HP1g_ChIP","mC5","hmC5")
counts <- counts[,c("chr","start","end","gene_id",
                    "EUseq","H3K4me1_ChIP","H3K4me3_ChIP","H3K9ac_ChIP","H3K27ac_ChIP",
                    "H3K27me3_ChIP","H3K36me3_ChIP","H3K9me3_ChIP","PolII_ChIP",
                    "EZH2_ChIP","HDAC1_ChIP","HDAC2_ChIP","HP1g_ChIP","mC5","hmC5")]

####计算CPM####
####Calculate CPM value####
counts$EUseq <- counts$EUseq/sum(counts$EUseq)*1000000
counts$H3K4me1_ChIP <- counts$H3K4me1_ChIP/sum(counts$H3K4me1_ChIP)*1000000
counts$H3K4me3_ChIP <- counts$H3K4me3_ChIP/sum(counts$H3K4me3_ChIP)*1000000
counts$H3K9ac_ChIP <- counts$H3K9ac_ChIP/sum(counts$H3K9ac_ChIP)*1000000
counts$H3K27ac_ChIP <- counts$H3K27ac_ChIP/sum(counts$H3K27ac_ChIP)*1000000
counts$H3K27me3_ChIP <- counts$H3K27me3_ChIP/sum(counts$H3K27me3_ChIP)*1000000
counts$H3K36me3_ChIP <- counts$H3K36me3_ChIP/sum(counts$H3K36me3_ChIP)*1000000
counts$H3K9me3_ChIP <- counts$H3K9me3_ChIP/sum(counts$H3K9me3_ChIP)*1000000
counts$PolII_ChIP <- counts$PolII_ChIP/sum(counts$PolII_ChIP)*1000000
counts$EZH2_ChIP <- counts$EZH2_ChIP/sum(counts$EZH2_ChIP)*1000000
counts$HDAC1_ChIP <- counts$HDAC1_ChIP/sum(counts$HDAC1_ChIP)*1000000
counts$HDAC2_ChIP <- counts$HDAC2_ChIP/sum(counts$HDAC2_ChIP)*1000000
counts$HP1g_ChIP <- counts$HP1g_ChIP/sum(counts$HP1g_ChIP)*1000000
counts$mC5 <- counts$mC5/sum(counts$mC5)*1000000
counts$hmC5 <- counts$hmC5/sum(counts$hmC5)*1000000

# CPM转化为Z score
# Convert CPM (Counts Per Million) to Z-score.
counts_clear <- counts[counts$EUseq<100000,]
counts_cleared <- counts_clear %>%
  group_by(gene_id) %>%
  filter(PolII_ChIP == max(PolII_ChIP))
table(duplicated(counts_cleared$gene_id))
counts_cleared <- counts_cleared[!duplicated(counts_cleared$gene_id),]
rownames(counts_cleared) <- counts_cleared$gene_id

z_scores <- counts_cleared[,5:ncol(counts_cleared)]
z_scores <- log2(z_scores+1)
z_scores <- scale(z_scores, center = TRUE, scale = TRUE)
z_scores_df <- as.data.frame(z_scores)
rownames(z_scores_df) <- counts_cleared$gene_id

####热图####
####pheatmap####
n <- data.frame()
i=1
repeat{
  n[i,1] <- median(z_scores_df[rownames(z_scores_df) %in% EtoI.up.gene,i])
  n[i,2] <- median(z_scores_df[rownames(z_scores_df) %in% EtoA.up.gene,i])
  n[i,3] <- median(z_scores_df[rownames(z_scores_df) %in% others.gene,i])
  rownames(n)[i] <- colnames(z_scores)[i]
  i=i+1
  if(i>(ncol(z_scores)-1)){
    break
  }
}
colnames(n) <- c("IRAK1i_related","ATMi_related","others")
n <- n[order(n$ATMi_related,decreasing = T),]

n2 <- n[rownames(n) %in% c("H3K4me3_ChIP","PolII_ChIP",
                           "H3K36me3_ChIP","H3K27me3_ChIP","H3K9me3_ChIP"),]
my_palette <- colorRampPalette(c("#4575B4", "white", "#D73027"))

p1 <- pheatmap::pheatmap(n2,show_colnames =T,show_rownames = T,
                         fontsize=7,
                         #scale = "row",
                         angle_col=90,
                         color = my_palette(100),
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         breaks = seq(-0.3,0.5,0.008),
                         fontsize_col = 12,             
                         fontsize_row = 12)
p1
ggsave(p1,filename = "./results/genes where mutations increase after treatment with IRAK1 inhibitors or ATMi are more likely to be highly transcribed-Heatmap.pdf",
       width=4,height=5)
