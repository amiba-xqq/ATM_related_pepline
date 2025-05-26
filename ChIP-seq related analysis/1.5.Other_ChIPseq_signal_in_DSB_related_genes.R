rm(list=ls())
library(dplyr)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
load(file="DSB.Rdata")
IRAK1_symbol <- IRAK1_table[IRAK1_table$ENSEMBL %in% IRAK1_gene,]$SYMBOL
others_symbol <- others_table[others_table$ENSEMBL %in% others_gene,]$SYMBOL

####导入其它ChIPseq在基因TSS±5kb范围的counts矩阵####
####Import ChIP-seq counts matrices within ±5 kb around gene TSS####
counts <- read.table("../TSS_ALL_seq_counts.txt",header=F)
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
z_scores2 <- z_scores_df[rownames(z_scores_df) %in% IRAK1_symbol | rownames(z_scores_df) %in% others_symbol, ]
z_scores2$group <- ifelse(rownames(z_scores2)%in% IRAK1_symbol, "IRAK1_regulated",
                        ifelse(rownames(z_scores2) %in% others_symbol,"others","NA"))
table(z_scores2$group)
colnames(z_scores2)
z_scores2 <- arrange(z_scores2,group,desc(EUseq))

median(z_scores2[z_scores2$group=="IRAK1_regulated",]$PolII_ChIP)
median(z_scores2[z_scores2$group=="others",]$PolII_ChIP)
wilcox.test(z_scores2[z_scores2$group=="IRAK1_regulated",]$PolII_ChIP,
            z_scores2[z_scores2$group=="others",]$PolII_ChIP)

n <- data.frame()
i=1
repeat{
  n[i,1] <- median(z_scores2[z_scores2$group=="IRAK1_regulated",i])
  n[i,2] <- median(z_scores2[z_scores2$group=="others",i])
  rownames(n)[i] <- colnames(z_scores2)[i]
  i=i+1
  if(i>(ncol(z_scores2)-1)){
    break
  }
}
colnames(n) <- c("IRAK1_regulated","others")
n <- n[order(n$IRAK1_regulated,decreasing = T),]

n2 <- n[c("PolII_ChIP","H3K4me1_ChIP","EUseq","H3K36me3_ChIP",
          "H3K27ac_ChIP","HP1g_ChIP","EZH2_ChIP","H3K27me3_ChIP","H3K9me3_ChIP"),]
my_palette <- colorRampPalette(c("#4575B4", "white", "#D73027"))

p1 <- pheatmap::pheatmap(n2,show_colnames =T,show_rownames = T,
                   fontsize=7,
                   #scale = "row",
                   angle_col=90,
                   color = my_palette(100),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   breaks = seq(-0.2,0.6,0.008),
                   fontsize_col = 12,             
                   fontsize_row = 12)
ggsave(p1,filename = "IRAK1调控的DSB倾向于在转录活跃区-热图.pdf",
       width=5,height=8)
