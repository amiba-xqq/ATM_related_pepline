rm(list=ls())
library(dplyr)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

####导入counts矩阵####
####Load counts matrix####
counts1 <- read.table("../RAPseq_and_pATM_ChIPseq_in_TSS_10kb.txt",header=F)
colnames(counts1) <- c("chromatin","upstrand","downstrand","gene_id","none","strand_direction",
                      "input1","DMSO_U1","DMSO_U2","DMSO_U6",
                      "input2","pOHT_U1","pOHT_U2","pOHT_U6",
                      "input3","pOHTIRAK1i_U1","pOHTIRAK1i_U2","pOHTIRAK1i_U6")
counts2 <- read.table("../TSS_histone_motification_counts.txt",header=F)
colnames(counts2) <- c("chromatin","upstrand","downstrand","gene_id","none","strand_direction",
                       "H3K4me1","H3K4me3","H3K9ac","H3K27ac","H3K27me3","H3K36me3","H3K9me3")
counts <- cbind(counts1,counts2[,7:ncol(counts2)])
counts <- counts[counts$input1 < 10000,]

####计算CPM####
####Calculate CPM value####
counts$input1 <- counts$input1/17502194*1e6
counts$DMSO_U1 <- counts$DMSO_U1/56648459*1e6
counts$DMSO_U2 <- counts$DMSO_U2/35649682*1e6
counts$DMSO_U6 <- counts$DMSO_U6/69741853*1e6
counts$input2 <- counts$input2/18890395*1e6
counts$pOHT_U1 <- counts$pOHT_U1/68145662*1e6
counts$pOHT_U2 <- counts$pOHT_U2/64166610*1e6
counts$pOHT_U6 <- counts$pOHT_U6/57188244*1e6
counts$input3 <- counts$input3/22112882*1e6
counts$pOHTIRAK1i_U1 <- counts$pOHTIRAK1i_U1/78395372*1e6
counts$pOHTIRAK1i_U2 <- counts$pOHTIRAK1i_U2/35607175*1e6
counts$pOHTIRAK1i_U6 <- counts$pOHTIRAK1i_U6/54407715*1e6

####去除U1/U2/U6 IP组CPM都为0的行####
####Remove rows where CPM of U1/U2/U6 IP group is 0####
table(rowSums(counts[,7:18]>0)>1)
counts_filt <- counts[rowSums(counts[,7:18]>0)>1,]

a <- counts_filt[,7:ncol(counts_filt)]
a <- min(a[a != 0])/2
counts_filt[,7:ncol(counts_filt)][counts_filt[,7:ncol(counts_filt)]==0] <- a

####转化成Z score####
####Transfer to Z score####
z_scores <- counts_filt[,c("DMSO_U1","DMSO_U2","DMSO_U6","pOHT_U1","pOHT_U2","pOHT_U6",
                           "pOHTIRAK1i_U1","pOHTIRAK1i_U2","pOHTIRAK1i_U6","H3K4me1",
                           "H3K4me3","H3K9ac","H3K27ac","H3K27me3","H3K36me3","H3K9me3")]
z_scores <- log2(z_scores)
z_scores <- scale(z_scores, center = TRUE, scale = TRUE)
z_scores_df <- as.data.frame(z_scores)

####计算相关系数,并做热图####
#### Calculate correlation coefficients and create heatmap####
cor_matrix <- cor(z_scores_df, method = "pearson")
cor_matrix_df <- as.data.frame(cor_matrix)

n2 <- cor_matrix_df[c("H3K4me1","H3K4me3","H3K9ac","H3K27ac","H3K27me3","H3K36me3","H3K9me3"),
                    c("DMSO_U1","DMSO_U2","DMSO_U6",
                      "pOHT_U1","pOHT_U2","pOHT_U6",
                      "pOHTIRAK1i_U1","pOHTIRAK1i_U2","pOHTIRAK1i_U6")]
n2 <- n2[order(n2$pOHT_U1,decreasing = T),]
my_palette <- colorRampPalette(c("#4575B4", "white", "#D73027"))
p1 <- pheatmap::pheatmap(n2,show_colnames =T,show_rownames = T,
                         fontsize=10,
                         #scale = "row",
                         angle_col=90,
                         color = my_palette(100),
                         cluster_rows = F,
                         cluster_cols = F,
                         breaks = seq(0.15,0.35,0.002))
ggsave(p1,filename = './results/heatmap_correlation between RAPseq and histone motification related ChIPseq.pdf',width = 8,height =6)
