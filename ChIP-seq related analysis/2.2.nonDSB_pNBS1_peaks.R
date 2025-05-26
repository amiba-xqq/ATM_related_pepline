rm(list=ls())
library(tidyverse)
pNBS1 <- read.table("pNBS1_counts_in_total_pNBS1_unDSB_peaks_2kb.txt")
colnames(pNBS1) <- c("chr","start","end","DMSO","DSB","DSBIRAK1i")
# 将counts数归一化为CPM值，10949187是pNBS1-DMSO组的ChIP-seq的成功比对到基因组的总counts数，其它以此类推
# Normalize the count values to CPM value, where 10,949,187 is the total ChIP-seq counts for the pNBS1-DMSO group that successfully mapped to the hg38 genome, and similar calculations are done for other groups.
pNBS1$DMSO <- pNBS1$DMSO/10949187*1e6
pNBS1$DSB <- pNBS1$DSB/24466496*1e6
pNBS1$DSBIRAK1i <- pNBS1$DSBIRAK1i/31576054*1e6

pNBS1[pNBS1$DMSO==0,]$DMSO <- min(pNBS1[pNBS1$DMSO!=0,]$DMSO)/2
pNBS1[pNBS1$DSBIRAK1i==0,]$DSBIRAK1i <- min(pNBS1[pNBS1$DSBIRAK1i!=0,]$DSBIRAK1i)/2

pNBS1$log2FC_DSBvsDMSO <- log2(pNBS1$DSB/pNBS1$DMSO)
pNBS1$log2FC_DSBIRAK1ivsDSB <- log2(pNBS1$DSBIRAK1i/pNBS1$DSB)

a <- pNBS1[pNBS1$log2FC_DSBvsDMSO>log2(1.5),]
table(a$log2FC_DSBIRAK1ivsDSB< -log2(1.5))

a1 <- a[a$log2FC_DSBIRAK1ivsDSB < -log2(1.5),]
a2 <- a[a$log2FC_DSBIRAK1ivsDSB >= 0,]

# pNBS1_noDSB_peaks.bed即记录了non-DSB pNBS1 peaks的BED文件
# 'pNBS1_noDSB_peaks.bed' file records the positions of non-DSB pNBS1 peaks.
write.table(a[,1:3],"pNBS1_noDSB_peaks.bed",col.names = F,row.names = F,sep="\t",quote = F)
write.table(a1[,1:3],"IRAK1_regulated_pNBS1_noDSB_peaks.bed",col.names = F,row.names = F,sep="\t",quote = F)
write.table(a2[,1:3],"others_pNBS1_noDSB_peaks.bed",col.names = F,row.names = F,sep="\t",quote = F)
save(pNBS1,file="pNBS_noDSB_peaks.Rdata")





