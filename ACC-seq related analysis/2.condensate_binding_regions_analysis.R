rm(list=ls())
library(tidyverse)
####分析DMSO组的condensate binding regions####
####Analyze condensate binding regions in the DMSO group####
DMSO <- read.table("../DMSO_groups_in_peaks_5kb.txt")
colnames(DMSO) <- c("chr","start","end","DMSO_native","DMSO_fix","DMSO_fixHex")

# 将counts数归一化为CPM值，188515059是DMSO_native组的ACC-seq的成功比对到hg38基因组的总counts数，其它以此类推
# Normalize the count values to CPM value, where 188,515,059 is the total ACC-seq counts for the DMSO_native group that successfully mapped to the hg38 genome, and similar calculations are done for other groups.
DMSO$DMSO_native <- DMSO$DMSO_native/188515059*1e6
DMSO$DMSO_fix <- DMSO$DMSO_fix/234117802*1e6
DMSO$DMSO_fixHex <- DMSO$DMSO_fixHex/172421224*1e6

DMSO$log2FC_nativevsfix <- log2(DMSO$DMSO_native/DMSO$DMSO_fix)
DMSO$log2FC_fixHexvsfix <- log2(DMSO$DMSO_fixHex/DMSO$DMSO_fix)
table(DMSO$log2FC_nativevsfix>log2(1.5)&
        DMSO$log2FC_fixHexvsfix>log2(1.5))

DMSO_U <- DMSO[DMSO$log2FC_nativevsfix>log2(1.5)&
                 DMSO$log2FC_fixHexvsfix>log2(1.5),]

DMSO_U_bed <- DMSO_U[,1:3]
DMSO_U_bed$start <- DMSO_U_bed$start+2500
DMSO_U_bed$end <- DMSO_U_bed$end-2500
write.table(DMSO_U_bed,"DMSO_Ushape.bed",col.names = F,row.names = F,quote = F)

####分析pOHT组的condensate binding regions####
####Analyze condensate binding regions in the pOHT group####
pOHT <- read.table("../pOHT_groups_in_peaks_5kb.txt")
colnames(pOHT) <- c("chr","start","end","pOHT_native","pOHT_fix","pOHT_fixHex")
pOHT$pOHT_native <- pOHT$pOHT_native/206564070*1e6
pOHT$pOHT_fix <- pOHT$pOHT_fix/173504623*1e6
pOHT$pOHT_fixHex <- pOHT$pOHT_fixHex/167656085*1e6

pOHT$log2FC_nativevsfix <- log2(pOHT$pOHT_native/pOHT$pOHT_fix)
pOHT$log2FC_fixHexvsfix <- log2(pOHT$pOHT_fixHex/pOHT$pOHT_fix)
table(pOHT$log2FC_nativevsfix>log2(1.5)&
        pOHT$log2FC_fixHexvsfix>log2(1.5))

pOHT_U <- pOHT[pOHT$log2FC_nativevsfix>log2(1.5)&
                 pOHT$log2FC_fixHexvsfix>log2(1.5),]

pOHT_U_bed <- pOHT_U[,1:3]
pOHT_U_bed$start <- pOHT_U_bed$start+2500
pOHT_U_bed$end <- pOHT_U_bed$end-2500
write.table(pOHT_U_bed,"pOHT_Ushape.bed",col.names = F,row.names = F,quote = F)

####分析pOHTIRAK1i组的condensate binding regions####
####Analyze condensate binding regions in the pOHTIRAK1i group####
pOHTIRAK1i <- read.table("../pOHTIRAK1i_groups_in_peaks_5kb.txt")
colnames(pOHTIRAK1i) <- c("chr","start","end","pOHTIRAK1i_native","pOHTIRAK1i_fix","pOHTIRAK1i_fixHex")
pOHTIRAK1i$pOHTIRAK1i_native <- pOHTIRAK1i$pOHTIRAK1i_native/175787469*1e6
pOHTIRAK1i$pOHTIRAK1i_fix <- pOHTIRAK1i$pOHTIRAK1i_fix/165602668*1e6
pOHTIRAK1i$pOHTIRAK1i_fixHex <- pOHTIRAK1i$pOHTIRAK1i_fixHex/139949413*1e6

pOHTIRAK1i$log2FC_nativevsfix <- log2(pOHTIRAK1i$pOHTIRAK1i_native/pOHTIRAK1i$pOHTIRAK1i_fix)
pOHTIRAK1i$log2FC_fixHexvsfix <- log2(pOHTIRAK1i$pOHTIRAK1i_fixHex/pOHTIRAK1i$pOHTIRAK1i_fix)
table(pOHTIRAK1i$log2FC_nativevsfix>log2(1.5)&
        pOHTIRAK1i$log2FC_fixHexvsfix>log2(1.5))

pOHTIRAK1i_U <- pOHTIRAK1i[pOHTIRAK1i$log2FC_nativevsfix>log2(1.5)&
                             pOHTIRAK1i$log2FC_fixHexvsfix>log2(1.5),]

pOHTIRAK1i_U_bed <- pOHTIRAK1i_U[,1:3]
pOHTIRAK1i_U_bed$start <- pOHTIRAK1i_U_bed$start+2500
pOHTIRAK1i_U_bed$end <- pOHTIRAK1i_U_bed$end-2500
write.table(pOHTIRAK1i_U_bed,"pOHTIRAK1i_Ushape.bed",col.names = F,row.names = F,quote = F)

# 最后得到的DMSO_Ushape.bed、pOHT_Ushape.bed、pOHTIRAK1i_Ushape.bed就是各组的condensate binding regions对应的BED文件
# Finally obtained DMSO_Ushape.bed, pOHT_Ushape.bed, and pOHTIRAK1i_Ushape.bed are the BED files corresponding to the condensate binding regions for each group.
