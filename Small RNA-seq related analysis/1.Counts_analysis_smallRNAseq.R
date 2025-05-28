rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyverse)
library(DESeq2)
library(FactoMineR)
library(factoextra)
library(edgeR)
library(statmod)
####读取PartII library组的all gene的counts，并转化成CPM####
####Read the counts of all genes in the PartII library group and convert them to CPM####
smallRNAB <- read.table("smallRNA_big_clear_all_gene.txt")
colnames(smallRNAB) <- c("chr","start","end","ensembl","strand","DMSO_r1",
                       "DMSO_r2","DMSO_r3","IRAK1i_r1","IRAK1i_r2",
                       "IRAK1i_r3","IRDMSO_r1","IRDMSO_r2",
                       "IRIRAK1i_r1","IRIRAK1i_r2","IRIRAK1i_r3")
smallRNAB$DMSO_r1 <- smallRNAB$DMSO_r1/68634559*1000000
smallRNAB$DMSO_r2 <- smallRNAB$DMSO_r2/56050541*1000000
smallRNAB$DMSO_r3 <- smallRNAB$DMSO_r3/45633381*1000000
smallRNAB$IRAK1i_r1 <- smallRNAB$IRAK1i_r1/55175288*1000000
smallRNAB$IRAK1i_r2 <- smallRNAB$IRAK1i_r2/67834650*1000000
smallRNAB$IRAK1i_r3 <- smallRNAB$IRAK1i_r3/57328017*1000000
smallRNAB$IRDMSO_r1 <- smallRNAB$IRDMSO_r1/47349462*1000000
smallRNAB$IRDMSO_r2 <- smallRNAB$IRDMSO_r2/31372143*1000000
smallRNAB$IRIRAK1i_r1 <- smallRNAB$IRIRAK1i_r1/43455415*1000000
smallRNAB$IRIRAK1i_r2 <- smallRNAB$IRIRAK1i_r2/46216991*1000000
smallRNAB$IRIRAK1i_r3 <- smallRNAB$IRIRAK1i_r3/33410037*1000000

####读取PartI library组的all gene的counts，并转化成CPM####
####Read the counts of all genes in the PartI library group and convert them to CPM####
smallRNAG <- read.table("smallRNA_great_clear_all_gene.txt")
colnames(smallRNAG) <- c("chr","start","end","ensembl","strand","DMSO_r1",
                         "DMSO_r2","DMSO_r3","IRAK1i_r1","IRAK1i_r2",
                         "IRAK1i_r3","IRDMSO_r1","IRDMSO_r2",
                         "IRIRAK1i_r1","IRIRAK1i_r2","IRIRAK1i_r3")
smallRNAG$DMSO_r1 <- smallRNAG$DMSO_r1/27406531*1000000
smallRNAG$DMSO_r2 <- smallRNAG$DMSO_r2/27199452*1000000
smallRNAG$DMSO_r3 <- smallRNAG$DMSO_r3/26956778*1000000
smallRNAG$IRAK1i_r1 <- smallRNAG$IRAK1i_r1/28463463*1000000
smallRNAG$IRAK1i_r2 <- smallRNAG$IRAK1i_r2/33018732*1000000
smallRNAG$IRAK1i_r3 <- smallRNAG$IRAK1i_r3/25054229*1000000
smallRNAG$IRDMSO_r1 <- smallRNAG$IRDMSO_r1/31158649*1000000
smallRNAG$IRDMSO_r2 <- smallRNAG$IRDMSO_r2/24682784*1000000
smallRNAG$IRIRAK1i_r1 <- smallRNAG$IRIRAK1i_r1/19452707*1000000
smallRNAG$IRIRAK1i_r2 <- smallRNAG$IRIRAK1i_r2/20318195*1000000
smallRNAG$IRIRAK1i_r3 <- smallRNAG$IRIRAK1i_r3/17284124*1000000

####合并两组的CPM值####
####Merge the CPM values of two groups####
a <- (smallRNAB[,6:16]+smallRNAG[,6:16])/2
b <- smallRNAB[,1:5]
smallRNA_CPM <- cbind(b,a)
table(rowSums(smallRNA_CPM[,6:16]>0)>=2)

# 删去所有组CPM都为0的组
# Remove genes where all groups have CPM = 0.
smallRNA_CPM_clear <- smallRNA_CPM[rowSums(smallRNA_CPM[,6:16]>0)>=2,] 

# 删去IRDMSO/IRIRAK1i组CPM都为0的组
# Remove genes where both IRDMSO and IRIRAK1i groups have CPM = 0.
smallRNA_CPM_clear2 <- smallRNA_CPM[rowSums(smallRNA_CPM[,12:16]>0)>=2,c(1:5,12:16)]

####计算IR+IRAK1i vs IR+DMSO的Fold change####
####Calculate the Fold change of IR+IRAK1i vs IR+DMSO####
a <- min(smallRNA_CPM_clear2[smallRNA_CPM_clear2$IRDMSO_r1!=0,"IRDMSO_r1"])
smallRNA_CPM_clear2[smallRNA_CPM_clear2$IRDMSO_r1==0,"IRDMSO_r1"] <- a/2
a <- min(smallRNA_CPM_clear2[smallRNA_CPM_clear2$IRDMSO_r2!=0,"IRDMSO_r2"])
smallRNA_CPM_clear2[smallRNA_CPM_clear2$IRDMSO_r2==0,"IRDMSO_r2"] <- a/2
a <- min(smallRNA_CPM_clear2[smallRNA_CPM_clear2$IRIRAK1i_r1!=0,"IRIRAK1i_r1"])
smallRNA_CPM_clear2[smallRNA_CPM_clear2$IRIRAK1i_r1==0,"IRIRAK1i_r1"] <- a/2
a <- min(smallRNA_CPM_clear2[smallRNA_CPM_clear2$IRIRAK1i_r2!=0,"IRIRAK1i_r2"])
smallRNA_CPM_clear2[smallRNA_CPM_clear2$IRIRAK1i_r2==0,"IRIRAK1i_r2"] <- a/2
a <- min(smallRNA_CPM_clear2[smallRNA_CPM_clear2$IRIRAK1i_r3!=0,"IRIRAK1i_r3"])
smallRNA_CPM_clear2[smallRNA_CPM_clear2$IRIRAK1i_r3==0,"IRIRAK1i_r3"] <- a/2

smallRNA_CPM_clear2$FC_IRIRAK1ivsIRDMSO <- rowMeans(smallRNA_CPM_clear2[,8:10])/rowMeans(smallRNA_CPM_clear2[,6:7])
smallRNA_CPM_clear2$log2FC_IRIRAK1ivsIRDMSO <- log2(smallRNA_CPM_clear2$FC_IRIRAK1ivsIRDMSO)

####使用edgeR计算差异表达基因####
####Use edgeR to calculate differentially expressed genes(DEGs)####
# 设定 实验组exp / 对照组ctr
# Set experimental group (exp) and control group (ctr).
exp="IR_IRAK1i"
ctr="IR_DMSO"
name_list <- c("IR_DMSO","IR_DMSO",
               "IR_IRAK1i","IR_IRAK1i","IR_IRAK1i")
nlgl <- data.frame(row.names=colnames(smallRNA_CPM_clear2[,6:10]),
                   name_list=name_list,
                   group_list=name_list)
group_list <- nlgl$group_list
CPM <- smallRNA_CPM_clear2[,6:10]
rownames(CPM) <- smallRNA_CPM_clear2$ensembl

group <- factor(group_list)
group <- relevel(group,ctr)     
design <- model.matrix(~0+group)
rownames(design) <- rownames(nlgl)
colnames(design) <- levels(group)

# 表达矩阵DGEList构建与过滤低表达基因
# Construction of DGEList and filtering out lowly expressed genes.
dge <- DGEList(counts=CPM, group=group)
#keep.exprs <- filterByExpr(dge) 
#dge <- dge[keep.exprs,,keep.lib.sizes=FALSE] 
dge <- calcNormFactors(dge, method = 'TMM') 
dge <- estimateDisp(dge, design, robust=T) 

# To perform quasi-likelihood(QL) F-test tests:  bulk RNA-seq 
fit <- glmQLFit(dge, design, robust=T)  
lt <- glmQLFTest(fit, contrast=c(-1,1))

tempDEG <- topTags(lt, n = Inf) #sort by PValue, n is the number of genes/tags to return
tempDEG <- as.data.frame(tempDEG)
DEG_edgeR <- na.omit(tempDEG)
DEG_edgeR$ensembl <- rownames(DEG_edgeR)

#### 导出IR+IRAK1i vs IR下调基因的gmt文件####
####Export GMT file of downregulated genes in IR+IRAK1i vs IR####
smallRNA_IRIRAK1ivsIR.gmx <- as.data.frame(DEG_edgeR[DEG_edgeR$logFC < -log2(1.5)&
                                   DEG_edgeR$PValue<0.05,"ensembl"])
colnames(smallRNA_IRIRAK1ivsIR.gmx) <- "SYMBOL"
smallRNA_IRIRAK1ivsIR.gmx <- rbind(data.frame(SYMBOL = "smallRNA_IRIRAK1ivsIR"), smallRNA_IRIRAK1ivsIR.gmx)
smallRNA_IRIRAK1ivsIR.gmx <- rbind(data.frame(SYMBOL = "smallRNA_IRIRAK1ivsIR"), smallRNA_IRIRAK1ivsIR.gmx)
smallRNA_IRIRAK1ivsIR.gmt <- t(smallRNA_IRIRAK1ivsIR.gmx)
write.table(smallRNA_IRIRAK1ivsIR.gmt,"smallRNA_IRIRAK1ivsIR.gmt",sep="\t",col.names = F,row.names = F,quote = F)
