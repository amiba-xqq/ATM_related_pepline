rm(list=ls())
library(ChIPseeker)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(pathview)
library(DOSE)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
####读取各组DSB bed文件####
####Read bed files of DSBs in different groups####
total <- readPeakFile("./bed/DSB.bed")
IRAK1 <- readPeakFile("./bed/IRAK1_regulated_DSB.bed")
others <- readPeakFile("./bed/other_DSB.bed")

a <- list(total=total,
          IRAK1_regulated=IRAK1,
          others=othersI)
peakAnnoList <- lapply(a, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 1000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")

####提取注释信息####
####Extract Annotation Information####
total_table <- as.data.frame(peakAnnoList[[1]])
IRAK1_table <- as.data.frame(peakAnnoList[[2]])
others_table <- as.data.frame(peakAnnoList[[3]])

# IRAK1_gene是IRAK1调控的ATM激活的DSB位点所在的基因；others_gene是IRAK1不依赖的ATM激活的DSB位点所在的基因
# IRAK1_gene: Genes where ATM activation of DSBs is regulated by IRAK1 (IRAK1-dependent DSB sites);others_gene: Genes where ATM activation of DSBs is independent of IRAK1 (IRAK1-independent DSB sites).
total_gene <- total_table[total_table$annotation!="Distal Intergenic"&
                             !duplicated(total_table$ENSEMBL)&
                             !is.na(total_table$ENSEMBL),"ENSEMBL"]
IRAK1_gene <- IRAK1_table[IRAK1_table$annotation!="Distal Intergenic"&
                            !duplicated(IRAK1_table$ENSEMBL)&
                            !is.na(IRAK1_table$ENSEMBL),"ENSEMBL"]
others_gene <- others_table[others_table$annotation!="Distal Intergenic"&
                            !duplicated(others_table$ENSEMBL)&
                            !is.na(others_table$ENSEMBL),"ENSEMBL"]
pNBS1 <- cbind(pNBS1,total_table[,c("ENSEMBL","annotation")])
pNBS1 <- pNBS1[!is.na(pNBS1$ENSEMBL),]
save(pNBS1,total_table,IRAK1_table,others_table,total_gene,IRAK1_gene,others_gene,
     file="DSB.Rdata")






