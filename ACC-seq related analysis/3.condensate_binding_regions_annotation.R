rm(list=ls())
library(ChIPseeker)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
library(dplyr)
library(Vennerable)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(enrichplot)
library(pathview)
library(DOSE)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
dir.create("results")

####读取各组condensate binding regions的bed文件####
####Read each group's BED file of condensate binding regions####
DMSO_U <- readPeakFile("DMSO_Ushape.bed")
pOHT_U <- readPeakFile("pOHT_Ushape.bed")
pOHTIRAK1i_U <- readPeakFile("pOHTIRAK1i_Ushape.bed")

Ushape <- list(DMSO_U=DMSO_U,
               pOHT_U=pOHT_U,
               pOHTIRAK1i_U=pOHTIRAK1i_U)
peakAnnoList <- lapply(Ushape, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 1000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")

####基因region注释####
####Gene region annotation####
pdf("./results/ACCseq U_shape regions peakAnnoList.pdf",width=8,height=4)
plotAnnoBar(peakAnnoList)
dev.off()

####转化为table####
####Convert to table####
DMSO_U_table <- as.data.frame(peakAnnoList[[1]])
pOHT_U_table <- as.data.frame(peakAnnoList[[2]])
pOHTIRAK1i_U_table <- as.data.frame(peakAnnoList[[3]])
write.csv(DMSO_U_table,"./results/DMSO U shape peaks大表.csv")
write.csv(pOHT_U_table,"./results/pOHT U shape peaks大表.csv")
write.csv(pOHTIRAK1i_U_table,"./results/pOHTIRAK1i U shape peaks大表.csv")

DMSO_U_gene <- DMSO_U_table[DMSO_U_table$annotation!="Distal Intergenic"&
                              !duplicated(DMSO_U_table$ENSEMBL)&
                              !is.na(DMSO_U_table$ENSEMBL),"ENSEMBL"]

pOHT_U_gene <- pOHT_U_table[pOHT_U_table$annotation!="Distal Intergenic"&
                              !duplicated(pOHT_U_table$ENSEMBL)&
                              !is.na(pOHT_U_table$ENSEMBL),"ENSEMBL"]

pOHTIRAK1i_U_gene <- pOHTIRAK1i_U_table[pOHTIRAK1i_U_table$annotation!="Distal Intergenic"&
                              !duplicated(pOHTIRAK1i_U_table$ENSEMBL)&
                              !is.na(pOHTIRAK1i_U_table$ENSEMBL),"ENSEMBL"]
save(DMSO_U_table,pOHT_U_table,pOHTIRAK1i_U_table,
     DMSO_U_gene,pOHT_U_gene,pOHTIRAK1i_U_gene,file="ACCseq_Ushape_regions_gene.Rdata")

####gene群的venn图####
####Venn plot of gene clusters####
genes <- list(DMSO_U_gene=DMSO_U_gene,
              pOHT_U_gene=pOHT_U_gene,
              pOHTIRAK1i_U_gene=pOHTIRAK1i_U_gene)
pdf("./results/ACCseq U shape peaks注释gene的venn图(DMSO and pOHT).pdf",width=7,height=5)
vennplot(genes[1:2], by='Vennerable')
dev.off()

####通路富集####
####Pathway enrichment####
DMSO_U_geneid <-  DMSO_U_table[DMSO_U_table$annotation!="Distal Intergenic"&
                                 !duplicated(DMSO_U_table$ENSEMBL)&
                                 !is.na(DMSO_U_table$ENSEMBL),"geneId"]

pOHT_U_geneid <-  pOHT_U_table[pOHT_U_table$annotation!="Distal Intergenic"&
                                 !duplicated(pOHT_U_table$ENSEMBL)&
                                 !is.na(pOHT_U_table$ENSEMBL),"geneId"]

pOHTIRAK1i_U_geneid <-  pOHTIRAK1i_U_table[pOHTIRAK1i_U_table$annotation!="Distal Intergenic"&
                                 !duplicated(pOHTIRAK1i_U_table$ENSEMBL)&
                                 !is.na(pOHTIRAK1i_U_table$ENSEMBL),"geneId"]

geneIDs <- list(DMSO_U=DMSO_U_geneid,
                pOHT_U=pOHT_U_geneid,
                pOHTIRAK1i_U=pOHTIRAK1i_U_geneid)

cc = compareCluster(geneCluster = geneIDs, 
                    fun="enrichGO", 
                    OrgDb="org.Hs.eg.db", 
                    ont= "BP",
                    pvalueCutoff=0.2,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.2)

pdf("./results/ACCseq U shape peaks注释gene的GO-BP通路富集.pdf",width=5,height=5)
enrichplot::dotplot(cc,showCategory=7,includeAll=TRUE,font.size=7)
dev.off()

kegg <- compareCluster(geneCluster = geneIDs, 
                       fun="enrichKEGG", 
                       organism="hsa",
                       pvalueCutoff=0.2,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2)

pdf("./results/ACCseq U shape peaks注释gene的KEGG通路富集.pdf",width=5,height=5)
enrichplot::dotplot(kegg,showCategory=10,includeAll=TRUE,font.size=10)
dev.off()

