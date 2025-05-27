rm(list=ls())
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(rtracklayer)

####读取GTF文件的信息####
#####Read information from the GTF file####
gtf_data <- import('Homo_sapiens.GRCh38.108.gtf') # GTF file is downloaded from ENSEMBL
gtf_df <- as.data.frame(gtf_data)
gtf_id <- gtf_df[,c("gene_id","gene_name","gene_biotype","width")]
colnames(gtf_id)[1] <- "ensembl"
gtf_id <- gtf_id[!duplicated(gtf_id$ensembl),]

####读取NBS1 RIP的all gene的counts，并转化成CPM####
####Read the gene expression counts for all genes in NBS1 RIP and convert them to CPM####
NBS1RIP <- read.table("NBS1_RIP_all_gene_merged.txt")
colnames(NBS1RIP) <- c("chr","start","end","SYMBOL","strand",
                       "mock","noIR","IR","IRIRAK1i")
colSums(NBS1RIP[,6:9])
NBS1RIP$mock <- NBS1RIP$mock/49948596*1e6
NBS1RIP$noIR <- NBS1RIP$noIR/82615987*1e6
NBS1RIP$IR <- NBS1RIP$IR/157551543*1e6
NBS1RIP$IRIRAK1i <- NBS1RIP$IRIRAK1i/94012891*1e6

####计算NBS1 RIP的IRIRAK1i vs IR的Fold change####
####Calculate the Fold change of NBS1 RIP (IRIRAK1i versus IR)####
NBS1RIP_clear <- NBS1RIP[,c(1:4,8:9)]
NBS1RIP_clear <- NBS1RIP_clear[NBS1RIP_clear$IR>0 | NBS1RIP_clear$IRIRAK1i>0,]

a <- min(NBS1RIP_clear[NBS1RIP_clear$IR!=0,"IR"])
NBS1RIP_clear[NBS1RIP_clear$IR==0,"IR"] <- a/2
a <- min(NBS1RIP_clear[NBS1RIP_clear$IRIRAK1i!=0,"IRIRAK1i"])
NBS1RIP_clear[NBS1RIP_clear$IRIRAK1i==0,"IRIRAK1i"] <- a/2

NBS1RIP_clear$FC_IRIRAK1ivsIR_RIP <- NBS1RIP_clear$IRIRAK1i/NBS1RIP_clear$IR
NBS1RIP_clear$log2_FC_IRIRAK1ivsIR_RIP <- log2(NBS1RIP_clear$FC_IRIRAK1ivsIR_RIP)

####编写导出gmt文件的function####
####Write a function to export the GMT file####
gmt <- function(a,k){
  b <- as.data.frame(a[,"ensembl"])
  colnames(b) <- "SYMBOL"
  b <- rbind(data.frame(SYMBOL = k), b)
  b <- rbind(data.frame(SYMBOL = k), b)
  c <- t(b)
  write.table(c,paste0(k,".gmt"),sep="\t",col.names = F,row.names = F,quote = F)
  return(0)
}

####分析NBS1-RIPseq的IRIRAK1ivsIR的差异基因，导出GSEA分析的rnk和gmt，用于分析各类型RNA的相关性####
####Analyze NBS1-RIPseq data comparing IRIRAK1 vs IR; export GSEA rnk and gmt files; analyze correlations between RNA types and IRAK1i induced NBS1-binding RNA changes####
NBS1RIP_IRvsIRIRAK1i.rnk <- NBS1RIP_clear[,c(4,8)]
# 取|FC|>1.5作为IRvsIR+IRAK1i富集差异的Gene
# Filter genes with |FC| > 1.5 for IR vs IR+IRAK1i enrichment differences
NBS1RIP_IRvsIRIRAK1i.rnk <- NBS1RIP_IRvsIRIRAK1i.rnk[NBS1RIP_IRvsIRIRAK1i.rnk$log2_FC_IRIRAK1ivsIR_RIP< -1,]
NBS1RIP_IRvsIRIRAK1i.rnk$log2_FC_IRIRAK1ivsIR_RIP <- -NBS1RIP_IRvsIRIRAK1i.rnk$log2_FC_IRIRAK1ivsIR_RIP
colnames(NBS1RIP_IRvsIRIRAK1i.rnk)[2] <- "log2_FC_IRvsIRIRAK1i_RIP"
NBS1RIP_IRvsIRIRAK1i.rnk <- NBS1RIP_IRvsIRIRAK1i.rnk[order(NBS1RIP_IRvsIRIRAK1i.rnk$log2_FC_IRvsIRIRAK1i_RIP,decreasing = T),]
write.table(NBS1RIP_IRvsIRIRAK1i.rnk,"NBS1RIP_IRvsIRIRAK1i.rnk",sep = "\t",col.names = F,row.names = F,quote = F)

colnames(NBS1RIP_clear)[4] <- "ensembl"
NBS1RIP_gtf <- inner_join(NBS1RIP_clear,gtf_id,by="ensembl") 
table(gtf_id$gene_biotype)
NBS1RIP_gtf2 <- NBS1RIP_gtf[NBS1RIP_gtf$log2_FC_IRIRAK1ivsIR_RIP< -1,]
a <- NBS1RIP_gtf2[NBS1RIP_gtf2$gene_biotype %in% c("snoRNA","scaRNA"),]
gmt(a,"NBS1RIP_snoRNA")

a <- NBS1RIP_gtf2[NBS1RIP_gtf2$gene_biotype %in% c("snRNA"),]
gmt(a,"NBS1RIP_snRNA")

a <- NBS1RIP_gtf2[NBS1RIP_gtf2$gene_biotype %in% c("miRNA"),]
gmt(a,"NBS1RIP_miRNA")

a <- NBS1RIP_gtf2[NBS1RIP_gtf2$gene_biotype %in% c("protein_coding"),]
gmt(a,"NBS1RIP_mRNA")

a <- NBS1RIP_gtf2[NBS1RIP_gtf2$gene_biotype %in% c("lncRNA"),]
gmt(a,"NBS1RIP_lncRNA")

a <- NBS1RIP_gtf2[grepl("pseudogene", NBS1RIP_gtf2$gene_biotype),]
gmt(a,"NBS1RIP_pseudogene")

