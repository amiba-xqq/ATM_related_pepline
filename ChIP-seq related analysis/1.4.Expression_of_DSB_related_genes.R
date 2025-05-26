rm(list=ls())
library(tidyverse)
library(ggpubr)
library("rtracklayer")
library("SummarizedExperiment")

####读取U2OS-counts矩阵####
####Read U2OS-counts matirx#####
counts_U2OS <- read.table("counts_U2OS.txt",header = T)
counts_U2OS <- counts_U2OS[,c(1,6:9)]
# "U2OS1","U2OS2","U2OS3"指3个重复组
# "U2OS1","U2OS2","U2OS3" represent three replicate groups
colnames(counts_U2OS)[3:5] <- c("U2OS1","U2OS2","U2OS3")
counts_U2OS$U2OS <- rowMeans(counts_U2OS[3:5])
counts_all <- counts_U2OS[,c(1,2,6)]
geneid_efflen <- subset(counts_all,select = c("Geneid","Length"))
colnames(geneid_efflen) <- c("geneid","efflen")  
counts <- data.frame(U2OS=counts_all$U2OS)
rownames(counts) <- counts_all[,1]
efflen <- geneid_efflen[match(rownames(counts),
                              geneid_efflen$geneid),
                        "efflen"]

####计算TPM####
####Calculate TPM(Transcripts Per Kilobase Million) value####
counts2TPM <- function(count=count, efflength=efflen){
  RPK <- count/(efflength/1000)   
  PMSC_rpk <- sum(RPK)/1e6        
  RPK/PMSC_rpk              
}  
tpm <- as.data.frame(apply(counts,2,counts2TPM))
colSums(tpm) 

####筛选出至少在重复样本数量内的表达量counts大于1的行（基因）####
####Filter out rows (genes) where the expression count is greater than 1 in at least one replicate sample####
keep_feature <- rowSums(counts>1)==1
table(keep_feature) 
counts_filt <- data.frame(U2OS=counts[keep_feature, ]) 
rownames(counts_filt) <- rownames(counts)[keep_feature]
tpm_filt <- data.frame(U2OS=tpm[keep_feature, ])
rownames(tpm_filt) <- rownames(tpm)[keep_feature]
tpm_filt$gene_id <- rownames(tpm_filt)

####分析各组DSB相关基因在正常U2OS细胞中的表达水平####
####Analyze the expression levels of DSB-related genes in normal U2OS cells across different groups####
load(file="DSB.Rdata")
a1 <- tpm_filt[rownames(tpm_filt) %in% total_gene,c("U2OS","gene_id")]
a1 <- data.frame(U2OS=a1$U2OS,
                 group="Total")
a2 <- tpm_filt[rownames(tpm_filt) %in% IRAK1_gene,c("U2OS","gene_id")]
a2 <- data.frame(U2OS=a2$U2OS,
                 group="IRAK1_regulated")
a3 <- tpm_filt[rownames(tpm_filt) %in% others_gene & !(rownames(tpm_filt) %in% IRAK1_gene)
               ,c("U2OS","gene_id")]
a3 <- data.frame(U2OS=a3$U2OS,
                 group="others")

a <- rbind(a3,a2)
a$U2OS <- log2(a$U2OS)

my_comparisons <- list( c("others", "IRAK1_regulated"))
ggboxplot(a,
          x="group",
          y="U2OS",
          fill="group",
          xlab= "",
          ylab= "TPM(log2) in U2OS cell",
          legend.title="",
          palette = c("#00B9E3","#F8766D"))+
  theme(axis.text.x = element_text(size = 14, angle = 270, vjust = 1, hjust = 0, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14),
        legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
ggsave("IRAK1调控DSB位点倾向于转录活跃基因.pdf",width=2.5,height=6)

