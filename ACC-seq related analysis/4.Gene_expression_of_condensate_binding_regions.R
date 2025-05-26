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
  RPK <- count/(efflength/1000)   #每千碱基reads (Reads Per Kilobase) 长度标准化
  PMSC_rpk <- sum(RPK)/1e6        #RPK的每百万缩放因子 (“per million” scaling factor ) 深度标准化
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

####分析各组condensate binding regions相关基因在正常U2OS细胞中的表达水平####
####Analyze the expression levels of condensate binding regions-related genes in normal U2OS cells across different groups####
load(file="ACCseq_Ushape_regions_gene.Rdata")
a1 <- tpm_filt[rownames(tpm_filt) %in% DMSO_U_gene,c("U2OS","gene_id")]
a1 <- data.frame(U2OS=a1$U2OS,
                 group="DMSO")
a2 <- tpm_filt[rownames(tpm_filt) %in% pOHT_U_gene,c("U2OS","gene_id")]
a2 <- data.frame(U2OS=a2$U2OS,
                 group="pOHT")
a3 <- tpm_filt[rownames(tpm_filt) %in% pOHTIRAK1i_U_gene,c("U2OS","gene_id")]
a3 <- data.frame(U2OS=a3$U2OS,
                 group="pOHTIRAK1i")

a <- rbind(a1,a2,a3)
a$U2OS <- log2(a$U2OS)

ggboxplot(a,
          x="group",
          y="U2OS",
          fill="group",
          xlab= "",
          ylab= "TPM(log2) in U2OS cell",
          legend.title="",
          palette = c("#00B9E3","#F8766D","orange"))+
  theme(axis.text.x = element_text(size = 14, angle = 270, vjust = 1, hjust = 0, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14),
        legend.position = "none")
ggsave(filename="./results/tpm_boxplot_U2OS_merged.pdf",
       width = 6,height = 8)

test <- wilcox.test(a[a$group=="DMSO","U2OS"],
                    a[a$group=="pOHT","U2OS"])
test$p.value

test <- wilcox.test(a[a$group=="pOHTIRAK1i","U2OS"],
                    a[a$group=="pOHT","U2OS"])
test$p.value


