rm(list = ls())
library(tidyverse)
library(TCGAbiolinks)
library(maftools)
library(TCGAmutations)
library(ggpubr)
load("TCGA_BRCA_TPM.Rdata")
gexp <- as.data.frame(t(tpm["ATM",]))
gexp$Tumor_Sample_Barcode <- gsub("\\.", "-", rownames(gexp))
ATM_High_ID <- gexp[gexp$ATM > quantile(gexp$ATM,0.75),]$Tumor_Sample_Barcode
ATM_Low_ID <- gexp[gexp$ATM < quantile(gexp$ATM,0.25),]$Tumor_Sample_Barcode
####读取ATM高表达肿瘤突变MAF文件####
query_SNV <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  barcode = ATM_High_ID,
  access = "open",
)
GDCdownload(query_SNV)
GDCprepare(query_SNV, save = T,save.filename = "TCGA-BRCA_ATM_high_SNP.Rdata") #运行时间非常长
load(file = "TCGA-BRCA_ATM_high_SNP.Rdata")
mafs <- data
maf_raw_high <- read.maf(mafs)
maf_table_high <- as.data.frame(maf_raw_high@data)
mafSummary_high <- mafSummary(maf_raw_high)

####读取ATM低表达肿瘤突变MAF文件####
query_SNV <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  barcode = ATM_Low_ID,
  access = "open",
)
GDCdownload(query_SNV)
GDCprepare(query_SNV, save = T,save.filename = "TCGA-BRCA_ATM_low_SNP.Rdata") #运行时间非常长
load(file = "TCGA-BRCA_ATM_low_SNP.Rdata")
mafs <- data
maf_raw_low <- read.maf(mafs)
maf_table_low <- as.data.frame(maf_raw_low@data)
mafSummary_low <- mafSummary(maf_raw_low)

#计算突变率
high_gene_mutation <- mafSummary_high[["gene.summary"]]
low_gene_mutation <- mafSummary_low[["gene.summary"]]
high_gene_mutation$freq_high <- high_gene_mutation$MutatedSamples/nrow(mafSummary_high[["variants.per.sample"]])*-100
low_gene_mutation$freq_low <- low_gene_mutation$MutatedSamples/nrow(mafSummary_low[["variants.per.sample"]])*100

gene_mutation <- full_join(high_gene_mutation[,c("Hugo_Symbol","freq_high")],
                           low_gene_mutation[,c("Hugo_Symbol","freq_low")],by="Hugo_Symbol")
gene_mutation[is.na(gene_mutation)] <- 0

####挑选ATMi诱导高突变gene####
load("3.各组突变相关gene分析.Rdata")
total.gene$log2FC_EtoIvsEto <- log2((total.gene$EtoI+1)/(total.gene$Eto+1))
total.gene$log2FC_EtoAvsEto <- log2((total.gene$EtoA+1)/(total.gene$Eto+1))
EtoI.up.gene <- total.gene[total.gene$log2FC_EtoIvsEto>=log2(4),]$SYMBOL
EtoA.up.gene <- total.gene[total.gene$log2FC_EtoAvsEto>=log2(4),]$SYMBOL

gene_mutation2 <- gene_mutation[gene_mutation$Hugo_Symbol %in% EtoA.up.gene,]
#作箱型图
plot_data <- data.frame(
  Gene = gene_mutation2$Hugo_Symbol,
  mutant_ratio = c(abs(gene_mutation2$freq_high), abs(gene_mutation2$freq_low)),
  Group = rep(c("ATM_High", "ATM_Low"), each = nrow(gene_mutation2))
)

plot_data$Group <- factor(plot_data$Group, levels = c("ATM_High", "ATM_Low")) #设置好箱型图x轴各组的分布顺序
my_comparisons <- list( c("ATM_High", "ATM_Low"))
ggviolin(plot_data,
         x="Group",
         y="mutant_ratio",
         #color="Group",
         xlab= "",
         ylab= "Mutation Rate (%)",
         legend.title="",
         palette = c("#00B9E3","#F8766D"),
         add = "mean_se",#增加散点
         shape = 21,
         fill="Group",
         size=0.75)+
  theme(axis.text.x = element_text(size = 16, angle = 270, vjust = 1, hjust = 0, margin = margin(t = 10)),
        axis.title.y = element_text(size = 16),
        legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
ggsave(filename="./results/ATM低表达肿瘤中ATMi诱导高发突变gene的突变率更高.pdf",
       width = 3,height = 6)

####挑选IRAK1i诱导高突变gene####
gene_mutation2 <- gene_mutation[gene_mutation$Hugo_Symbol %in% EtoI.up.gene &
                                  gene_mutation$Hugo_Symbol!="TTN",]
#作箱型图
plot_data <- data.frame(
  Gene = gene_mutation2$Hugo_Symbol,
  mutant_ratio = c(abs(gene_mutation2$freq_high), abs(gene_mutation2$freq_low)),
  Group = rep(c("ATM_High", "ATM_Low"), each = nrow(gene_mutation2))
)

plot_data$Group <- factor(plot_data$Group, levels = c("ATM_High", "ATM_Low")) #设置好箱型图x轴各组的分布顺序
my_comparisons <- list( c("ATM_High", "ATM_Low"))
ggviolin(plot_data,
         x="Group",
         y="mutant_ratio",
         #color="Group",
         xlab= "",
         ylab= "Mutation Rate (%)",
         legend.title="",
         palette = c("#00B9E3","#F8766D"),
         add = "mean_se",#增加散点
         shape = 21,
         fill="Group",
         size=0.75)+
  theme(axis.text.x = element_text(size = 16, angle = 270, vjust = 1, hjust = 0, margin = margin(t = 10)),
        axis.title.y = element_text(size = 16),
        legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
ggsave(filename="./results/ATM低表达肿瘤中IRAK1i诱导高发突变gene的突变率更高.pdf",
       width = 3,height = 6)
