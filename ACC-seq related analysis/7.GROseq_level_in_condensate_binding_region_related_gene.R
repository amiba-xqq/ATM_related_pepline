rm(list=ls())
library(tidyverse)
library(ggpubr)
####读取GRO-seq的tpm文件,GROseq_U2OS_tpm.Rdata来自于6.GROseq_U2OS_tpm.R####
####Read the TPM file of GRO-seq data (GROseq_U2OS_tpm.Rdata), which comes from `6.GROseq_U2OS_tpm.R`####
load("GROseq_U2OS_tpm.Rdata")

load(file="ACCseq_Ushape_regions_gene.Rdata")
a1 <- tpm_filt[tpm_filt$ENSEMBL %in% DMSO_U_gene,c("U2OS","gene_id")]
a1 <- data.frame(U2OS=a1$U2OS,
                 group="DMSO")
a2 <- tpm_filt[tpm_filt$ENSEMBL %in% pOHT_U_gene,c("U2OS","gene_id")]
a2 <- data.frame(U2OS=a2$U2OS,
                 group="pOHT")
a3 <- tpm_filt[tpm_filt$ENSEMBL %in% pOHTIRAK1i_U_gene,c("U2OS","gene_id")]
a3 <- data.frame(U2OS=a3$U2OS,
                 group="pOHTIRAK1i")

a <- rbind(a1,a2,a3)
a$U2OS <- log2(a$U2OS)

ggboxplot(a,
          x="group",
          y="U2OS",
          fill="group",
          xlab= "",
          ylab= "TPM(log2) of GRO-seq per gene",
          legend.title="",
          palette = c("#00B9E3","#F8766D","orange"))+
  theme(axis.text.x = element_text(size = 14, angle = 270, vjust = 1, hjust = 0, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14),
        legend.position = "none")
ggsave(filename="./results/tpm_boxplot_U2OS_merged_GROseq.pdf",
       width = 6,height = 8)

test <- wilcox.test(a[a$group=="DMSO","U2OS"],
                    a[a$group=="pOHT","U2OS"])
test$p.value

test <- wilcox.test(a[a$group=="pOHTIRAK1i","U2OS"],
                    a[a$group=="pOHT","U2OS"])
test$p.value
