rm(list = ls())
library(tidyverse)
library(pheatmap)
library(ggpubr)
library("rtracklayer")
library("SummarizedExperiment")
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggpubr)

####读取HeLa-counts矩阵####
####Read HeLa-counts matirx#####
counts_HeLa <- read.table("counts_HeLa.txt",header = T)
counts_HeLa <- counts_HeLa[,c(1,6:9)]
colnames(counts_HeLa)[3:5] <- c("HeLa1","HeLa2","HeLa3")
counts_HeLa$HeLa <- rowMeans(counts_HeLa[3:5])
counts_all <- counts_HeLa[,c(1,2,6)]
geneid_efflen <- subset(counts_all,select = c("Geneid","Length"))
colnames(geneid_efflen) <- c("geneid","efflen")  
counts <- data.frame(HeLa=counts_all$HeLa)
rownames(counts) <- counts_all[,1]
efflen <- geneid_efflen[match(rownames(counts),
                              geneid_efflen$geneid),
                        "efflen"]

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
counts_filt <- data.frame(HeLa=counts[keep_feature, ]) 
rownames(counts_filt) <- rownames(counts)[keep_feature]
tpm_filt <- data.frame(HeLa=tpm[keep_feature, ])
rownames(tpm_filt) <- rownames(tpm)[keep_feature]
tpm_filt$gene_id <- rownames(tpm_filt)
save(tpm_filt,file = "HeLa_tpm_filt.Rdata")

####计算各组gene mutant的Fold change####
####Calculate the Fold Change of Gene Mutants Across Groups####
rm(list = ls())
load("2.mutant_gene.Rdata")
load("HeLa_tpm_filt.Rdata")
total.gene$log2FC_EtoIvsEto <- log2((total.gene$EtoI+1)/(total.gene$Eto+1))
total.gene$log2FC_EtoAvsEto <- log2((total.gene$EtoA+1)/(total.gene$Eto+1))
df <- bitr(total.gene$SYMBOL, 
           fromType = "SYMBOL",
           toType =  "ENSEMBL",
           OrgDb = org.Hs.eg.db) 
total.gene <- inner_join(total.gene,df,by="SYMBOL")

####分析IRAK1i/ATMi处理后突变增加的gene的表达水平####
####Analyze the expression levels of genes where mutations have increased following treatment with IRAK1i/ATMi####
EtoI.up.gene <- total.gene[total.gene$log2FC_EtoIvsEto>=log2(4),]$ENSEMBL
EtoA.up.gene <- total.gene[total.gene$log2FC_EtoAvsEto>=log2(4),]$ENSEMBL
Eto.top.gene <- total.gene[total.gene$Eto>quantile(total.gene$Eto,0.9),]$ENSEMBL
EtoI.top.gene <- total.gene[total.gene$EtoI>quantile(total.gene$EtoI,0.9),]$ENSEMBL
EtoA.top.gene <- total.gene[total.gene$EtoA>quantile(total.gene$EtoA,0.9),]$ENSEMBL
others.gene <- total.gene[total.gene$log2FC_EtoIvsEto<0 & 
                            total.gene$log2FC_EtoAvsEto<0,]$ENSEMBL

a1 <- tpm_filt[rownames(tpm_filt) %in% EtoI.up.gene,c("HeLa","gene_id")]
a1 <- data.frame(HeLa=a1$HeLa,
                 group="Eto+IRAK1i")
a2 <- tpm_filt[rownames(tpm_filt) %in% EtoA.up.gene,c("HeLa","gene_id")]
a2 <- data.frame(HeLa=a2$HeLa,
                 group="Eto+ATMi")
a3 <- tpm_filt[rownames(tpm_filt) %in% Eto.top.gene,c("HeLa","gene_id")]
a3 <- data.frame(HeLa=a3$HeLa,
                 group="Eto.top10%")
a4 <- tpm_filt[rownames(tpm_filt) %in% EtoI.top.gene,c("HeLa","gene_id")]
a4 <- data.frame(HeLa=a4$HeLa,
                 group="Eto+IRAK1i.top10%")
a5 <- tpm_filt[rownames(tpm_filt) %in% EtoA.top.gene,c("HeLa","gene_id")]
a5 <- data.frame(HeLa=a5$HeLa,
                 group="Eto+ATMi.top10%")
a6 <- tpm_filt[rownames(tpm_filt) %in% total.gene$ENSEMBL,c("HeLa","gene_id")]
a6 <- data.frame(HeLa=a6$HeLa,
                 group="total")
a7 <- tpm_filt[rownames(tpm_filt) %in% others.gene,c("HeLa","gene_id")]
a7 <- data.frame(HeLa=a6$HeLa,
                 group="others")

a <- rbind(a1,a2,a7)
a$HeLa <- log2(a$HeLa)

my_comparisons <- list( c("others", "Eto+IRAK1i"),
                        c("others","Eto+ATMi"))
a$group <- factor(a$group ,level=c("Eto+IRAK1i","Eto+ATMi","others"))
ggboxplot(a,
          x="group",
          y="HeLa",
          fill="group",
          xlab= "",
          ylab= "TPM(log2) in HeLa cell",
          legend.title="",
          palette = c("#00B9E3","#F8766D","orange","lightblue","lightgreen","grey"))+
  theme(axis.text.x = element_text(size = 14, angle = 270, vjust = 1, hjust = 0, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14),
        legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
ggsave("./results/Genes where mutations have increased due to either IRAK1i or ATMi tend to be those with higher transcriptional activity.pdf",width=3,height=6)

mean(a[a$group=="Eto+IRAK1i",]$HeLa)
mean(a[a$group=="Eto+ATMi",]$HeLa)
mean(a[a$group=="others",]$HeLa)

a <- rbind(a3,a4,a5,a6)
a$HeLa <- log2(a$HeLa)

my_comparisons <- list( c("total","Eto.top10%"),
                        c("total","Eto+IRAK1i.top10%"),
                        c("total","Eto+ATMi.top10%"))
a$group <- factor(a$group ,level=c("Eto.top10%","Eto+IRAK1i.top10%",
                                   "Eto+ATMi.top10%","total"))
ggboxplot(a,
          x="group",
          y="HeLa",
          fill="group",
          xlab= "",
          ylab= "TPM(log2) in HeLa cell",
          legend.title="",
          palette = c("#00B9E3","#F8766D","orange","grey"))+
  theme(axis.text.x = element_text(size = 14, angle = 270, vjust = 1, hjust = 0, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14),
        legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
ggsave("./results/High-mutation genes induced by Etoposide tend to be genes with low transcriptional activity.pdf",width=3,height=6)

####将Etoposide组按照突变总数分类，再分析各组表达水平####
####Classify the Etoposide group based on total mutation counts and then analyze the expression levels of each category####
Eto <- total.gene[,c("Eto","ENSEMBL")]
colnames(Eto)[2] <- "gene_id"
Eto <- inner_join(Eto,tpm_filt,by="gene_id")
Eto <- Eto[order(Eto$Eto,decreasing = T),]
k <- round(nrow(Eto)/5)
Eto$group <- c(rep("Q1",k),rep("Q2",k),rep("Q3",k),rep("Q4",k),rep("Q5",(nrow(Eto)-4*k)))
Eto$HeLa <- log2(Eto$HeLa)

my_comparisons <- list( c("Q1","Q2"),
                        c("Q1","Q3"),
                        c("Q1","Q4"),
                        c("Q1","Q5"))
ggboxplot(Eto,
          x="group",
          y="HeLa",
          fill="group",
          xlab= "",
          ylab= "TPM(log2) in HeLa cell",
          legend.title="",
          palette = c("#0072B2", "#009E79", "#D55E00", "#CC79A7", "#FA8072"))+
  theme(axis.text.x = element_text(size = 14, angle = 270, vjust = 1, hjust = 0, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14),
        legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
ggsave("./results/High-mutation genes induced by Etoposide tend to be genes with low transcriptional activity_2.pdf",width=3,height=6)

####先挑选出转录活跃基因，然后分析各组中这些转录活跃基因的突变个数####
####First select the highly transcribed genes, then analyze the number of mutations in these genes across groups####
rm(list = ls())
load("2.mutant_gene.Rdata")
load("HeLa_tpm_filt.Rdata")
colnames(tpm_filt)[2] <- "ENSEMBL"
total.gene$log2FC_EtoIvsEto <- log2((total.gene$EtoI+1)/(total.gene$Eto+1))
total.gene$log2FC_EtoAvsEto <- log2((total.gene$EtoA+1)/(total.gene$Eto+1))
df <- bitr(total.gene$SYMBOL, 
           fromType = "SYMBOL",
           toType =  "ENSEMBL",
           OrgDb = org.Hs.eg.db) 
total.gene <- inner_join(total.gene,df,by="SYMBOL")
total.gene <- inner_join(total.gene,tpm_filt,by="ENSEMBL")
total.gene <- total.gene[order(total.gene$HeLa,decreasing = T),]
k <- round(nrow(total.gene)/5)
total.gene$group <- c(rep("Q1",k),rep("Q2",k),rep("Q3",k),
                      rep("Q4",k),rep("Q5",nrow(total.gene)-4*k))

total.gene$Eto <- log2(total.gene$Eto+1)
my_comparisons <- list( c("Q1","Q2"),
                        c("Q1","Q3"),
                        c("Q1","Q4"),
                        c("Q1","Q5"))
ggboxplot(total.gene,
          x="group",
          y="Eto",
          fill="group",
          xlab= "",
          ylab= "Number of mutant per gene (log2)",
          legend.title="",
          palette = c("#0072B2", "#009E79", "#D55E00", "#CC79A7", "#FA8072"))+
  theme(axis.text.x = element_text(size = 14, angle = 270, vjust = 1, hjust = 0, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14),
        legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
ggsave("./results/highly transcribed genes are less likely to undergo mutations.pdf",width=3,height=6)


