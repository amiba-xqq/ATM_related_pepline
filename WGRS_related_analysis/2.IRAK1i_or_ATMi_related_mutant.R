rm(list = ls())
library(tidyverse)
library(FactoMineR)
library(factoextra)  
####读取vcf文件####
####Read VCF file####
Eto <- read.csv("./vcf_clean/1Eto_vcf_diff.csv")
EtoI <- read.csv("./vcf_clean/3EtoI_vcf_diff.csv")
EtoA <- read.csv("./vcf_clean/2EtoA_vcf_diff.csv")

####分析各个gene的突变总数，只保留genic的突变####
####Analyze the total number of mutations for each gene and retain only genic mutations####
Eto.gene <- Eto[Eto$Func.refGene!="intergenic",]
Eto.gene <- Eto.gene %>%
  mutate(Gene.refGene = str_replace_all(str_extract(Gene.refGene, "^(.*?)(?:;|$)"), ";", "")) 
Eto.gene <- as.data.frame(table(Eto.gene$Gene.refGene))
colnames(Eto.gene) <- c("SYMBOL","Eto")

EtoA.gene <- EtoA[EtoA$Func.refGene!="intergenic",]
EtoA.gene <- EtoA.gene %>%
  mutate(Gene.refGene = str_replace_all(str_extract(Gene.refGene, "^(.*?)(?:;|$)"), ";", "")) 
EtoA.gene <- as.data.frame(table(EtoA.gene$Gene.refGene))
colnames(EtoA.gene) <- c("SYMBOL","EtoA")

EtoI.gene <- EtoI[EtoI$Func.refGene!="intergenic",]
EtoI.gene <- EtoI.gene %>%
  mutate(Gene.refGene = str_replace_all(str_extract(Gene.refGene, "^(.*?)(?:;|$)"), ";", "")) 
EtoI.gene <- as.data.frame(table(EtoI.gene$Gene.refGene))
colnames(EtoI.gene) <- c("SYMBOL","EtoI")

####合并所有突变####
####Merge all mutants####
total.gene <- full_join(Eto.gene,EtoI.gene,by="SYMBOL")
total.gene <- full_join(total.gene,EtoA.gene,by="SYMBOL")
total.gene[is.na(total.gene)] <- 0

total.gene$log2FC_EtoIvsEto <- log2((total.gene$EtoI+1)/(total.gene$Eto+1))
total.gene$log2FC_EtoAvsEto <- log2((total.gene$EtoA+1)/(total.gene$Eto+1))
table(total.gene$log2FC_EtoIvsEto>=log2(2))
table(total.gene$log2FC_EtoIvsEto<= -log2(2))
table(total.gene$log2FC_EtoAvsEto>=log2(2))
table(total.gene$log2FC_EtoAvsEto<= -log2(2))

####分析IRAK1i与ATMi导致的基因突变的变化是否相关####
####Analyze whether changes in gene mutations caused by IRAK1i and ATMi are correlated####
total.gene.filt <- total.gene[(total.gene$Eto+total.gene$EtoI+total.gene$EtoA)>=5,]
test <- cor.test(total.gene.filt$log2FC_EtoIvsEto,
                 total.gene.filt$log2FC_EtoAvsEto,method="pearson",use="complete.obs")
p_value <- test$p.value
correlation_coef <- test$estimate 

ggplot(total.gene.filt,aes(x=log2FC_EtoIvsEto,y=log2FC_EtoAvsEto))+
  geom_point(color = "#00B9E3",alpha=0.4,size=2)+
  geom_smooth(method = "lm", se = T, color = "orange") + 
  theme_classic(base_size = 16)+
  annotate("text", x = min(total.gene.filt$log2FC_EtoIvsEto) + 0.5, y = max(total.gene.filt$log2FC_EtoAvsEto) + 0.5,
           label = paste("Pearson's r =", round(correlation_coef, 2), "\n", 
                         "p =", signif(p_value, digits = 3)), 
           hjust = 0, vjust = 1, size = 6, color = "black")+
  xlab("Mutant number change (Eto+IRAK1i vs Eto)")+
  ylab("Mutant number change (Eto+ATMi vs Eto)")
ggsave('./results/IRAK1i和ATMi带来的突变变化是类似的.pdf',width = 6,height =6)































