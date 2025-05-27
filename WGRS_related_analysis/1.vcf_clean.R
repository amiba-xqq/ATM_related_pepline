rm(list = ls())
library(tidyverse)
####读取注释后的vcf文件####
####Read the annotated VCF file####
ctrl <- read.csv("./vcf_annotation/1WT.hg38_multianno.csv")
ctrl$type <- paste0(ctrl$Chr,"_",ctrl$Start,"_",ctrl$End,"_",ctrl$Ref,"_",ctrl$Alt)
#2Eto
Eto <- read.csv("./vcf_annotation/1Eto.hg38_multianno.csv")
Eto$type <- paste0(Eto$Chr,"_",Eto$Start,"_",Eto$End,"_",Eto$Ref,"_",Eto$Alt)
Eto_diff <- Eto[(!(Eto$type %in% ctrl$type))&(Eto$Chr != "chrY")&(Eto$Chr != "chrY_KI270740v1_random"),]
Eto_diff <- Eto_diff[,-61]
write.csv(Eto_diff,"./vcf_clean/1Eto_vcf_diff.csv")
rm(Eto,Eto_diff)

#3EtoA
EtoA <- read.csv("./vcf_annotation/2EtoA.hg38_multianno.csv")
EtoA$type <- paste0(EtoA$Chr,"_",EtoA$Start,"_",EtoA$End,"_",EtoA$Ref,"_",EtoA$Alt)
EtoA_diff <- EtoA[(!(EtoA$type %in% ctrl$type))&(EtoA$Chr != "chrY")&(EtoA$Chr != "chrY_KI270740v1_random"),]
EtoA_diff <- EtoA_diff[,-61]
write.csv(EtoA_diff,"./vcf_clean/2EtoA_vcf_diff.csv")
rm(EtoA,EtoA_diff)

#4CPTI
EtoI <- read.csv("./vcf_annotation/3EtoI.hg38_multianno.csv")
EtoI$type <- paste0(EtoI$Chr,"_",EtoI$Start,"_",EtoI$End,"_",EtoI$Ref,"_",EtoI$Alt)
EtoI_diff <- EtoI[(!(EtoI$type %in% ctrl$type))&(EtoI$Chr != "chrY")&(EtoI$Chr != "chrY_KI270740v1_random"),]
EtoI_diff <- EtoI_diff[,-61]
write.csv(EtoI_diff,"./vcf_clean/3EtoI_vcf_diff.csv")
rm(EtoI,EtoI_diff)

# 以ctrl组作为对照；对于实验组，去除对照组中出现的突变，最后分析得到实验组特异的突变位点，保存在vcf_clean文件夹中
# Use the ctrl group as the control; for the experimental group, remove mutations that appear in the control group, and save the resulting experiment-specific mutation sites in the vcf_clean folder.

