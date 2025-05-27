rm(list = ls())
library(tidyverse)
library(maftools)
library(TCGAbiolinks)
library(TCGAmutations)
library(SummarizedExperiment)
library("IOBR")

####下载TCGA数据库的RNA-seq表达矩阵####
####Download the RNA-seq expression matrix from the TCGA database####
query_rnaseq <- GDCquery(project = "TCGA-BRCA", 
                         experimental.strategy = "RNA-Seq",
                         data.category = "Transcriptome Profiling", 
                         data.type = "Gene Expression Quantification",
                         workflow.type = "STAR - Counts")

GDCdownload(query_rnaseq)
Rnaseq <- GDCprepare(query_rnaseq)
# 得到样本的基因counts矩阵，每行为一个基因，每列为一个样本
# Obtain the gene count matrix for the samples, where each row represents a gene and each column represents a sample.
count_matrix <- assay(Rnaseq)
# selection of tumor samples "TP"
#TP: PRIMARY SOLID TUMOR
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(count_matrix),typesample = c("TP"))
total_sample_ID <- sapply(strsplit(samplesTP,'-'),function(x) paste0(x[1:3],collapse="-"))
colnames(count_matrix)[1:3]

count_matrix <- count_matrix[,samplesTP]
count_matrix[1:5,1:5]
rownames(count_matrix) <- substr(rownames(count_matrix), 1, 15)
nrow(count_matrix)

# 使用IOBR包的count2tpm函数将counts矩阵转化成TPM
# Use the count2tpm function from the IOBR package to convert counts into TPM.
tpm <- count2tpm(countMat = count_matrix, source = "local", idType = "ensembl") 
colSums(tpm)
colnames(tpm) <- sapply(strsplit(colnames(tpm),'-'),function(x) paste0(x[1:3],collapse="-"))
save(tpm,file="TCGA_BRCA_TPM.Rdata")
