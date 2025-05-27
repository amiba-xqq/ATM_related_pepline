# WGRS_related_analysis

`0.pepline_WGRS.sh`用于对WGRS进行上游分析，GATK之前的步骤建立构建conda环境，依赖Python3(>=3.8)，需要安装下面的软件：
* fastqc
* multiqc
* trim-galore
* samtools
* bwa

此外还需要使用GATK，从而得到VCF文件；最后使用ANNOVAR对VCF进行注释。

`1.vcf_clean.R`利用经过ANNOVAR注释的VCF文件，去除对照组中出现的突变，获得各组特异的突变位点信息。

`2.IRAK1i_or_ATMi_related_mutant.R`分析了IRAK1i和ATMi导致的全局基因突变变化的相关性。

`3.Gene_expression_of_mutants.R`分析了IRAK1i或ATMi处理后突变增加的基因在正常HeLa细胞的表达水平；此外还分析了在Etoposide单处理组，基因转录水平与基因突变水平的相关性。

`4.Other_ChIPseq_signal_in_mutanted_genes.R`分析了其它ChIP-seq在IRAK1i或ATMi处理后突变增加的基因区域的信号分布水平。

`5.Get_TPM_matrix_from_TCGA_BRCA.R`下载TCGA-BRCA的表达矩阵，并转化为TPM。

`6.TCGA_BRCA_related_analysis.R`利用`5.Get_TPM_matrix_from_TCGA_BRCA.R`得到的TPM矩阵，并将BRCA肿瘤分为ATM高表达/低表达肿瘤，最后分析了IRAK1i或ATMi处理后突变增加的基因在这两群肿瘤的突变率。


`0.pepline_WGRS.sh` is used for upstream analysis of WGRS. Before GATK steps, set up a Conda environment. Requires Python 3 (>=3.8). The following software needs to be installed:
* multiqc
* trim-galore
* samtools
* bwa

Additionally, GATK is needed to generate a VCF file; finally, ANNOVAR is used to annotate the VCF.

`1.vcf_clean.R`: Uses the VCF file annotated with ANNOVAR to remove mutations present in the control group, thereby obtaining specific mutation site information for each experimental group.

`2.IRAK1i_or_ATMi_related_mutant.R`: Analyzes the correlation between global gene mutation changes caused by IRAK1 inhibitors(IRAK1i) or ATMi inhibitors(ATMi).

`3.Gene_expression_of_mutants.R`: Examines expression levels of genes with increased mutations after IRAK1i or ATMi treatment in normal HeLa cells. Additionally, it assesses the correlation between gene transcription levels and mutation levels in Etoposide-treated groups.
 
`4.Other_ChIPseq_signal_in_mutated_genes.R`: Investigates the distribution of other ChIP-seq signals across genes with increased mutations post-IRAK1i or ATMi treatment.

`5.Get_TPM_matrix_from_TCGA_BRCA.R`: Downloads and processes gene expression data from TCGA-BRCA to create a TPM (Transcripts Per Million) matrix.

`6.TCGA_BRCA_related_analysis.R`: Uses the TPM matrix generated in `5.Get_TPM_matrix_from_TCGA_BRCA.R`. It classifies BRCA tumors into ATM high-expression and low-expression groups and analyzes mutation rates of genes with increased mutations after IRAK1i or ATMi treatment across these tumor groups.

