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
