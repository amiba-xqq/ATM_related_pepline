# WGRS_related_analysis

`0.pepline_WGRS.sh`用于对WGRS进行上游分析，GATK之前的步骤建立构建conda环境，依赖Python3(>=3.8)，需要安装下面的软件：
* fastqc
* multiqc
* trim-galore
* samtools
* bwa

此外还需要使用GATK，从而得到VCF文件；最后使用ANNOVAR对VCF进行注释。
