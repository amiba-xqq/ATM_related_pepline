# Small RNA-seq related analysis

`0.pepline_smallRNAseq.sh`用于对small RNA-seq进行上游分析，建议构建conda环境，依赖Python3(>=3.8)，并需要安装下面的软件：
* fastqc
* multiqc
* fastp
* bowtie2
* samtools
* bedtools

`1.Counts_analysis_smallRNAseq.R`利用`0.pepline_smallRNAseq.sh`最后得到的`smallRNA_great_clear_all_gene.txt`和`smallRNA_big_clear_all_gene.txt`counts矩阵，分析IR+IRAK1i vs IR+DMSO组差异表达基因，并提取下调的差异表达基因，导出为gmt文件，同`RIP-seq related analysis`部分得到的RIP-seq相关的rnk文件一起，用于后续GSEA分析。


`0.pipeline_smallRNAseq.sh` is used for upstream analysis of small RNA-seq data. It is suggested to build a conda environment, depends on Python3 (>=3.8), and requires the installation of the following software:
* fastqc
* multiqc
* fastp
* bowtie2
* samtools
* bedtools

`1.Counts_analysis_smallRNAseq.R` uses the counts matrices obtained from `0.pipeline_smallRNAseq.sh`, specifically `smallRNA_great_clear_all_gene.txt` and `smallRNA_big_clear_all_gene.txt`, to analyze differentially expressed genes (IR+IRAK1i versus IR+DMSO). It extracts downregulated differentially expressed genes and exports them as a GMT file, which will be used together with the rnk files generated from the `RIP-seq related analysis` part for subsequent GSEA analysis.
