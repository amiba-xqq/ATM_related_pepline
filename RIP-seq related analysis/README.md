# RIP-seq related analysis

`0.pepline_RIPseq.sh`用于RIP-seq的上游分析，建议构建conda环境，依赖Python3(>=3.8)，并需要下载下面的软件：
* fastqc
* multiqc
* fastp
* STAR
* samtools
* bedtools

`Homo_sapiens.GRCh38.108.bed`与`black_bamCoverage.bed`用于`0.pepline_RIPseq.sh`中。`Homo_sapiens.GRCh38.108.bed`记录了hg38基因组每个基因的坐标位置，`black_bamCoverage.bed`记录了hg38基因组中21号染色体上rDNA附近区域。

`1.prepare_rnk_for_GSEAPreranked.R`利用`0.pepline_RIPseq.sh`最后获得的文件`NBS1_RIP_all_gene_merged.txt`,并分析IR vs IR+IRAK1i时同NBS1结合能力下调的RNA，并导出根据IR vs IR+IRAK1i的Fold change排序的rnk文件；根据GTF文件信息，导出各种类型RNA的gmt文件；rnk和gmt文件后续可以使用GSEA的Preranked模式进行分析，从而判断IRAK1主要促进NBS1结合哪些类型的RNA。
