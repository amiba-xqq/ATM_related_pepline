# RIP-seq related analysis

`0.pepline_RIPseq.sh`用于RNA immunoprecipitation and high-throughput sequencing (RIP-seq)的上游分析，建议构建conda环境，依赖Python3(>=3.8)，并需要下载下面的软件：
* fastqc
* multiqc
* fastp
* STAR
* samtools
* bedtools

`Homo_sapiens.GRCh38.108.bed`与`black_bamCoverage.bed`用于`0.pepline_RIPseq.sh`中。`Homo_sapiens.GRCh38.108.bed`记录了hg38基因组每个基因的坐标位置，`black_bamCoverage.bed`记录了hg38基因组中21号染色体上rDNA附近区域。

`1.prepare_rnk_and_gmt_for_GSEAPreranked.R`利用`0.pepline_RIPseq.sh`最后获得的文件`NBS1_RIP_all_gene_merged.txt`,并分析IR+IRAK1i vs IR时同NBS1结合能力下调的RNA，并导出根据IR vs IR+IRAK1i的Fold change排序的rnk文件；根据GTF文件信息，导出各种类型RNA的gmt文件；rnk和gmt文件后续可以使用GSEA的Preranked模式进行分析，从而判断IRAK1主要促进NBS1结合哪些类型的RNA。

`0.pipeline_RIPseq.sh` is used for upstream analysis of RNA immunoprecipitation and high-throughput sequencing (RIP-seq) and it is recommended to build a conda environment, which requires Python3 (version >=3.8) and the following software needs to be downloaded:
* fastqc
* multiqc
* fastp
* STAR
* samtools
* bedtools

`Homo_sapiens.GRCh38.108.bed` and `black_bamCoverage.bed` are used in the script `0.pipeline_RIPseq.sh`. The file `Homo_sapiens.GRCh38.108.bed` records the coordinate positions of each gene in the hg38 genome, while `black_bamCoverage.bed` records the regions near rDNA on chromosome 21 in the hg38 genome.

`1.prepare_rnk_and_gmt_for_GSEAPreranked.R` utilizes the file `NBS1_RIP_all_gene_merged.txt`, which is obtained from the final step of `0.pipeline_RIPseq.sh`. It analyzes RNA with reduced binding ability to NBS1 when comparing IR+IRAK1i versus IR and exports an rnk file sorted by Fold change based on this comparison. Additionally, it generates a gmt file for various types of RNA using information from the GTF file. These rnk and gmt files can subsequently be used in GSEA's Preranked mode to determine which types of RNA are primarily associated with NBS1 that regulated by IRAK1.




