# RAP-seq related analysis

`0.pepline_RAPseq.sh`用于对RNA antisense purification coupled sequencing (RAP-seq)进行上游分析，建议构建conda环境，依赖Python3(>=3.8)，并且需要安装下面的软件：
* fastqc
* multiqc
* trim-galore
* STAR
* sambamba
* samtools
* deeptools

`1.Relationship_between_other_ChIPseq_and_RAPseq.R`用于分析各组RAP-seq与其它ChIP-seq在全基因的信号分布特征是否具有相关性。


`0.pipeline_RAPseq.sh` is used for upstream analysis of RNA antisense purification coupled sequencing (RAP-seq). It is recommended to create a conda environment, requires Python3 (≥3.8), and install the software below:
* fastqc
* multiqc
* trim-galore
* STAR
* sambamba
* samtools
* deeptools

`1.Relationship_between_other_ChIPseq_and_RAPseq.R` is used to analyze whether the signal distribution characteristics of other ChIP-seq across all genes are correlated with those of the RAP-seq groups.
