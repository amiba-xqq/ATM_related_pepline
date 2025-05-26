# ACC-seq related analysis

`0.pepline_ACCseq.sh`用于对ACC-seq原始数据进行上游分析。建议构建conda环境，依赖Python3(>=3.8)，并已安装下面的软件：
* fastqc
* multiqc
* fastp
* STAR
* sambamba
* samtools
* deeptools

`1.condensate_binding_regions_analysis.sh`与 `2.condensate_binding_regions_analysis.R`用于得到ACC-seq各组的condendate binding regions。该shell脚本需要使用MACS2，建议构建conda环境，依赖Python2(>=2.7)。
`hg38.chrom.sizes`用于`1.condensate_binding_regions_analysis.sh`脚本。

`3.condensate_binding_regions_annotation.R`用于对condendate binding regions进行注释，分析在基因组分布的特征，以及进行相关基因的通路富集分析。

`4.Gene_expression_of_condensate_binding_regions.R`用于分析各组condendate binding regions相关基因在正常U2OS细胞中的表达水平。

`5.pepline_GROseq_U2OS.sh`用于对正常U2OS细胞的GRO-seq原始数据进行上游分析。建议构建conda环境，依赖Python3(>=3.8)，并已安装下面的软件：
* sra-tools
* rename
* fastqc
* multiqc
* bowtie2
* samtools
* deeptools
* featureCounts

`6.GROseq_U2OS_tpm.R`用于正常U2OS细胞的GRO-seq原始数据的下游分析，从而得到TPM表达矩阵，用于后续`7.GROseq_level_in_condensate_binding_region_related_gene.R`的分析。

`7.GROseq_level_in_condensate_binding_region_related_gene.R`用于分析各组condendate binding regions相关基因在GRO-seq中的TPM值。

`8.distance_between_condensate_binding_regions_and_other_genomic_sites.R`用于分析condensate_binding_region与`ChIP-seq related analysis`与最临近的ATM激活DSB位点或non-DSB pNBS1 peaks之间的基因组距离。


`0.pepline_ACCseq.sh` is used for upstream analysis of ACC-seq raw data. It is recommended to set up a conda environment, which requires Python3 (version 3.8 or higher), and the following software should be installed:
* fastqc
* multiqc
* fastp
* STAR
* sambamba
* samtools
* deeptools

`1.condensate_binding_regions_analysis.sh` and `2.condensate_binding_regions_analysis.R` are used to obtain condensate binding regions for each group of ACC-seq data. This shell script requires MACS2, so it's recommended to set up a conda environment depending on Python2 (version >=2.7). The `hg38.chrom.sizes` file is needed by the `1.condensate_binding_regions_analysis.sh` script.

`3.condensate_binding_regions_annotation.R` annotates condensate binding regions, analyzes their distribution across the genome, and performs pathway enrichment analysis for genes associated with these regions.

`4.Gene_expression_of_condensate_binding_regions.R` analyzes the expression levels of genes related to condensate binding regions in normal U2OS cells.

`5.pepline_GROseq_U2OS.sh` is designed for upstream analysis of raw GRO-seq data from normal U2OS cells. It is recommended to create a conda environment that depends on Python3 (version 3.8 or higher) and install the following software tools:
* sra-tools
* rename
* fastqc
* multiqc
* bowtie2
* samtools
* deeptools
* featureCounts

`6.GROseq_U2OS_tpm.R` performs downstream analysis of GRO-seq raw data to generate TPM (Transcripts Per Million) expression matrices, which are then used for subsequent analyses by the script `7.GROseq_level_in_condensate_binding_region_related_gene.R`.

`7.GROseq_level_in_condensate_binding_region_related_gene.R` evaluates the TPM values of genes associated with condensate binding regions within the GRO-seq data.

`8.distance_between_condensate_binding_regions_and_other_genomic_sites.R` calculates the genomic distances between condensate binding regions and other significant genomic sites, including those related to ChIP-seq analysis and the nearest ATM-activated double-strand break (DSB) sites or non-DSB pNBS1 peaks.

