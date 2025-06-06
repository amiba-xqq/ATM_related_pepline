# ChIP-seq related analysis

`0.pepline_ChIPseq.sh`用于对ChIP-seq的raw data进行上游分析。建议构建conda环境，依赖Python3(>=3.8)，并已安装下面的软件：
* fastqc
* multiqc
* trim-galore
* hisat2
* sambamba
* samtools
* deeptools

`1.1.ATM_activated_DSB_DIvA_cells.sh`与`1.2.ATM_activated_DSB_DIvA_cells.R`用于获取DIvA细胞中ATM激活的DSB位点。`AsiSI.total.bed`是`1.1.ATM_activated_DSB_DIvA_cells.sh`使用的文件，标注了hg38基因组所有的AsiSI酶切位点。

`1.3.DSB_annotation.R`用于对ATM激活的DSB位点进行注释，并分析其中IRAK1依赖与不依赖的DSB位点。

`1.4.Expression_of_DSB_related_genes.R`用于分析IRAK1依赖与不依赖的DSB位点相关基因在正常的U2OS细胞中的表达水平。

`1.5.Other_ChIPseq_signal_in_DSB_related_genes.R`用于分析IRAK1依赖与不依赖的DSB位点相关基因附近组蛋白修饰相关ChIPseq信号的相对水平。`TSS_ALL_seq_counts.txt`是该脚本使用的文件，是各组ChIPseq在所有基因的TSS±5kb范围内的counts矩阵。

`2.1.nonDSB_pNBS1_peaks.sh`与`2.2.nonDSB_pNBS1_peaks.R`用于分析non-DSB pNBS1 peak位置。该shell脚本需要使用MACS2，建议构建conda环境，依赖Python2(>=2.7)。


`0.ChIPseq_pipeline.sh` is used for upstream analysis of raw data in ChIP-seq experiments. It is recommended to build a conda environment with Python3 (version >=3.8) and install 
the following software:
* fastqc
* multiqc
* trim-galore
* hisat2
* sambamba
* samtools
* deeptools

`1.1.ATM_activated_DSB_DIvA_cells.sh` and `1.2.ATM_activated_DSB_DIvA_cells.R` are used to obtain ATM-activated DSB sites in DIvA cells. `AsiSI.total.bed` is the file used by `1.1.ATM_activated_DSB_DIvA_cells.sh`, marking all AsiSI enzyme cutting sites in the hg38 genome.

`1.3.DSB_annotation.R` is used to annotate the ATM-activated DSB sites and analyze the DSB sites that are IRAK1-dependent or IRAK1-independent.

`1.4.Expression_of_DSB_related_genes.R` is used to analyze the expression levels of genes related to IRAK1-dependent and IRAK1-independent DSB sites in normal U2OS cells.

`1.5.Other_ChIPseq_signal_in_DSB_related_genes.R` is used to analyze the relative levels of histone modification-related ChIPseq signals near genes associated with IRAK1-dependent and IRAK1-independent DSB sites. `TSS_ALL_seq_counts.txt` is the file used by this script, containing a matrix of counts for all groups of ChIPseq data within ±5kb of the transcription start site (TSS) of all genes.

`2.1.nonDSB_pNBS1_peaks.sh` and `2.2.nonDSB_pNBS1_peaks.R` are used to analyze pNBS1 peak positions at non-DSB regions. This shell script requires MACS2, so it is recommended to build a conda environment with Python 2 (version ≥2.7).
