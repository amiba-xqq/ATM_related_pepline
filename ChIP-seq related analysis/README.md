# ChIP-seq related analysis

0.pepline_ChIPseq.sh用于对ChIP-seq的raw data进行上游分析。建议构建conda环境，依赖Python3(>=3.8)，并已安装下面的软件：fastqc;multiqc;trim-galore;hisat2;sambamba;samtools;deeptools。

1.1.ATM_activated_DSB_DIvA_cells.sh与1.2.ATM_activated_DSB_DIvA_cells.R用于获取DIvA细胞中ATM激活的DSB位点。AsiSI.total.bed是1.1.ATM_activated_DSB_DIvA_cells.sh使用的文件，标注了hg38基因组所有的AsiSI酶切位点。

1.3.DSB_annotation.R用于对ATM激活的DSB位点进行注释，并分析其中IRAK1依赖与不依赖的DSB位点。

1.4.Expression_of_DSB_related_genes.R用于分析IRAK1依赖与不依赖的DSB位点相关基因在正常U2OS细胞中的表达水平。

1.5.Other_ChIPseq_signal_in_DSB_related_genes.R用于分析IRAK1依赖与不依赖的DSB位点相关基因附近组蛋白修饰相关ChIPseq信号的相对水平。TSS_ALL_seq_counts.txt是该脚本使用的文件，是各组ChIPseq在所有基因的TSS±5kb范围内的counts矩阵。

2.1.nonDSB_pNBS1_peaks.sh与2.2.nonDSB_pNBS1_peaks.R用于分析non-DSB pNBS1 peak位置。该shell脚本需要使用MACS2，建议构建conda环境，依赖Python2(>=2.7)。

######
0.ChIPseq_pipeline.sh is used for upstream analysis of raw data in ChIP-seq experiments. It is recommended to build a conda environment with Python3 (version >=3.8) and install 
the following software: fastqc, multiqc, trim-galore, hisat2, sambamba, samtools, deeptools.

1.1.ATM_activated_DSB_DIvA_cells.sh and 1.2.ATM_activated_DSB_DIvA_cells.R are used to obtain ATM-activated DSB sites in DIvA cells. AsiSI.total.bed is the file used by 1.1.ATM_activated_DSB_DIvA_cells.sh, marking all AsiSI enzyme cutting sites in the hg38 genome.



