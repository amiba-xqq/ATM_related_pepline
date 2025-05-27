# HiChIP and HiC related analysis

`0.pepline_HiChIP.sh`用于对HiChIP和HiC进行上游分析，建议使用conda环境，依赖Python3(>=3.8)，且需要安装下面的软件：
* fastqc
* multiqc
* fastp
* bwa
* samtools
* pairtools

此外还需要安装HiChiP(https://github.com/dovetail-genomics/HiChiP) 。

`hg38.genome`用于`0.pepline_HiChIP.sh`。

`1.FitHiChIP_get_pNBS1_HiChIP_loops.sh`用于分析得到pNBS1 HiChIP loops，并进行可视化，需要提前安装`FitHiChIP`(https://github.com/ay-lab/FitHiChIP) 和`HiCExplorer`(https://github.com/deeptools/HiCExplorer) 。

`configfile_HiChIP_4OHT_DSB`和`configfile_HiChIP_DMSO_DSB`是`1.FitHiChIP_get_pNBS1_HiChIP_loops.sh`中使用`FitHiChIP_HiCPro.sh`取pNBS1 HiChIP loops所需要的设置参数的文档。

`FitHiChIP_DSB_link.table`是`1.FitHiChIP_get_pNBS1_HiChIP_loops.sh`中使用`hicPlotTADs`作图需要的设置参数的文档。

`2.APA_analysis.sh`用于获取DIvA细胞中同一染色体内ATM激活DSB位点与相邻的condensate binding regions/non-DSB pNBS1 peaks配对的loops，以及随机配对loops,最后进行APA分析并作图。需要提前安装`fanc`(https://github.com/vaquerizaslab/fanc) 。
