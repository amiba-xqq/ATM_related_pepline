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

`1.FitHiChIP_get_pNBS1_HiChIP_loops.sh`用于分析得到pNBS1 HiChIP loops，并进行可视化，需要提前安装FitHiChIP(https://github.com/ay-lab/FitHiChIP) 和HiCExplorer(https://github.com/deeptools/HiCExplorer) 。


