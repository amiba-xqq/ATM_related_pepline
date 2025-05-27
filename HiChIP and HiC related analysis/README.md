# HiChIP and HiC related analysis

`0.pepline_HiChIP.sh`用于对HiChIP和HiC进行上游分析，建议使用conda环境，依赖Python3(>=3.8)，且需要安装下面的软件：
* fastqc
* multiqc
* fastp
* bwa
* samtools
* pairtools

此外还需要安装[HiChIP](https://github.com/dovetail-genomics/HiChiP)。

`hg38.genome`用于`0.pepline_HiChIP.sh`。

`1.FitHiChIP_get_pNBS1_HiChIP_loops.sh`用于分析得到pNBS1 HiChIP loops，并进行可视化，需要提前安装`FitHiChIP`(https://github.com/ay-lab/FitHiChIP) 和`HiCExplorer`(https://github.com/deeptools/HiCExplorer) 。

`configfile_HiChIP_4OHT_DSB`和`configfile_HiChIP_DMSO_DSB`是`1.FitHiChIP_get_pNBS1_HiChIP_loops.sh`中使用`FitHiChIP_HiCPro.sh`取pNBS1 HiChIP loops所需要的设置参数的文档。

`FitHiChIP_DSB_link.table`是`1.FitHiChIP_get_pNBS1_HiChIP_loops.sh`中使用`hicPlotTADs`作图需要的设置参数的文档。

`2.APA_analysis.sh`用于获取DIvA细胞中同一染色体内ATM激活DSB位点与相邻的condensate binding regions/non-DSB pNBS1 peaks配对的loops，以及随机配对loops,最后进行APA分析并作图。需要提前安装`fanc`(https://github.com/vaquerizaslab/fanc) 。


`0.pepline_HiChIP.sh`:is used for upstream analysis of HiChIP and HiC data. It is recommended to use a conda environment with Python3 (>=3.8) installed, along with the following software:
* fastqc
* multiqc
* fastp
* bwa
* samtools
* pairtools

Additionally, you need to install HiChIP from the following GitHub repository: [HiChIP](https://github.com/dovetail-genomics/HiChiP).

`hg38.genome` file is used in `0.pepline_HiChIP.sh`.

`1.FitHiChIP_get_pNBS1_HiChIP_loops.sh`:analyze pNBS1 HiChIP loops and generate visualizations using this script, ensure you have the following tools installed:

- **FitHiChIP:** Available at [FitHiChIP](https://github.com/ay-lab/FitHiChIP).
- **HiCExplorer:** Obtainable from [HiCExplorer](https://github.com/deeptools/HiCExplorer).

`configfile_HiChIP_4OHT_DSB` and `configfile_HiChIP_DMSO_DSB` are essential for `FitHiChIP_HiCPro.sh` within `1.FitHiChIP_get_pNBS1_HiChIP_loops.sh`. These files specify parameters needed to process HiChIP 
data.

For creating plots using `hicPlotTADs`, ensure you have the configuration file `FitHiChIP_DSB_link.table` in your directory. This file is used within the script to set up plotting parameters correctly.

`2.APA_analysis.sh`: identifies loops formed by ATM-activated DSB sites within the same chromosome, pairing them with adjacent condensate binding regions or non-DSB pNBS1 peaks. Additionally, it generates random pairings for comparative analysis.Finally, performs APA analysis to analyze these loop structures. Before running `2.APA_analysis.sh`, install [fanc](https://github.com/vaquerizaslab/fanc).


