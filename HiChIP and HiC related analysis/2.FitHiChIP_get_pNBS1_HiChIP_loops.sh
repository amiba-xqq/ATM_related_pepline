#!/bin/bash
dir=/mnt/NC/HiChIP
cd ${dir}/FitHiChIP
mkdir DMSOHiChIP
# 需要提前准备好configfile_HiChIP文档
# The configfile_HiChIP document needs to be prepared in advance.
~/FitHiChIP/FitHiChIP_HiCPro.sh -C configfile_HiChIP_4OHT_DSB
~/FitHiChIP/FitHiChIP_HiCPro.sh -C configfile_HiChIP_DMSO_DSB

# 使用HicPlotTADs作图，展示各组中同DSB相互作用的loops
# Use HicPlotTADs to create plots and display the loops interacting with the DSB sites for each group.
cd /mnt/NC/HiChIP/FitHiChIP/
hicPlotTADs --tracks FitHiChIP_DSB_link4.table --width 20 --height 10 --region chr2:53,369,673-57,369,673 -o ./figures/FitHiChIP_4OHTvsDSB_links.pdf
