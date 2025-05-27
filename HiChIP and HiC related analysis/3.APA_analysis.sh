# 产生同一染色体内部DSB与相邻condensate binding regions/non-DSB pNBS1 peaks配对的loops
# Generate loops between intra-chromosomal DSBs and adjacent condensate binding regions or non-DSB pNBS1 peaks.
cd /mnt/NC/HiChIP/FitHiChIP/bed/
bedtools window -a DSBs.bed -b pOHT_Ushape.bed -l 2000000 -r 2000000 -sm > DSB_pOHT_Ushape_pairs.bedpe
bedtools window -a DSB.bed -b pNBS1_noDSB_peaks.bed -l 2000000 -r 2000000 -sm > DSB_pNBS1_noDSB_peaks.bedpe

# 产生同一染色体上随机配对的loops，作为对照
# Generate randomly intra-chromosomal paired loops as a control.
awk '{ print $1"\t"($2-2000000)}' hg38.genome > tmp0
bedtools random -l 1 -n 799 -g tmp0 > tmp
awk '$2 > 10000 {print $1"\t"$2"\t"$3"\t"int(rand() * (2000000 - 10000) ) + 10000 }' tmp > tmp2
awk '{print $1"\t"$2"\t"$3"\t"$1"\t"($2+$4)"\t"($3+$4)}' tmp2 > tmp3
sort -k1,1V -k2,2n -k3,3n tmp3 > random_pairs.bedpe
rm tmp*

# 将hic转化为cool
# Convert hic to cool.
cd /mnt/NC/HiChIP/hic/
hicConvertFormat --matrices 4OHTHiC.hic --inputFormat hic \
--outputFormat cool -o 4OHTHiC.cool --resolutions 50000

# APA分析作图
# Create plots for APA analysis.
cd /mnt/NC/HiChIP/APA/
fanc aggregate /mnt/NC/HiChIP/hic/4OHTHiC_50000.cool \
/mnt/NC/HiChIP/FitHiChIP/bed/DSB_pOHT_Ushape_pairs.bedpe \
4OHTHiC_in_DSB_pOHT_Ushape_loops.agg \
-p 4OHTHiC_in_DSB_pOHT_Ushape_loops.pdf  \
-m 4OHTHiC_in_DSB_pOHT_Ushape_loops.table \
--colormap coolwarm -e \
--vmax 1.6 --vmin 0.5

cd /mnt/NC/HiChIP/APA/
fanc aggregate /mnt/NC/HiChIP/hic/4OHTHiC_50000.cool \
/mnt/NC/HiChIP/FitHiChIP/bed/DSB_pNBS1_noDSB_peaks.bedpe \
4OHTHiC_in_DSB_pNBS1_noDSB_peaks_loops.agg \
-p 4OHTHiC_in_DSB_pNBS1_noDSB_peaks_loops.pdf  \
-m 4OHTHiC_in_DSB_pNBS1_noDSB_peaks_loops.table \
--colormap coolwarm -e  \
--vmax 1.6 --vmin 0.5  

cd /mnt/NC/HiChIP/APA/
fanc aggregate /mnt/NC/HiChIP/hic/4OHTHiC_50000.cool \
/mnt/NC/HiChIP/FitHiChIP/bed/random_pairs.bedpe \
4OHTHiC_in_random_pairs.agg \
-p 4OHTHiC_in_random_pairs.pdf  \
-m 4OHTHiC_in_random_pairs.table \
--colormap coolwarm -e  \
--vmax 1.6 --vmin 0.5
