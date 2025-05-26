# 检测pNBS1在AsiSI位点±5kb范围的分布，并使用kmeans算法进行分群
# Detect the distribution of pNBS1 within ±5kb around AsiSI sites and cluster them using the k-means algorithm.
cd /media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA/bw/
computeMatrix reference-point -S  pNBS1_DMSO_IP.rmdupsorted.bw pNBS1_4OHT_IP.merged.bw \
-R /media/niechen/niechen3/IRAK1_study/gH2AX_ChIPseq_DIvA/peaks_AsiSI/AsiSI.total.bed \
-p 20 -a 5000 -b 5000 --referencePoint center --skipZeros --missingDataAsZero --averageTypeBins mean \
-o ../deeptools/all_pNBS1_group_in_AsiSI_5kb.gz \
--outFileSortedRegions ../deeptools/all_pNBS1_group_in_AsiSI_5kb.bed;
plotProfile -m ../deeptools/all_pNBS1_group_in_AsiSI_5kb.gz \
-out ../deeptools/all_pNBS1_group_in_AsiSI_5kb.Profile.pdf \
--plotFileFormat pdf --perGroup --dpi 720  --samplesLabel pNBS1_DMSO pNBS1_DSB   \
--plotHeight 6 --plotWidth 8  --refPointLabel AsiSI --yAxisLabel CPM;
plotHeatmap -m ../deeptools/all_pNBS1_group_in_AsiSI_5kb.gz  \
-out ../deeptools/all_pNBS1_group_in_AsiSI_5kb.Heatmap.pdf \
--plotFileFormat pdf  --dpi 720 --colorMap RdBu_r --samplesLabel pNBS1_DMSO pNBS1_DSB  \
--heatmapHeight 15 --heatmapWidth 4 --refPointLabel AsiSI \
--kmeans 5 \
--outFileSortedRegions /media/niechen/niechen3/IRAK1_study/gH2AX_ChIPseq_DIvA/deeptools/all_pNBS1_group_in_AsiSI_5kb.2.bed

# 取clutster_4做后续分析
# Proceed with subsequent analysis using cluster_4.
cd /media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA/deeptools/
awk '($1 !~ /^#/ ) && (($13 == "cluster_4")) {print}' all_pNBS1_group_in_AsiSI_5kb.2.bed > all_pNBS1_group_in_AsiSI_5kb.cluster_4.bed

# 统计各组ChIP-seq在clutster_4 peaks中心±1kb范围内的counts
# Count the ChIP-seq reads within ±1kb around the center of each peak in cluster_4 for different groups
awk '{print $1"\t"$2"\t"$3}' all_pNBS1_group_in_AsiSI_5kb.cluster_4.bed > /media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA/peaks_AsiSI/tmp
cd /media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA/peaks_AsiSI/
sort -k1,1V -k2,2n tmp > AsiSI_pNBS1_top.bed
awk '{print $1"\t"($2-1000)"\t"($3+1000)}' AsiSI_pNBS1_top.bed > AsiSI_pNBS1_top_2kb.bed
cd /media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA/bam_sort
bedtools multicov -bams pNBS1_DMSO_IP.rmdupsorted.bam pNBS1_4OHT_IP.merged.bam pNBS1_4OHTIRAK1i_IP.merged.bam \
-bed /media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA/peaks_AsiSI/AsiSI_pNBS1_top_2kb.bed > ../pNBS1_analysis/pNBS1_counts_in_DSB_2kb.txt

# pNBS1_counts_in_DSB_2kb.txt文件后续使用R脚本分析得到ATM activated-DSB sites的信息
# Analyze the 'pNBS1_counts_in_DSB_2kb.txt' file using an R script to obtain information about ATM-activated DSB sites.
