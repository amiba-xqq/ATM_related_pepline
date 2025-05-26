# 分析得到所有组的pNBS1 peaks
# Analysis to obtain pNBS1 peak groups for all
cd /media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA/bam_sort/
macs2 callpeak -c DMSO_input.rmdupsorted.bam \
-t  pNBS1_DMSO_IP.rmdupsorted.bam \
-n pNBS1_DMSO -f BAM -g hs --outdir ../peaks/ -q 0.05 --keep-dup all

cd /media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA/bam_sort/
macs2 callpeak -c 4OHT_input.rmdupsorted.bam  \
-t pNBS1_4OHT_IP.merged.bam \
-n pNBS1_4OHT -f BAM -g hs --outdir ../peaks/ -q 0.05 --keep-dup all

cd /media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA/peaks/
awk -F '\t' '{print $1"\t"$2"\t"$3}' pNBS1_DMSO_summits.bed > tmp1
awk -F '\t' '{print $1"\t"$2"\t"$3}' pNBS1_4OHT_summits.bed > tmp2
cat tmp1 tmp2 > tmp
sort -k1,1V -k2,2n tmp > tmp.bed
bedtools merge -i tmp.bed > total_pNBS1_peaks.bed
rm tmp*

# 分析得到DSB以外的pNBS1 peaks的位点,取Peaks的±1kb范围的regions,并分析各组ChIP-seq在这些region的counts数
# DSB.bed记录了DIvA细胞中90个ATM激活的DSB(Double Strand Break)位点所在位置
# Analyze pNBS1 peaks outside of DSB regions, take ±1kb regions around these peaks, and analyze the ChIP-seq counts within these regions for each group.
# The DSB.bed file records the positions of 90 ATM-activated DSB (Double Strand Break) sites in DIvA cells.
cd /media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA/peaks/
bedtools intersect -v -a total_pNBS1_peaks.bed -b /media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA/pNBS1_analysis/R_analysis/bed/DSB.bed > total_pNBS1_unDSB_peaks.bed
cd /media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA/peaks/
awk '{print $1"\t"($2-1000)"\t"($3+1000)}' total_pNBS1_unDSB_peaks.bed > total_pNBS1_unDSB_peaks_2kb.bed
cd /media/niechen/niechen3/IRAK1_study/gH2AX_ChIPseq_DIvA/bam_sort
bedtools multicov -bams pNBS1_DMSO_IP.rmdupsorted.bam pNBS1_4OHT_IP.merged.bam pNBS1_4OHTIRAK1i_IP.merged.bam \
-bed /media/niechen/niechen3/IRAK1_study/gH2AX_ChIPseq_DIvA/peaks/total_pNBS1_unDSB_peaks_2kb.bed > ../pNBS1_analysis/pNBS1_counts_in_total_pNBS1_unDSB_peaks_2kb.txt

# pNBS1_counts_in_total_pNBS1_unDSB_peaks_2kb.txt文件后续使用R脚本分析得到non-DSB pNBS1 peaks的信息
# 'pNBS1_counts_in_total_pNBS1_unDSB_peaks_2kb.txt' will be analyzed using an R script to obtain information about non-DSB pNBS1 peaks.

