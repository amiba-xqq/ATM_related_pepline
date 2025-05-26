# 使用MACS2取各组的peaks
# Use MACS2 to extract peaks for each group
cd /media/niechen/niechen3/IRAK1_study/ACCseq/bam_sort_merged/
file_bamsort3=*bam
ls $file_bamsort3 | while read id
do
macs2 callpeak -t $id --nomodel --nolambda -n $(basename -s .bam $id) -f BAM -g hs --outdir ../peaks/ -q 0.05 --keep-dup all
done

# 合并DMSO组的narrowpeaks
# Merge narrowpeaks of DMSO group
cd /media/niechen/niechen3/IRAK1_study/ACCseq/peaks/
awk -F '\t' '{print $1"\t"$2"\t"$3}' DMSO_native_peaks.narrowPeak > tmp1
awk -F '\t' '{print $1"\t"$2"\t"$3}' DMSO_fix_peaks.narrowPeak > tmp2
awk -F '\t' '{print $1"\t"$2"\t"$3}' DMSO_fixHex_peaks.narrowPeak > tmp3
cat tmp1 tmp2 tmp3 > tmp
sort -k1,1V -k2,2n tmp > tmp.bed
bedtools merge -i tmp.bed > total_DMSO_ACCseq_peaks.bed
rm tmp*

# 合并pOHT组的narrowpeaks
# Merge narrowpeaks of pOHT group
cd /media/niechen/niechen3/IRAK1_study/ACCseq/peaks/
awk -F '\t' '{print $1"\t"$2"\t"$3}' pOHT_native_peaks.narrowPeak > tmp1
awk -F '\t' '{print $1"\t"$2"\t"$3}' pOHT_fix_peaks.narrowPeak > tmp2
awk -F '\t' '{print $1"\t"$2"\t"$3}' pOHT_fixHex_peaks.narrowPeak > tmp3
cat tmp1 tmp2 tmp3 > tmp
sort -k1,1V -k2,2n tmp > tmp.bed
bedtools merge -i tmp.bed > total_pOHT_ACCseq_peaks.bed
rm tmp*

# 合并pOHTIRAK1i组的narrowpeaks
# Merge narrowpeaks of pOHTIRAK1i group
cd /media/niechen/niechen3/IRAK1_study/ACCseq/peaks/
awk -F '\t' '{print $1"\t"$2"\t"$3}' pOHTIRAK1i_native_peaks.narrowPeak > tmp1
awk -F '\t' '{print $1"\t"$2"\t"$3}' pOHTIRAK1i_fix_peaks.narrowPeak > tmp2
awk -F '\t' '{print $1"\t"$2"\t"$3}' pOHTIRAK1i_fixHex_peaks.narrowPeak > tmp3
cat tmp1 tmp2 tmp3 > tmp
sort -k1,1V -k2,2n tmp > tmp.bed
bedtools merge -i tmp.bed > total_pOHTIRAK1i_ACCseq_peaks.bed
rm tmp*

# 统计各组的native/fix/fixHex的ACC-seq，在narrowpeaks ±2.5kb的counts
# Perform statistical analysis on the native/fix/fixHex ACC-seq counts within 2.5kb of narrowpeaks for each group
cd /media/niechen/niechen3/IRAK1_study/ACCseq/peaks/
awk 'BEGIN{FS=OFS="\t"} {mean=int(($2 + $3) / 2); print $1, mean, (mean+1)}' total_DMSO_ACCseq_peaks.bed > tmp
awk '$2>2500 {print $1"\t"($2-2500)"\t"($3+2500)}' tmp > tmp2
awk 'NR==FNR {chrom_sizes[$1]=$2; next} $1 in chrom_sizes && $3 <= chrom_sizes[$1]'  hg38.chrom.sizes tmp2 > total_DMSO_ACCseq_peaks_5kb.bed

cd /media/niechen/niechen3/IRAK1_study/ACCseq/peaks/
awk 'BEGIN{FS=OFS="\t"} {mean=int(($2 + $3) / 2); print $1, mean, (mean+1)}' total_pOHT_ACCseq_peaks.bed > tmp
awk '$2>2500 {print $1"\t"($2-2500)"\t"($3+2500)}' tmp > tmp2
awk 'NR==FNR {chrom_sizes[$1]=$2; next} $1 in chrom_sizes && $3 <= chrom_sizes[$1]'  hg38.chrom.sizes tmp2 > total_pOHT_ACCseq_peaks_5kb.bed

cd /media/niechen/niechen3/IRAK1_study/ACCseq/peaks/
awk 'BEGIN{FS=OFS="\t"} {mean=int(($2 + $3) / 2); print $1, mean, (mean+1)}' total_pOHTIRAK1i_ACCseq_peaks.bed > tmp
awk '$2>2500 {print $1"\t"($2-2500)"\t"($3+2500)}' tmp > tmp2
awk 'NR==FNR {chrom_sizes[$1]=$2; next} $1 in chrom_sizes && $3 <= chrom_sizes[$1]'  hg38.chrom.sizes tmp2 > total_pOHTIRAK1i_ACCseq_peaks_5kb.bed
rm tmp*

cd /media/niechen/niechen3/IRAK1_study/ACCseq/bam_sort_merged/
bedtools multicov -bams DMSO_native.bam DMSO_fix.bam DMSO_fixHex.bam \
-bed /media/niechen/niechen3/IRAK1_study/ACCseq/peaks/total_DMSO_ACCseq_peaks_5kb.bed > /media/niechen/niechen3/IRAK1_study/ACCseq/ACCseq_in_peaks_regions/DMSO_groups_in_peaks_5kb.txt

bedtools multicov -bams pOHT_native.bam pOHT_fix.bam pOHT_fixHex.bam \
-bed /media/niechen/niechen3/IRAK1_study/ACCseq/peaks/total_pOHT_ACCseq_peaks_5kb.bed > /media/niechen/niechen3/IRAK1_study/ACCseq/ACCseq_in_peaks_regions/pOHT_groups_in_peaks_5kb.txt

bedtools multicov -bams pOHTIRAK1i_native.bam pOHTIRAK1i_fix.bam pOHTIRAK1i_fixHex.bam \
-bed /media/niechen/niechen3/IRAK1_study/ACCseq/peaks/total_pOHTIRAK1i_ACCseq_peaks_5kb.bed > /media/niechen/niechen3/IRAK1_study/ACCseq/ACCseq_in_peaks_regions/pOHTIRAK1i_groups_in_peaks_5kb.txt

