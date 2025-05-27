#!/bin/bash
# 下面两行路径根据实际情况更改，dir是分析的最上级文件夹路径,gemome是基因组比对的索引文件的所在路径
# The following two lines of paths should be adjusted according to your actual setup: dir is the top-level folder path for analysis, and genome is the path where the genome alignment index file is located.
dir=/media/niechen/niechen3/IRAK1_study/RAPseq
genome=/media/niechen/niechen3/STAR_index/
cd ${dir}/
mkdir cleandata bam bam_sort bw bam_rmdup peaks;

cd  ${dir}/raw/
mkdir fastq_results;
ls *.gz | xargs fastqc -t 12 -o  ./fastq_results/;
multiqc ./fastq_results/ -n multiqc_raw -o ./multiqcresults/;

ls *_1*  >1
ls *_2*  >2
paste 1 2 >config
cat config | while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    trim_galore -j 4 -q 25  --phred33 --length 20 --stringency 4 --paired --gzip -o ${dir}/cleandata/ $fq1 $fq2 
done

cd ${dir}/cleandata/;
mkdir trimed;
ls *.gz | xargs fastqc -t 12 -o ./trimed;
multiqc ./trimed -n multi -o ./multiqcresults/;
mv *.txt ./multiqcresults/;

ls *_1*  >1
ls *_2*  >2
paste 1 2 >config
cat config | while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    STAR \
    --runMode alignReads \
    --runThreadN 18 \
    --genomeDir ${genome} \
    --readFilesCommand zcat \
    --readFilesIn $fq1 $fq2 \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped None \
    --outFileNamePrefix ${dir}/bam/$(basename -s _val_1.fq.gz $fq1).bam \
    --outFilterMismatchNmax 10 \
    --outTmpDir /home/niechen/niechen \
    --outFilterMatchNmin 10 \
    --alignEndsType Local \
    --seedSearchLmax 10 \
    --outFilterScoreMinOverLread 0.1 \
    --outFilterMatchNminOverLread 0.1
done

cd ${dir}/bam/;
files_bam=*.bam
ls $files_bam | while read id 
do
echo $id
sambamba markdup -r -p -t 6 $id ${dir}/bam_rmdup/$(basename -s .bamAligned.out.bam $id).rmdup.bam
done

cd ${dir}/bam_rmdup/;
files=*.bam
ls $files | while read id
do
 samtools sort -@ 12 -O bam -o ${dir}/bam_sort/$(basename -s .bam $id)sorted.bam  ${id}
done
cd ${dir}/bam_sort/;
ls *.bam | xargs  -i  samtools index {};

cd ${dir}/bam_sort/;
sample=*.bam
ls $sample | while read id
do
echo $id
bamCoverage --normalizeUsing CPM -b $id -o ${dir}/bw/$(basename -s .bam $id).bw
done 

cd ${dir}/bam_sort/;
file_bamsort=*bam
ls $file_bamsort | while read id
do
samtools flagstat $id > $(basename -s .bam $id).stat
done

cd ${dir}/bam_sort/
samtools merge -o 1U1.merged.bam 1U1-r1.rmdupsorted.bam 1U1-r2.rmdupsorted.bam
samtools merge -o 1U6.merged.bam 1U6-r1.rmdupsorted.bam 1U6-r2.rmdupsorted.bam
samtools merge -o 2U1.merged.bam 2U1-r1.rmdupsorted.bam 2U1-r2.rmdupsorted.bam
samtools merge -o 2U2.merged.bam 2U2-r1.rmdupsorted.bam 2U2-r2.rmdupsorted.bam
samtools merge -o 2U6.merged.bam 2U6-r1.rmdupsorted.bam 2U6-r2.rmdupsorted.bam
samtools merge -o 3U1.merged.bam 3U1-r1.rmdupsorted.bam 3U1-r2.rmdupsorted.bam
samtools merge -o 3U6.merged.bam 3U6-r1.rmdupsorted.bam 3U6-r2.rmdupsorted.bam
ls *merged.bam | xargs  -i  samtools index {};
file_bamsort2=*merged.bam
ls $file_bamsort2 | while read id
do
samtools flagstat $id > $(basename -s .bam $id).stat
done
ls $file_bamsort2 | while read id
do
echo $id
bamCoverage --normalizeUsing CPM -b $id -o ${dir}/bw/$(basename -s .bam $id).bw
done 

# 计算各组RAPseq在TSS±5kb范围分布的counts
# Calculate distribution of RAPseq counts within TSS ±5kb range for each group.
# `ucsc.refseq.tss.sort.bed` is downloaded from UCSC website, marks the TSS sites of all gene in hg38 genome.
cd /media/niechen/niechen3/IRAK1_study/RAPseq/bam_sort/
bedtools multicov -bams 1input.rmdupsorted.bam 1U1.merged.bam 1U2-r1.rmdupsorted.bam 1U6.merged.bam \
2input.rmdupsorted.bam 2U1.merged.bam 2U2.merged.bam 2U6.merged.bam \
3input.rmdupsorted.bam 3U1.merged.bam 3U2-r1.rmdupsorted.bam 3U6.merged.bam \
-bed ucsc.refseq.tss.sort.bed > /media/niechen/niechen3/IRAK1_study/RAPseq/RAPseq_and_pATM_ChIPseq_relationship/RAPseq_and_pATM_ChIPseq_in_TSS_10kb.txt
