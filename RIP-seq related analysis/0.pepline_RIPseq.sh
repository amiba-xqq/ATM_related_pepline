#!/bin/bash
# 下面两行路径根据实际情况更改，dir是分析的最上级文件夹路径,gemome是基因组比对的索引文件的所在路径
# The following two lines of paths should be adjusted according to your actual setup: dir is the top-level folder path for analysis, and genome is the path where the genome alignment index file is located.
dir=/media/niechen/niechen2/NBS1_RIPseq/
genome=/media/niechen/niechen3/STAR_index/
cd ${dir}/
mkdir cleandata bam bam_sort bw bam_rmdup peaks

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
    time fastp \
    --in1 $fq1 \
    --in2 $fq2 \
    --out1 ../cleandata/$(basename -s .fq.gz $fq1).fq.gz \
    --out2 ../cleandata/$(basename -s .fq.gz $fq2).fq.gz \
    --json ../cleandata/$(basename -s _1.fq.gz $fq1).json \
    --html ../cleandata/$(basename -s _1.fq.gz $fq1).html \
    --trim_poly_g --poly_g_min_len 6 \
    --trim_poly_x --poly_x_min_len 6 \
    --cut_front --cut_tail --cut_window_size 4 \
    --qualified_quality_phred 15 \
    --low_complexity_filter \
    --complexity_threshold 30 \
    --length_required 4 \
    --thread 4
done

cd ${dir}/cleandata/
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
    --runThreadN 14 \
    --genomeDir ${genome} \
    --readFilesCommand zcat \
    --readFilesIn $fq1 $fq2 \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped None \
    --outFileNamePrefix ${dir}/bam_STAR/$(basename -s _1.fq.gz $fq1).bam \
    --outFilterMismatchNmax 10 \
    --outTmpDir /home/niechen/niechen \
    --outFilterMatchNmin 10 \
    --alignEndsType Local \
    --seedSearchLmax 10 \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0
done

cd ${dir}/bam_STAR/;
files=*.bam
ls $files | while read id
do
 samtools sort -@ 18 -O bam -o ${dir}/bam_STAR_sort/$(basename -s .bamAligned.out.bam $id)sorted.bam  ${id}
done

cd ${dir}/bam_STAR_sort/;
ls *.bam | xargs  -i  samtools index {}

file_bamsort2=*bam
ls $file_bamsort2 | while read id
do
bedtools intersect -a ${id} -b ../black_bamCoverage.bed -v > ../bam_STAR_sort_clear/$(basename -s sorted.bam $id).filted.bam
done

cd ${dir}/bam_STAR_sort_clear/
ls *.bam | xargs  -i  samtools index {};
file_bamsort3=*bam
ls $file_bamsort3 | while read id
do
samtools flagstat $id > $(basename -s .bam $id).stat
done
