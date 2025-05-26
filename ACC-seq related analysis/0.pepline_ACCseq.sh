#!/bin/bash
# 下面两行路径根据实际情况更改，dir是分析的最上级文件夹路径,gemome是基因组比对的索引文件的所在路径
# The following two lines of paths should be adjusted according to your actual setup: dir is the top-level folder path for analysis, and genome is the path where the genome alignment index file is located.
dir=/media/niechen/niechen3/IRAK1_study/ACC_seq
genome=/media/niechen/niechen3/STAR_index/
cd ${dir}/
mkdir cleandata bam bam_sort bw bam_rmdup peaks

cd  ${dir}/raw/
mkdir fastq_results;
ls *.gz | xargs fastqc -t 12 -o  ./fastq_results/;
multiqc ./fastq_results/ -n multiqc_raw -o ./multiqcresults/;

cd ${dir}/raw_merge/;
paired
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
mkdir trimed;
ls *.gz | xargs fastqc -t 12 -o ./trimed;
multiqc ./trimed -n multi -o ./multiqcresults/;
mv *.txt ./multiqcresults/;
mv *.json ./multiqcresults/;
mv *.html ./multiqcresults/;

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
    --runThreadN 18 \
    --genomeDir ${genome} \
    --readFilesCommand zcat \
    --readFilesIn $fq1 $fq2 \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped None \
    --outFileNamePrefix ${dir}/bam/$(basename -s _1.fq.gz $fq1).bam \
    --outFilterMismatchNmax 10 \
    --outTmpDir /home/niechen/niechen \
    --outFilterMatchNmin 10 \
    --alignEndsType Local \
    --seedSearchLmax 10 \
    --outFilterScoreMinOverLread 0.05 \
    --outFilterMatchNminOverLread 0.05
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
samtools merge -o ../bam_sort_merged/DMSO_native.bam  1DMSO_native_R1.rmdupsorted.bam 2DMSO_native_R2.rmdupsorted.bam
samtools merge -o ../bam_sort_merged/pOHT_native.bam  3pOHT_native_R1.rmdupsorted.bam 4pOHT_native_R2.rmdupsorted.bam
samtools merge -o ../bam_sort_merged/pOHTIRAK1i_native.bam  5pOHTIRAK1i_native_R1.rmdupsorted.bam 6pOHTIRAK1i_native_R2.rmdupsorted.bam
samtools merge -o ../bam_sort_merged/DMSO_fix.bam  7DMSO_fix_R1.rmdupsorted.bam 8DMSO_fix_R2.rmdupsorted.bam
samtools merge -o ../bam_sort_merged/pOHT_fix.bam  9pOHT_fix_R1.rmdupsorted.bam 10pOHT_fix_R2.rmdupsorted.bam
samtools merge -o ../bam_sort_merged/pOHTIRAK1i_fix.bam  11pOHTIRAK1i_fix_R1.rmdupsorted.bam 12pOHTIRAK1i_fix_R2.rmdupsorted.bam
samtools merge -o ../bam_sort_merged/DMSO_fixHex.bam  13DMSO_fixHex_R1.rmdupsorted.bam 14DMSO_fixHex_R2.rmdupsorted.bam
samtools merge -o ../bam_sort_merged/pOHT_fixHex.bam  15pOHT_fixHex_R1.rmdupsorted.bam 16pOHT_fixHex_R2.rmdupsorted.bam
samtools merge -o ../bam_sort_merged/pOHTIRAK1i_fixHex.bam  17pOHTIRAK1i_fixHex_R1.rmdupsorted.bam 18pOHTIRAK1i_fixHex_R2.rmdupsorted.bam
cd ${dir}/bam_sort_merged/
ls *.bam | xargs  -i  samtools index {};
file_bamsort2=*bam
ls $file_bamsort2 | while read id
do
samtools flagstat $id > $(basename -s .bam $id).stat
done
sample2=*.bam
ls $sample2 | while read id
do
echo $id
bamCoverage --normalizeUsing CPM -b $id -o ${dir}/bw_merged/$(basename -s .bam $id).bw
done 

