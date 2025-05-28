#!/bin/bash
# 下面两行路径根据实际情况更改，dir是分析的最上级文件夹路径,index_hg38是基因组比对的索引文件的所在路径
# The following two lines of paths should be adjusted according to your actual setup: `dir` is the top-level folder path for analysis, and `index_hg38` is the path where the genome alignment index file is located.
dir=/media/niechen/niechen3/small_RNAseq
index_hg38=/media/niechen/niechen3/refs_smallRNA/Homo_sapiens_v86/Homo_sapiens_v86.chromosomes.fasta
cd ${dir}/;
mkdir cleandata bam bam_sort bam_clear counts_clear

cd ${dir}/raw/;
mkdir raw_trimed;
ls *.gz | xargs fastqc -t 6 -o ./raw_trimed;
multiqc ./raw_trimed -n multiraw -o ./multiqcresults_raw/;

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

cd ${dir}/cleandata/;
mkdir clean_trimed;
ls *.gz | xargs fastqc -t 6 -o ./clean_trimed;
multiqc ./clean_trimed -n multiclean -o ./multiqcresults_clean/;
mv *.txt ./multiqcresults_clean;
mv *.json ./multiqcresults_clean;
mv *.html ./multiqcresults_clean;

ls *_1*  >1
ls *_2*  >2
paste 1 2 >config
cat config | while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    bowtie2 -p 12 -x $index_hg38 -1 $fq1 -2 $fq2 -D 15 -R 2 -N 1 -k 1 -L 12 -i S,1,1.15 | samtools view -Sb  > ../bam/$(basename -s _1.fq.gz $fq1).bam
done

cd ${dir}/bam/;
files=*.bam
ls $files | while read id
do
 samtools sort -@ 12 -O bam -o ${dir}/bam_sort/$(basename -s .bam $id)sorted.bam  ${id}
done

cd ${dir}/bam_sort/;
ls *.bam | xargs  -i  samtools index {};

# 清除chr21的rDNA附近的reads
# clear reads near the rDNA region on chromosome 21
file_bamsort2=*bam
ls $file_bamsort2 | while read id
do
bedtools intersect -a ${id} -b ../black_bamCoverage.bed -v > ../bam_clear/$(basename -s sorted.bam $id).filted.bam
done

cd ${dir}/bam_clear/
ls *.bam | xargs  -i  samtools index {};
file_bamsort3=*bam
ls $file_bamsort3 | while read id
do
samtools flagstat $id > $(basename -s .bam $id).stat
done

# 计算PartI library的所有Gene的counts数
# Calculate the count numbers for all genes in PartI library
# `Homo_sapiens.GRCh38.108.bed` is listed in `RAP-seq related analysis` part.
cd /media/niechen/niechen3/small_RNAseq_great/bam_clear/
bedtools multicov -bams 1DMSOG.filted.bam 2DMSOG.filted.bam 3DMSOG.filted.bam \
4IRAK1iG.filted.bam 5IRAK1iG.filted.bam 6IRAK1iG.filted.bam \
7IRDMSOG.filted.bam 8IRDMSOG.filted.bam \
10IRIRAK1iG.filted.bam 11IRIRAK1iG.filted.bam 12IRIRAK1iG.filted.bam \
-bed Homo_sapiens.GRCh38.108.bed > ../counts_clear/smallRNA_great_clear_all_gene.txt

# 计算PartII library的所有Gene的counts数
# Calculate the count numbers for all genes in PartII library
cd ${dir}/bam_clear/
bedtools multicov -bams 1DMSOB.filted.bam 2DMSOB.filted.bam 3DMSOB.filted.bam \
4IRAK1iB.filted.bam 5IRAK1iB.filted.bam 6IRAK1iB.filted.bam \
7IRDMSOB.filted.bam 8IRDMSOB.filted.bam \
10IRIRAK1iB.filted.bam 11IRIRAK1iB.filted.bam 12IRIRAK1iB.filted.bam \
-bed Homo_sapiens.GRCh38.108.bed > ../counts_clear/smallRNA_big_clear_all_gene.txt

