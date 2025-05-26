#!/bin/bash
# 下面两行路径根据实际情况更改，dir是分析的最上级文件夹路径,gemome是基因组比对的索引文件的所在路径
# The following two lines of paths should be adjusted according to your actual setup: dir is the top-level folder path for analysis, and genome is the path where the genome alignment index file is located.
dir=/media/niechen/niechen3/IRAK1_study/GROseq_U2OS
genome=/media/niechen/niechen3/refs_smallRNA/Homo_sapiens.GRCh38.86.chr.gtf
#cd ${dir}/
#cat > idname
#SRR1916552
#SRR1916554

# 下载GRO-seq原始数据SRR1916552/SRR1916554，并转化为FASTQ格式
# Download raw data SRR1916552/SRR1916554 of GRO-seq and convert it to FASTQ format
echo -e "\n \n \n prefetch sra !!! \n \n \n "
date
mkdir raw;
cd ${dir}/raw/;
mkdir sra fq qc1;
cd ${dir}/raw/sra/
pwd

cat  ${dir}/idname | while read id ; \
do
      ( prefetch -X 30G -p -O ./ $id & )
done

date
echo  -e "\n \n \n  111  move files !!! \n \n \n  "
cd ${dir}/raw/sra/
cat ${dir}/idname | while read id
do
mv $id/*  ./
rm -rf $id/
done
date

echo  -e "\n \n \n  111  sra>>>fq !!! \n \n \n  "
cd ${dir}/raw/fq/
pwd
ls ../sra/*.sra |while read id 
do
echo " PROCESS $(basename $id) "
fasterq-dump -3 -e 12 -O ./ $id
pigz -p 12 ${dir}/raw/fq/*q
done
date

cd ${dir}/;
mkdir cleandata bam bam_sort bw counts;

cd ${dir}/raw/fq/;
rename 's/fastq/fq/' *
rename 's/SRR1916552/GROseq_r1/' *
rename 's/SRR1916554/GROseq_r2/' *
mkdir raw_trimed;
ls *.gz | xargs fastqc -t 6 -o ./raw_trimed;
multiqc ./raw_trimed -n multiraw -o ./multiqcresults_raw/;
rawdata=*.fq.gz
ls $rawdata | while read id
do
    time fastp \
    --in1 $id \
    --out1 ../../cleandata/$(basename -s .fq.gz $id).fq.gz \
    --json ../../cleandata/$(basename -s .fq.gz $id).json \
    --html ../../cleandata/$(basename -s .fq.gz $id).html \
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

cd ${dir}/cleandata/;
mv *.txt ./multiqcresults_clean;
mv *.json ./multiqcresults_clean;
mv *.html ./multiqcresults_clean;
index_hg38=/media/niechen/niechen3/refs_smallRNA/Homo_sapiens_v86/Homo_sapiens_v86.chromosomes.fasta
clean=*.fq.gz
ls $clean | while read id
do
    bowtie2 -p 12 -x $index_hg38 -U $id -D 15 -R 2 -N 1 -k 1 -L 12 -i S,1,1.15 | samtools view -Sb  > ../bam/$(basename -s .fq.gz $id).bam
done

cd ${dir}/bam/;
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

cd ${dir}/bam_sort/
file_bamsort=*bam
ls $file_bamsort | while read id
do
samtools flagstat $id > $(basename -s .bam $id).stat
done
featureCounts -T 8 -a ${genome} -o ../counts/counts_GROseq_U2OS.txt  *.bam
