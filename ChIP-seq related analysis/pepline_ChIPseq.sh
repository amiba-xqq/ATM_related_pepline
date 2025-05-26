#!/bin/bash
# 下面两行路径根据实际情况更改，dir是分析的最上级文件夹路径,gemome是基因组比对的索引文件的所在路径
# The following two lines of paths should be adjusted according to your actual setup: dir is the top-level folder path for analysis, and genome is the path where the genome alignment index file is located.
dir=/media/niechen/niechen3/IRAK1_study/pNBS1_ChIPseq_DIvA
genome=/home/DDR/genome/hisat2_GRCh38/hisat2_GRCh38
# 创建各个文件夹存放ChIP-seq各步骤产生文件
# Create directories to store files generated at each step of the ChIP-seq process.
cd ${dir}/
mkdir cleandata bam bam_sort bw bam_rmdup peaks

# 针对raw reads进行质检
# Quality control for raw reads.
cd  ${dir}/raw/
mkdir fastq_results
ls *.gz | xargs fastqc -t 12 -o  ./fastq_results/
multiqc ./fastq_results/ -n multiqc_raw -o ./multiqcresults/

# 清洗reads,只针对pair-end reads
# Clean reads, only for paired-end reads.
cd ${dir}/raw/
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

# 对clean reads做质检
# Quality control on cleaned reads.
cd ${dir}/cleandata/
mkdir trimed
ls *.gz | xargs fastqc -t 18 -o ./trimed
multiqc ./trimed -n multi -o ./multiqcresults/
mv *.txt ./multiqcresults/

# reads比对基因组，得到bam文件,只针对pair-end reads
# Align reads to the genome to obtain BAM files, only for paired-end reads.
ls *_1*  >1
ls *_2*  >2
paste 1 2 >config
cat config | while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    hisat2 -p 18 --dta-cufflinks --un-conc -x ${genome} -1 $fq1 -2 $fq2 | samtools view -Sb > ${dir}/bam/$(basename -s _1_val_1.fq.gz $fq1).bam
done

# 去除重复的reads
# Remove duplicate reads.
cd ${dir}/bam/
files_bam=*.bam
ls $files_bam | while read id 
do
echo $id
sambamba markdup -r -p -t 18 $id ${dir}/bam_rmdup/$(basename -s .bam $id).rmdup.bam
done

# 对bam文件进行排序
# Sort the BAM files.
cd ${dir}/bam_rmdup/
files=*.bam
ls $files | while read id
do
 samtools sort -@ 18 -O bam -o ${dir}/bam_sort/$(basename -s .bam $id)sorted.bam  ${id}
done

# 根据情况决定是否需要合并生物学重复组（需要是sorted好的bam）,合并得到的merged的bam文件再导出为bigwig
# Determine whether to merge biological replicates based on the situation (sorted BAM files required), and export the merged BAM file as a BigWig file.
cd ${dir}/bam_sort/
samtools merge -o pNBS1_4OHTIRAK1i_IP.merged.bam pNBS1_4OHTIRAK1i_IP_rep1.rmdupsorted.bam pNBS1_4OHTIRAK1i_IP_rep2.rmdupsorted.bam
samtools merge -o pNBS1_4OHT_IP.merged.bam pNBS1_4OHT_IP_rep1.rmdupsorted.bam pNBS1_4OHT_IP_rep2.rmdupsorted.bam

# bam文件添加索引
# Add index to BAM files.
cd ${dir}/bam_sort/
ls *.bam | xargs  -i  samtools index {}
# 计算bam文件比对率
# Calculate the alignment rate of BAM files.
cd ${dir}/bam_sort/
file_bamsort=*bam
ls $file_bamsort | while read id
do
samtools flagstat $id > $(basename -s .bam $id).stat
done

# 导出bigwig文件
# Export BigWig files.
cd ${dir}/bam_sort/
sample=*.bam
ls $sample | while read id
do
echo $id
bamCoverage -p 18 --normalizeUsing CPM -b $id -o ${dir}/bw/$(basename -s .bam $id).bw
done 
