#!/bin/bash
# 下面两行路径根据实际情况更改，dir是分析的最上级文件夹路径,gemome是基因组比对的索引文件的所在路径
# The following two lines of paths should be adjusted according to your actual setup: dir is the top-level folder path for analysis, and genome is the path where the genome alignment index file is located.
dir=/mnt/NC/HiChIP
genome=/home/DDR/genome/bwa_hg38/hg38.fa
hg38_genome=hg38.genome
cd ${dir}/
mkdir cleandata bam hic FitHiChIP

# 针对raw reads进行质检
# Perform quality control on raw reads.
cd  ${dir}/rawe/
mkdir fastq_results
ls *.gz | xargs fastqc -t 12 -o  ./fastq_results/
multiqc ./fastq_results/ -n multiqc_raw -o ./multiqcresults/

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
    --qualified_quality_phred 20 \
    --low_complexity_filter \
    --complexity_threshold 30 \
    --length_required 4 \
    --thread 4
done

cd  ${dir}/cleandata/
mkdir multiqcresults
mv *.json ./multiqcresults/;
mv *.html ./multiqcresults/;

# reads比对基因组，得到bam文件,只针对pair-end reads
# Align reads to the genome and generate a BAM file, but only for paired-end reads.
ls *_1*  >1
ls *_2*  >2
paste 1 2 >config
cat config | while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    bwa mem -5SP -T 0 -t 18 ${genome} $fq1 $fq2 | pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 12 --chroms-path ${hg38_genome} | pairtools sort --tmpdir=/mnt/NC/ --nproc 18 | pairtools dedup --nproc-in 8 --nproc-out 12 --mark-dups --output-stats ../bam/$(basename -s _1.fq.gz $fq1).stats.txt | pairtools split --nproc-in 8 --nproc-out 12 --output-pairs ../bam/$(basename -s _1.fq.gz $fq1).pairs --output-sam - | samtools view -bS -@ 18 | samtools sort -@ 18 -o ../bam/$(basename -s _1.fq.gz $fq1).bam
done

cd ${dir}/bam/
ls *.bam | xargs  -i  samtools index {}

#下面是qc过程
#The follows are the process of quality control:Proximity-ligation assessment
#Suggestion：nondupe pairs cis > 1,000 bp is greater than 20% of the total mapped No-Dup pairs
txt=*.stats.txt
ls $txt | while read id
do
python3 ~/HiChiP/get_qc.py -p ${id} > $(basename -s .stats.txt $id)_qc.table
done

#将pairs文件转化成hic文件，从而可视化
#将pairs文件转化成pairs.gz文件，用于后续FitHiChIP分析
# Convert the pairs file to a hic file for visualization.  
# Compress the pairs file into a pairs.gz file for subsequent FitHiChIP analysis.
pairs=*.pairs
ls $pairs | while read id
do
java -Xmx48000m -jar ~/HiChiP/juicertools.jar pre --threads 18 ${id} ../hic/$(basename -s .pairs $id).hic ${hg38_genome}
grep -v '#' ${id} | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > $(basename -s .pairs $id).pairs.gz
done

