#!/bin/bash
# 下面两行路径根据实际情况更改，dir是分析的最上级文件夹路径,index是基因组比对的索引文件的所在路径
# The following two lines of paths should be adjusted according to your actual setup: `dir` is the top-level folder path for analysis, and `index` is the path where the genome alignment index file is located.
dir=/mnt/NC/WGRS_202503
index=/home/DDR/genome/bwa_hg38/hg38.fa
cd ${dir}/
mkdir cleandata bam bam_markdup vcf vcf_annotation
cd ${dir}/raw/
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
    trim_galore -j 4 -q 25  --phred33 --length 35 --stringency 3 --paired --gzip -o ${dir}/cleandata/ $fq1 $fq2 
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
	  sample=$(basename -s _1_val_1.fq.gz $fq1)
    bwa mem -t 20 -R "@RG\tID:$sample\tSM:$sample\tLB:WGRS\tPL:Illumina" $index $fq1 $fq2 | samtools view -Sb | samtools sort -@ 20 -O bam  > ${dir}/bam/$(basename -s _1_val_1.fq.gz $fq1).bam
done

cd ${dir}/bam/
ls *.bam | xargs  -i  samtools index {}
bam_stat=*bam
ls $bam_stat | while read id
do
samtools flagstat $id > $(basename -s .bam $id).stat
done

####gatk获取VCF####
####Get VCF by gatk####
cd ${dir}/bam/;
bam=*.bam
ls $bam | while read id
do
  gatk --java-options "-XX:ParallelGCThreads=18 -Xmx20G -Djava.io.tmpdir=./" CollectWgsMetrics -I ${id} -R $index -O $(basename -s .bam $id).metrics
  gatk --java-options "-XX:ParallelGCThreads=18" MarkDuplicates -I ${id} -M ${dir}/bam_markdup/$(basename -s .bam $id).markdup_metrics.txt -O ${dir}/bam_markdup/$(basename -s .bam $id).rmdup.bam
done

cd ${dir}/bam_markdup/;
ls *.rmdup.bam | xargs  -i  samtools index {};
bam_markdups=*.rmdup.bam
ls $bam_markdups | while read id
do
  gatk --java-options "-XX:ParallelGCThreads=18 -Xmx20G -Djava.io.tmpdir=./" BaseRecalibrator \
  -R $index -I ${id} \
  --known-sites /home/DDR/genome/bwa_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-sites /home/DDR/genome/bwa_hg38/dbsnp_146.hg38.vcf.gz \
  --sequence-dictionary /home/DDR/genome/bwa_hg38/hg38.dict \
  -O $(basename -s .rmdup.bam $id).recal_data.table  
  gatk --java-options "-XX:ParallelGCThreads=18 -Xmx20G -Djava.io.tmpdir=./" ApplyBQSR \
  -R $index -I ${id} \
  -bqsr $(basename -s .rmdup.bam $id).recal_data.table \
  -O $(basename -s .rmdup.bam $id).BQSR.bam
  gatk GetPileupSummaries \
  -I $(basename -s .rmdup.bam $id).BQSR.bam \
  -V /home/DDR/genome/bwa_hg38/af_only_gnomad.hg38.INDEL_SNP_SYMBOLIC.biallelic.vcf.gz \
  -L /home/DDR/genome/bwa_hg38/hg38.interval_list \
  -O $(basename -s .rmdup.bam $id).pileups.table 
  gatk --java-options "-XX:ParallelGCThreads=18" CalculateContamination \
  -I $(basename -s .rmdup.bam $id).pileups.table \
  -O ${dir}/vcf/$(basename -s .rmdup.bam $id).calculatecontamination.table 
  gatk --java-options "-XX:ParallelGCThreads=18 -Xmx20G -Djava.io.tmpdir=./" Mutect2 \
  -R $index \
  -I $(basename -s .rmdup.bam $id).BQSR.bam \
  -tumor $(basename -s .rmdup.bam $id) \
  --germline-resource /home/DDR/genome/bwa_hg38/af_only_gnomad.hg38.vcf.gz \  #download from GATK database
  --af-of-alleles-not-in-resource 0.0000025 \
  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
  --bam-output ${dir}/vcf/$(basename -s .rmdup.bam $id).Mutect2.bam \
  --native-pair-hmm-threads 15 \
  -O ${dir}/vcf/$(basename -s .rmdup.bam $id).Mutect2.vcf
  gatk FilterMutectCalls \
  -V ${dir}/vcf/$(basename -s .rmdup.bam $id).Mutect2.vcf \
  -R $index \
  --contamination-table ${dir}/vcf/$(basename -s .rmdup.bam $id).calculatecontamination.table \
  -O ${dir}/vcf/$(basename -s .Mutect2.vcf $id).filtered.Mutect2.vcf
done

cd ${dir}/vcf/
VCF=*.filtered.Mutect2.vcf
ls $VCF | while read id
do
awk '$1 ~ /^#/ || $7 == "PASS" {print}' ${id} > $(basename -s .filtered.Mutect2.vcf $id).PASS.vcf
done

cd /home/DDR/bin/annovar/humandb/
annotate_variation.pl -downdb -buildver hg38 -webfrom annovar avsnp150 ./
annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene ./
annotate_variation.pl -downdb -buildver hg38 -webfrom annovar exac03 ./
annotate_variation.pl -downdb -buildver hg38 -webfrom annovar esp6500siv2_all ./
annotate_variation.pl -downdb -buildver hg38 -webfrom annovar clinvar_20220320 ./
annotate_variation.pl -downdb -buildver hg38 -webfrom annovar dbnsfp30a ./
annotate_variation.pl -downdb -buildver hg38 -webfrom annovar cosmic70 ./

cd ${dir}/vcf/
VCF_filtered=*.PASS.vcf
ls $VCF_filtered | while read id
do
  convert2annovar.pl -format vcf4 -allsample -withfreq ${id} > $(basename -s .PASS.vcf $id).avinput
  table_annovar.pl $(basename -s .PASS.vcf $id).avinput /home/DDR/bin/annovar/humandb/ \
  -buildver hg38 -out ../vcf_annotation/$(basename -s .PASS.vcf $id) \
  -remove -protocol refGene,avsnp150,exac03,esp6500siv2_all,clinvar_20220320,dbnsfp30a,cosmic70 \
  -operation g,f,f,f,f,f,f -nastring NA -csvout
done

  
