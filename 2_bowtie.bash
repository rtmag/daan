################################################################################################
bowtie /root/resources/hg38_bowtie/hg38 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/bowtie/OB1_multimap.fastq \
<( zcat /root/daan/fastq_trim/Ageing.OB.RNA.Rep1.20170522.R1.OK_trimmed.fq.gz ) \
/root/daan/bowtie/OB1_uniq.sam

bowtie /root/resources/hg38_bowtie/hg38 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/bowtie/OB2_multimap.fastq \
<( zcat /root/daan/fastq_trim/Ageing.OB.RNA.Rep2.20170522.R1.OK_trimmed.fq.gz ) \
/root/daan/bowtie/OB2_uniq.sam

bowtie /root/resources/hg38_bowtie/hg38 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/bowtie/WE1_multimap.fastq \
<( zcat /root/daan/fastq_trim/Ageing.WEHI.RNA.Rep1.20160606.R1.OK_trimmed.fq.gz ) \
/root/daan/bowtie/WE1_uniq.sam

bowtie /root/resources/hg38_bowtie/hg38 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/bowtie/WE2_multimap.fastq \
<( zcat /root/daan/fastq_trim/Ageing.WEHI.RNA.Rep2.20160606.R1.OK_trimmed.fq.gz ) \
/root/daan/bowtie/WE2_uniq.sam

bowtie /root/resources/hg38_bowtie/hg38 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/bowtie/YB1_multimap.fastq \
<( zcat /root/daan/fastq_trim/Ageing.YB.RNA.Rep1.20160606.R1.OK_trimmed.fq.gz ) \
/root/daan/bowtie/YB1_uniq.sam

bowtie /root/resources/hg38_bowtie/hg38 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/bowtie/YB2_multimap.fastq \
<( zcat /root/daan/fastq_trim/Ageing.YB.RNA.Rep2.20160606.R1.OK_trimmed.fq.gz ) \
/root/daan/bowtie/YB2_uniq.sam
################################################################################################


samtools view -bS HCT_siControl_HWN2YCCXX_L5_uniq.sam | samtools sort - HCT_siControl_HWN2YCCXX_L5_uniq
samtools index HCT_siControl_HWN2YCCXX_L5_uniq.bam


python /home/roberto/myPrograms/RepEnrich/RepEnrich.py \
/home/roberto/references/hg19_repeatmasker_clean.txt \
/home/roberto/deepa/novogene/repenrich HCT116_siC_DMSO \
/home/roberto/references/RepEnrich_hg19/ \
/home/roberto/deepa/novogene/repenrich_bowtie/HCT116_siC_DMSO_multimap_1.fastq \
--fastqfile2 /home/roberto/deepa/novogene/repenrich_bowtie/HCT116_siC_DMSO_multimap_2.fastq \
/home/roberto/deepa/novogene/repenrich_bowtie/HCT116_siC_DMSO_uniq.bam \
--cpus 30 --pairedend TRUE
