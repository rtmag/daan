################################################################################################
bowtie-build /root/resources/tmp/mm10.fasta /root/resources/mm10_bowtie/mm10
################################################################################################
bowtie /root/resources/mm10_bowtie/mm10 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/bowtie/OB1_multimap.fastq \
<( zcat /root/daan/fastq_trim/Ageing.OB.RNA.Rep1.20170522.R1.OK_trimmed.fq.gz ) \
/root/daan/bowtie/OB1_uniq.sam

bowtie /root/resources/mm10_bowtie/mm10 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/bowtie/OB2_multimap.fastq \
<( zcat /root/daan/fastq_trim/Ageing.OB.RNA.Rep2.20170522.R1.OK_trimmed.fq.gz ) \
/root/daan/bowtie/OB2_uniq.sam

bowtie /root/resources/mm10_bowtie/mm10 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/bowtie/WE1_multimap.fastq \
<( zcat /root/daan/fastq_trim/Ageing.WEHI.RNA.Rep1.20160606.R1.OK_trimmed.fq.gz ) \
/root/daan/bowtie/WE1_uniq.sam

bowtie /root/resources/mm10_bowtie/mm10 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/bowtie/WE2_multimap.fastq \
<( zcat /root/daan/fastq_trim/Ageing.WEHI.RNA.Rep2.20160606.R1.OK_trimmed.fq.gz ) \
/root/daan/bowtie/WE2_uniq.sam

bowtie /root/resources/mm10_bowtie/mm10 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/bowtie/YB1_multimap.fastq \
<( zcat /root/daan/fastq_trim/Ageing.YB.RNA.Rep1.20160606.R1.OK_trimmed.fq.gz ) \
/root/daan/bowtie/YB1_uniq.sam

bowtie /root/resources/mm10_bowtie/mm10 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/bowtie/YB2_multimap.fastq \
<( zcat /root/daan/fastq_trim/Ageing.YB.RNA.Rep2.20160606.R1.OK_trimmed.fq.gz ) \
/root/daan/bowtie/YB2_uniq.sam
################################################################################################
samtools view -bS OB1_uniq.sam | samtools sort - -o OB1_uniq.bam
samtools index OB1_uniq.bam

samtools view -bS OB2_uniq.sam | samtools sort - -o OB2_uniq.bam
samtools index OB2_uniq.bam

samtools view -bS WE1_uniq.sam | samtools sort - -o WE1_uniq.bam
samtools index WE1_uniq.bam

samtools view -bS WE2_uniq.sam | samtools sort - -o WE2_uniq.bam
samtools index WE2_uniq.bam

samtools view -bS YB1_uniq.sam | samtools sort - -o YB1_uniq.bam
samtools index YB1_uniq.bam

samtools view -bS YB2_uniq.sam | samtools sort - -o YB2_uniq.bam
samtools index YB2_uniq.bam
################################################################################################

