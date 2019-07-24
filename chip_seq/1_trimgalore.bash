trim_galore --illumina -q 20 --fastqc -o /root/daan/fastq_trim/ \
Ageing.OB.RNA.Rep1.20170522.R1.OK.fastq.gz &

##################
#ALL Variant NO NORMAL
for bamfile in /home/rtm/vivek/navi/wes/bam/*_recalibrated.bam;
do ls -lh $bamfile; 
name=${bamfile//\/home\/rtm\/vivek\/navi\/wes\/bam\/} ;
sample=${name//\_recalibrated\.bam} ;
java -Xmx200G -jar /home/rtm/myprograms/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar Mutect2 \
-R /home/references/broadhg38/broad_hg38/Homo_sapiens_assembly38.fasta \
-I $bamfile \
-tumor $sample \
-O /home/rtm/vivek/navi/wes/test_vcf/all_$sample.vcf.gz
done
##################
