##################
#TRIM
for bamfile in /root/daan/chip-seq/fastq/*fastq.gz ;
do echo $bamfile; 
name=${bamfile//\/root\/daan\/chip-seq\/fastq\/Ageing\.} ;
sample=$(echo $name|cut -f1,2,3 -d'.') ;
trim_galore --illumina -q 20 --fastqc -o /root/daan/chip-seq/fastq_trim/ $bamfile ;
done
##################
##################
#Bowtie 
for fastq in /root/daan/chip-seq/fastq_trim/*_trimmed.fq.gz ;
do echo $fastq; 
name=${fastq//\/root\/daan\/chip-seq\/fastq_trim\/Ageing\.} ;
name=${name//\.R1\_trimmed\.fq\.gz} ;
bowtie /root/resources/mm10_bowtie/mm10 -p 30 -t -m 1 -S --chunkmbs 4000 \
--max /root/daan/chip-seq/bowtie/$name\_multimap.fastq \
<( zcat $fastq ) \
/root/daan/chip-seq/bowtie/$name\_uniq.sam &> /root/daan/chip-seq/bowtie/$name.bowtie.log ;
done
##################
##################
#SAMTOOLS
for sam in /root/daan/chip-seq/bowtie/*_uniq.sam ;
do echo $sam; 
samtools view -bS $sam | samtools sort - -o sam.bam ;
samtools index sam.bam ;
done
##################
#RepEnrich
for rep in /root/daan/chip-seq/bowtie/*_uniq.sam ;
do echo "####################################"; 
name=${rep//\/root\/daan\/chip-seq\/bowtie\/} ;
name=${name//\_uniq\.sam} ;
echo $name;
python /root/myPrograms/RepEnrich/RepEnrich.py  \
/root/resources/mm10_repeat/mm10_repeatmasker_clean.txt \
/root/daan/chip-seq/repenrich/ $name \
/root/resources/mm10_repeat/setup_folder_mm10/ \
/root/daan/chip-seq/bowtie/$name\_multimap.fastq \
/root/daan/chip-seq/bowtie/$name\_uniq.sam.bam \
--cpus 30 --pairedend FALSE ;
done
##################
##################
