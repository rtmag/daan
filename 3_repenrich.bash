more mm10_repeatmasker.txt|grep -P -v "Simple_repeat|Low_complexity" > mm10_repeatmasker_clean.txt

python /root/myPrograms/RepEnrich/RepEnrich_setup.py /root/resources/mm10_repeat/mm10_repeatmasker_clean.txt \
/root/resources/mm10_repeat/mm10.fasta /root/resources/mm10_repeat/setup_folder_mm10

########################################################################################################################

python /root/myPrograms/RepEnrich/RepEnrich.py  \
/root/resources/mm10_repeat/mm10_repeatmasker_clean.txt \
/root/daan/repenrich/ OB1 \
/root/resources/mm10_repeat/setup_folder_mm10/ \
/root/daan/bowtie/OB1_multimap.fastq \
/root/daan/bowtie/OB1_uniq.bam \
--cpus 30 --pairedend FALSE

python /root/myPrograms/RepEnrich/RepEnrich.py  \
/root/resources/mm10_repeat/mm10_repeatmasker_clean.txt \
/root/daan/repenrich/ OB2 \
/root/resources/mm10_repeat/setup_folder_mm10/ \
/root/daan/bowtie/OB2_multimap.fastq \
/root/daan/bowtie/OB2_uniq.bam \
--cpus 30 --pairedend FALSE

python /root/myPrograms/RepEnrich/RepEnrich.py  \
/root/resources/mm10_repeat/mm10_repeatmasker_clean.txt \
/root/daan/repenrich/ WE1 \
/root/resources/mm10_repeat/setup_folder_mm10/ \
/root/daan/bowtie/WE1_multimap.fastq \
/root/daan/bowtie/WE1_uniq.bam \
--cpus 30 --pairedend FALSE

python /root/myPrograms/RepEnrich/RepEnrich.py  \
/root/resources/mm10_repeat/mm10_repeatmasker_clean.txt \
/root/daan/repenrich/ WE2 \
/root/resources/mm10_repeat/setup_folder_mm10/ \
/root/daan/bowtie/WE2_multimap.fastq \
/root/daan/bowtie/WE2_uniq.bam \
--cpus 30 --pairedend FALSE

python /root/myPrograms/RepEnrich/RepEnrich.py  \
/root/resources/mm10_repeat/mm10_repeatmasker_clean.txt \
/root/daan/repenrich/ YB1 \
/root/resources/mm10_repeat/setup_folder_mm10/ \
/root/daan/bowtie/YB1_multimap.fastq \
/root/daan/bowtie/YB1_uniq.bam \
--cpus 30 --pairedend FALSE

python /root/myPrograms/RepEnrich/RepEnrich.py  \
/root/resources/mm10_repeat/mm10_repeatmasker_clean.txt \
/root/daan/repenrich/ YB2 \
/root/resources/mm10_repeat/setup_folder_mm10/ \
/root/daan/bowtie/YB2_multimap.fastq \
/root/daan/bowtie/YB2_uniq.bam \
--cpus 30 --pairedend FALSE

###################################################################
##### FOR ASIA; NOTES FOR RUNNING REPENRICH

 python /root/myPrograms/RepEnrich/RepEnrich.py  \         Path to repenrich
 /root/resources/mm10_repeat/mm10_repeatmasker_clean.txt \ Path to the repeatmasker text file
 /root/daan/repenrich/ \                                   Path to directory where you want to save the results  
 YB2 \                                                     Label of your sample
 /root/resources/mm10_repeat/setup_folder_mm10/ \          Path to rep_enrich directory (the tar file you downloaded)
 /root/daan/bowtie/YB2_multimap.fastq \                    Path to multimap fastq
 /root/daan/bowtie/YB2_uniq.bam \                          Path to uniquely mapped reads bam file
 --cpus 16 --pairedend TRUE                               Number of CPUs and whether your samples are PAIREDEND
 
 



