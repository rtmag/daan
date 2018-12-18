more mm10_repeatmasker.txt|grep -P -v "Simple_repeat|Low_complexity" > mm10_repeatmasker_clean.txt

python /root/myPrograms/RepEnrich/RepEnrich_setup.py /root/resources/mm10_repeat/mm10_repeatmasker_clean.txt \
/root/resources/mm10_repeat/mm10.fasta /root/resources/mm10_repeat/setup_folder_mm10

########################################################################################################################
python /root/myPrograms/RepEnrich/RepEnrich.py  \
/home/roberto/references/hg19_repeatmasker_clean.txt \
/home/roberto/deepa/novogene/repenrich HCT116_siC_DMSO \
/home/roberto/references/RepEnrich_hg19/ \
/home/roberto/deepa/novogene/repenrich_bowtie/HCT116_siC_DMSO_multimap_1.fastq \
--fastqfile2 /home/roberto/deepa/novogene/repenrich_bowtie/HCT116_siC_DMSO_multimap_2.fastq \
/home/roberto/deepa/novogene/repenrich_bowtie/HCT116_siC_DMSO_uniq.bam \
--cpus 30 --pairedend TRUE
