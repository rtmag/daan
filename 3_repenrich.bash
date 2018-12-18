
python /home/roberto/myPrograms/RepEnrich/RepEnrich.py \
/home/roberto/references/hg19_repeatmasker_clean.txt \
/home/roberto/deepa/novogene/repenrich HCT116_siC_DMSO \
/home/roberto/references/RepEnrich_hg19/ \
/home/roberto/deepa/novogene/repenrich_bowtie/HCT116_siC_DMSO_multimap_1.fastq \
--fastqfile2 /home/roberto/deepa/novogene/repenrich_bowtie/HCT116_siC_DMSO_multimap_2.fastq \
/home/roberto/deepa/novogene/repenrich_bowtie/HCT116_siC_DMSO_uniq.bam \
--cpus 30 --pairedend TRUE
