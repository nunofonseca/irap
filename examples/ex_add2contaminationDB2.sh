#!/bin/bash
# Create a ecoli contamination DB to be used with bowtie2
set -ea
wget ftp://ftp.ensemblgenomes.org/pub/release-37/bacteria/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa.gz -O ecoli.fa.gz
gunzip ecoli.fa.gz
irap_add2contaminationDB -n e_coli -f ecoli.fa -d ../data/contamination/ -m bowtie2
echo "IRAP: cont_index=e_coli"
exit 0
