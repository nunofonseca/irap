#!/bin/bash
# Create a fungi and microbial contamination DB
set -e
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/fungi.16.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archea/microbial.69.1.genomic.fna.gz

gunzip fungi.16.1.genomic.fna.gz
gunzip microbial.69.1.genomic.fna.gz

cat fungi.16.1.genomic.fna microbial.69.1.genomic.fna > combined.fa
irap_add2contaminationDB -n fungi_16_1_microbial_99_1 -f combined.fa
rm -f fungi.16.1.genomic.fna microbial.69.1.genomic.fna combined.fa

echo "IRAP: cont_index=fungi_16_1_microbial_99_1"
exit 0
