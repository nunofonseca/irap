#!/bin/bash
# Create a ecoli contamination DB to be used with bowtie2
set -e
wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_536_uid58531/NC_008253.fna
irap_add2contaminationDB -n e_coli -f NC_008253.fna -d ../data/contamination/ -m bowtie2
echo "IRAP: cont_index=e_coli"
exit 0
