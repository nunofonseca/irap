#!/bin/bash
# Create a fungi and microbial contamination DB
set -e
mkdir -p $IRAP_DIR/data/reference/ecoli_k12
mkdir -p $IRAP_DIR/data/raw_data/ecoli_k12

pushd $IRAP_DIR/data/reference/ecoli_k12
wget -c ftp://ftp.ensemblgenomes.org/pub/release-19/bacteria//fasta/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.19.dna.toplevel.fa.gz 
wget -c ftp://ftp.ensemblgenomes.org/pub/release-19/bacteria//gtf/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.19.gtf.gz
popd


pushd $IRAP_DIR/data/raw_data/ecoli_k12
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR933/SRR933983/SRR933983.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR933/SRR933984/SRR933984.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR933/SRR933985/SRR933985.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR933/SRR933989/SRR933989.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR933/SRR933990/SRR933990.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR933/SRR933991/SRR933991.fastq.gz
popd

cat <<EOF > ecoli_ex.conf

# minimal configuration file for DE
# experiment name
name=ecoli_ex
# species
species=ecoli_k12
# reference genome
reference=Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.19.dna.toplevel.fa.gz
# gtf file
gtf_file=Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.19.gtf.gz
# Enable filtering based on quality
qual_filtering=on
# Use a contamination data set to filter out reads
cont_index=no
# Toplevel directory with the data
# data_dir=

# some contrasts...
# GA=Group A
contrasts=GAvsGB GBvsGA
GAvsGB=GA GB
GBvsGA=GB GA
GA=FA FB FC
GB=FD FE FF

se=FA FB FC FD FE FF

FA=SRR933983.fastq.gz
FA_rs=50
FA_qual=33

FB=SRR933984.fastq.gz
FB_rs=50
FB_qual=33

FC=SRR933985.fastq.gz
FC_rs=50
FC_qual=33

FD=SRR933989.fastq.gz
FD_rs=50
FD_qual=33

FE=SRR933990.fastq.gz
FE_rs=50
FE_qual=33

FF=SRR933990.fastq.gz
FF_rs=50
FF_qual=33
EOF
