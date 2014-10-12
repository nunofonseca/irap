#!/bin/bash

DATA_DIR=$IRAP_DIR/data
set -e

echo Installing data to $DATA_DIR
echo Edit the script to install to a different folder

mkdir -p $DATA_DIR/reference/homo_sapiens
mkdir -p $DATA_DIR/raw_data/homo_sapiens


echo "Downloading GTF file..."
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz && echo "Downloading complete."

echo "Downloading genome fasta file (chromosomes only)"
wget http://www.ebi.ac.uk/~nf/GRCh37.p13.genome.chr_only.fa.gz && echo "Download complete"

echo "Moving files to their location and uncompressing them"
mv GRCh37.p13.genome.chr_only.fa.gz $DATA_DIR/reference/homo_sapiens
mv gencode.v19.annotation.gtf.gz $DATA_DIR/reference/homo_sapiens
gunzip $DATA_DIR/reference/homo_sapiens/GRCh37.p13.genome.chr_only.fa.gz
gunzip $DATA_DIR/reference/homo_sapiens/gencode.v19.annotation.gtf.gz

# template irap configuration file
cat <<EOF 
# irap configuration template
# Mapping SOP
sop=pawg3_th2_mapping
# uncomment the following line to use Star for mapping
#sop=pawg3_star_mapping


# Folder where all results will be placed
name=pcawg_example

######################
# Add the samples 
se=e100831_UNC2-RDR300275_00022_FC_62ERGAAXX.6
pe=e110302_UNC11-SN627_0067_BB04EPABXX.6

e100831_UNC2-RDR300275_00022_FC_62ERGAAXX.6=100831_UNC2-RDR300275_00022_FC_62ERGAAXX.6.fastq
e100831_UNC2-RDR300275_00022_FC_62ERGAAXX.6_rs=76
e100831_UNC2-RDR300275_00022_FC_62ERGAAXX.6_qual=33


e110302_UNC11-SN627_0067_BB04EPABXX.6=110302_UNC11-SN627_0067_BB04EPABXX.6_1.fastq 110302_UNC11-SN627_0067_BB04EPABXX.6_2.fastq
e110302_UNC11-SN627_0067_BB04EPABXX.6_rs=50
e110302_UNC11-SN627_0067_BB04EPABXX.6_qual=33
e110302_UNC11-SN627_0067_BB04EPABXX.6_ins=350
e110302_UNC11-SN627_0067_BB04EPABXX.6_sd=300



EOF
