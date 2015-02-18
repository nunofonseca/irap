#!/bin/bash

if [ "$1-" == "-" ]; then
    DATA_DIR=$IRAP_DIR/data
else
    DATA_DIR=$1
fi

set -e

echo Installing data to $DATA_DIR
echo Edit the script to install to a different folder

mkdir -p $DATA_DIR/reference/homo_sapiens
mkdir -p $DATA_DIR/raw_data/homo_sapiens


echo "Downloading GTF file..."
wget -c ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz && echo "Downloading complete."
# TODO: check md5
# md5sum: fdc985b094ef823e5a51164cd109f4e8  gencode.v19.annotation.gtf.gz

echo "Downloading genome fasta file (chromosomes only)"
wget -c http://www.ebi.ac.uk/~nf/hs37d5.genome.chr_only.fa.gz && echo "Download complete"
#wget ftp://ftp.sanger.ac.uk/pub/project/PanCancer/genome.fa.gz && echo "Download complete"

# 
echo "Moving files to their location and uncompressing them"
mv hs37d5.genome.chr_only.fa.gz $DATA_DIR/reference/homo_sapiens
mv gencode.v19.annotation.gtf.gz $DATA_DIR/reference/homo_sapiens
gunzip $DATA_DIR/reference/homo_sapiens/hs37d5.genome.chr_only.fa.gz
# Fix the gtf
zcat $DATA_DIR/reference/homo_sapiens/gencode.v19.annotation.gtf.gz | tail -n +6 | sed -e "s/^chrM/MT/g;s/^chr//g" > $DATA_DIR/reference/homo_sapiens/gencode.v19.annotation.hs37d5_chr.gtf

# metadata.tsv
wget http://www.ebi.ac.uk/~nf/pcawg/metadata.tsv -O $IRAP_DIR/metadata.tsv


# template irap configuration file
cat <<EOF > $IRAP_DIR/pcawg.conf
# irap configuration template
# Mapping SOP
sop=pawg3_th2_mapping
# uncomment the following line to use Star for mapping
# or override in the command line
#sop=pawg3_star_mapping
data_dir=$IRAP_DIR/data

# Folder where all results will be placed
name=pcawg

######################
# Add the samples 
# se=e100831_UNC2-RDR300275_00022_FC_62ERGAAXX.6
# pe=e110302_UNC11-SN627_0067_BB04EPABXX.6

# e100831_UNC2-RDR300275_00022_FC_62ERGAAXX.6=100831_UNC2-RDR300275_00022_FC_62ERGAAXX.6.fastq
# e100831_UNC2-RDR300275_00022_FC_62ERGAAXX.6_rs=76
# e100831_UNC2-RDR300275_00022_FC_62ERGAAXX.6_qual=33


# e110302_UNC11-SN627_0067_BB04EPABXX.6=110302_UNC11-SN627_0067_BB04EPABXX.6_1.fastq 110302_UNC11-SN627_0067_BB04EPABXX.6_2.fastq
# e110302_UNC11-SN627_0067_BB04EPABXX.6_rs=50
# e110302_UNC11-SN627_0067_BB04EPABXX.6_qual=33
# e110302_UNC11-SN627_0067_BB04EPABXX.6_ins=350
# e110302_UNC11-SN627_0067_BB04EPABXX.6_sd=300


EOF

# Generate the index for TH2 and Star
pushd $IRAP_DIR
#irap conf=$IRAP_DIR/pcawg.conf sop=pawg3_th2_mapping stage0
#irap conf=$IRAP_DIR/pcawg.conf sop=pawg3_star_mapping stage0
