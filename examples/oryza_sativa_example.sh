#!/bin/bash
#################################
# how many cores/threads to use
CORES=2
#################################

if [ "$IRAP_DIR-" == "-" ]; then
    echo "ERROR: irap is not properly configured"
    exit 1
fi
echo "IRAP=$IRAP_DIR"


SPECIES=oryza

set -e

echo "* Creating directory structure..."
mkdir -p data/reference/oryza
mkdir -p data/raw_data/oryza
echo "* Creating directory structure...done."

echo "* Downloading the reference genome and annotation..."
pushd data/reference/oryza
# Download the genome
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-31/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.31.dna.toplevel.fa.gz
# Download the annotation 
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-31/gtf/oryza_sativa/Oryza_sativa.IRGSP-1.0.31.gtf.gz
# Remove comments from the GTF
zcat Oryza_sativa.IRGSP-1.0.31.gtf.gz | grep -v "^#" > tmp && mv tmp Oryza_sativa.IRGSP-1.0.31.gtf
# some gene names have ;
sed -i -E  "s/([a-zA-Z0-9])\;/\1/g" Oryza_sativa.IRGSP-1.0.31.gtf
popd
echo "* Downloading the reference genome and annotation...done."

# download a library 
echo "* Downloading a library to test..."
pushd data/raw_data/$SPECIES
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA209/ERA209245/fastq/nippon_control_1hr_rep1.fastq.gz
# pick a subset of the reads (this is just a test anyway)
zcat nippon_control_1hr_rep1.fastq.gz | head -n 40000 > test.fastq
popd
echo "* Downloading a library to test...complete."

# generate the contamination DB
echo "* Creating the contamination DB..."
mkdir -p data/contamination
pushd data/contamination
# Create a contamination DB - fungi genomes
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/fungi.16.1.genomic.fna.gz
gunzip -f fungi.16.1.genomic.fna.gz
mv fungi.16.1.genomic.fna fungi.16.1.genomic.fa
irap_add2contaminationDB -n fungi.16.1 -f fungi.16.1.genomic.fa -d .
rm -f fungi.16.1.genomic.fa
popd
echo "* Creating the contamination DB...done."

# create a control/configuration file for iRAP
echo "* Creating a configuration file for iRAP..."
cat <<EOF > oryza_sativa_ex.conf
# minimal configuration file
atlas_run=y
data_dir=data
species=$SPECIES
reference=Oryza_sativa.IRGSP-1.0.31.dna.toplevel.fa.gz
gtf_file=Oryza_sativa.IRGSP-1.0.31.gtf
cont_index=data/contamination/fungi.16.1
EOF
echo "* Creating a configuration file for iRAP...done."


# Now everything is ready...
# Index the genome and perform some other initial preprocessing
# This needs to be executed only once for each genome/annotation
echo "* Initializing..."
irap_single_lib.sh -c oryza_sativa_ex.conf -1 test.fastq -0 -A -t $CORES -o out

#
irap_single_lib.sh -c oryza_sativa_ex.conf -1 data/raw_data/oryza/test.fastq  -A -t $CORES -o out
# 
exit 0

