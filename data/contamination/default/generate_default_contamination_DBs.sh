#!/bin/bash
# Create the default set of cont. databases
set -e

## Use refseq
## Fungi
if [ ! -e fungi.fa.gz ]; then
    wget -c ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/fungi.53.rna.fna.gz -O fungi.fa.gz
fi

## viruses
if [ ! -e viral.fa.gz ]; then
    wget -c ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz -O viral.fa.gz
fi

## phi-X174
if [ ! -e phiX174.fa.gz ]; then
    wget -c 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=9626372&maxplex=1' -O phiX174.fa
    gzip phiX174.fa
fi

## bacterial
##if [ ! -e bacterial.fa.gz ]; then
##    wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.1075.rna.fna.gz -O bacterial.fa.gz
##fi

## ecoli
if [ ! -e ecoli.fa.gz ]; then
    wget -c ftp://ftp.ensemblgenomes.org/pub/release-37/bacteria/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa.gz -O ecoli.fa.gz
fi


function gen_multi_species_DB {
    name=$1
    shift 1
    fasta_files=$*
    rm -f combined.fa
    for f in $fasta_files; do
	no_suffix=`echo $f|sed "s/.fa.gz$//"`;
	## add the group (fasta file prefix) to the seq name
	zcat $f | sed "s/^>/>$no_suffix:/" >> combined.fa
    done
    irap_add2contaminationDB -n $name -f combined.fa -d .. -m bowtie2
    rm -f combined.fa
    echo "IRAP: cont_index=$name"
}
## DB1 - bact/ecoli
gen_multi_species_DB  e_coli ecoli.fa.gz phiX174.fa.gz

## DB2 - bact_fungi_viral
gen_multi_species_DB ecoli_fungi_viral ecoli.fa.gz viral.fa.gz fungi.fa.gz phiX174.fa.gz

## DB3 - bact_viral
gen_multi_species_DB ecoli_viral ecoli.fa.gz viral.fa.gz phiX174.fa.gz

rm -f *.fa
echo "All done!"
exit 0
