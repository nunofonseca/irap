#!/bin/env bash
set -e
folder=$1
shift
index_params=$*
# 1 - folder with the fasta files
files=$(ls -1 $folder/*.fa)

for f in $files; do
    echo "Indexing $f..."
    irap_map.sh mapsplice bowtie-build  $index_params $f $(dirname $f)/$(basename $f .fa)
done

exit 0
