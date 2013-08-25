#!/bin/bash
if [ "-$1" = "-" ]; then
    echo "ERROR! usage: num_reads.sh fastq_filename" >&2
    exit 1
fi

if [ -e $1 ]; then 
    grep -c "^+$" $1
else 
if [ -e $1.gz ]; then 
    zgrep -c "^+$" $1.gz
else
    echo "ERROR! File $1(.gz) not found." >&2
    exit 1
fi
fi
exit 0
