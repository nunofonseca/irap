#!/bin/bash
if [ "-$1" = "-" ]; then
    echo "ERROR! usage: num_reads.sh fastq_filename" >&2
    exit 1
fi

if [ -e $1 ]; then 
    fastq_num_reads $1
else 
    echo "ERROR! $1 file not found" >&2
    exit 1
if [ -e $1.gz ]; then 
    echo "ERROR! .gz files are not supported" >&2
    exit 1
#    zgrep -c "^+$" $1.gz
else
    echo "ERROR! File $1(.gz) not found." >&2
    exit 1
fi
fi
exit 0
