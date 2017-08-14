#!/bin/bash
if [ "-$1" = "-" ]; then
    echo "ERROR! usage: num_reads.sh fastq_filename" >&2
    exit 1
fi

if [ -e $1 ]; then 
    fastq_num_reads $1
    exit 0
else 
    echo "ERROR! $1 file not found" >&2
    echo NA
    exit 0
fi

