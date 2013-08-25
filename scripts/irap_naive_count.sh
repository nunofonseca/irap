#!/bin/env bash

BAM_FILE=$1
GTF_FILE=$2

if [ "$GTF_FILE-" = "-" ]; then
    echo "Usage: irap_naive_count.sh BAM_FILE GTF_FILE"  >&2
    exit 1
fi

if [ "$BAM_FILE-" = "-" ]; then
    echo "Usage: irap_naive_count.sh BAM_FILE GTF_FILE"  >&2
    exit 1
fi

if [ ! -e "$GTF_FILE" ]; then
    echo "$GTF_FILE not found"  >&2
    exit 1
fi

if [ ! -e "$BAM_FILE" ]; then
    echo "$BAM_FILE not found"  >&2
    exit 1
fi

# gene count
# by default we get the count by transcript
bedtools coverage -abam $BAM_FILE -b $GTF_FILE -split | grep -w CDS filt.bed.cov | cut -f 9,11|sed "s/.*transcript_id \"\([^\"]*\)\".*;/\1 /" 
