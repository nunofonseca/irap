#!/bin/bash
# $Id: 0.1.3 Nuno Fonseca Fri Dec 21 14:59:29 2012$
# Merge TSV files using the first column 
FILES=$*
set -e
if [ "$1-" == "-stdin-" ]; then
    cat | irap_merge_tsv_NA.R -stdin > $1.tmp
else
    irap_merge_tsv_NA.R $FILES  > $1.tmp
fi
sed  -e "s/\.[^\.]*\.tsv//g;s/\(.raw\|.rpkms\|.genes\|.nlib\|.exons\|.transcripts\)//g;s/\(\.pe\|\.se\)//g" $1.tmp
rm -f $1.tmp
exit 0
