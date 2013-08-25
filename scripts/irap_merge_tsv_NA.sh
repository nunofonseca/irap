#!/bin/bash
# $Id: 0.1.3 Nuno Fonseca Fri Dec 21 14:59:29 2012$
# Merge TSV files by the first column 
FILES=$*
set -e
irap_merge_tsv_NA.R $FILES  > $1.tmp
sed  -e "s/\.[^\.]*\.tsv//g;s/\(.raw\|.rpkms\|.genes\|.nlib\|.exons\|.transcripts\)//g;s/\(\.pe\|\.se\)//g" $1.tmp
rm -f $1.tmp
exit 0
