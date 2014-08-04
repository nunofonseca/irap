#!/bin/bash
# $Id: 0.1.3 Nuno Fonseca Fri Dec 21 14:59:29 2012$
# Merge TSV files using the first column 
FILES=$*
set -e
if [ "$1-" == "-" ]; then
    echo "Missing arguments" > /dev/stderr
    exit 2
fi
OFILE=$1
if [ "$1-" == "-stdin-" ]; then
    OFILE=`mktemp -p .`
    cat - | irap_merge_tsv_NA.R -stdin > $OFILE.tmp
else
    irap_merge_tsv_NA.R $FILES  > $OFILE.tmp
fi
sed  -e "s/\.[^\.]*\.tsv//g;s/\(.raw\|.rpkms\|.genes\|.nlib\|.exons\|.transcripts\)//g;s/\(\.pe\|\.se\)//g" $OFILE.tmp 
rm -f $OFILE.tmp
exit 0
