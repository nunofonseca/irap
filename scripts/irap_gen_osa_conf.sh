#!/bin/sh

file1=$1
threads=$2
expression=$3
outfile=$4
file2=$5
ins=$6
sd=$7

if [ "$outfile-" = "-" ]; then
    echo "Usage: irap_gen_osa_conf.sh file1 threads expression outfile [file2 ins sd]" >&2
    exit 1
fi

PAIRED=False
if [ "$sd-" == "-" ]; then
    ins=30
    sd=20
else
    PAIRED=True
fi
outdir=`dirname $outfile`
outfile=`basename $outfile`

cat <<EOF
<Files>
$file1
$file2

<Options>
ThreadNumber=$threads
PairedEnd=$PAIRED
FileFormat=FASTQ
Gzip=False
AutoPenalty=True 
ExpressionMeasurement=$expression
SearchNovelExonJunction=True
InsertSizeStandardDeviation=$sd //Default value=40
ExpectedInsertSize=$ins //Default value=300
InsertOnSameStrand=False // Possible values: True, False. Default value=False

<Output>
OutputName=$outfile
OutputPath=$outdir 

EOF
exit 0
