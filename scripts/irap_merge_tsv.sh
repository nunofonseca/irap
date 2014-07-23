#!/bin/bash
# $Id: 0.1.3 Nuno Fonseca Fri Dec 21 14:59:29 2012$
FILES=$*

NAMES=$FILES
if [ "$*-" == "-stdin-" ]; then
    set -e
    read -t 60 NAMES
    set -- $NAMES
    set +e
else
    if [ "$1-" == "-" ]; then
	echo "ERROR: no file names provided"
	exit 1
    fi
fi
#echo $*
f1=$1
shift 1
if [ ! -e $f1 ]; then
    echo "ERROR: file $f1 not found"
    exit 1
fi

# check if order is ok
# exclude entries added by CUFFlinks 
#cut -f 1 $f1 > $f1.tmp
#cut -f 2 $f1 > $f1.tmp2
grep -v CUFF $f1| cut -f 1 > $f1.tmp
grep -v CUFF $f1| cut -f 2 > $f1.tmp2
echo $NAMES >&2
for f in $*; do
    echo Comparing $f1 $f >&2
    grep -v CUFF $f|cut -f 1 > b.tmp
    #cut -f 1 $f > b.tmp
    DIFF=`diff -q $f1.tmp b.tmp`
    if [ "$DIFF-" = "-" ]; then
       echo "Gene order OK" >&2
    else
       echo "ERROR." >&2
       exit 1
    fi
    grep -v CUFF $f|cut -f 2  > $f.tmp2
    #cut -f 2  $f > $f.tmp2
done
N=""
for n in $NAMES; do
    N="$N `basename $n`"
done
N=`echo $N|sed -e "s/\.[^\.]*\.tsv//g;s/\(.raw\|.rpkms\|.genes\|.nlib\|.exons\|.transcripts\)//g;s/\(\.pe\|\.se\)//g"`
echo "Gene $N" |tr " " "\t"  
FILES2=`echo $NAMES|sed  "s/.tsv/.tsv.tmp2/g"`
#echo $FILES2
paste $f1.tmp $FILES2
rm $FILES2
exit 0
