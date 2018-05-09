#!/bin/bash

STDERR=/dev/stderr

if [ "$*-" == "-"  ]; then
    echo "fastq_info.sh .fastq [.fastq]" > $STDERR
    exit 1
fi

A="`fastq_info $* 2>/dev/stdout | tee $STDERR |tail -n 5`"  
ret=$?
if [ $ret != 0 ]; then
    echo "ERROR" > $STDERR
    exit $ret
fi
if [ `echo "$A"|grep -c -i "Error" ` != 0 ]; then
    exit 1
fi

f1=$1
let nfiles=`echo $* | wc -w`
if [ "$1" == "-f" ] || [ "$2" == "-f" ] || [ "$3" == "-f" ]; then
    f1=$2
    let nfiles=nfiles-1
fi

#
# note: some SE file names may end with _[12].fastq but it will not work with the following code
libname=`basename $1|sed -E "s/((_[12])?\.(fastq|fq)[^ ]*)//"`
dir=`dirname $1` 

nreads=`echo $A| sed -E "s/.*Number of reads: ([0-9]+).?([0-9]*).*/\1 \2/" ` 
qual_range=`echo $A| sed -E "s/.*Quality encoding range: ([0-9]+) ([0-9]+).*/\1 \2/" ` 
qual=`echo $A| sed -E "s/.*Quality encoding: ([0-9]+|solexa).*/\1/" ` 


if [ "$qual-"  == "solexa-" ]; then
    echo "ERROR: Quality encoding not supported by iRAP -  $qual" > /dev/stderr
    exit 1
fi

# use the median read size
rs=`echo $A| sed -E "s/.*Read length: ([0-9]+) ([0-9]+) ([0-9]+).*/\3/" ` 
rs_range=`echo $A| sed -E "s/.*Read length: ([0-9]+) ([0-9]+) ([0-9]+).*/\1 \2 \3/" ` 
# validate libname
# if libname end with _{1,2} then change it
# 2018-05-09: no longer necessary...
#libname=`echo $libname|sed -E "s/_([12])$/.\1/"`

echo "#nreads=$nreads" 
echo $libname=`basename $1` `basename $2 2>/dev/null` 
echo ${libname}_dir=$dir
echo "#rs=$rs_range" 
echo ${libname}_rs=$rs
echo "#range: $qual_range" 
echo ${libname}_qual=$qual
echo ${libname}_strand=both
# PE
if [ $nfiles -gt 1 ]; then
    echo ${libname}_ins=400
    echo ${libname}_sd=350
    echo "#pe=$libname"
else
    echo "#se=$libname"
fi


exit 0
