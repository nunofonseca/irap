#!/bin/bash

if [ "$*-" == "-"  ]; then
    echo "fastq_info.sh .fastq [.fastq]" 
    exit 1
fi

A="`fastq_info $*`" 
ret=$?
if [ $ret != 0 ]; then
    echo "ERROR" > /dev/stderr
    exit $ret
fi

#echo 
libname=l`basename $1|sed -E "s/((_[12])?\.fastq[^ ]*)//" ` 
dir=`dirname $1` 

nreads=`echo $A| sed -E "s/.*Number of reads: ([0-9]+)\s?([0-9]*).*/\1 \2/" ` 
qual_range=`echo $A| sed -E "s/.*Quality encoding range: ([0-9]+) ([0-9]+).*/\1 \2/" ` 
qual=`echo $A| sed -E "s/.*Quality encoding: ([0-9]+).*/\1/" ` 
rs=`echo $A| sed -E "s/.*Read length: ([0-9]+).*/\1/" ` 
rs_range=`echo $A| sed -E "s/.*Read length: ([0-9]+) ([0-9]+).*/\1 \2/" ` 


# validate libname
# if libname end with _{1,2} then change it
libname=`echo $libname|sed -E "s/_([12])$/.\1/"`

echo "#nreads=$nreads" 
echo $libname=`basename $1` `basename $2 2>/dev/null` 
echo ${libname}_dir=$dir
echo "#rs=$rs_range" 
echo ${libname}_rs=$rs
echo "#range: $qual_range" 
echo ${libname}_qual=$qual
echo ${libname}_strand=both
# PE
if [ "$2-" != "-" ]; then
    echo ${libname}_ins=400
    echo ${libname}_sd=350
    echo "#pe=$libname"
else
    echo "#se=$libname"
fi


exit 0
