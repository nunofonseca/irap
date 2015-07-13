#!/bin/bash

FILES=$*

NAMES=$FILES
if [ "$*-" == "-stdin-" ]; then
    set -e
    read -e -t  60 NAMES
    set -- $NAMES
    set +e
else
    if [ "$1-" == "-" ]; then
	echo "ERROR: no file names provided" 
	exit 1
    fi
fi
#echo $* >&2
#set -e
f1=$1
shift 1

if [ ! -e $f1 ]; then
    echo "ERROR: file $f1 not found"
    exit 1
fi

#echo $f1 >&2
FEATURE=`head -n 1 $f1| cut -f 1|grep -i -E "(Gene|Exon|Transcript)"`
HEADER=yes
if [ "$FEATURE-" == "-" ]; then
    FEATURE="Gene"
    HEADER=no
fi
#echo $FEATURE >&2
if [ $HEADER == "no" ]; then
    filter_header=cat
else
    filter_header="tail -n +2"
fi

#
if [ "-$USE_BASENAME" == "-" ]; then
    USE_BASENAME=yes
fi

# check if order is ok
# exclude entries added by CUFFlinks 
grep -v CUFF $f1| $filter_header | cut -f 1 > $f1.tmp
grep -v CUFF $f1| $filter_header | cut -f 2 > $f1.tmp2

#echo $NAMES >&2
lfiles_name=`mktemp`
lfiles_name_m=`mktemp`
# merged file
paste $f1.tmp $f1.tmp2 > $lfiles_name_m
#echo $lfiles_name >&2
echo -n "`basename $f1.tmp2`" > $lfiles_name
for f in $*; do
    echo Comparing $f1 $f >&2    
    #cut -f 1 $f > b.tmp
    DIFF=`grep -v CUFF $f| $filter_header | cut -f 1 | diff -q $f1.tmp -`
    if [ "$DIFF-" = "-" ]; then
       echo "$FEATURE order OK" >&2
    else
       echo "ERROR." >&2
       exit 1
    fi
    rm -f b.tmp
    grep -v CUFF $f| $filter_header | cut -f 2  | paste $lfiles_name_m - > $lfiles_name_m.tmp
    mv $lfiles_name_m.tmp $lfiles_name_m
    if [ $USE_BASENAME==yes ]; then
	echo -n " `basename $f.tmp2`" >> $lfiles_name
    else
	echo -n " $f" >> $lfiles_name
    fi
done
echo  >> $lfiles_name
# labels/header
sed  -e "s/\.[^\.]*\.tsv//g;s/\(.tmp2\|.raw\|.rpkms\|.rpkm\|.genes\|.nlib\|.exons\|.transcripts\)//g;s/\(\.pe\|\.se\)//g;s/ /\t/g;s/^/$FEATURE\t/;" $lfiles_name 


cat $lfiles_name_m
rm -f $FILES2 $f1.tmp* $lfiles_name $lfiles_name_m
exit 0
