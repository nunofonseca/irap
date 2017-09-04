#!/bin/bash
# Script to merge a set of quantification matrices
FILES=$*

## Read the list of files from stdin
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

#set -e
f1=$1
shift 1

if [ ! -e $f1 ]; then
    echo "ERROR: file $f1 not found"
    exit 1
fi

# 0 - no, 1 - yes
GZIP_FILE=`file $f1|grep -c gzip`
if [ $GZIP_FILE == 1 ]; then
    CAT=zcat
    GREP=zgrep
else
    CAT=cat
    GREP=grep
fi
###################################
# Does the matrix have an header?
FEATURE=`$CAT $f1 | head -n 1 | cut -f 1|grep -i -E "(Gene|Exon|Transcript|Feature)"`
HEADER=yes
if [ "$FEATURE-" == "-" ]; then
    FEATURE="Gene"
    HEADER=no
fi
#echo $FEATURE >&2
if [ $HEADER == "no" ]; then
    filter_header=
else
    filter_header="tail -n +2 |"
fi
##################################
# Number of columns
NCOLS=`$CAT $f1 | head -n 1| tr "\t" "\n" | wc -l`
#echo "NCOLS=$NCOLS|"

if [ "$NCOLS-" != "2" ]; then
    KEEP_COL_NAME=yes
else
    KEEP_COL_NAME=no
    ## replace the column name by the file name
    ##################################
    if [ "-$USE_BASENAME" == "-" ]; then
	USE_BASENAME=yes    
    fi
    echo BASENAME=$USE_BASENAME >&2    
fi


lfiles_name=`mktemp`
lfiles_name_m=`mktemp`
    
if [ "$KEEP_COL_NAME" == "yes" ]; then

    $CAT $f1 |  cut -f 1 > $f1.tmp

    ## merged file
    $CAT  $f1 >  $lfiles_name_m

    for f in $*; do
	echo Comparing $f1 $f >&2    
	#cut -f 1 $f > b.tmp
	DIFF=`$GREP -v CUFF $f| $filter_header  cut -f 1 | diff -q $f1.tmp -`
	if [ "$DIFF-" = "-" ]; then
	    echo "$FEATURE order OK" >&2
	else
	    echo "ERROR: order of the rows is different - $f1 $f" >&2
	    exit 1
	fi
	rm -f b.tmp
	$GREP -v CUFF $f| $filter_header cut -f 2-  | paste $lfiles_name_m - > $lfiles_name_m.tmp
	mv $lfiles_name_m.tmp $lfiles_name_m
    done

else
    
    # check if order is ok
    # exclude entries added by CUFFlinks 
    $GREP -v CUFF $f1| $filter_header  cut -f 1 > $f1.tmp
    $GREP -v CUFF $f1| $filter_header  cut -f 2 > $f1.tmp2

    # merged file
    paste $f1.tmp $f1.tmp2 > $lfiles_name_m

    if [ "$USE_BASENAME" == "yes" ]; then
	echo -n "`basename $f1.tmp2`" > $lfiles_name
    else
	echo -n " $f1" > $lfiles_name
    fi

    for f in $*; do
	echo Comparing $f1 $f >&2    
	#cut -f 1 $f > b.tmp
	DIFF=`$GREP -v CUFF $f| $filter_header  cut -f 1 | diff -q $f1.tmp -`
	if [ "$DIFF-" = "-" ]; then
	    echo "$FEATURE order OK" >&2
	else
	    echo "ERROR: order of the rows is different - $f1 $f" >&2
	    exit 1
	fi
	rm -f b.tmp
	$GREP -v CUFF $f| $filter_header  cut -f 2  | paste $lfiles_name_m - > $lfiles_name_m.tmp
	mv $lfiles_name_m.tmp $lfiles_name_m
	if [ "$USE_BASENAME" == "yes" ]; then
	    echo -n " `basename $f.tmp2`" >> $lfiles_name
	else
	    echo -n " $f" >> $lfiles_name
	fi
    done
    echo  >> $lfiles_name
fi

# labels/header
sed  -e "1s/\.riu\.[^\. \t]*//g;1s/\.[^\.]*\.tsv//g;1s/\(.tmp2\|.raw\|.fpkm\|.rpkms\|.rpkm\|.genes\|.nlib\|.exons\|.transcripts\)//g;1s/\(\.pe\|\.se\)//g;s/ /\t/g;1s/^/$FEATURE\t/;1s/.gz$//g" $lfiles_name 

cat $lfiles_name_m
rm -f $FILES2 $f1.tmp* $lfiles_name $lfiles_name_m
exit 0
