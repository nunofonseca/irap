#!/bin/bash
CONF=$1
MAX_FASTQ=$2
data_dir=$3

species=test_species

if [ "$1-" == "-" ] ; then
    echo "Usage: irap_distrib_fastq.sh irap_conf_file max_files_per_folder [data_dir]" > /dev/stderr
    exit 1
fi

if [ "$CONF-" == "-" ] || [ ! -e $CONF ]; then
    echo "Unable to open configuration file" > /dev/stderr 
    exit 1
fi
if [ "$data_dir-" == "-" ]; then
    data_dir=`grep "^data_dir=" $CONF|cut -f 2 -d=`
fi

species=`grep "^species=" $CONF|cut -f 2 -d=`
echo data_dir=$data_dir
echo species=$species

# libs used
libs_se=`grep "^se=" $CONF|cut -f 2 -d=` 
libs_pe=`grep "^pe=" $CONF|cut -f 2 -d=`

libs="$libs_se $libs_pe"
nlibs=`echo $libs|wc -w`

echo $libs
echo $nlibs

CUR_NFILES=0
CUR_DIR=0
NEW_CONF=`echo $CONF|sed "s/.conf//"`
echo $NEW_CONF
grep -v "_dir=" $CONF > $NEW_CONF.new.conf
echo data_dir=$data_dir >> $NEW_CONF.new.conf

for l in $libs; do

   if [ $CUR_NFILES -gt $MAX_FASTQ ]; then
       CUR_NFILES=0
   fi
   if [ $CUR_NFILES -eq 0 ]; then
       # next folder
       CUR_DIR=`expr $CUR_DIR + 1`
       CUR_DIR_NAME=d$CUR_DIR
       PATH2DIR=$data_dir/raw_data/$species/$CUR_DIR_NAME
       echo "New folder: $PATH2DIR"
       mkdir -p $PATH2DIR
       CUR_NFILES=1
   fi
   echo "lib=$l"
   fastq_files=`grep "^$l=" $CONF|cut -f 2 -d=`
   echo $fastq_files
   
   echo "${l}_dir=$CUR_DIR_NAME" >> $NEW_CONF.new.conf
   echo "${l}_dir=$CUR_DIR_NAME"
   for f in $fastq_files; do
       if [ -e $data_dir/raw_data/$species/$f ]; then
	   mv $data_dir/raw_data/$species/$f $PATH2DIR
       else
	   if [ -e $PATH2DIR/$f ]; then
	       echo "file $f already moved"
	   else
	       echo "Missing file $f"
	       #exit 2
	   fi
       fi
       CUR_NFILES=`expr $CUR_NFILES + 1`
   done
done

exit 0


