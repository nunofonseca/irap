#!/bin/bash
# wrapper to fastq_info

FILES=$*

function file_extension {
    echo $*|sed -E "s/([^ \.]+)\.//g" 
}
set -e
# check extension
ext=`file_extension $1`
if [ "$ext-" == "bam-" ]; then
    #hmm, this now validates bams...kind of
    # samtools version should be 1 or above
    f=$1    
    echo "BAM file"
    # check if the BAM contains unaligned reads
    echo "Checking for unmapped reads"
    UN=`samtools view -c -f 4 $f`
    if [ "$UN-" == "0-" ]; then
	echo "ERROR: No unaligned reads found in $f."  1>&2
	exit 1
    fi
    named_pipe=.`basename .$f`.pipe.fastq
    mkfifo $named_pipe
    echo "Converting BAM to fastq"
    samtools bam2fq $f > $named_pipe &
    
    FILES2PROCESS=$named_pipe
    FILES2DELETE=$named_pipe
else
    FILES2PROCESS=
    FILES2DELETE=
    for f in $FILES; do
	ext=`file_extension $f`
	if [ "-$ext" == "-bz2" ] || [ "-$ext" == "-bzip2" ] ; then
	    echo BZIP file
	    named_pipe=.`basename .$f`.pipe.fastq
	    mkfifo $named_pipe
	    bunzip2 -k  -c $f > $named_pipe  &
	    FILES2PROCESS="$FILES2PROCESS $named_pipe"
	    FILES2DELETE="$FILES2DELETE $named_pipe"
	else
	    FILES2PROCESS="$FILES2PROCESS $f"
	fi
    done
fi
fastq_info $FILES2PROCESS

if [ "-$FILES2DELETE" != "-" ]; then
    #echo -n "Removing named pipes..."
    rm -f $FILES2DELETE
    #echo "done."
fi
exit 0
