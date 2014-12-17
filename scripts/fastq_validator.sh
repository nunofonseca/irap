#!/bin/bash
# wrapper to fastq_info

FILES=$*

function file_extension {
    echo $*|sed -E "s/([^ \.]+)\.//g" 
}
set -e
# check extension
FILES2PROCESS=
FILES2DELETE=
for f in $FILES; do
    ext=`file_extension $f`
    if [ "-$ext" == "-bz2" ] || [ "-$ext" == "-bzip2" ] ; then
	echo BZIP file
	named_pipe=.$f.pipe.fastq
	mkfifo $named_pipe
	bunzip2 -k  -c $f > $named_pipe  &
	FILES2PROCESS="$FILES2PROCESS $named_pipe"
	FILES2DELETE="$FILES2DELETE $named_pipe"
    else
	FILES2PROCESS="$FILES2PROCESS $f"
    fi
done

fastq_info $FILES2PROCESS

if [ "-$FILES2DELETE" != "-" ]; then
    #echo -n "Removing named pipes..."
    rm -f $FILES2DELETE
    #echo "done."
fi
exit 0
