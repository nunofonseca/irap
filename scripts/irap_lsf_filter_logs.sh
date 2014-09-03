#!/bin/bash

if [ "$1-" == "-" ]; then
    echo "Usage: irap_lsf_filter_logs log_dir" > /dev/stderr
    exit 1
fi

INSUF_MEM=0
# Provide a bit more information
function check_jobs_status {
  dir=$1
  let success_jobs=0
  let failed_jobs=0

  function check_jobs_in_dir {    
    local f
    local out_files
    out_files=`ls -1 $1/*.out 2>/dev/null`
    for f in $out_files; do
       OK=`grep -c -E  "^Successfully completed." $f`
       #echo "OK=$OK" > /dev/stderr
       if [ "$OK" == "1" ]; then
         let success_jobs=$success_jobs+1
       else 
         # error
         LAST_LINES=`tail -n 50 $f|head -n 10`
         let failed_jobs=$failed_jobs+1
         EXIT_STATUS=`grep "Exited with exit code" $f|tail -n 1 |cut -f 5 -d\ |sed "s/\.$//"`
         CMD=`grep "LSBATCH: User input" -A 1 $f|tail -n 1`
	 if [ "$EXIT_STATUS" == "1" ]; then 
	     INSUF_MEM=1
	 fi
	 echo "------------------------------------------" 
	 echo "Out file: $f" 
         echo "Exit status: $EXIT_STATUS" 
	 echo "$CMD" 
	 echo "$LAST_LINES" 
    fi
    done
    # checking subfolders
    local DIRS=`ls -d $1/*/ 2> /dev/null`
    local d
    for d in $DIRS; do
      check_jobs_in_dir $d
    done 
   }
   check_jobs_in_dir $dir
   
   echo JOBS OK=$success_jobs
   if [ $failed_jobs \> 0 ]; then
       echo JOBS FAILED=$failed_jobs
       if [ $INSUF_MEM == 1 ]; then 
	   EXIT_STATUS=1
       else
	   EXIT_STATUS=2
       fi
   else
       EXIT_STATUS=0
   fi
   
}
check_jobs_status $1
exit $EXIT_STATUS
