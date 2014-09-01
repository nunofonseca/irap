#!/bin/bash
# $Id: 0.1.1$
# TO BE executed by irap
jobname=$1
JOB_ID=$2
MEM=$3
LOG_DIR=$4
shift 4
JOB_CMD=$*
if [ "$JOB_CMD-" = "-" ]; then
    echo "ERROR! Usage: irap_lsf_job_status.sh jobname jobid MEM log_dir job_cmd" >&2
    exit 1
fi
echo $JOB_CMD
# sleep for a while otherwise you get the dependent job as being in the RUN state
sleep 60
STATUS=`bjobs  -a -J $jobname| cut -c 17-22|tail -n 1| sed -E "s/([a-zA-Z]*).*/\1/"`

if [ "$STATUS-" == "-" ]; then
    echo "ERROR: Job $job information not found" > /dev/stderr
    exit 1
fi

if [ "$STATUS-" == "DONE-" ]; then
 EXIT_CODE=0
 echo "Job finished successfully ($STATUS)"
else
 EXIT_CODE=1
 echo "Job failed ($STATUS)"
 irap_lsf_filter_logs.sh $LOG_DIR
 EXIT_CODE=$?
fi
# print some info about the job
bjobs  -a -l -J $jobname
cat<<EOF


----------------------------------------------------------------------------------
EOF

if [  $EXIT_CODE == 1  ] &&  [ "$JOB_MAX_MEM-" != "-"  ]; then 
    if [ "$JOB_MEM_INCR-" == "-" ]; then
	echo JOB_MEM_INCR undefined
	exit 1
  fi
  # get the jobs name prefix to capture potential failures
  jobname_prefix=`echo $jobname|sed "s/..$//"`
  MEM_ERROR=`bjobs -a -l -J $jobname_prefix* | grep -c TERM_MEMLIMIT`
  echo MEM_ERROR=$MEM_ERROR
  if [ $MEM_ERROR != 0 ]; then
     CUR_MEM=$MEM
     NEW_MAX_MEM=`expr $CUR_MEM + $JOB_MEM_INCR`
     #echo "$NEW_MAX_MEM=$CUR_MEM + $JOB_MEM_INCR"
     if [ $NEW_MAX_MEM -lt $JOB_MAX_MEM ]; then
       # resubmit
       echo "Job(s) failed due to lack of memory...resubmitting job(s) with more memory ($NEW_MAX_MEM)."
       echo $JOB_CMD
       MEM=$NEW_MAX_MEM $JOB_CMD
     else
       echo "Job(s) failed due to lack of memory...not resubmitting because max. memory limit was reached."
     fi
  fi
fi

exit $EXIT_CODE
