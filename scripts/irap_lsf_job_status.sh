#!/bin/bash
# $Id: 0.1.1$
# TO BE executed by irap
jobname=$1
JOB_ID=$2
shift 2
JOB_CMD=$*
if [ "$JOB_CMD-" = "-" ]; then
    echo "ERROR! Usage: irap_lsf_job_status.sh jobname jobid job_cmd" >&2
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
     CUR_MEM=`bjobs -a -l -J  $jobname_prefix* | grep "MAX MEM" | cut -f 2 -d: | sed "s/Mbytes//"|head -n 1`
     MEM=`expr $CUR_MEM + $JOB_MEM_INCR`
     if [ $MEM -lt $JOB_MAX_MEM ]; then
       # resubmit
       echo "Job(s) failed due to lack of memory...resubmitting job(s) with more memory ($MEM)."
       echo $IRAP_JOB_CMD
       MEM=$MEM $IRAP_JOB_CMD
     fi
  fi
fi

exit $EXIT_CODE
