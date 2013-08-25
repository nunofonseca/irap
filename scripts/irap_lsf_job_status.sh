#!/bin/bash
# $Id: 0.1.1$
# TO BE executed by irap
jobname=$1
JOB_ID=$2
out_dir=$3
shift 3
extra=$*
if [ "$out_dir-" = "-" ]; then
    echo "ERROR! Usage: irap_lsf_job_status.sh jobname jobid out_dir [ignored]" >&2
    exit 1
fi
echo $extra
# sleep for a while otherwise you get the dependent job as being in the RUN state
sleep 60
STATUS=`bjobs  -a -J $jobname| cut -c 17-22|tail -n 1| sed -E "s/([a-zA-Z]*).*/\1/"`

if [ "$STATUS-" = "DONE-" ]; then
EXIT_CODE=0
echo "Job finished successfully ($STATUS)"
else
EXIT_CODE=1
echo "Job failed ($STATUS)"
fi
echo $extra

# print some info about the job
bjobs  -a -l -J $jobname
echo OUTFILE=$out_dir/$jobname-$JOB_ID.out
cat<<EOF


----------------------------------------------------------------------------------
EOF

exit $EXIT_CODE
