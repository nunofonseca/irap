# =========================================================
# Copyright 2012-2014,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
# This file is part of iRAP.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
#
#
#    $Id: 0.1.3 Nuno Fonseca Wed Dec 26 16:20:16 2012$
# =========================================================


#source $IRAP_DIR/aux/sh/irap_fun.sh

if [ "$QUEUE-" = "-" ]; then
    #QUEUE="research-rh6"
    echo "ERROR in $LSF_CMD: LSF queue not defined. Please check the  $IRAP_DIR/irap_setup.sh file."
    exit 1
fi

if [ "$IRAP_LSF_PARAMS-" = "-" ]; then 
    IRAP_LSF_PARAMS=
fi

export QUEUE
export IRAP_LSF_GROUP
export IRAP_LSF_PARAMS
export MEM
export THREADS

# directory to keep the logs
if [ "$log_dir-" == "-" ]; then
    log_dir=`pwd`/$name/lsf_logs
fi
LOG_DIR=$log_dir/$jobname_prefix
mkdir -p $LOG_DIR

# Check JOB_MAX_MEM
if [ "$JOB_MAX_MEM-" != "-" ]; then 
  if [ "$JOB_MEM_INCR-" == "-" ]; then
    echo ERROR: JOB_MEM_INCR should be defined when JOB_MAX_MEM is defined! 
    exit 1
  fi
fi

function check_dependency {
    jobname=$1
    sleep 1
    FOR_RUNNING=`bjobs -a -J $jobname|grep JOBID| wc -l`
    if [ $FOR_RUNNING == 0 ]; then
	echo ''
    else
	echo "-w ended(\"$jobname\")"
    fi
    
}

function stop_job {
    if [ "$DEBUG-" == "1-" ]; then
	echo "stop/suspend job $1"
    else
	bstop -J $1
    fi
}

function resume_job {
    if [ "$DEBUG-" == "1-" ]; then
	echo "resume job $1"
    else
	bresume -J $1
    fi
}

function submit_job_status {
    jobname=$1
    WAITFOR=
    if [ "$DEBUG-" = "1-" ]; then
        ECHO=echo
    else
	JOB_ID=`bjobs -a -J $jobname| tail -n 1|sed -E "s/([0-9]*).*/\1/"`
	WAITFOR=`check_dependency $jobname`
	ECHO= 
    fi
    #p_info "WAITFOR (id)=$JOB_ID $jobname   $WAITFOR"
    # in spite of the checks, the job may have finished before lunching the new one, hence catch the error and submit a new one if an error occurs
    $ECHO bsub $IRAP_LSF_PARAMS -M 1000 -R "select[mem>=1000]  rusage[mem=1000]" -q $QUEUE  -J "${jobname}n" $WAITFOR  irap_lsf_job_status.sh $jobname $JOB_ID  `get_maxmem $MEM` $LOG_DIR  $IRAP_PAR_CMD
    #2> /dev/null
    if [ $? != 0 ]; then
	p_info "$jobname not  found...probably it has already finished"
	WAITFOR=
	$ECHO bsub $IRAP_LSF_PARAMS -M 1000 -R "select[mem>=1000]  rusage[mem=1000]" -q $QUEUE  -J "${jobname}n" $WAITFOR  irap_lsf_job_status.sh $jobname $JOB_ID `get_maxmem $MEM` $LOG_DIR $IRAP_PAR_CMD 
    fi
}


################
## Job functions (computer farm)
# length of jobname needs to be small otherwise lsf dependencies will not work
function submit_job {
    jobname=$1
    shift
    cmd2e=$*
    #echo "$jobname: $* max_threads=$THREADS  data_dir=$DATA_DIR/data" 
    if [ "$DEBUG-" = "1-" ]; then
        ECHO=echo
    else
	ECHO=
    fi
    #########################################################
    # limit the number of parallel jobs by using lsf's groups
    
    # default group: irap
    # TODO: define groups by stage: 
    #    irap_qc, ...
    GROUP=
    if [ "$IRAP_LSF_GROUP-" != "-" ]; then
	GROUP="-g $IRAP_LSF_GROUP"
    fi
    #########################################################
    #-R  "span[ptile=$THREADS]"
    MAX_MEM=`get_maxmem $MEM`
    if [ "$WAIT_FOR_IDS-" != "-" ]; then
	$ECHO bsub $IRAP_LSF_PARAMS -q $QUEUE -n $THREADS  -M $MAX_MEM -R "select[mem>=$MEM] rusage[mem=$MEM]"  -w "$WAIT_FOR_IDS"  -cwd `pwd` -o "`get_path2logfile`/$jobname-%J.out" -e "`get_path2logfile`/$jobname-%J.err" -J $jobname  $cmd2e max_threads=$THREADS  data_dir=$DATA_DIR max_mem=$MEM
    else
	$ECHO bsub $IRAP_LSF_PARAMS -q $QUEUE  $GROUP -n $THREADS  -M $MAX_MEM -R "select[mem>=$MEM]  rusage[mem=$MEM]"   -cwd `pwd` -o "`get_path2logfile`/$jobname-%J.out" -e "`get_path2logfile`/$jobname-%J.err" -J $jobname  $cmd2e max_threads=$THREADS  data_dir=$DATA_DIR max_mem=$MEM	
    fi
}
