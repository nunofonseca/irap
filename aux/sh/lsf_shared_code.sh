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

# Shared code by all *_lsf scripts
#cmd="bsub_wrapper.sh "
cmd="irap "
RAND=`perl -e "print int(rand()*10);"`
#DATE=`date "+%w%H%M"`
DATE=`date "+%w%H%M%S"`
###################################
# length of jobname needs to be 
# small otherwise lsf dependencies 
# will not work
jobname_prefix="$RAND${DATE}e"
###################################
## Check computer farm requirements
#Number of threads
if [ "$THREADS-" = "-" ]; then 
    THREADS=8
fi
#Memory
if [ "$MEM-" = "-" ]; then 
    MEM=16000
fi
#Debug status
if [ "$DEBUG-" = "1-" ]; then
    cmd="irap"
fi
#Queue
if [ "$QUEUE-" = "-" ]; then
    #QUEUE="research-rh6"
    echo "ERROR in irap_lsf: LSF queue not defined. Please check the  $IRAP_DIR/irap_setup.sh file."
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

#########################
## Check input parameters 
#Configuration file
if [ "$1-" = "-" ]; then
	echo "Missing conf. file"
	echo "Usage: $LSF_CMD conf=my_irap.conf [IRAP options]"
	exit 1
fi
conf=`echo $1|cut -f 2 -d=`
if [ ! -e $conf ]; then
    echo "$conf file not found"
    exit 1
fi
conf_var=`echo $1|cut -f 1 -d=`
if [ "$conf_var-" != "conf-" ]; then
	echo "Usage: $LSF_CMD conf=my_irap.conf [IRAP options]"
	exit 1
fi

conf=`echo $1|cut -f 2 -d=`
if [ ! -e "$conf" ]; then
	echo "Conf. file $conf does not exist"
	exit 1
fi
#DATA_DIR=$2
data_dir=`grep "^data_dir=" $conf|cut -f 2 -d=`

shift 1
IRAP_PARAMS=$*
#Data directory

IRAP_PAR_CMD="$0 conf=$conf $IRAP_PARAMS"
###################
## Load config file
echo " * Trying to load configuration file $conf..."
#name
name=`grep "^name=" $conf|cut -f 2 -d=`
se="`grep "^se=" $conf|cut -f 2 -d=`"
pe="`grep "^pe=" $conf|cut -f 2 -d=`"
mapper="`grep "^pe=" $conf|cut -f 2 -d=`"
contrasts="`grep "^contrasts=" $conf|cut -f 2 -d=`"
# log_dir
log_dir=`grep "^log_dir=" $conf|cut -f 2 -d=`
report_dir=$name/report
echo " * Configuration loaded."

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

# override irap_params
for p in $IRAP_PARAMS; do
    pv=`echo $p|sed "s/.*=.*/=/g"`
    if [ "$pv-" = "=-" ]; then
	var=`echo $p|cut -f 1 -d\=`
	value=`echo $p|cut -f 2 -d\=`
	export $var=$value
	echo "$var=$value"
    fi
done

DATA_DIR=$data_dir
if [ "$DATA_DIR-" = "-" ]; then
	echo "Missing DATA_DIR"
	exit 1
fi
if [ ! -e $DATA_DIR ]; then
	echo "DATA_DIR $DATA_DIR does not exist"
	exit 1
fi

#####################
## Set default values
# TODO: get default values from IRAP
if [ "$mapper-" = "-" ]; then
    mapper=tophat1
fi
if [ "$name-" = "-" ]; then
 echo "!missing argument name!"
	exit 1
fi
#echo PE=$pe
#echo SE=$se 

# cur_stage
CUR_STAGE=
function get_path2logfile {
    mkdir -p $LOG_DIR/$CUR_STAGE
    echo $LOG_DIR/$CUR_STAGE
}

function p_info {
    echo "[INFO] $*" >&2
}

function get_param_value {
    param_name=$1
    conf_file=$2
    val=`grep "^$param_name=" $conf_file|cut -f 2 -d=`
    echo $val
}

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

function get_maxmem {
    let MAX_MEM=($1/4000+1)*4000
    echo $1
}

function get_cached_value {

    if [ ! -e $name/cached_vars.mk ]; then
	echo ERROR: File $name/cached_vars.mk not found > /dev/stderr
	exit 1
    fi
    val=`grep "^$1=" $name/cached_vars.mk|cut -f 2 -d=`
    echo $val
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
