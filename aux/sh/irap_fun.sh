# =========================================================
# Copyright 2012-2018,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
# =========================================================
# Shared code by all *_lsf scripts

irap_path=`which irap`
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
data_dir=`grep "^data_dir=" $conf|cut -f 2 -d=|tail -n 1`

shift 1
IRAP_PARAMS=$*
#Data directory

IRAP_PAR_CMD="$0 conf=$conf $IRAP_PARAMS"
###################
## Load config file
echo " * Trying to load configuration file $conf..."
#name
name=`grep "^name=" $conf|cut -f 2 -d= |tail -n 1`
se="`grep "^se=" $conf|cut -f 2 -d= |tail -n 1`"
pe="`grep "^pe=" $conf|cut -f 2 -d= |tail -n 1`"
mapper="`grep "^mapper=" $conf|cut -f 2 -d= |tail -n 1`"
contrasts="`grep "^contrasts=" $conf|cut -f 2 -d= | tail -n 1`"
# log_dir
log_dir=`grep "^log_dir=" $conf|cut -f 2 -d= | tail -n 1`
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

# cur_targets
CUR_TARGETS=
# cur_stage
CUR_STAGE=

declare -i s1=0
declare -i s2=0

#######################################
# Functions
function p_info {
    echo "[INFO] $*" >&2
}

function get_maxmem {
    let MAX_MEM=($1/4000+1)*4000
    echo $1
}

# cur_stage
CUR_STAGE=
function get_path2logfile {
    mkdir -p $LOG_DIR/$CUR_STAGE
    echo $LOG_DIR/$CUR_STAGE
}

function get_param_value {
    param_name=$1
    conf_file=$2
    val=`grep "^$param_name=" $conf_file|cut -f 2 -d=`
    echo $val
}
function get_cached_value {

    if [ ! -e $name/cached_vars.mk ]; then
	echo ERROR: File $name/cached_vars.mk not found > /dev/stderr
	exit 1
    fi
    val=`grep "^$1=" $name/cached_vars.mk|cut -f 2 -d=`
    echo $val
}

function get_targets_4_level {
    level=$1
    if [ "$level-" == "-" ]; then
	echo "ERROR: Internal error in get_targets_4_level $level"
	exit 1
    fi
    WAVE_FILE=.$$.wave
    targets=`$cmd conf=$conf $IRAP_PARAMS print_wave_${level}_targets wave_file=$WAVE_FILE`
    #echo $targets
    cat $WAVE_FILE
}

_INIT_SUSP_JOB=
function irap_init_job {

    # #############
    ## Submit jobs
    p_info " * Initialization..."
    # save variables
    $cmd conf=$conf $IRAP_PARAMS save_cache
    ##load the file with the variables already computed
    export IRAP_PARAMS="$IRAP_PARAMS use_cached_vars=y"
    p_info " * Starting initial job"
    CUR_STAGE=stage0
    submit_job "${jobname_prefix}i" $cmd conf=$conf bootstrap $IRAP_PARAMS -j $THREADS
    stop_job  "${jobname_prefix}i"
    _INIT_SUSP_JOB=${jobname_prefix}i
    CUR_TARGETS=bootstrap
    p_info " * First job suspended until all jobs are submitted."
    echo "${jobname_prefix}i"
}


function submit_jobs4stage {
    waitfor=$1
    level=$2
    targets=`get_targets_4_level $level`
    if [ "$targets-" == "-" ]; then
	DEPS=$1
    else
	p_info "* Checking current status..."
	$cmd conf=$conf $IRAP_PARAMS run_wave_$level -n -q 2> /dev/null
	let ret=$?
	if [ $ret -eq 0 ]; then
	    p_info "Skipping submission of jobs in level $level - all done"
	    DEPS=$1	    
	else
	    p_info "* Stage ${level}..."
	    DEPS="${jobname_prefix}${level}.*"
	    local nx=($targets)
	    local NJOBS=${#nx[@]}
	    p_info "Submiting $NJOBS jobs"
	    let JOBSPERBATCH=800
	    let BATCH=1
	    BATCHGROUP=
	    if [ $NJOBS -gt $JOBSPERBATCH ]; then
		BATCHGROUP=b${BATCH}.
		p_info "Too many jobs...grouping them in batches of $JOBSPERBATCH"
	    else
		BATCHGROUP=.
	    fi
	    
	    let i=1	   
	    for f in $targets; do    
		submit_job "${jobname_prefix}${level}${BATCHGROUP}$i"  -w "ended($waitfor)" $cmd conf=$conf  $IRAP_PARAMS $f
		let i=$i+1
		if [ $i -gt $JOBSPERBATCH ]; then
		    p_info "Batch $BATCH submitted"
		    THREADS=1 MEM=4000 submit_job "${jobname_prefix}${level}.$BATCH"  -w "ended(${jobname_prefix}${level}${BATCHGROUP}*)" echo nop
		    let i=1
		    let BATCH=$BATCH+1
		    BATCHGROUP=b${BATCH}.
		fi
	    done
	    if [ $i -gt 1 ] && [ "$BATCHGROUP-" != ".-" ]; then
		p_info "Batch $BATCH submitted"
		THREADS=1 MEM=4000 submit_job "${jobname_prefix}${level}.$BATCH"  -w "ended(${jobname_prefix}${level}${BATCHGROUP}*)" echo nop
	    fi

	    if [ "$BATCHGROUP-" != ".-" ]; then
		p_info "* Stage ${level}...$i jobs/$BATCH"
	    else
		p_info "* Stage ${level}...$i jobs"
	    fi
	    CUR_TARGETS=run_wave_$level
	fi
    fi
    echo $DEPS
}


function final_job {
    CUR_STAGE=
    waitfor=$1
    susp_job=$2
    level=$3
    #########################################
    # To finalize, run the whole pipeline
    # If everything went ok then nothing should be done
    # otherwise it should fail and an email will be sent
    all_levels=run_wave_b
    let cur_level=0
    while [ $cur_level -le $level ]; do
	all_levels="$all_levels run_wave_${cur_level}"
	let cur_level=cur_level+1
    done
    CUR_STAGE=  submit_job "${jobname_prefix}f" `check_dependency $waitfor`   "$cmd conf=$conf $all_levels $IRAP_PARAMS -n -q"
    
    # send an email to the user
    CUR_STAGE=  submit_job_status "${jobname_prefix}f"

    resume_job $susp_job
    # TODO: report
    echo JOBS=$LOG_DIR/${jobname_prefix}*.out
    echo JOBNAME=${jobname_prefix}f
}


##########################################################
# Report
declare -i counter=0



function qc_report {
    waitfor=$1
    echo " * QC Report "
    let counter=counter+1
    submit_job "${jobname_prefix}f[$counter]" "-w ended(\"$waitfor\")"  "$cmd conf=$conf qc_report $IRAP_PARAMS browsing=n"
}

function init_page_report {
    waitfor=$1    
    echo " * Information page "
    let counter=counter+1
    submit_job "${jobname_prefix}f[$counter]" "-w ended(\"${waitfor}\")"  "$cmd conf=$conf info_report $IRAP_PARAMS browsing=n"
}


function mapping_report {
    waitfor=$1
    mappers_dirs=`get_cached_value MAPPING_DIRS`
    report_qc_only=`get_param_value report_qc_only $conf`
    if [ "$report_qc_only-" != "y-" ] && [ "$mappers_dirs-" != "-" ] ; then
# one report for each mapper
	WAIT_FOR_IDS_=$WAIT_FOR_IDS
	declare -i mapper_id=1
	for d in $mappers_dirs; do
	    declare -i mapper_counter=0
	    mapper=`basename $d`
	    p_info "Bam files generated using $mapper..."
	    let c=1
	    FILES=$d/*.hits.bam
            #echo $FILES
	    if [ "$FILES" != '$d/*.hits.bam' ]; then
		for p in $pe; do 
		    let mapper_counter=mapper_counter+1
		    submit_job "${jobname_prefix}m${mapper_id}j[$mapper_counter]" "-w  ended($waitfor)"  $cmd conf=$conf $IRAP_PARAMS browsing=n se= pe=$p $report_dir/mapping/$mapper.html_doreq
		done
		for f in $se ; do
		    let mapper_counter=mapper_counter+1
		    submit_job "${jobname_prefix}m${mapper_id}j[$mapper_counter]" "-w  ended($waitfor)"  $cmd conf=$conf $IRAP_PARAMS browsing=n se=$f pe= $report_dir/mapping/$mapper.html_doreq
		done
		# when all libs are processed then merge the stats
		submit_job "${jobname_prefix}m${mapper_id}f[1]" "-w  ended(\"${jobname_prefix}m${mapper_id}j*\")"  $cmd conf=$conf $IRAP_PARAMS browsing=n $name/$mapper/featstats_raw.tsv
		submit_job "${jobname_prefix}m${mapper_id}f[2]" "-w  ended(\"${jobname_prefix}m${mapper_id}j*\")"  $cmd conf=$conf $IRAP_PARAMS browsing=n $name/$mapper/genestats_raw.tsv
		submit_job "${jobname_prefix}m${mapper_id}f[3]" "-w  ended(\"${jobname_prefix}m${mapper_id}j*\")"  $cmd conf=$conf $IRAP_PARAMS browsing=n $name/$mapper/stats_raw.tsv
	    fi
            # wait for the generation of all tsv files
	    submit_job "${jobname_prefix}mr[$mapper_id]"  "-w ended(\"${jobname_prefix}m${mapper_id}f*\")"  "$cmd conf=$conf $IRAP_PARAMS browsing=n $report_dir/mapping/$mapper.html"
	    let mapper_id=mapper_id+1
	done
        # wait for the reports of all mappers"22164809em2f*: N
	let counter=counter+1
	submit_job "${jobname_prefix}f[$counter]"  "-w ended(\"${jobname_prefix}mr*\")"  "$cmd conf=$conf $IRAP_PARAMS browsing=n $report_dir/mapping/comparison.html"
	echo ${jobname_prefix}f
    else
	echo $waitfor

    fi
}


function quant_report {
    waitfor=$1
    WAIT_FOR_IDS=

# GE
    p_info "GE: quant_report"
    quant_ofiles=`get_cached_value QUANT_HTML_FILES`
    for f in $quant_ofiles; do
	let counter=counter+1
	submit_job "${jobname_prefix}f[$counter]"  "-w ended(\"$waitfor\")"  "$cmd conf=$conf $IRAP_PARAMS browsing=n $f"
    done
}

function DE_report {
    waitfor=$1
####
# DE
    de_ofiles=`get_cached_value DE_HTML_FILES`
    p_info "DE"
    for f in $de_ofiles; do
	let counter=counter+1
	p_info "DE: $f"
	submit_job "${jobname_prefix}f[$counter]"  "-w ended(\"${waitfor}\")"  "$cmd conf=$conf  $IRAP_PARAMS browsing=n $f"
    done
}

function gsa_report {
    waitfor=$1
  #####
  # GSA
    p_info "GSA"
    gse_ofiles=`get_cached_value GSE_HTML_FILES`
    for f in $gse_ofiles; do
	let counter=counter+1
	p_info "GSA: $f"
	submit_job "${jobname_prefix}f[$counter]"  "-w ended(\"$waitfor\")"  "$cmd conf=$conf  $IRAP_PARAMS browsing=n $f"
    done
}

function finish_report {
# update the menu 
   submit_job "${jobname_prefix}u" -w "ended(\"${jobname_prefix}f*\")"  "$cmd conf=$conf report  $IRAP_PARAMS browsing=n"
   echo ${jobname_prefix}u
}

#
function submit_jobs4libs {
    waitfor=$1
    target=$2
    level=$3
    libs_ids=$all_libs
    if [ "$libs_ids-" == "-" ]; then
	DEPS=$1
    else
	OLDDEPS=$1
	DEPS="${jobname_prefix}${level}.*"
	local nx=($libs_ids)
	local NJOBS=${#nx[@]}
	p_info "Submiting a maximum of $NJOBS jobs"
	let JOBSPERBATCH=1000
	let BATCH=1
	BATCHGROUP=
	if [ $NJOBS -gt $JOBSPERBATCH ]; then
	    BATCHGROUP=b${BATCH}.
	    p_info "Too many jobs...grouping them in batches of $JOBSPERBATCH"
	else
	    BATCHGROUP=.
	fi
	let i=1	   
	for lib in $libs_ids; do
	    p_info "* Checking current status for $lib..."
	    $cmd conf=$conf $IRAP_PARAMS $lib $target -n -q 2> /dev/null
	    let ret=$?
	    if [ $ret -eq 0 ]; then
		p_info "Skipping submission of job for $lib"
	    else
		p_info "* Submitting jobs for $lib..."	    
		submit_job "${jobname_prefix}${level}${BATCHGROUP}$i"  -w "ended($waitfor)" $cmd conf=$conf  $IRAP_PARAMS $lib $target
		let i=$i+1
		if [ $i -gt $JOBSPERBATCH ]; then
		    p_info "Batch $BATCH submitted"
		    THREADS=1 MEM=4000 submit_job "${jobname_prefix}${level}.$BATCH"  -w "ended(${jobname_prefix}${level}${BATCHGROUP}*)" echo nop
		    let i=1
		    let BATCH=$BATCH+1
		    BATCHGROUP=b${BATCH}.
		fi
	    fi
	done	
	if [ $i -gt 1 ] && [ "$BATCHGROUP-" != ".-" ]; then
	    p_info "Batch $BATCH submitted"
	    THREADS=1 MEM=4000 submit_job "${jobname_prefix}${level}.$BATCH"  -w "ended(${jobname_prefix}${level}${BATCHGROUP}*)" echo nop
	fi
	if [ $i -eq 1  ]; then
	    ## no jobs where submitted
	    DEPS=$OLDDEPS
	fi
	if [ "$BATCHGROUP-" != ".-" ]; then
	    p_info "* Stage ${level}...$i jobs/$BATCH"
	else
	    p_info "* Stage ${level}...$i jobs"
	fi
	CUR_TARGETS=$target
    fi
    echo $DEPS
}

function submit_jobs4target {
    waitfor=$1
    level=$2
    shift 2
    target=$*

    DEPS="${jobname_prefix}${level}"
    p_info "* Checking current status for $target ($waitfor)..."
    $cmd conf=$conf $IRAP_PARAMS $target -n -q 2> /dev/null
    let ret=$?
    if [ $ret -eq 0 ]; then
	p_info "Skipping submission of job for $target"
	DEPS=$waitfor
    else
	p_info "* Submitting jobs for $target ($IRAP_PARAMS)..."
	submit_job "${jobname_prefix}${level}"  -w "ended($waitfor)" $cmd conf=$conf  $IRAP_PARAMS $target
    fi
    echo $DEPS
}
		
function get_pe_libs {
    #echo $cmd conf=$conf $IRAP_PARAMS print_pe_libs se= > /dev/stderr
    $cmd conf=$conf $IRAP_PARAMS print_pe_libs se= |tail -n 1
}
function get_se_libs {
    #echo $cmd conf=$conf $IRAP_PARAMS print_se_libs pe= > /dev/stderr
    $cmd conf=$conf $IRAP_PARAMS print_se_libs pe= | tail -n 1
}
