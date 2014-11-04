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

_INIT_SUSP_JOB=
function irap_init_job {
 #############
 # Submit jobs
 p_info " * Initialization "
 # save variables
 $cmd conf=$conf $IRAP_PARAMS save_cache

 #load the file with the variables already computed
 export IRAP_PARAMS="$IRAP_PARAMS use_cached_vars=y"
 p_info " * Starting initial job"
 CUR_STAGE=stage0
 submit_job "${jobname_prefix}i" $cmd conf=$conf setup_dirs $IRAP_PARAMS
 stop_job  "${jobname_prefix}i"
 _INIT_SUSP_JOB=${jobname_prefix}i
 p_info " * First job suspended until all jobs are submitted."
 echo "${jobname_prefix}i"
}

function stage0_jobs {
    waitfor=$1
    stage0_targets=`irap conf=$conf $IRAP_PARAMS print_stage0_files|tail -n 1`
    let i=2
    for f in $stage0_targets; do    
	submit_job "${jobname_prefix}0[$i]"  -w "ended($waitfor)" $cmd conf=$conf setup_dirs $IRAP_PARAMS $f
	let i=$i+1
    done
}

function stage1_jobs {    
    waitfor=$1

    declare -i s1=0
    declare -i s2=0

    # Paired end files
    CUR_STAGE=stage1
    p_info "* Step&1 (PE)"
    for p in $pe ; do 
	let s1=s1+1
	p_info "Lib (PE): $p"
	CUR_STAGE=stage1 submit_job "${jobname_prefix}s1[${s1}]" -w "ended($waitfor)" "$cmd conf=$conf  pe=$p se=  $IRAP_PARAMS stage1"
    done
    # Single end files
    echo "***********************************************"
    echo "*Step1&2 (SE)"
    for f in $se ; do
	let s1=s1+1
	p_info "Lib (SE): $f"
	CUR_STAGE=stage1 submit_job "${jobname_prefix}s1[$s1]" -w  "ended($waitfor)"  "$cmd conf=$conf stage1 pe=  se=$f  $IRAP_PARAMS"
    done
    echo "${jobname_prefix}s1*"
}

# run stage1&2
function stage12_jobs {    
    waitfor=$1
    declare -i s1=0
    declare -i s2=0
    # Paired end files
    p_info "* Step&21 (PE)"
    for p in $pe ; do 
	let s1=s1+1
	let s2=s2+1
	p_info "Lib (PE): $p"
	CUR_STAGE=stage12 submit_job "${jobname_prefix}s12[${s1}]" -w "ended($waitfor)" "$cmd conf=$conf  pe=$p se=  $IRAP_PARAMS stage1 stage2"
    done

    # Single end files
    echo "***********************************************"
    echo "*Step1&2 (SE)"
    for f in $se ; do
	let s2=s2+1
	let s1=s1+1
	p_info "Lib (SE): $f"
	CUR_STAGE=stage12 submit_job "${jobname_prefix}s12[${s1}]" -w "ended($waitfor)" "$cmd conf=$conf  pe= se=$f  $IRAP_PARAMS stage1"
    done
    echo "${jobname_prefix}s12*"
}

# run stage1&2&3
function stage123_jobs {    
    waitfor=$1
    declare -i s1=0
    declare -i s2=0
    # Paired end files
    p_info "* Step&21 (PE)"
    for p in $pe ; do 
	let s1=s1+1
	let s2=s2+1
	p_info "Lib (PE): $p"
	CUR_STAGE=stage123 submit_job "${jobname_prefix}s123[${s1}]" -w "ended($waitfor)" "$cmd conf=$conf  pe=$p se=  $IRAP_PARAMS stage1 stage2 stage3as"
    done

    # Single end files
    echo "***********************************************"
    CUR_STAGE=stage123
    echo "*Step1&2 (SE)"
    for f in $se ; do
	let s2=s2+1
	let s1=s1+1
	p_info "Lib (SE): $f"
	CUR_STAGE=stage123 submit_job "${jobname_prefix}s123[${s1}]" -w "ended($waitfor)" "$cmd conf=$conf  pe= se=$f  $IRAP_PARAMS stage1 stage2 stage3as"
    done
    echo "${jobname_prefix}s123*"
}

# run stage3 - Normalization/merging
function stage3_jobs {    
    waitfor=$1
    WAIT_FOR_IDS=  CUR_STAGE=stage3 

    p_info "Stage3"
    submit_job "${jobname_prefix}q"  -w "ended($waitfor)"  "$cmd conf=$conf  $IRAP_PARAMS stage3"
    echo "${jobname_prefix}q"
}

function DE_jobs {
    #####
    # DE
    waitfor=$1
    CUR_STAGE=DE
    de_ofiles=`irap conf=$conf de_files $IRAP_PARAMS|tail -n 1`
    let s2=1
    for f in $de_ofiles; do
	let s2=s2+1
	p_info "DE: $f"
	submit_job "${jobname_prefix}de[$s2]"  `check_dependency $waitfor`  "irap conf=$conf  $IRAP_PARAMS $f"
    done
    if [ "$de_ofiles-" = "-" ]; then
	DEP=`check_dependency ${jobname_prefix}q`
    else
	DEP=`check_dependency ${jobname_prefix}de`
    fi
    echo $DEP
}

function GSA_jobs {
    DEP=$1
    ######
    # GSA
    CUR_STAGE=GSA
    if [ "$de_ofiles-" != "-" ]; then
        #####
        # GSA
	gse_ofiles=`irap conf=$conf $IRAP_PARAMS GSE_files|tail -n 1`
	let s2=1
	for f in $gse_ofiles; do
	    let s2=s2+1
	    p_info "GSE: $f"
	    submit_job "${jobname_prefix}gse[$s2]" $DEP  "irap conf=$conf  $IRAP_PARAMS $f"
	done
	if [ "$gse_ofiles-" != "-" ]; then
	    DEP=`check_dependency ${jobname_prefix}gse`
	fi    
    fi
    echo $DEP
}

function final_job {
    CUR_STAGE=
    DEP=$1
    #########################################
    # To finalize, run the whole pipeline
    # If everything went ok then nothing should be done
    # otherwise it should fail and an email will be sent
    CUR_STAGE= submit_job "${jobname_prefix}f" $DEP   "irap conf=$conf  $IRAP_PARAMS -n -q"
    
    # send an email to the user
    CUR_STAGE= submit_job_status "${jobname_prefix}f"
    
    resume_job $_INIT_SUSP_JOB
    # TODO: report
    echo JOBS=$LOG_DIR/${jobname_prefix}*.out
    echo JOBNAME=${jobname_prefix}f
}

##########################################################
# Report
declare -i counter=0


# Report jobs: with the exception of the initial jobs, most jobs have the same prefix
function qc_report {
    waitfor=$1
    echo " * QC Report "
    let counter=counter+1
    submit_job "${jobname_prefix}f[$counter]" "-w ended(\"$waitfor\")"  "$cmd conf=$conf qc_report $IRAP_PARAMS"
}

function init_page_report {
    waitfor=$1    
    echo " * Information page "
    let counter=counter+1
    submit_job "${jobname_prefix}f[$counter]" "-w ended(\"${waitfor}\")"  "$cmd conf=$conf info_report $IRAP_PARAMS"
}


function mapping_report {
    waitfor=$1
    mappers_dirs=`get_cached_value MAPPING_DIRS`
    report_qc_only=`get_param_value report_qc_only $conf`
    if [ "$report_qc_only-" != "y-" ] && [ "$mappers_dirs-" != "-" ] ; then
	declare -i mapper_counter=0
# one report for each mapper
	WAIT_FOR_IDS_=$WAIT_FOR_IDS
	for d in $mappers_dirs; do
	    mapper=`basename $d`
	    p_info "Bam files generated using $mapper..."
	    let c=1
	    FILES=$d/*.hits.bam
    #echo $FILES
	    if [ "$FILES" != '$d/*.hits.bam' ]; then
		let mapper_counter=mapper_counter+1
		if [ "$DEBUG-" = "1-" ]; then	    
		    p_info "running irap_bam_report_lsf $mapper"
		else
		    NEW_JOB=`WAIT_FOR_IDS=${waitfor} irap_bam_report_lsf conf=$conf $IRAP_PARAMS mapper=$mapper | grep "JOB=" | cut -f 2 -d\=`	
		    p_info "running irap_bam_report: $NEW_JOB"
		    if [ "$NEW_JOB-"  != "-" ] ; then		
			WAIT_FOR_IDS=
			submit_job "${jobname_prefix}m[$mapper_counter]"  "-w  ended($NEW_JOB)"  "irap conf=$conf $IRAP_PARAMS mapper=$mapper quant_method=none quant_norm_method=none de_method=none gse_tool=none $name/report/mapping/$mapper.html"
			WAIT_FOR_IDS=$WAIT_FOR_IDS_
		    else 
			p_info "Failed to start irap_bam_report_lsf"
		    fi
		fi
	    fi
	done
        # wait for the reports of all mappers
	let counter=counter+1
	submit_job "${jobname_prefix}f[$counter]"  "-w ended(\"${jobname_prefix}m*\")"  "irap conf=$conf $IRAP_PARAMS $report_dir/mapping/comparison.html"
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
	submit_job "${jobname_prefix}f[$counter]"  "-w ended(\"$waitfor*\")"  "irap conf=$conf $IRAP_PARAMS $f"
    done
}

function DE_report {
    waitfor=
####
# DE
    de_ofiles=`get_cached_value DE_HTML_FILES`
    p_info "DE"
    for f in $de_ofiles; do
	let counter=counter+1
	p_info "DE: $f"
	submit_job "${jobname_prefix}f[$counter]"  "-w ended(\"${waitfor}\")"  "irap conf=$conf  $IRAP_PARAMS $f"
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
	submit_job "${jobname_prefix}f[$counter]"  "-w ended(\"$waitfor\")"  "irap conf=$conf  $IRAP_PARAMS $f"
    done
}

function finish_report {
# update the menu 
   submit_job "${jobname_prefix}u" -w "ended(\"${jobname_prefix}f*\")"  "$cmd conf=$conf report  $IRAP_PARAMS"
   echo ${jobname_prefix}u
}
