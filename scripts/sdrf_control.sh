#!/bin/env bash
# Script to be run as cron

SDRF_FOLDER=$1
TOPLEVEL_FOLDER=$2

QUEUE=production-rh7
LSF_GROUP=/download
control_folder=$TOPLEVEL_FOLDER/.control
jobs_folder=$TOPLEVEL_FOLDER/.control/jobs
EMAIL_CC=nf@ebi.ac.uk
MEM=8000
has_changes=n


if [ "$TOPLEVEL_FOLDER-" == "-" ]; then
    echo "ERROR: usage: sdrf_control.sh SDRF_FOLDER TOPLEVEL_FOLDER"
    exit 1
fi

if [ ! -e $SDRF_FOLDER ]; then
    echo "$SDRF_FOLDER  file not found"
    exit 1
fi

#set -eux
set -e
function filepath2sdrf_id {
    x=`basename $1|sed "s/.txt$//;s/.sdrf//"`
    echo $x
}

function cached_sdrf_file {
    echo $control_folder/`basename $1`    
}

function quiet_echo {
    if [ "$nchanges-" != "0-" ]; then
	echo $*
    fi
}

############################################################
# 
lock_file=$control_folder/lock
## lock
if [ -e $lock_file ]; then
    echo "ERROR: Unable to continue - found lock file $lock_file"
    exit 2
fi

function exec_on_exit {
    rm -f $lock_file
}
trap exec_on_exit EXIT
touch $lock_file


pushd $TOPLEVEL_FOLDER &> /dev/null
mkdir -p $control_folder
mkdir -p $jobs_folder

SDRF_FILES=`ls --color=never $SDRF_FOLDER/*.sdrf.txt`
files=( $SDRF_FILES )
set +e

## keep a copy of the sdrf file
let new_sdrfs=0
let mod_sdrfs=0
let new_idfs=0
let mod_idfs=0
for f in $SDRF_FILES; do
    bf=`basename $f`
    cbf=`cached_sdrf_file $f`
    cidf=$(echo $cbf|sed "s/.sdrf./.idf./")
    idf=$(echo $f|sed "s/.sdrf./.idf./")
    if [ ! -e $cbf ]; then 
	## new 
	cp $f $cbf
	echo "New SDRF: $bf"
	let new_sdrfs=$new_sdrfs+1
    else
	if [  $f -nt $cbf ]; then
	    cp $f $cbf
	    echo "Updated SDRF: $bf"
	    let mod_sdrfs=$mod_sdrfs+1
	fi
    fi
    if [ ! -e $cidf ]; then 
	## new 
	cp $idf $cidf
	## force update
	touch $cbf
	echo "New IDF: $bf"
	let new_idfs=$new_idfs+1
    else
	if [  $idf -nt $cidf ]; then
	    cp $idf $cidf
	    ## force update
	    touch $cbf
	    echo "Updated SDRF: $bf"
	    let mod_idfs=$mod_idfs+1
	fi
    fi
done

set -e
nchanges=$(echo $new_sdrfs + $mod_sdrfs + $new_idfs + $mod_idfs |bc )
quiet_echo "Found ${#files[@]} SDRF files"
quiet_echo New sdrfs: $new_sdrfs
quiet_echo Modified sdrfs: $mod_sdrfs
quiet_echo New idfs: $new_idfs
quiet_echo Modified idfs: $mod_idfs

########################################################
##
## status=new,mod,run,conf (generated and data downloaded)
function rerun_wrapper {
   sdrf=$1
   status=`get_sdrf_status $sdrf|cut -f 1 -d\ `
   ## kill the previous job if it was already running
   if [ "$status" == "run" ]; then
       stop_wrapper $sdrf
   fi
   run_wrapper $*
}
function run_wrapper {
    sdrf=$1
    shift 1
    set -o pipefail
    rm -f $jobs_folder/$sdrf.log
    id=`bsub  -o $jobs_folder/$sdrf.log -M $MEM -g $LSF_GROUP $* | tail -n 1 |cut -f 1 -d\>|cut -f 2 -d\<`
    if [ $? -ne 0 ]; then
	echo "Unable to submit job"
	exit 1
    fi
    if [ "$id-" == "-" ]; then
	echo "Unable to get job id"
	exit 1
    fi
    ## change status
    set_sdrf_status $sdrf run $id
 }
function stop_wrapper {
    sdrf=$1
    jobid=`get_sdrf_status $sdrf|cut -f 2 -d\ `
    bkill $jobid
    ## change status
    set_sdrf_status $sdrf mod
}
function status_wrapper {
    sdrf=$1
    jobid=`get_sdrf_status $sdrf|cut -f 2 -d\ `
    jobstatus=`get_sdrf_status $sdrf|cut -f 1 -d\ `
    if [ "$jobstatus-" == "run-" ]; then
	status=`bjobs -o 'stat'  $jobid | tail -n 1`
	if [ "$status-" == "Job <$jobid> is not found" ]; then
	    x=`tail -n 1 $jobs_folder/$sdrf.log | grep -c "Files placed"`
	    if [ $x -eq 1 ]; then
		status=DONE
	    fi
	fi
	echo $status
    fi    
}
 
## was the sdrf modified,created?
function get_file_status {
    f1=$1
    bf=`filepath2sdrf_id $f`
    if [ ! -e "$control_folder/$bf.sdrf.status" ]; then 
	echo new
    else
	if [  "$control_folder/$bf.sdrf.txt" -nt "$control_folder/$bf.sdrf.status" ]; then
	    echo mod
	    set_sdrf_status $bf mod
	else
	    echo nochange
	fi
    fi    
}
## status: new/mod/unchanged
function set_sdrf_status {
    id=$1
    shift 1
    echo $* > $control_folder/$id.sdrf.status
}
function get_sdrf_status {
    id=$1
    if [ -e  $control_folder/$id.sdrf.status ]; then
	cat $control_folder/$id.sdrf.status
    else
	set_sdrf_status $id new
	get_sdrf_status $id
    fi
}

function notify_by_email {
    suser=$1
    sdrf=$2
    status=$3
    owner=$EMAIL_CC

    echo "mail -b $owner -s 'SDRF processing: $sdrf status=$status' $suser@ebi.ac.uk <  $jobs_folder/$sdrf.log"
    gzip -c $jobs_folder/$sdrf.log > $jobs_folder/$sdrf.log.txt.gz
    tail -n 80 $jobs_folder/$sdrf.log | mail -b $owner -a  $jobs_folder/$sdrf.log.txt.gz -s "SDRF processing: $sdrf status=$status" $suser
    rm -f $jobs_folder/$sdrf.log.txt.gz
}
########################################################
## Now, handle the new/modified files
for f in $SDRF_FILES; do
    bf=`filepath2sdrf_id $f`
    cbf=`cached_sdrf_file $f`
    status=`get_file_status $bf|cut -f 1 -d\ `
    case $status in
	new) run_wrapper $bf irap_AE_setup.sh -t $TOPLEVEL_FOLDER -c -i $cbf -a;;
	mod) rerun_wrapper $bf irap_AE_setup.sh -t $TOPLEVEL_FOLDER -c -i $cbf -a;;
	nochange);; ## do nothing
    esac
done

quiet_echo Checking jobs...
## Now, handle the jobs
for f in $SDRF_FILES; do
    bf=`filepath2sdrf_id $f`
    cbf=`cached_sdrf_file $f`
    status=`get_sdrf_status $bf|cut -f 1 -d\ `
    suser=`stat -c "%U" $f`
    case $status in
	run) s=`status_wrapper $bf`
	     if [ "$s-" == "DONE-" ]; then
		 set_sdrf_status $bf conf
		 notify_by_email $suser $bf Success
		 rm -f $jobs_folder/$bf.log
		 echo $bf complete/conf
	     else
		 if [ "$s-" != "RUN-" ] &&  [ "$s-" != "PEND-" ];  then
		     set_sdrf_status $bf failed
		     notify_by_email $suser $bf Failed
		     rm -f $jobs_folder/$bf.log
		     echo "$bf conversion failed: $s"
		 else
		     quiet_echo $bf $s
		 fi
	     fi
	     ;;
	mod) rerun_wrapper $bf irap_AE_setup.sh -t $TOPLEVEL_FOLDER -c -i $cbf -a
	     echo $bf RUN
	     ;;
	failed) quiet_echo $bf failed;;
	nochange) quiet_echo $bf nochange;;
	conf) quiet_echo $bf conf;; ## do nothing
    esac
done
exec_on_exit
exit 0
./sdrf_control.sh /nfs/production3/ma/home/atlas3-production/singlecell/experiment /homes/nf/storage3/EA/sc

