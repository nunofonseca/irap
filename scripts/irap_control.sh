#!/bin/env bash
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
# =========================================================

# Script to be executed as cron

TOPLEVEL_FOLDER=$1

if [ "$TOPLEVEL_FOLDER-" == "-" ]; then
    echo "ERROR: usage: irap_control.sh TOPLEVEL_FOLDER [id]"
    exit 1
fi

# optional
ID=$2


if [ -e $TOPLEVEL_FOLDER/irap_control.config ]; then
    source $TOPLEVEL_FOLDER/irap_control.config
fi

if [ "$THREADS-" == "-" ]; then
    THREADS=4
fi
if [ "$QUEUE-" == "-" ]; then    
    QUEUE=production-rh7
fi
if [ "$LSG_GROUP-" == "-" ]; then
    LSF_GROUP=/irap_sc1
fi
if [ "$MEM1-" == "-" ]; then
    MEM1=12000
fi
if [ "$MEM2-" == "-" ]; then
    MEM2=20000
fi
if [ "$MEM3-" == "-" ]; then
    MEM3=30000
fi
if [ "$MEM4-" == "-" ]; then
    MEM4=45000
fi
if [ "$MEM5-" == "-" ]; then
    MEM5=55000
fi
if [ "$MEM6-" == "-" ]; then
    MEM6=65000
fi
if [ "$MEM7-" == "-" ]; then
    MEM7=85000
fi
max_mem_level=7
control_folder=$TOPLEVEL_FOLDER/.control

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
###########################################################

set -eu
set -e
function filepath2sdrf_id {
    x=`basename $1|sed "s/.txt$//;s/.sdrf.*//"`
    echo $x
}

function cached_sdrf_file {
    echo $control_folder/`basename $1`    
}

pushd $TOPLEVEL_FOLDER
mkdir -p $control_folder


if [ "$ID-" == "-" ]; then
    SDRF_FILES=`grep conf .control/*.sdrf.status|cut -f 1 -d:`
    files=( $SDRF_FILES )
    echo "Found ${#files[@]} SDRF files"
else
    SDRF_FILES=".control/$ID.sdrf.status"
    if [ ! -e $SDRF_FILES ]; then
	echo $SDRF_FILES not found
	exit 1
    fi    
fi
set +e

function kill_id_jobs {
    set +e    
    lsf_list|grep $LSF_GROUP|grep "=$id.conf"|cut -f 1 -d,|while read n; do bkill $n; done
    set -e
}

function id2runstatus {
    echo $control_folder/$1.run.status
}

function id2conf_fp {
    local id=$1
    cat $(id2runinfo $id)
}

function id2runinfo {
    echo $control_folder/$1.run.info
}

function id2sdrf {
    echo $control_folder/$1.sdrf.txt
}
function id2idf {
    echo $control_folder/$1.idf.txt
}
function rundir2bundle_dir {
    echo $1/$(basename $1)/sc_bundle
}
##
set -e
for f in $SDRF_FILES; do
    id=`filepath2sdrf_id $f`
    rs=`id2runstatus $id`
    ri=`id2runinfo $id`
    csf=`cached_sdrf_file $f`
    echo "rs=$rs ri=$ri id=$id"
    if [ ! -e $ri ]; then 
	## new 
	mfolder=`find . -maxdepth 2 -name "$id"`	
	if [ "$mfolder-" == "-" ]; then
	    echo "unable to find $id.conf"
	    exit 1
	fi
	echo mfolder=$mfolder
	pf="$mfolder/$id.conf"
	echo $pf > $ri
	echo new > $rs
    else
	pf=`cat $ri`	
    fi
    if [ ! -e $pf ]; then
	echo "Oops...unable to find $pf"
	exit 1
    else
	echo Found $pf
    fi
    ## conf more recent than the status file
    ## and there is a difference in the conf file
    if [ $pf -nt $rs ]; then
	#diff -q TODO
	#res=$?
	#if [ $res != 0 ]; then
	    # modified
	    echo Modified $id
	    ## killing jobs if necessary
	    kill_id_jobs $id
	    ##
	    echo mod > $rs
	#fi
    fi
done

## double check theat the bundle is created
## and copy the sdrf/idf files to the bundle folder
function is_all_done {
    local id=$1
    local run_dir=$2
    
    conf=$(id2conf_fp $id)

    pushd $run_dir >/dev/null
    irap_sc conf=$conf atlas_bundle > /dev/stderr
    let ret=$?
    ##echo ret=$ret > /dev/stderr
    popd >/dev/null
    if [ "$ret-" != "0-" ]; then
	echo "n"
    else
	for f in `id2sdrf $id` `id2idf $id`; do
	    cp -a $f `rundir2bundle_dir $run_dir`
	done
	echo y
    fi
}
set -e
########################################################
##
## status=new,mod,run,conf (generated and data downloaded)
function rerun_wrapper {
   sdrf=$1
   status=`get_sdrf_status $sdrf|cut -f 1 -d\ `
   if [ "$status" == "run" ]; then
       stop_wrapper $sdrf
   fi
   run_wrapper $*
}

function run_wrapper {
    local id=$1
    local status=$2
    local wd=$3
    shift 3
    
    pushd $wd
    set -e
    if [ $status == "runa" ]; then
	jid=`bash -c "set -o pipefail; $* | tail -n 1 |cut -f 2 -d="`
	rets=$?
    else
	jid=`bash -c "set -o pipefail;$* | tail -n 1 |cut -f 1 -d\>|cut -f 2 -d\<"`
	rets=$?
    fi
    popd
    echo jid=$jid
    if [ $rets -ne 0 ]; then
	echo "Unable to submit job"
	exit 1
    fi
    ##
    if [ "$jid" == "All done - no need to submit jobs" ]; then
	jid="DONE"
    fi
    ## change status
    set_run_status $id $status $jid
}
function stop_wrapper {
    local id=$1
    local jobid=`get_run_status $id|cut -f 2 -d\ `
    kill_id_jobs $id
    ## change status
    set_run_status $id mod
}
function status_wrapper {
    local jobid=`get_run_status $1|cut -f 2 -d\ `
    jobstatus=`get_run_status $1|cut -f 1 -d\ `
    status=
    set +e
    if [ "$jobid" == "DONE" ]; then
	status="DONE"
    else
	case $jobstatus in
	    runa)
		status=`bjobs -o 'stat' -J $jobid | tail -n 1`
		;;
	    runb)
		status=`bjobs -o 'stat' $jobid | tail -n 1`
		;;
	esac
    fi
    set -e
    echo $status
}

## status: new/mod/unchanged
function set_run_status {
    local id=$1
    shift 1
    echo $* > $control_folder/$id.run.status
}

function get_run_status {
    local id=$1
    if [ -e  $control_folder/$id.run.status ]; then
	cat $control_folder/$id.run.status
    else
	set_run_status $id new
	get_run_status $id
    fi
}

function set_mem_level {
    local id=$1
    shift 1
    echo $* > $control_folder/$id.run.mem
}

function get_mem_level {
    local id=$1
    if [ -e  $control_folder/$id.run.mem ]; then
	cat $control_folder/$id.run.mem
    else
	set_mem_level $id 1
	get_mem_level $id
    fi
}



########################################################
## 
echo Handling jobs...
## Now, handle the jobs
for f in $SDRF_FILES; do
    id=`filepath2sdrf_id $f`
    rs=`id2runstatus $id`
    ri=`id2runinfo $id`
    conf=`cat $ri`
    run_dir=`dirname $conf`
    conf2=`basename $conf`
    scprot=`grep sc_protocol= $conf | cut -f 2 -d=`
    status=`get_run_status $id|cut -f 1 -d\ `
    echo id=$id $status
    case $status in
	new)
	    run_wrapper $id runa $run_dir IRAP_LSF_GROUP=$LSF_GROUP THREADS=$THREADS MEM=$MEM1 QUEUE=$QUEUE irap_lsf2 -s conf=$conf2
	    ;;
	mod)
	    run_wrapper $id runa $run_dir IRAP_LSF_GROUP=$LSF_GROUP THREADS=$THREADS MEM=$MEM1 QUEUE=$QUEUE irap_lsf2 -s conf=$conf2
	    ;;
	reruna)
	    rs=`get_mem_level $id`
	    v="MEM$rs"
	    echo $rs
	    if [ "$rs-" != "$max_mem_level-" ]; then
		MEMX=${!v}
		run_wrapper $id runa $run_dir IRAP_LSF_GROUP=$LSF_GROUP THREADS=$THREADS MEM=$MEMX QUEUE=$QUEUE irap_lsf2 -s conf=$conf2
	    fi
	    ;;
	rerunb)
	    run_wrapper $id runb $run_dir bsub -M $MEM irap_sc conf=$conf2 atlas_bundle
	    ;;
	runa) s=`status_wrapper $id`
	      if [ "$s-" == "DONE-" ]; then
		  set_run_status $id rerunb
	      else
		  if [ "$s-" != "RUN-" ] &&  [ "$s-" != "PEND-" ];  then
		      ml=`get_mem_level $id`
		      if [ $ml -ge $max_mem_level ]; then
			  set_run_status $id onhold
		      else
			  set_mem_level $id `expr $ml + 1`
			  set_run_status $id reruna
			  echo "$id processing failed: $s $ml"
		      fi
		  else
		      echo "$id $s"
		  fi
	      fi
	      ;;
	runb) s=`status_wrapper $id`
	     if [ "$s-" == "DONE-" ]; then
		 set_run_status $id run_complete
		 echo $id $run_dir run_complete
	     else
		 if [ "$s-" != "RUN-" ] &&  [ "$s-" != "PEND-" ];  then
		     rs=`get_mem_level $id`
		     if [ "$rs-" == "6-" ]; then
			 set_run_status $id onhold $run_dir
		     else
			 set_run_status $id rerunb $run_dir
		     fi
		     echo "$id analysis failed"
		 else
		     echo $id $s
		 fi
	     fi
	     ;;
#	complete)
	run_complete)
	    s=`is_all_done $id $run_dir`
	    echo "s=$s" > /dev/stderr
	    if [ "$s" == "y" ]; then
		# set as all done
		ddd=$(readlink -f `rundir2bundle_dir $run_dir`)
		set_run_status $id all_done $ddd
		echo $id $run_dir all_done $ddd
		touch all.done.txt
		echo $id $ddd >> all.done.txt
		sort -u all.done.txt > all.tmp && mv all.tmp all.done.txt
		
	    else
		set_run_status $id complete_onhold
		echo $id $run_dir complete_onhold
	    fi
	    ;;
	all_done) echo $id all_done;;
	onhold) echo $id onhold;;
	complete_onhold) echo $id complete_onhold;;
    esac
done
exec_on_exit
exit 0
./sdrf_control.sh /nfs/production3/ma/home/atlas3-production/singlecell/experiment /homes/nf/storage3/EA/sc

