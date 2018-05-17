#!/usr/bin/env bash
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

## requires bash 4.4
#
# Download data from ArrayExpress and setup the files necessary
# to run iRAP
function usage {
    echo "irap_AE_setup.sh [ -s  -t root_dir  -p prot -d data_dir -c -h] -i id"
    echo <<EOF
 -t root_dir - toplevel directory where a folder <id> will be created and the fastq files will be downloaded
 -p prot     - protocol (blk,10xV1,drop-seq,smart-seq2,...)
 -d data_dir - path to iRAP's data directory
 -i id       - array express experiment ID (e.g., E-MTAB-4411)
 -c       - skip downloading if files are already present
 -C       - skip downloading or checking if files are already present
 -S       - skip checking IDF
 -s       - single cell
 -a       - try to download private data (only works at EBI)
 -T threads - maximum number of threads that may be used to download the files
 -P program_name - name of the script/program to be used to download the data (instead of wget). Shoud take two arguments: 1) uri of the file to download; 2) local filename.
 -b batch_number 
 -B batch_size   - number of cells/samples per batch
EOF
}

function pinfo {
    echo "== $*"
}

function perror {
    echo "ERROR: $*" > /dev/stderr
}

## fix the filename by replacing some characters
## that may complicate things
function fix_filename {
    echo ${1//\#/_};
}

function create_dir {
    FOLDER=$1
    if [ ! -e $FOLDER ]; then
	pinfo creating $FOLDER
	set -e
	mkdir -p $FOLDER
	set +e
    fi
}

function get_column {
    sdrf_file=$1
    colname=$2
    exit_on_error=$3
    COL=""
    set -e
    COL=$(head -n 1 $sdrf_file | tr "\t" "\n" | grep -i -n -F "$colname" | cut -f 1 -d:)
    set +e
    if [ "$COL-" == "-" ] && [ "$exit_on_error-" != "-" ]; then
	perror column $colname not found in $sdrf_file
	exit 1
    fi
    echo $COL
}
####################################
# Look for the species column in the SDRF
function get_species {
    COL=""
    COL=$(get_column $1 "[organism]")
    if [ "$COL-" == "-" ]; then
	perror organism not found in $1
	exit 1
    fi
    if [ "$(echo $COL | wc -w)" != "1" ]; then
	cols="$(echo $COL| tr ' ' ',')"
	perror found multiple organism columns in $1: $(head -n 1 $1|cut -f $cols)
	exit 1
    fi
#    NSPECIES=$(tail -n +2 $SDRF_FILE|cut -f $COL |sort -u |wc -l)
#    if [ $NSPECIES -ne 1 ]; then
#	perror "SDRF files with data from different ($NSPECIES) organisms not supported yet"
#	exit 1
    #    fi
    tail -n +2 $1|cut -f $COL |sort -u
}

###################################
ROOT_DIR=$PWD
DATA_DIR=$ROOT_DIR/data
PROTOCOL=
ID=
DEBUG=0
USE_CACHE=0
AT_EBI=0
DOWNLOAD_SDRF=n
download_private_data=n
skip_idf=n
be_quiet=n
THREADS=1
batch_number=0
batch_size=0
skip_file_checking=0
download_program=

while getopts "b:B:T:t:p:d:i:P:CahceDdsSq"  Option; do
    case $Option in
	P ) download_program=$OPTARG;;
	a ) download_private_data=y;;
	e ) AT_EBI=1;;
	c ) USE_CACHE=1;;
	C ) skip_file_checking=1;;
	D ) DOWNLOAD_SDRF="y";;
	s ) irap_cmd=irap_sc;;
        t ) ROOT_DIR=$OPTARG;;
	p ) PROTOCOL=$OPTARG;;
	d ) DATA_DIR=$OPTARG;;
	S ) skip_idf="y";;
	i ) ID=$OPTARG;;
	q ) be_quiet="y";;
	T ) THREADS=$OPTARG;;
	B ) batch_size=$OPTARG;;
	b ) batch_number=$OPTARG;;
	h) usage; exit;;
    esac
done

if [ "$ID-" == "-" ]; then
    perror "no value passed to -i"
    usage
    exit 1
fi

habemus_batches=0
if [ "$batch_size-" != "0-" ]; then
    habemus_batches=1
fi
set +e
set -o pipefail
# folder where previously downloaded data may be found
pinfo "root_dir=$ROOT_DIR"
pinfo "data_dir=$DATA_DIR"
pinfo "protocol=$PROTOCOL"

if [ $habemus_batches != 0 ]; then
    pinfo "batch=$batch_number"
    pinfo "batch size=$batch_size"
fi
###################################

if [ "$DOWNLOAD_SDRF" == "y" ]; then
    pinfo "DOWNLOADING SDRF "
    ###################################
    ## Download sdrf
    SDRF_ROOT_URL="https://www.ebi.ac.uk/arrayexpress/files"

    SDRF_FILE=$ID.sdrf.txt
    SDRF_FILE_FP=$PWD/$ID.sdrf.txt
    SDRF_URL="$SDRF_ROOT_URL/$ID/$SDRF_FILE"
    IDF_FILE=$ID.idf.txt
    IDF_FILE_FP=$PWD/$ID.idf.txt
    IDF_URL="$SDRF_ROOT_URL/$ID/$IDF_FILE"
    
    if [ $USE_CACHE -eq 1 ]  &&  [ -e $SDRF_FILE_FP ]; then
	pinfo "skipping downloading $SDRF_FILE"
    else
	pinfo "downloading $SDRF_URL"
	set -e
	wget  --show-progress -c $SDRF_URL -O $SDRF_FILE_FP
	set +e
    fi
    if [ $USE_CACHE -eq 1 ]  &&  [ -e $IDF_FILE_FP ]; then
	pinfo "skipping downloading $IDF_FILE"
    else
	pinfo "downloading $IDF_URL"
	set -e
	wget  --show-progress -c $IDF_URL -O $IDF_FILE_FP
	set +e
    fi
else
    if [ ! -e $ID ]; then
	perror "File $ID not found"
	exit 1
    fi
    SDRF_FILE=$(basename $ID)
    SDRF_FILE_FP=$(readlink -f $ID)
    ID=$(basename $ID|sed 's/.sdrf.txt.*//')
    IDF_FILE=$(dirname $SDRF_FILE_FP)/$ID.idf.txt
    if [ ! -e $IDF_FILE ]; then
	perror "File $ID_FILE not found"
	exit 1
    fi
    IDF_FILE_FP=$(readlink -f $IDF_FILE)    
fi

## to keep backwards compatibility do not distribute fastq files
## if the number of entries/rows in the SDRF is below 5100
distribute_by_subfolders=0
N_FQ=$(cut -f 1 $SDRF_FILE_FP|wc -l )
if [ $N_FQ -gt 5090 ]; then
    distribute_by_subfolders=1
fi

## avoid keeping the local download folder under ROOT_DIR
local_download_folder=$(readlink -f $ROOT_DIR/../ManuallyDownloaded/$ID)
pinfo "previously downloaded data folder=$local_download_folder"
###################################
set -e
#set -ex
pinfo "id=$ID"
TOP_FOLDER=$ROOT_DIR/tmp/$ID
SPECIES_CONF_DIR=$ROOT_DIR/species_conf
CONF_FILE_PREF=$TOP_FOLDER/$ID
CONF_FILE=$CONF_FILE_PREF.conf
create_dir $TOP_FOLDER
create_dir $SPECIES_CONF_DIR
pushd $TOP_FOLDER 2>/dev/null
SPECIES=`get_species $SDRF_FILE_FP|tr ' ' '_'|tr 'A-Z' 'a-z'`

pinfo species=$SPECIES
if [ "$SPECIES-" == "-" ]; then
    exit 1
fi

###################################
##
NL=$(wc -l $SDRF_FILE_FP|cut -f 1 -d\ )
pinfo $SDRF_FILE with $NL lines
####################################
# batches
if [ $habemus_batches != 0 ]; then
    ##
    pwd
    cd $ROOT_DIR
    tmp_file_pref=$ROOT_DIR/$SPECIES/$ID/.$(basename $SDRF_FILE_FP)
    tmp_file_lock=$tmp_file_pref.lock
    tmp_file_split=$tmp_file_pref.split.$batch_size
    mkdir -p $ROOT_DIR/$SPECIES/$ID/
    if [ -e $tmp_file_split.1 ] && [ $tmp_file_split.1 -nt $SDRF_FILE_FP ]; then
	echo "Using cached $tmp_file_split.1 "
    else
	pinfo "Splitting sdrf into $batch_size long chunks"
	if [ ! -e $tmp_file_lock ]; then
	    touch $tmp_file_lock
	    tail -n +2 $SDRF_FILE_FP|cut -f 1 |sort -u > ${tmp_file_pref}.libs
	    echo "Libs=$(wc -l ${tmp_file_pref}.libs)"
	    split -a 4 --lines $batch_size  --numeric-suffixes=1  ${tmp_file_pref}.libs $tmp_file_split.
	    rename split.$batch_size.0 split.$batch_size. $tmp_file_split.0*
	    rename split.$batch_size.0 split.$batch_size. $tmp_file_split.0*
	    rename split.$batch_size.0 split.$batch_size. $tmp_file_split.0*
	    rm -f $tmp_file_lock
	else
	    echo "ERROR: unable to proceed - please delete lock file $tmp_file_lock"
	    exit
	fi
    fi
    set -e
    if [ "$batch_number-" == "0-" ]; then
	set +e
	all_batches=$(shopt -s nullglob && echo $tmp_file_split.*{0,1,2,3,4,5,6,7,8,9}*|sed "s|$tmp_file_split.||g"|tr ' ' '\n' |grep -v .txt|sort -n -u)
	set -e
	echo $all_batches
	for batch_num in $all_batches; do
	    BCONF_FILE_PREF=$ROOT_DIR/$SPECIES/$ID/$ID.B$batch_num.$batch_size 
	    BCONF_FILE=$BCONF_FILE_PREF.conf
	    echo $BCONF_FILE
	    echo $SDRF_FILE_FP
	    if [ -e $BCONF_FILE ] && [ $BCONF_FILE -nt $SDRF_FILE_FP ]; then
		echo Skipping generation of $BCONF_FILE
	    else
		set +e
		let jobs_run=0
		let jobs_run=$(jobs -r | wc -l)
		set -e		
		if [ $jobs_run -gt $THREADS ]; then
		    echo -n "waiting for slot..."
		    builtin wait -n
		    ret=$?
		    if [ "$ret-" != "0-" ] && [ "$ret-" != "127-" ] ; then
			echo "Failed to process batch."
			exit 1
		    fi
		    echo "done."
		fi
		irap_AE_setup.sh $* -b $batch_num -d $DATA_DIR > $BCONF_FILE.log &
	    fi
	done
	builtin wait
	if [ "$?-" != "0-" ]; then
	    echo "Failed to process batch."
	    exit 1
	fi
	# finally...generate a single conf file
	irap_AE_setup.sh $* -b 0 -B 0 -d $DATA_DIR
	exit
    fi
    ## split the sdrf
    splitted_sdrf=$tmp_file_split.$batch_number.sdrf.txt
    if [ -e $splitted_sdrf ] && [ $splitted_sdrf -nt $SDRF_FILE_FP ]; then
	echo Using cached  $splitted_sdrf
    else
	echo "Generating $splitted_sdrf"
	head -n 1 $SDRF_FILE_FP > $splitted_sdrf.tmp
	grep -F -w -f $tmp_file_split.$batch_number $SDRF_FILE_FP>> $splitted_sdrf.tmp
	mv $splitted_sdrf.tmp $splitted_sdrf
    fi
    SDRF_FILE=$splitted_sdrf
    SDRF_FILE_FP=$(readlink -f $splitted_sdrf)
    CONF_FILE_PREF=$TOP_FOLDER/$ID.B$batch_number.$batch_size
    CONF_FILE=$CONF_FILE_PREF.conf
    echo habemus_batches
    pinfo "CONF_FILE=$CONF_FILE"
    pinfo "SDRF_FILE=$SDRF_FILE_FP"
fi
####################################
# IDF
expected_clusters_params=""
if [ $skip_idf != "y" ]; then
    pinfo "Processing IDF..."
    set +e
    is_single_cell=$(grep AEExperimentType $IDF_FILE_FP|grep -c "single cells")
    if [ $is_single_cell == 1 ]; then
	pinfo "single cell RNA-seq"
	sc_params="--sop atlas_sc --sc"
	irap_cmd=irap_sc
	## expected number of clusters
	if [ $(grep -c EAExpectedClusters $IDF_FILE_FP) -eq 0 ]; then
	    echo ERROR: EAExpectedClusters not found in IDF > /dev/stderr
	    exit 1
	fi
	expected_clusters=$(grep EAExpectedClusters $IDF_FILE_FP|cut -f 2)
	if [ "$expected_clusters-" != "-" ]; then
	    expected_clusters_params="--nc $expected_clusters"
	else
	    expected_clusters_params=""
	fi
    else
	## assume, by default, that is bulk RNA-seq
	pinfo "bulk RNA-seq"
	sc_params=
    fi
    pinfo "Processing IDF...done"
    set -e
fi
####################################
# SDRF
SOURCE_NAME_COL=$(get_column $SDRF_FILE_FP "SOURCE NAME" 1)

set +e
#ENA_SAMPLE for tech replicates
ENA_SAMPLE_COL=$(get_column $SDRF_FILE_FP "[ENA_SAMPLE]")
pinfo "ENA_SAMPLE_COL: $ENA_SAMPLE_COL"

TECH_REPL_COL=$(get_column $SDRF_FILE_FP "[technical replicate group]")
pinfo "TECH_REPL_COL: $TECH_REPL_COL"


ENA_RUN_COL=$(get_column $SDRF_FILE_FP "[ENA_RUN]")
pinfo "ENA_RUN_COL: $ENA_RUN_COL"
set -e
if [ "$ENA_RUN_COL-" == "-" ]; then
    ENA_RUN_COL=$(get_column $SDRF_FILE_FP "[RUN]" 1)
    pinfo "RUN_COL: $ENA_RUN_COL"
fi

STRAND_COL=$(get_column $SDRF_FILE_FP "[LIBRARY_STRAND]")
pinfo "LIBRARY_STRAND_COL: $STRAND_COL"

LAYOUT_COL=$(get_column $SDRF_FILE_FP "[LIBRARY_LAYOUT]" 1)
pinfo "LIBRARY_LAYOUT_COL: $LAYOUT_COL"

LAYOUT_COL=$(get_column $SDRF_FILE_FP "[LIBRARY_LAYOUT]" 1)
pinfo "LIBRARY_LAYOUT_COL: $LAYOUT_COL"

# optional
# [single cell quality]
#SC_QUAL_COL=$(get_column $SDRF_FILE_FP "[single cell quality]" 1)
#pinfo "single cell quality: $SC_QUAL_COL"

# [end bias]
#END_BIAS_COL=$(get_column $SDRF_FILE_FP "[end bias]" 1)
#pinfo "END_BIAS_COL: $END_BIAS_COL"

# not_applicable => ignore

#[LIBRARY_STRATEGY]==RNA_SEQ
LIB_STRATEGY_COL=$(get_column $SDRF_FILE_FP "[LIBRARY_STRATEGY]" 1)
pinfo "LIBRARY_STRATEGY_COL: $LIB_STRATEGY_COL"
#[LIBRARY_SOURCE]==TRANSCRIPTOMIC
LIBRARY_SOURCE_COL=$(get_column $SDRF_FILE_FP "[LIBRARY_SOURCE]" 1)
pinfo "LIBRARY_SOURCE_COL: $LIBRARY_SOURCE_COL"
#[LIBRARY_SELECTION]==other => ignore otherwise fail

# quickly validate the sdrf before continuing
set -e
irap_sdrf2conf --name $ID --sdrf $SDRF_FILE_FP --out_conf $CONF_FILE_PREF --species=$SPECIES --data_dir=$DATA_DIR --raw_dir=$SPECIES/$ID/fastq $sc_params $expected_clusters_params -c --atlas
set +e


# get unique runs
RUNS=$(cut -f $ENA_RUN_COL $SDRF_FILE_FP|tail -n +2|sort -u)

################################################
# create folders
FASTQ_FOLDER=$DATA_DIR/raw_data/$SPECIES/$ID/fastq
create_dir $FASTQ_FOLDER
pushd $FASTQ_FOLDER

################################################
# Download fastq files
# by default assumes that data is in ENA
#FTP_URL_PREF=ftp.sra.ebi.ac.uk/vol1/fastq/
FTP_DIR_PREF=

ALT_NAME_COL=`get_column $SDRF_FILE_FP "[SUBMITTED_FILE_NAME]" ''`
echo "SUBMITTED_FILE_NAME=$ALT_NAME_COL"
# files to download
FASTQ_COL=
FASTQ_COL=`get_column $SDRF_FILE_FP "[FASTQ_URI]" ''`
echo "FASTQ_URI=$FASTQ_COL"
if [ "$FASTQ_COL-" == "-" ]; then
    ##Comment[ArrayExpress FTP file]
    echo "FASTQ_URI not found...Looking for alternative column"    
    FASTQ_COL=`get_column $SDRF_FILE_FP '[ArrayExpress FTP file]' ''`
    if [ "$FASTQ_COL-" == "-" ]; then
	perror "FASTQ_URI and ArrayExpress FTP file columns not found"
	exit 1
    fi
    echo "FASTQ_LOCATION=$FASTQ_COL"
fi

FASTQ_FILES=$(cut -f $FASTQ_COL $SDRF_FILE_FP | tail -n +2 )
N_FQ=$(echo $FASTQ_FILES|wc -w)

# TODO: handle the cases where the fastq_uri  may point to non-fastq files (bam/cram) or be empty
function DOWNLOAD_PUBLIC {
    local url=$1
    local ofile=$2
    local download_program=$3
    rm -f $ofile.tmp
    ##echo "<<<<<<<<<<<<<<<<<"
    ##echo wget  -nv -c $url -O $ofile.tmp
    if [ "-$download_program" != "-" ]; then
	echo "Running: $download_program $url $ofile.tmp"
	$download_program $url $ofile.tmp && mv $ofile.tmp $ofile
    else
	wget  -nv -c $url -O $ofile.tmp && mv $ofile.tmp $ofile
    fi
    [ -e $ofile ]
}
function DOWNLOAD_PRIVATE {
    url=$1
    ofile=$2
    id=$3
    sdrf_file=$4
    ## submitted file
    if [ "$ALT_NAME_COL-" != "-" ]; then
	set -e
	sub_file=`grep "$url" $sdrf_file|cut -f $ALT_NAME_COL`
	echo $sub_file
	echo `pwd`
	chmod 777 .
	rm -f $ofile.tmp
	set +e
	sudo -u fg_cur scp sra-login-2:/fire/staging/aexpress/$id-*/$sub_file $ofile.tmp
	sudo -u fg_cur chmod 766 $ofile.tmp
	mv $ofile.tmp $ofile
    fi
    [ -e $ofile ]
}

function DOWNLOAD_PRIVATE2 {
    local url=$1
    local ofile=$2
    local id=$3
    local sdrf_file=$4
    ## submitted file
    set +e
    rm -f $ofile.tmp
    FOLDER_AT_EBI="/ebi/ftp/pub/databases/microarray/data/experiment/MTAB"
    file_loc="$FOLDER_AT_EBI/$id/$ofile"
    if [ -e $file_loc ]; then
	echo "File found @ $file_loc"
	cp -uva $file_loc $ofile.tmp
	mv $ofile.tmp $ofile
    else
	echo "file $file_loc not found"
    fi
    [ -e $ofile ]
}

function download {
    file=$1
    fn=$2
    local download_program=$3
    ##
    if [ -e $local_download_folder ] && [ "-$local_download_folder" != "-" ]; then
	pinfo "Skipping download - looking for $fn in $local_download_folder"
	if [ -e $local_download_folder/$file ]; then
	    ## rename file if necessary
	    ## note: the subfolder, if used, is based on the $file name
	    fn=$(fix_filename $fn)
	    # avoid duplicating
	    if [ ! -h $(readlink -f $fn) ]; then
		ln -s $(readlink -f  $local_download_folder/$file) $(readlink -f $fn)
	    fi
	else
	    echo "File $local_download_folder/$file  not found"
	    exit 1
	fi
    else
	## 
	pinfo "downloading $fn"
	if [ "$dl_fun" == "DOWNLOAD_PUBLIC" ]; then
	    set +e
	    DOWNLOAD_PUBLIC $file $fn $download_program
	    if [ $? -ne 0 ]; then
	    echo "Failed to download from ENA  $file ..." > /dev/stderr
	    dl_fun="DOWNLOAD_PRIVATE"
	fi
	    set -e
	fi
	
	if [ "$dl_fun" == "DOWNLOAD_PRIVATE" ]; then
	    stdbuf -o0 echo -n
	    set +e
	    DOWNLOAD_PRIVATE $file $fn $ID $SDRF_FILE_FP
	    if [ $? -ne 0 ]; then
		echo "Failed to download from private location $file ..." > /dev/stderr
		stdbuf -o0 echo -n
		DOWNLOAD_PRIVATE2 $file $fn $ID $SDRF_FILE_FP
		if [ $? -ne 0 ]; then
		    echo "Failed to download from private location2 $file ...given up" > /dev/stderr
		    exit 1
		fi
	    fi
	    set -e		
	fi
    fi
}

if [ $skip_file_checking == 1 ]; then
    pinfo Skipping downloading and checking $N_FQ fastq files
else
    dl_fun=DOWNLOAD_PUBLIC
    dest_dir=.
    pinfo Downloading $N_FQ fastq files
    for file in $FASTQ_FILES; do
	fn=$(basename $file)
	if [ $distribute_by_subfolders == 1 ]; then
	    ## avoid relying on the download path since
	    ## not all files come from ENA
	    dest_dir=$(get_lib_folder $file)
	    mkdir -p $dest_dir
	fi
	fn=$dest_dir/$fn
	if ( ( [ -e $fn ] || [ -h $fn ] ) && [ $USE_CACHE -eq 1 ] ) || ( [ $DEBUG -eq 1 ] ); then
	pinfo "skipped downloading $fn"
	else
	    set +e	
	    let jobs_run=$(jobs -r | wc -l)
	    set -e
	    if [ $jobs_run -ge $THREADS ]; then
		echo -n "."
		builtin wait -n
		ret=$?
		if [ "$ret-" != "0-" ] && [ "$ret-" != "127-" ] ; then
		    echo "Failed to download."
		    exit 1
		fi
	    fi
	    if [ $AT_EBI -eq 1 ]; then
		echo file exists?
		echo creating symlink?
		echo TODO
	    else
		download $file $fn $download_program &
	    fi
	fi
    done
    builtin wait
    if [ "$?-" != "0-" ]; then
	echo "Failed to download."
	exit 1
    fi

fi



   
popd 2>/dev/null
## for each species...
##
extra_params=
if [ $distribute_by_subfolders == 1 ]; then
    extra_params=--subfolder
fi
echo irap_sdrf2conf --name $ID --sdrf $SDRF_FILE_FP --out_conf $CONF_FILE_PREF  --species=$SPECIES --data_dir=$DATA_DIR --raw_dir=$SPECIES/$ID/fastq --sop atlas_sc $sc_params $expected_clusters_params $extra_params --atlas

set -e
irap_sdrf2conf --name $ID --sdrf $SDRF_FILE_FP --out_conf $CONF_FILE_PREF  --species=$SPECIES --data_dir=$DATA_DIR --raw_dir=$SPECIES/$ID/fastq $sc_params $expected_clusters_params --sc --atlas $extra_params --atlas
## comment reference/gtf lines
## add an include to the species.conf file
pinfo "Adding last details to the conf. file"
sed -i "s/^reference/#reference/;s|^gtf_file.*|include $SPECIES_CONF_DIR/$SPECIES.conf|" $CONF_FILE
#sed -i "s/^reference/#reference/;s/.*_dir=.*fastq//;s|^gtf_file.*|include $SPECIES_CONF_DIR/$SPECIES.conf|" $CONF_FILE
## Move to the destination dir
##
DEST_DIR=$ROOT_DIR/$SPECIES/$ID
mkdir -p $DEST_DIR   
if [ -e $CONF_FILE_PREF.conf ]; then
    mv $CONF_FILE_PREF.* $DEST_DIR
fi
echo "Files placed in $DEST_DIR"
if [ $batch_number == 0 ]; then
    rm -rf $TOP_FOLDER
fi
exit 0
####################################


