# Shell functions used in iRAP scripts
STDERR=/dev/stderr
################################################################
# IO
function pinfo { 
    echo "[INFO: $*]"
}

function perror { 
    echo "[ERROR: $*]" >> $STDERR
}

function pwarning { 
    echo "[WARNING: $*]" >> $STDERR
}

###############################################################
# iRAP configuration file functions

# conf_get_rs conf_file
function conf_get_rs {
    conf_file=$1
    # lookup in the conf file
    # 
    d=`grep -E "^\s*.*_rs=" $conf_file | cut -f 2 -d\=`
    echo $d
}

# conf_get_data_dir conf_file [irap_cmdline_options]
function conf_get_data_dir {
    local conf=$1
    local irap_options=$2
    # lookup in the conf file
    # 
    d=`echo $irap_options|grep "data_dir="`
    if [ "$d-" != "-" ]; then
	d=`echo $irap_options|sed -E "s/.*\s?data_dir=([^\s]+).*/\1/"`
    else
	d=`grep -E "^\s*data_dir=" $conf | cut -f 2 -d\=`
	if [ "$d-" == "-" ]; then
	    # check if it is already defined
	    if [ "$data_dir-" == "-" ]; then
		perror "Unable to find iRAP's data directory"
		exit 1
	    fi
	fi
    fi
    echo $d
}

function conf_get_mapper {
    conf=$1
    cmdline_options=$2
    echo `conf_get_var_value mapper $conf "$cmdline_options"`
}

function conf_get_quant_method {
    conf=$1
    cmdline_options=$2
    echo `conf_get_var_value quant_method $conf "$cmdline_options"`
}

function conf_get_name {
    conf=$1
    cmdline_options=$2
    echo `conf_get_var_value name $conf "$cmdline_options"`
}

function conf_get_reference {
    conf=$1
    cmdline_options=$2
    echo `conf_get_var_value reference $conf "$cmdline_options"`
}

function conf_get_gtf {
    conf=$1
    cmdline_options=$2
    echo `conf_get_var_value gtf_file $conf "$cmdline_options"`
}

function conf_get_species {
    conf=$1
    cmdline_options=$2
    echo `conf_get_var_value species $conf "$cmdline_options"`
}

# 1 - conf file
# 2 - cmdline options
function conf_ok {
    local conf=$1
    local irap_options=$2
    echo  "irap conf=$conf $irap_options -n" > /dev/stderr
    set +e
    irap conf=$conf $irap_options -n 2>/dev/null    
    if [ $? == 2 ]; then
	irap conf=$conf $irap_options -n
	echo "ERROR: error in iRAP options." > /dev/stderr	
	exit 2
    fi
    set -e
}

# conf_get_var_value var conf_file [irap_cmdline_options]
function conf_get_var_value {
    local conf_var=$1
    local conf=$2
    local irap_options=$3
    # lookup in the conf file
    #
    echo  "irap conf=$conf $irap_options -n | grep $conf_var=" > /dev/stderr
    d=`irap conf=$conf $irap_options pe= se= -n 2>/dev/null|grep -E "*\s*$conf_var\s*="|tail -n 1`
    if [ "$d-" == "-" ]; then
	echo "ERROR: unable to get $conf_var value" > /dev/stderr	
	exit 1
    else
	d=`echo $d|cut -f 2 -d=`
    fi
    echo "debug:$conf_var=$d" > /dev/stderr
    echo $d
}



###############################################################

# get_fastq_fileprefix filename
# removes .(fq|fastq).* from the file name
function get_fastq_fileprefix {
    filename=`basename $1`
    filepref=`echo $filename|sed -E "s/.(fastq|fq).*//"`
    echo $filepref
}


# is_relative path - paths to the libraries should be relative although it will work
# with full paths
function is_relative_path {
    path=$1
    var_name=$2
    if [ `echo "$1"|cut -b 1` == "/" ]; then
	pwarning "Relative path expected but a full path was  provided ($var_name)."
    fi
}

# two levels - 
# L1: ~1.7K folders
# get_new_folder path
function get_new_folder {
    filename=`basename $1`
    filepref=`get_fastq_fileprefix $1`
    md5=`echo $filename|md5sum`
    fl=`echo $md5|cut -b 1,2`
    echo $fl/$filepref
}

##########################################################################
# Command execution

#
#
cmd_debug=0
function enable_cmd_debug {
    cmd_debug=1
}
function disable_cmd_debug {
    cmd_debug=0
}
function run_cmd {

    if [ $cmd_debug -eq 1 ]; then
	echo $*
    else 
	$*
    fi
}

##############################
# run and time cmds
declare -i classified_error=0

function run_AND_timeIt {
    label=$1
    logfile=$2
    shift 2
    datetime=`date "+%F %R"`
    d=`dirname $logfile`
    mkdir -p $d/logs
    sout=$d/logs/$label.out
    serr=$d/logs/$label.err
    #echo `pwd`
    echo "CMD: $*" 
    # Redirect stderr and stdout to a file
# `W'
#     Number of times the process was swapped out of main memory (KB)
#`I'
#     Number of file system inputs by the process.
#`O'
#     Number of file system outputs by the process.
    # label\tTime elapsed\tTime leapsed(Seconds)\taximum resident set size\tcommand\t exit status
    # label |Time elapsed |Time elapsed(Seconds)| maximum resident memory |date|command\t | exit status | ....
    /usr/bin/time -o $logfile -a --format "$label\t%E\t%e\t%M\t$datetime\t$*\t%x\t%I\t%O\t%W" bash -c "$*" 2> $serr 
    EXIT_STATUS=$?
    # output stderr
    #cat $sout >/dev/stdout
    cat $serr > /dev/stderr

    if [ $EXIT_STATUS -ne 0 ]; then
	classify_error $label $serr
	exit $EXIT_STATUS
    else
	rm -f $sout $serr
    fi
    #pinfo $EXIT_STATUS
    return $EXIT_STATUS
}

########################################################
# Try to classify errors on each stage
function classify_error {
    perror "LOG files: $2"    
    ${1}_errors $2
    #echo "!!!!!!!!!!!!!!$classified_error<<<<<<<<<<<<<<<<<<<<" > /dev/stderr
    print_classified_error
}

function print_classified_error {

    if [ "$classified_error" == "1" ]; then
	pinfo "Classified error...cleaning up before exiting."
	clean_up_all
    fi   
    perror $errmsg
}

function set_error {
    readonly errmsg="$*"
}

function set_classified_error {
    classified_error=1
    set_error $*
}

function io_error {
    errf=$1
    E=`grep -E "(disk I/O error|Stale|IOError:|write error|error 521|Cant't exec)" $errf`
    if [ $? -eq 0 ]; then
	echo 1
    else
	echo 0
    fi    
}

function stage0_errors {
    errf=$1
    is_io_error=`io_error $errf`
    #E=`grep -E "(disk I/O error|Stale|IOError:|write error|error 521)" $errf`
    if [ $is_io_error -eq 1 ]; then
	set_classified_error "iRAP stage0: I/O error"
    else
	set_error "iRAP Stage0: unclassified error"
    fi

}


function fastqInfo_errors {
    errf=$1

    is_io_error=`io_error $errf`
    if [ $is_io_error -eq 1 ]; then
	set_classified_error "FastqInfo: I/O error"
    else
	e=`grep -E "(ERROR:|error:|Error)" $errf|tail -n 1`
	# Error in the file
	if [ $? -eq 0 ]; then
	    msg=`echo "$e" | grep -E "(unable|duplicated|header|truncated|identifier|character|length|implemented|unpaired|encoding|invalid)"|tail  -n 1|cut -f 2- -d: `
	    if [ $? -eq 0 ]; then
		set_classified_error "FastqInfo: $msg"
	    else
		set_classified_error "FastqInfo: unclassified error - $e"
	    fi
	else 
	    set_error "FastqInfo: unclassified error - $e"
	fi
    fi
}

function iRAP-QC_errors {
    errf=$1

    is_io_error=`io_error $errf`
    if [ $is_io_error -eq 1 ]; then
	set_classified_error "QC: I/O error"
    else
	E=`grep "mv: cannot stat.*filter1.stats"  $errf `
	if [ $? -eq 0 ]; then
	    set_classified_error "QC: 1 too short reads or below quality threshold"
	else
	    E=`grep "mv: cannot stat.*filter2.stats"  $errf `
	    if [ $? -eq 0 ]; then
		set_classified_error "QC: 2 contamination"
	    else
		E=`grep -E "reads with at least one reported alignment: [0-9]* \(100.00%\)"  $errf `
		if [ $? -eq 0 ]; then
		    set_classified_error "QC: 2 contamination"
		else
		    E=`grep "mv: cannot stat.*filter3.stats"  $errf `
		    if [ $? -eq 0 ]; then
			set_classified_error "QC: 3 reads with uncalled bases discarded"
		    else
			E=`grep "mv: cannot stat.*filter4.stats"  $errf `
			if [ $? -eq 0 ]; then
			    set_classified_error "QC: 4 reads without mates"
			else
			    E=`grep "Aborted.* bowtie" $errf`
			    if [ $? -eq 0 ]; then
				set_classified_error "QC: IO error?"
			    else
				E=`grep "Premature End-Of-File" $errf`
				if [ $? -eq 0 ]; then
				    set_classified_error "QC: 1 no reads pass the quality threshold"
				else
				    E=`grep "gzip: .*.gz: unexpected end of file" $errf`
				    if [ $? -eq 0 ]; then
					set_classified_error "gzip: unexpected end of file"			
				    else
					E=`grep -E "(disk I/O error|Stale|IOError:)" $errf`
					if [ $? -eq 0 ]; then
					    set_classified_error "QC: I/O error"
					else
					    E=`grep "gzip: .*.gz: unexpected end of file" $errf`
					    if [ $? -eq 0 ]; then
						set_classified_error "gzip: unexpected end of file"			
					    else
						E=`grep -E "0 paired reads! are the headers ok?" $errf`
						if [ $? -eq 0 ]; then
						    set_classified_error "QC: 5 reads without mates"			
						else					
						    set_error "QC: unclassified error"
						fi
					    fi
					fi
				    fi
				fi
			    fi
			fi
		    fi
		fi
	    fi
	fi
    fi
}

function iRAP-Mapping_errors {
    errf=$1
    is_io_error=`io_error $errf`
    if [ $is_io_error -eq 1 ]; then
	set_classified_error "iRAP Mapping: I/O error"
    else
	e=`grep "Error occured when reading beginning of SAM/BAM file." $errf`
	if [ $? -eq 0 ]; then
	    set_classified_error "iRAP Mapping: no aligned reads in BAM"
	else
	    e=`grep -E "bowtie2-align died with signal .* (core dumped)" $errf`
	    if [ $? -eq 0 ]; then
		echo hostname=`hostname`
		set_classified_error "iRAP Mapping: bowtie2-align  crashed: `hostname`"
	    else
		e=`grep -E "Error running.*tophat_reports" $errf`
		if [ $? -eq 0 ]; then
		    set_classified_error "iRAP Mapping: internal tophat2 error"
		else
		    E=`grep -E "(disk I/O error|Stale|IOError:)" $errf`
		    if [ $? -eq 0 ]; then
			set_classified_error "iRAP Mapping: I/O error"
		    else
			set_error "iRAP Mapping: unclassified error"
		    fi
		fi
	    fi
	    
	fi
    fi
}

function iRAP-Mapping-QC_errors {
    errf=$1
    is_io_error=`io_error $errf`
    if [ $is_io_error -eq 1 ]; then
	set_classified_error "iRAP Mapping QC: I/O error"
    else
	#E=`grep -E "(disk I/O error|Stale|IOError:)" $errf`
	E=`grep -E "database is locked" $errf`
	if [ $? -eq 0 ]; then
	    set_classified_error "iRAP Mapping QC: sqlite - database is locked"
	else	
	    set_error "iRAP Mapping QC: unclassified error"
	fi
    fi
}

function iRAP-Quant_errors {
    errf=$1
    is_io_error=`io_error $errf`
    if [ $is_io_error -eq 1 ]; then
	set_classified_error "iRAP Quant: I/O error"
    else
	e=`grep "Error occured when reading beginning of SAM/BAM file." $errf`
	if [ $? -eq 0 ]; then
	    set_classified_error "iRAP Quant: no aligned reads in BAM?"
	else
	    E=`grep -E "(disk I/O error|Stale|IOError:)" $errf`
	    if [ $? -eq 0 ]; then
		set_classified_error "iRAP Quant: I/O error"
	    else	    
		set_error "iRAP Quant: unclassified error"
	    fi
	fi
    fi
}

function iRAP-CRAM_errors {
    errf=$1
    is_io_error=`io_error $errf`
    if [ $is_io_error -eq 1 ]; then
	set_classified_error "iRAP CRAM: I/O error"
    else
	E=`grep -E "(disk I/O error|Stale|IOError:)" $errf`
	if [ $? -eq 0 ]; then
	    set_classified_error "iRAP CRAM: I/O error"
	else
	    set_error "iRAP CRAM: unclassified error"
	fi
    fi
}

function iRAP-tidyup_errors {
    errf=$1
    set_error "iRAP Tidyup: unclassified error"
}
