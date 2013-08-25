#!/bin/bash
set -e 
function pinfo {
    echo "[INFO] $*"
}

function usage {
    echo "Usage: track_add.sh  -o data_dir -d data -t type -i input_file -l label -m metainfo [-h]";
    echo "                     -d data     : gff3|bam|bam2|cov|covn|wig(bed)";
    echo "                     -m metainfo : var=value2; var1=value2";
}


function not_empty {
    var=$1
    value=$2
    if [ "$value-" = "-" ] ; then
	echo "Missing value for $var"
	usage
	exit 1
    fi
}

function add_cov {
    outdir=$1  
    file=$2
    label=$3
    key=$4;#=4
    mkdir -p $outdir/raw/wig/    
    cp $file $outdir/raw/wig/`basename $file`
    file=raw/wig/`basename $file`
    metadata=`get_metadata`
    cat <<EOF | add-track-json.pl $outdir/trackList.json
{ 
         "label"         : "$label",
         "key"           : "$key",         
         "storeClass"    : "BigWig",
         $metadata
         "urlTemplate"   : "$file",
         "type"          : "JBrowse/View/Track/Wiggle/XYPlot",
         "variance_band" : false,
         "min_score"     : 0,
         "style": {
             "pos_color"         : "black",
             "neg_color"         : "black",
             "clip_marker_color" : "red",
             "height"            : 100
          }
}
EOF
}

function add_covd {
    outdir=$1  
    file=$2
    label=$3
    key=$4;#=4
    mkdir -p $outdir/raw/wig/        
    cp $file $outdir/raw/wig/`basename $file`
    file=raw/wig/`basename $file`
    metadata=`get_metadata`
#
    cat <<EOF  | add-track-json.pl $outdir/trackList.json
{
         "label"         : "$label",
         "key"           : "$label",
         "storeClass"    : "BigWig",
         "urlTemplate"   : "$file",
         "type"          : "JBrowse/View/Track/Wiggle/Density",
         "bicolor_pivot" : 0,
         $metadata
         "style" : {
             "pos_color": "black",
             "neg_color": "red",
             "height"   : 25
         }
}
EOF
}


function add_foldchange {
    outdir=$1  
    file=$2
    label=$3
    key=$4;#=4
    # makes more sense for normalized values
    expression_level=0
    mkdir -p $outdir/raw/wig/    
    cp $file $outdir/raw/wig/`basename $file`
    file=raw/wig/`basename $file`
    metadata=`get_metadata`
#
    cat <<EOF  | add-track-json.pl $outdir/trackList.json
{ 
         "label"         : "$label",
         "key"           : "$key",         
         "storeClass"    : "BigWig",
         "urlTemplate"   : "$file",
         "type"          : "JBrowse/View/Track/Wiggle/XYPlot",
         "variance_band" : true,
         "bicolor_pivot" : $expression_level,
         $metadata
         "style": {
             "pos_color"         : "orange",
             "neg_color"         : "yellow",
             "clip_marker_color" : "red",
             "height"            : 80
          }
}
EOF
}


function add_quant {
    outdir=$1  
    file=$2
    label=$3
    key=$4;#=4
    # makes more sense for normalized values
    expression_level=1
    mkdir -p $outdir/raw/wig/    
    cp $file $outdir/raw/wig/`basename $file`
    file=raw/wig/`basename $file`
    metadata=`get_metadata`
#
    cat <<EOF  | add-track-json.pl $outdir/trackList.json
{ 
         "label"         : "$label",
         "key"           : "$key",         
         "storeClass"    : "BigWig",
         "urlTemplate"   : "$file",
         "type"          : "JBrowse/View/Track/Wiggle/XYPlot",
         "variance_band" : true,
         "bicolor_pivot" : $expression_level,
         "min_score"     : 0,
         "max_score"     : 20,
         $metadata
         "style": {
             "pos_color"         : "purple",
             "neg_color"         : "green",
             "clip_marker_color" : "red",
             "height"            : 80
          }
}
EOF
    cat <<EOF | add-track-json.pl $outdir/trackList.json
{
         "label"         : "$label-D",
         "key"           : "$label-D",
         "storeClass"    : "BigWig",
         "urlTemplate"   : "$file",
         "type"          : "JBrowse/View/Track/Wiggle/Density",
         "bicolor_pivot" : $expression_level,
         $metadata,
         "style": {
             "height"   : 25,
             "pos_color": "purple",
             "neg_color": "green"
         }
}
EOF
}

# style (rna cov)
#             "pos_color": "purple",
#             "neg_color": "green"
#         "max_score" : 100,

# Not working :(
function add_rna_bedgraph {
    outdir=$1  
    file=$2
    label=$3
    key=$4;#=4
    mkdir -p $outdir/raw/wig/
    cp $file $outdir/raw/wig/
    file=raw/wig/`basename $file`
    cat <<EOF | add-track-json.pl $outdir/trackList.json
{ 
         "label"         : "$label",
         "key"           : "$key",         
         $metadata
         "storeClass"    : "BigWig",
         "urlTemplate"   : "$file",
         "type"          : "JBrowse/View/Track/Wiggle/XYPlot",
         "variance_band" : true,
         "min_score"     : 0,
         "max_score"     : 2000,
         "style": {
             "pos_color"         : "#FFA600",
             "neg_color"         : "#005EFF",
             "clip_marker_color" : "red",
             "height"            : 100
          }
}
EOF
}

function add_gff {
    outdir=$1
    gff=$2
    label=$3
    type=$4
    key=$5;#=4
    extra=
    if [ "$type" = "mRNA" ]; then	
	extra="    --getSubs  \
               --getSubfeatures \
               --subfeatureClasses '{\"CDS\": \"transcript-CDS\", \"exon\": \"transcript-exon\"}' \
               --urlTemplate 'http://www.ensembl.org/Multi/Search/Results?species=all;idx=;q={name}' \
               --arrowheadClass arrowhead"
	#echo $extra
    fi
    #    --getSubs  \
    #--getSubfeatures               
    if [ "$type" = "gene" ]; then	
	extra="\
               --getSubfeatures \
               --class transcript \
               --subfeatureClasses '{\"CDS\": \"transcript-CDS\", \"UTR\": \"transcript-UTR\"}' \
               --arrowheadClass arrowhead"
	type='gene'
    fi
    cmd="flatfile-to-json.pl \
    --out $out \
    --gff $gff \
    --type $type \
    --trackLabel $label  \
    --getLabel \
    --autocomplete all \
    --key $key \
    --clientConfig ' { $metadata } '\
	$extra \
	;"
    echo $cmd
    bash -c "$cmd"
}

#
function add_bam2_track {
    outdir=$1
    bam=$2
    label=$3
    bam-to-json.pl        \
	--out $outdir     \
	--bam $bam        \
	--tracklabel "$label"   \
	--key "$label"     \
	--cssClass alignment        \
	--clientConfig  '{  "featureCss": "background-color: #66F; height: 8px","histCss": "background-color: #88F"  }'
}

function add_bam_track {
    outdir=$1
    bam=$2
    label=$3
    mkdir -p $outdir/raw/bam/
    # full path name
    if [ ! -e $bam ]; then
	echo "ERROR: File $bam not found!" >&2
	exit 1
    fi
    # check if file is already in target directory
    # FILENAME should be the label (that should be unique)
    #FILENAME=`basename $bam`
    FILENAME=$label.bam
    TDIR=$outdir/raw/bam
    metadata=`get_metadata`
    if [ -e $TDIR/$FILENAME ] ; then
	echo "Warning: File $bam already in $TDIR. Skipping move." >&2
    else
	set -e    
	# full path
	bam=`readlink -f $bam`
	# 
	TDIRF=`readlink -f $TDIR`
	cp -a $bam $TDIR/$FILENAME
	cp -a $bam.bai $TDIR/$FILENAME.bai
	# for now copy the files
	# create a symbolic link
	#pushd `dirname $bam` > /dev/null
	#ln -s $TDIRF/$FILENAME .
	#ln -s $TDIRF/$FILENAME.bai .
	#popd > /dev/null
    fi
	    
    cat <<EOF | add-track-json.pl $outdir/trackList.json
      {
         "storeClass"  : "JBrowse/Store/SeqFeature/BAM",
         "urlTemplate" : "raw/bam/$FILENAME",
         $metadata
         "label"       : "$label",
         "type"        : "JBrowse/View/Track/Alignments",
         "key"         : "$label",
         "style" : {
            "className"      : "alignment",
            "arrowheadClass" : "arrowhead",
            "labelScale"     : 500
         }
}
EOF
}

function add_wig_track {
    outdir=$1
    bed=$2
    label=$3
   
    cmd="wig-to-json.pl  --wig $bed  \
	--out $outdir           \
	--tracklabel $label   \
        --key '$label'"
    echo $cmd
    $cmd
}

function get_metadata {
    
    if [ "-" == "-$metainfo" ] ; then
	    echo "";
    fi
    echo " \"metadata\" : { $metainfo }, " 
}

############################
OPTERR=0
data=
out=
type=
files=
label=
metainfo=
while getopts "m:d:o:t:i:l:h"  Option
do
    case $Option in
	d ) data=$OPTARG;;
	o ) out=$OPTARG;;
	t ) type=$OPTARG;;
	i ) files=$OPTARG;;
	l ) label=$OPTARG;;
	m ) metainfo=$OPTARG;;
        h ) usage; exit;;
    esac
done

#echo `get_metadata`
metainfo=`echo $metainfo|sed "s/;/,\n/g"`
not_empty "-d" $data
not_empty "-o" $out
not_empty "-i" $files
not_empty "-l" $label

if [ "$data" = "gff3" ] ; then
    not_empty "-t" $type
    add_gff $out $files $label $type $type
    exit 0
fi

if [ "$data" = "bam" ] ; then
    add_bam_track $out $files $label
    exit 0
fi

if [ "$data" = "bam2" ] ; then
    add_bam2_track $out $files $label
    exit 0
fi

if [ "$data" = "wig" ] ; then
    add_wig_track $out $files $label
    exit 0
fi

if [ "$data" = "rna_bedgraph" ] ; then
    add_rna_bedgraph $out $files $label $label
    exit 0
fi

if [ "$data" = "cov" ] ; then
    add_cov $out $files $label $label
    exit 0
fi

if [ "$data" = "covd" ] ; then
    add_covd $out $files $label $label
    exit 0
fi

if [ "$data" = "foldchange" ] ; then
    add_foldchange $out $files $label $label
    exit 0
fi

echo "Invalid value for -d. Valid values are gff3,bam,wig and rna_bedgraph"

exit 1
