#!/bin/env bash
# =========================================================
# Copyright 2012-2017,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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

command -v tagBam   >/dev/null 2>&1 || { echo "ERROR: tagBAM from bedtools2 is required but was not found.  Aborting." >&2; exit 1; }
##command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools is required but was not found.  Aborting." >&2; exit 1; }

##########################################################
#
function pinfo {
    echo "[INFO] $*" > /dev/stderr
}

function perror {
    echo "[ERROR] $*" > /dev/stderr
}

function usage() {
    echo "Usage: bam_annotate.sh  -b bam_file [ -e exon_bed -i intron_bed -g gene_bed -t transcript_bed ]";
}

function check_file {
    file=$1
    name=$2
    param_name=$3
    if [ ! -e $file ] ; then
	perror "$name file $file passed with $param_name option not found"
	exit 1
    fi
}
##########################################################
##
# -b bam_file
# -e exon_bed
# -i intron_bed
# -g gene_bed
# -t transcript_bed
while getopts "b:e:i:g:t:h"  Option
do
    case $Option in
# update/reinstall
        b ) bam_file=$OPTARG;;
	e ) exon_bed_file=$OPTARG;;
	i ) intron_bed_file=$OPTARG;;
	g ) gene_bed_file=$OPTARG;;
	t ) trans_bed_file=$OPTARG;;
        h ) usage; exit 0;;
    esac
done

if [ "$*-" == "-" ]; then
    usage
    exit 1
fi

## mandatory
if [ ! -e  $bam_file ]; then
    perror "BAM file not found $bam_file"
    exit 1
fi
pinfo BAM=$bam_file

# 
annot_files=
annot_tags=
# gene|transcript
feat_annot_files=
feat_annot_tags=
let num_ann_files=0
let num_feat_ann_files=0
if [ "$exon_bed_file-" != "-" ]; then
    let num_ann_files=num_ann_files+1
    annot_files="$annot_files $exon_bed_file"
    annot_tags="$annot_tags exonic"
    check_file $exon_bed_file "exon bed" "-e"
    pinfo exons=$exon_bed_file	
fi

if [ "$intron_bed_file-" != "-" ]; then
    let num_ann_files=num_ann_files+1
    annot_files="$annot_files $intron_bed_file"
    annot_tags="$annot_tags intronic"
    check_file $intron_bed_file "intron bed" "-i"
    pinfo intron=$intron_bed_file	
fi

## gene or transcript is mandatory
if [ "$gene_bed_file-" != "-" ]; then
    let num_feat_ann_files=num_feat_ann_files+1
    feat_annot_files=$gene_bed_file
    feat_annot_tags="GX"
    check_file $gene_bed_file "gene bed" "-g"
    pinfo genes=$gene_bed_file	
fi

if [ "$trans_bed_file-" != "-" ]; then
    let num_feat_ann_files=num_feat_ann_files+1
    feat_annot_files="$feat_annot_files $trans_bed_file"
    feat_annot_tags="$feat_annot_tags TX"
    check_file $gene_bed_file "transcript bed" "-t"
    pinfo transcripts=$trans_bed_file       
fi


if [ $num_feat_ann_files -eq 0  ]; then
    perror "-g or -t option need to be provided"
    exit 1
fi


# bam_annotate.sh -b test_files/bam/Brain60R3.se.hits.bam -e test_files/bed/chr19.fa.exons.bed -g test_files/bed/chr19.fa.genes.bed
# TODO: overlap using the  strand info (-s option)
if [ "$annot_files-" != "-" ]; then
    cmd1="tagBam -i $bam_file  -files  $annot_files   -tag  YB  -labels $annot_tags | "
    input_file="stdin"
else
    cmd1=""
    input_file=$bam_file
fi


if  [ $num_feat_ann_files -eq 1 ]; then
    cmd="$cmd1 tagBam -i $input_file  -names   -tag $feat_annot_tags  -files $feat_annot_files"
else
    cmd="$cmd1 tagBam -i $input_file  -names  -tag GX  -files $gene_bed_file | tagBam -i stdin  -names  -tag  TX  -files $trans_bed_file"
fi
echo $cmd > /dev/stderr
exec bash -c "$cmd"
#$cmd1 tagBam -i $input_file  -names   -tag $feat_annot_tags  -files $feat_annot_files"
exit 0
