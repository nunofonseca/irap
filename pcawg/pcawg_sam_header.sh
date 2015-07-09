#!/bin/bash


user_analysis_id=$1
file_prefix=$2
meta_data_file=$3
SOP_VERSION=1

if [  "-$user_analysis_id" = "-" ]; then
    echo "ERROR: analysis_id not given" > /dev/stderr
    exit 1
fi
if [  "-$file_prefix" = "-" ]; then
    echo "ERROR: file_prefix not given" > /dev/stderr
    exit 1
fi
if [ ! -e "$meta_data_file" ]; then
    echo "ERROR: $meta_data_file not found"  > /dev/stderr
   exit 1
fi

if [  "$file_prefix" == "icgc" ]; then                                            
     file_prefix=$user_analysis_id                                                 
fi   
if [ ! -e $meta_data_file ]; then
    echo "ERROR: file $meta_data_file not found"  > /dev/stderr
    exit 1
fi

# 
sample_metadata=`grep -F "$file_prefix" $meta_data_file|tr '\t' '|'`
if [ "-$sample_metadata" = "-" ]; then
    echo "ERROR: Unable to find sample $file_prefix in $meta_data_file"  > /dev/stderr
    exit 1
fi

#study   project_code    icgc_donor_id   icgc_specimen_id        icgc_sample_id  submitted_donor_id      submitted_specimen_id   submitted_sample_id     specimen_type   donor_id        specimen_id     sample_id       date    aliquot_id      fastq_file_tarball      filesize        md5_checksum    analysis_id     uri     center_name     platform        platform_model  lib_id  rg_label        sample_uuid     fastq_files
analysis_id_col=18
fastq_files_col=26
center_name_col=20
platform_col=21
platform_model_col=22
sample_id_col=11
submitter_sample_id_col=8
lib_id_col=23
read_group_label_col=24
#sample_uuid_col=25
# removed 
fastq_files_col=25

#echo $sample_metadata
for var in analysis_id fastq_files center_name platform platform_model sample_id submitter_sample_id  lib_id read_group_label fastq_files; do
    #echo $var
    col_var=${var}_col
    #echo ${!col_var}
    export $var="`echo $sample_metadata|cut -f ${!col_var} -d\|`"
done

# print the variables
for var in analysis_id fastq_files center_name platform platform_model sample_id submitter_sample_id lib_id read_group_label fastq_files ; do
    #echo $var=${!var}
    if [ "${!var}-" == "-" ]; then
	echo "ERROR: Unable to get value for $var"  > /dev/stderr
	exit 1
    fi
done
# analysis id should match
if [ "$analysis_id" != "$user_analysis_id" ]; then
    echo "ERROR: Analysis id mismatch - expected $user_analysis_id and got $analysis_id"  > /dev/stderr
    exit 1
fi

                                                                                                                                                                                             
# ICGC                                                                                                                                                               
if [  "$file_prefix" = "$user_analysis_id" ]; then 
    file_prefix=$fastq_files                                                      
fi 

#@PG already included
cat << EOF
@CO	SOP:$SOP_VERSION
@CO	submitter_sample_id:$submitter_sample_id
@RG	ID:$center_name:$analysis_id:$file_prefix	LB:RNA-Seq:$center_name:$lib_id	PL:$platform	PM:$platform_model	PU:$center_name:$read_group_label	SM:$sample_id
EOF


exit 0

# wrong pair 
../scripts/pcawg_sam_header.sh  5a372c2d-ce65-47c3-9bdc-54d2dbbf7213 111005_UNC12-SN629_0147_AB01MFABXX.2 metadata.tsv
../scripts/pcawg_sam_header.sh  74925915-0801-448a-abb3-7ccc8e8ab278  111005_UNC12-SN629_0147_AB01MFABXX.2 metadata.tsv
# example in the wiki: Specimen_ID is used as sample_id?
../scripts/pcawg_sam_header.sh fea6fa88-4ebb-4971-b614-8c1e2eaa24ae 120921_UNC16-SN851_0181_AD14JFACXX_GATCAG_L006  metadata.tsv
