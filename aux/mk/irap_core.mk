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
# =============================================================

# Main code

###############################################################
# disable deletion of temporary files
.SECONDARY: 
.ONESHELL:

SHELL=bash
################################################################################
# Auxiliary functions
################################################################################

# Information messages
define p_info=
$(info * $(shell date "+%d/%m/%Y %H:%M:%S") : $(1))
endef

# Error messages
define p_error=
$(info * $(shell date "+%d/%m/%Y %H:%M:%S") : ERROR: $(1)) && $(error Fatal error: $(1))
endef

# check if the parameter has a value - prints an error if not
ifdef verbose
define check_ok=
$(if $($(1)),$(call p_info, *	$(1)=$($(1))),$(p_error Missing $(1)))
endef
else
define check_ok=
$(if $($(1)),,$(p_error Missing $(1)))
endef
endif

empty:=
space:=$(empty) $(empty)#

# complain if a file does not exist and exit
file_exists=$(if  $(realpath $(1)),,$(call p_error,$(1) not found))
file_exists_ce=$(if  $(realpath $(1)),,$(call p_error,$(1) - $(2)))

must_exist=$(if  $(realpath $(1)),,$(1))

# 1- pat
define quiet_ls=
$(shell ls --color=never -1 $(1) 2>/dev/null)
endef


#  check if a variable  $(1) is defined - return the variable name if it is defined or empty otherwise
is_defined=$(if $(subst undefined,,$(origin $(1))),$(1),)

get_cmdline_params=$(sort $(foreach p,$(.VARIABLES),$(if $(filter-out $(origin $(p)),command line),,$(p)='$($(p))')))


###################################################################
# Variables
data_dir?=$(IRAP_DIR)/data

ifdef def_lib_dir
override def_lib_dir:=$(patsubst %/,%,$(def_lib_dir))/
endif
#

ifneq ($(origin raw_data_dir),undefined)
def_lib_dir:=
endif

#raw_data_dir?=
def_lib_dir?=

#$(info ---------$(def_lib_dir))
# ensure that the dir name starts and ends with a /
#check_libdir_ok=$(if $(call is_defined,$(1)_dir),$(patsubst %/,%,$($(1)_dir))/,$(def_lib_dir))
check_libdir_ok=$(if $(call is_defined,$(1)_dir),$(patsubst %/,%,$($(1)_dir))/,$(dir $(word  1,$($(1)))))

# *****************
# Stranded data
# return y if the strand info was provided
lib_strand_info=$(if $(call is_defined,$(1)_strand),y,)

# strand value
VALID_STRAND_OPTIONS=first second both
# return non empty if option is invalid
strand_value_ok=$(if $(filter $($(1)),$(VALID_STRAND_OPTIONS)),invalid option $($(1)) for $(1),)

# check if stranded information is ok
check_libstrand_ok=$(if $(call is_defined,$(1)_strand),$(if $(call strand_value_ok,$(1)_strand),$(call set_has_stranded_data,$(1)) $($(1)_strand),$(call p_error,Invalid strand option for lib $(1) - $($(1)_strand))))

# validate the read group id info
check_rgid_ok=$(call get_rgid,$(1))
get_rgid=$(if $(call is_defined,$(1)_rgid),$($(1)_rgid),)

# sam header
# add @rg if not defined
check_sheader_ok=$(if $(call is_defined,$(1)_shl),$($(1)_shl),$(if $($(1)_rgid),'@RG	ID:$($(1)_rgid)',))

has_stranded_data=no
set_has_stranded_data=$(if $(filter both,$($(1)_strand)),,$(eval has_stranded_data=yes))
##########################
# barcode values - lib, bc
check_bc_value_ok=$(if $(call is_defined,$(1)_$(2)),$($(1)_$(2)),$(def_$(bc)))

# the following parameters may be defined once - and are applied as default values
# for all libraries/assays - or defined per assay
def_umi_read?=undef
def_umi_offset?=0
def_umi_size?=0
def_cell_read?=undef
def_cell_offset?=0
def_cell_size?=0
def_sample_read?=undef
def_sample_offset?=0
def_sample_size?=0
def_sample_name?=
def_read1_offset?=0
def_read1_size?=-1
def_read2_offset?=0
def_read2_size?=-1
def_index1?=undef
def_index2?=undef
def_index3?=undef
def_known_umi_file?=
def_known_cells_file?=

barcode_min_qual?=10

#######################
# subfolder names
qc_label?=.
quant_label?=.
mapper_label?=.
de_label?=.

#######################
# only method for now
qc_method:=irap_qc

########################
# path functions - 1=lib
# get the path for a filtered fastq file
lib2filt_folder=$(qc_toplevel_folder)/$($(1)_dir)
# get the path for a bam file
lib2bam_folder=$(mapper_toplevel_folder)/$($(1)_dir)
#
lib2snp_folder=$(snp_toplevel_folder)/$($(1)_dir)
#
lib2quant_folder=$(quant_toplevel_folder)/$($(1)_dir)
#
lib2fusion_folder=$(fusion_toplevel_folder)/$($(1)_dir)

# check if a parameter is defined
check_param_ok_verbose=$(if $(call is_defined,$(strip $(1))),$(info *	$(1)=$($(1))),$(error * Missing -$(1)-))
check_param_ok=$(if $(call is_defined,$(strip $(1))),,$(error * Missing -$(1)-))

# check  contrast name
# 1 - param 2-param_name
check_name=$(if $(shell echo $(1)|sed "s/[0-9]*//g"),,$(error * Invalid $(2) $(1)))

# the library name should be different from the prefix of the fastq files
# ex. of an invalid library name: lib1=lib1.fastq
# 1=libname
#check_se_libname_ok=$(if $(subst $(1).fastq,,$($(1))),, $(warning * no longer invalid libname $(1)))
#check_pe_libname_ok=$(if $(strip $(subst $(1)_2.fastq,,$(subst $(1)_1.fastq,,$($(1))))), ,$(warning * no longer invalid libname $(1)))
check_se_libname_ok=
check_pe_libname_ok=

# 1 - lib
# 2 - bam
define convert_bam2sefastq=
$(call file_exists,$(2))
$(eval override $(1):=$(notdir $(basename $(2))).fastq)
endef
define convert_bam2pefastq=
$(call file_exists,$(2))
$(eval override $(1):=$(notdir $(basename $(2)))_1.fastq $(notdir $(basename $(2)))_2.fastq)
endef


# which program to use to merge the tsv files
# 1 = quant_method
# TODO: reduce memory footprint of irap_merge_tsv_NA.sh 
define merge_tsv=
$(if $(findstring cufflinks,$(1)),irap_merge_tsv_NA.sh,MAX_THREADS=$(max_threads) irap_merge_tsv_p.sh)
endef

################################################################################
# Variables
################################################################################

## default is bulk (blk)
rnaseq_type?=blk

## single_cell_protocol=none 
sc_protocol?=none

ifndef debug	
debug=0
p_debug=
else
override debug:=1
p_debug=$(info DEBUG:$(1))
endif

ifndef pe
 pe=
else
override pe:=$(sort $(pe))
endif

ifndef se
 se=
else
override se:=$(sort $(se))
endif


#************************
#Version and license info
pname=IRAP
version=1.0.0b
contact=Developed by Nuno Fonseca (authorname (at) acm.org)
license=This pipeline is distributed  under the terms of the GNU General Public License 3



###########
## Internal 
isl_mode?=n

share_files_in_reference?=y

################################################################################
# START!
$(info *****************************************************)
$(info * $(pname) $(version))
$(info * $(contact))
$(info * $(license))
$(info *)

$(info * Initializing...)


###############################################
# Load configuration (mandatory)
ifdef conf
 $(call file_exists,$(conf))
 $(info * Trying to load configuration file $(conf)...)
 include $(conf)
 $(info * Configuration loaded.)
else
 $(call p_error,Configuration file missing)
endif

# load library specific information
ifdef lib.info
lib_info=$(lib.info)
endif


ifdef lib_info
 $(call file_exists,$(lib_info))
 $(info * Trying to load information about the libraries - file $(lib_info)...)
 include $(lib_info)
 $(info * Information about the libraries loaded.)
endif
lib_info?=

get_lib_info_option=$(if $(lib_info),--lib_info $(lib_info),)

###############################################
# Load some definitions
include $(irap_path)/../aux/mk/irap_sc_defs.mk
include $(irap_path)/../aux/mk/irap_defs.mk

################################################################################
# Default values
################################################################################
transcript_de_method?=none
exon_de_method?=none

def_gse_tool?=none

# max. memory (in MB)
def_max_mem?=6000

#Contamination file
def_cont_index?=$(data_dir)/contamination/e_coli

#Default number of threads to run on a computer farm
def_max_threads?=1

#Minimal base quality accepted
def_min_read_quality?=10

# Trim poly-A/T? y|n
def_trim_poly_at?=n

# minimum poly-at length
# by default, if a read has at least 10 consecutive A or T in the edges then it will be trimmed. This option is only used if trim_poly_at is set to y
def_trim_poly_at_len?=10

# trim trailing N bases independtly of quality
# 0 = no trimming
trim_n_trailing_bases?=0

#Trim all reads to the minimum read size after quality trimming - y/n
def_trim_reads?=y

# none is being kept for backwards compatibility  (alias: qc)
# Quality filtering - on/off/report
def_qual_filtering?=on

# Mapper to use in QC contamination check 
def_cont_mapper?=bowtie2

#Software for mapping (reads -> genome/transcriptome)
def_mapper?=tophat2
def_mapper_splicing?=yes

#Software for differential expression
def_de_method?=none

#Is the experiment designed to compare the results of different mappers? yes|no => undef conditions|contrasts
def_mapper_comparison?=no

# default method to count reads mapped to features (genes, exons, ...)
quant_method?=htseq2

# default method to normalize the expression values
def_quant_norm_method?=none

# default tool to compute the normalized expression values
def_quant_norm_tool?=none

# count-specific normalisation method to be used by analysis tools that use
# normalised counts (rather than e.g. TPM)- e.g. SC3 in single-cell. May be the
# same as def_quant_method_norm.

def_quant_norm_count_method?=none
def_quant_norm_count_tool?=none

# produce quantification per exon? default is gene level or transcript level depending on the method used
def_exon_quant?=n
def_transcript_quant?=n

# maximum number of hits reported by the mapper
def_max_hits?=10

# fix NH flags in bam files
def_fix_NH?=y

# use only a subset of the genes in the DE analysis (y|n)
def_de_annot_genes_only?=n

# Htseq - produce a sam file with annotations?y/n
htseq_sam_output_def?=n

# Dominant transcript fold-change - no (n) to disables it
# only used if transcript quantification is y
dt_fc?=2

# Relative isoform usage?
# if df_fc is different than 0 then Riu files are generated
gen_riu?=n

# DE
def_de_pvalue_cutoff?=0.05
def_de_num_genes_per_table?=300

#def_annot_tsv=auto
def_annot_tsv?=off

CSS_FILE?=irap.css

gen_html_report?=n
GEN_REPORT_QC_ONLY=


###############################################################
# Check and validate the parameters values
# fastq files -> pair-end, single-end or both

$(info * )
$(info * Required Parameters:)

#***********************
# Name of the experiment
ifndef name
 $(call p_error, missing argument name!)
else
 $(info *	name=$(name))
endif

## Toplevel paths
auxdata_toplevel_folder:=$(name)/data
report_toplevel_folder:=$(name)

#************************
# data_dir Data directory (directory where the data is expected to be)
# Directory organization
# data_dir/reference/species/(fasta+(cdna+gtf))
# default: data_dir/raw_data/$(raw_folder)/

# check if data dir exists and structure is OK
# TODO
ifndef data_dir
 $(call p_error,missing argument data_dir)
endif


$(info *	data_dir=$(data_dir))
$(call file_exists,$(data_dir))


#********
# Species
ifndef species
 $(call p_error,missing argument species)
else
 $(info *	species=$(species))
endif

#****************
# raw data folder
ifndef raw_folder
 raw_folder=$(species)
endif

#***************
# raw_data_dir - folder where the raw data (fastq/bam) may be found
# 
ifdef raw_data_dir
ifneq ($(raw_data_dir)-,-)
# not empty
override raw_data_dir:=$(raw_data_dir)/
endif
endif
# always add a end /
raw_data_dir?=$(data_dir)/raw_data/$(raw_folder)/
$(info *	raw_data_dir=$(raw_data_dir))

#**********************
# Reference genome file
ifndef reference
 $(call p_error,missing argument reference)
endif

# 
reference_dir:=$(abspath $(data_dir)/reference/$(species))
reference_abspath:=$(abspath $(reference_dir)/$(subst .gz,,$(reference)))
reference_prefix:=$(reference_abspath)
reference_basename:=$(notdir $(reference_abspath))
# remove .gz if the file is gziped...the creation of the uncompress file is automatic

$(info *	reference=$(reference))
$(call file_exists,$(reference_dir)/$(reference))
oreference=$(reference)
#********* 
# GTF file
ifndef gtf_file
 gtf_file?=$(subst .fa,.gtf,$(reference))
endif

# remove .gz if the file is gziped...the creation of the uncompress file is automatic
gtf_file_dir:=$(abspath $(data_dir)/reference/$(species))
gtf_file_abspath:=$(abspath $(gtf_file_dir)/$(subst .gz,,$(gtf_file)))


$(info *       gtf_file  = $(gtf_file))
$(call file_exists,$(gtf_file_dir)/$(gtf_file))

DEXSEQ_GFF:=$(gtf_file_abspath).DEXSeq.gff
ogtf_file=$(gtf_file)
# irap's gtf file
# lgtf_file_dir:=


# cDNA file
ifndef cdna_file
 cdna_file?=$(subst .dna.fa,.cdna.all.fa,$(reference))
 ifeq ($(cdna_file),$(reference))
   cdna_file=$(subst .fa,.cdna.all.fa,$(reference))
 endif
endif
cdna_dir:=$(abspath $(data_dir)/reference/$(species))
cdna_file_abspath:=$(abspath $(cdna_dir)/$(cdna_file))
cdna_file_fasta_abspath:=$(abspath $(cdna_dir)/$(subst .gz,,$(cdna_file)))

##
user_trans_biotypes?=protein_coding|IG_([a-zA-Z0-9]+)_gene|TR_([a-zA-Z0-9]+)_gene
# the fasta file with the transcripts can be the cdna file or generated by iRAP (auto)
# default: backwards compatibility
user_trans?=cdna
ifeq ($(user_trans),cdna)
user_trans_abspath:=$(cdna_file_abspath)
else
## auto
user_trans_abspath:=$(patsubst %.gtf,%.trans.irap.fa,$(gtf_file_abspath))
endif
## $(info *  $(user_trans_abspath))

#$(call file_exists,$(trans_file))
ifdef spikein_fasta
ifeq ($(spikein_fasta),)
override undefine spikein_fasta
endif
endif
# ************
# spikein data
ifdef spikein_fasta
ifeq ($(spikein_fasta),ERCC)
override spikein_fasta:=$(data_dir)/ercc/ercc.fasta.gz
spikein_fasta_abspath:=$(abspath $(spikein_fasta))
else
spikein_fasta_abspath:=$(abspath $(spikein_fasta))
endif

endif

ifdef spikein_fasta
spikein_data=y
$(info *	spikein_data=y)
#spikein_fasta_abspath:=$(raw_data_dir)/$(spikein_fasta)
$(call file_exists,$(spikein_fasta_abspath))
# remove .gz
override spikein_fasta_abspath:=$(subst .gz,,$(spikein_fasta_abspath))
# gtf and fasta files
# new reference file cannot be shared across experiments
ifeq ($(share_files_in_reference),n)
reference_dir:=$(auxdata_toplevel_folder)
reference_abspath:=$(abspath $(reference_dir))
endif
user_reference_abspath:=$(reference_abspath)
# newref: ref_prefix.spikein_prefix.fasta
spikein_fasta_prefix:=$(patsubst %.fa,%,$(patsubst %.fasta,%,$(patsubst %.gz,%,$(notdir $(spikein_fasta)))))
new_spike_ref_prefix:=$(dir $(reference_abspath))/$(patsubst %.fasta,%.$(spikein_fasta_prefix),$(patsubst %.fa,%.fasta,$(subst .gz,,$(reference))))
override reference_abspath:=$(abspath $(new_spike_ref_prefix).fa)
reference_prefix=$(reference_abspath)

reference_basename=$(notdir $(reference_abspath))

## transcripts
override trans_abspath:=$(dir $(user_trans_abspath))/$(patsubst %.fasta,%.$(spikein_fasta_prefix),$(patsubst %.fa,%.fasta,$(subst .gz,,$(notdir $(user_trans_abspath))))).fa

user_gtf_abspath:=$(gtf_file_abspath)
ifeq ($(share_files_in_reference),n)
gtf_file_dir:=$(auxdata_toplevel_folder)
gtf_file_abspath:=$(abspath $(auxdata_toplevel_folder)/$(subst .fasta,,$(spikein_fasta_prefix)).$(subst .gz,,$(notdir $(user_gtf_abspath))))
else
gtf_file_abspath:=$(abspath $(dir $(reference_abspath))/$(subst .fasta,,$(spikein_fasta_prefix)).$(subst .gz,,$(notdir $(user_gtf_abspath))))
gtf_file_dir:=$(dir gtf_file_abspath)
endif

spikein_gtf_file:=$(patsubst %.fasta,%.gtf,$(spikein_fasta_abspath))
override gtf_file:=$(notdir $(gtf_file_abspath))


## concentration (TSV file)
ifndef spikein_concentration
spikein_concentration=
$(info * Warning:  spikein concentration not provided)
endif

#
$(info * Currently spikeins are not used while performing differential expression analysis)
endif

#trans_abspath?=$(cdna_file_abspath)
trans_abspath?=$(user_trans_abspath)
$(info *       Transcripts = $(trans_abspath))

gtf_file_basename:=$(notdir $(gtf_file_abspath))

# mapping (exon/transcript to gene) file obtained from the gtf file
trans2gene_mapping_file:=$(subst .gz,,$(subst .gtf,.mapping_trans.tsv,$(gtf_file_abspath)))
feat_mapping_files?=$(trans2gene_mapping_file) $(subst .gz,,$(subst .gtf,.mapping_exons.tsv,$(gtf_file_abspath)))

# ****************
# single cell
ifndef single_cell
single_cell=n
endif


# ****************
# file with the junctions (as used by tophat)
# offset=0
# includes 1 base of the exons
juncs_file_abspath:=$(subst .gz,,$(gtf_file_abspath)).juncs
juncs_file:=$(notdir $(juncs_file_abspath))



ifndef gff3_file
# GFF3 file obtained from the gtf file
gff3_file_abspath:=$(auxdata_toplevel_folder)/$(patsubst %.gtf,%.gff3,$(subst .gz,,$(gtf_file)))
else
gff3_file_abspath:=$(subst .gz,,$(gff3_file))
endif

gff3_file:=$(notdir $(gff3_file_abspath))
$(info *       gff3_file  = $(gff3_file_abspath))



#**********
refgeneclass_file=$(subst .gz,,$(gtf_file_abspath)).gene_class.txt
refgeneannot_file=$(subst .gtf,,$(subst .gz,,$(gtf_file_abspath))).gene_annot.tsv

#************
# FASTQ files
# SE reads
all_se_files=
all_pe_files=

# should be defined by the code that includes this file
ALL_TARGETS?=
# no need to check the se and pe parameters if
# the target is stage0 only
TARGETS=
ifdef MAKECMDGOALS
TARGETS=$(MAKECMDGOALS)
endif

# dry_run?
dry_run=$(findstring n,$(firstword $(MAKEFLAGS)))

# 
ifeq ($(TARGETS),all)
tr_validation=y
endif

# skip validation of SE or PE
ifeq ($(TARGETS),stage0)
skip_lib_validation=y
#tr_validation=y
endif

ifeq ($(TARGETS),stage0_fix)
skip_lib_validation=y
endif

ifeq ($(TARGETS),report)
override gen_html_report:=y
endif

ifeq ($(TARGETS),qc_report)
override gen_html_report:=y
endif

ifdef do_stage0_only
skip_lib_validation=y
endif

#
ifdef skip_lib_validation
$(info * Skipping validation of se and pe)
override se:=
override pe:=
else
$(info * Validating se and pe)
ifdef se
ifneq ($(se),)
 $(foreach l,$(se),$(call check_param_ok,$(l)))
 $(info *	se=$(se))
 all_se_files:=$(foreach l,$(se),$($(l)))
 #$(foreach l,$(se),$(info $(l)=$($(l)))) 
 # check if fastq file is in a different directory
 $(foreach l,$(se),$(eval $(l)_dir=$(call check_libdir_ok,$(l))))
 $(foreach l,$(se),$(eval $(l)_strand=$(call check_libstrand_ok,$(l)))) 
# read group id
 $(foreach l,$(se),$(eval $(l)_rgid=$(call check_rgid_ok,$(l)))) 
# sam header (may/should include @RG)
 $(foreach l,$(se),$(eval $(l)_shl=$(call check_sheader_ok,$(l)))) 
 #$(foreach l,$(se),$(info $(l)_dir=$($(l)_dir)))
 #$(foreach l,$(se),$(info $(l)_dir=$($(l)_dir)))
 #$(foreach l,$(se),$(info $(l)_dir=$(call check_libdir_ok,$(l))))
 $(foreach l,$(se),$(call check_se_libname_ok,$(l)))
 $(foreach l,$(se),$(call check_param_ok,$(l)_rs))
 $(foreach l,$(se),$(call check_param_ok,$(l)_qual))
 $(foreach l,$(se),$(foreach bc,known_umi_file known_cells_file index1 index2 index3 umi_read umi_offset umi_size cell_read cell_offset cell_size sample_read sample_offset sample_size read1_offset read2_offset read1_size read2_size sample_name,$(eval $(l)_$(bc)=$(call check_bc_value_ok,$(l),$(bc)))))
 ifile_given=1
endif
endif
# PE files (libraries)
# lib=A B C
# 
# A=""
ifdef pe
ifneq ($(pe),)
 $(foreach l,$(pe),$(call check_param_ok,$(l)))
 all_pe_files:=$(foreach l,$(pe),$($(l)))
 $(info *	pe=$(pe))
 # $(info * debug * fastq_files=$(fastq_files))
 # check the definition of 
 #  insert size
 #  sd
 #  read size
 #  read group id
 #  shl (sam header lines)
 # spike_fasta - single lib only
 # spike_concentration - single lib only
 # for each lib
 rs_list=
 set_rs_list=$(eval rs_list+= $($(1)))
 $(foreach l,$(pe),$(call check_pe_libname_ok,$(l)))
 $(foreach l,$(pe),$(eval $(l)_dir=$(call check_libdir_ok,$(l))))
 $(foreach l,$(pe),$(eval $(l)_strand=$(call check_libstrand_ok,$(l))))
# read group id
 $(foreach l,$(pe),$(eval $(l)_rgid=$(call check_rgid_ok,$(l)))) 
# sam header lines (may/should include @RG)
 $(foreach l,$(pe),$(eval $(l)_shl=$(call check_sheader_ok,$(l)))) 
 #$(foreach l,$(pe),$(info $(l)_dir=$($(l)_dir)))
 $(foreach l,$(pe),$(call check_param_ok,$(l)_sd))
 $(foreach l,$(pe),$(call check_param_ok,$(l)_ins))
 $(foreach l,$(pe),$(call check_param_ok,$(l)_rs))
 $(foreach l,$(pe),$(call check_param_ok,$(l)_qual))
 $(foreach l,$(pe),$(call set_rs_list,$(l)_rs))
 $(foreach l,$(pe),$(foreach bc,known_umi_file known_cells_file index1 index2 index3 umi_read umi_offset umi_size cell_read cell_offset cell_size sample_read sample_offset sample_size read1_offset read2_offset read1_size read2_size sample_name,$(eval $(l)_$(bc)=$(call check_bc_value_ok,$(l),$(bc)))))
 #$(foreach l,$(pe),$(call check_param_ok,$(l)_mp))
 ifile_given=1
endif
endif

all_fq_files:=$(all_pe_files) $(all_se_files)
nfqfiles:=$(words $(all_fq_files))
nfqfilesu:=$(words $(sort $(all_fq_files)))

ifndef skip_lib_validation
ifneq ($(nfqfiles),$(nfqfilesu))
$(call p_error,Some of the libraries provided in se and/or pe have the same filename)
endif
#$(call p_info,has_stranded_data=$(has_stranded_data))

ifndef ifile_given
$(info *	pe parameter or se parameter should be defined and non-empty)
endif
endif

endif

#***********
# Contrasts
#***********

# backward compatibility
#ifdef conditions
#endif


ifdef contrasts
 $(info *	contrasts=$(contrasts))
# check contrast/group names
 $(foreach l,$(contrasts),$(call check_name,$(l),contrast);$(foreach g,$($(l)),$(call check_name,$(g),group/condition)))
 $(foreach l,$(contrasts),$(info *      $(l)=$($(l)));$(foreach g,$($(l)),$(call check_param_ok_verbose,$(g))))
else
contrasts=
endif


# Contrasts
# Ex.
# contrasts=contrast1 contrast2
# contrast1=g1 g2
# g1=Lib1 Lib2
# g2=Lib3 Lib4

########
# groups
# list of groups defined (only used during html report generation)
# ex.
# groups=group1 group2 group3 
ifndef groups
$(info *	groups=  parameter not defined, this should be defined if you intend to generate an HTML report for the inferred gene/level/transcript quantification)
groups=
else
$(foreach g,$(groups),$(call check_param_ok_verbose,$(g)))
$(info *	groups=$(groups))
endif

## names of the factors
## by default onky 
factors?=
## like the groups
# factorA=...

#*********************
# Technical replicates
#*********************
ifdef technical_replicates_labels
ifeq ($(technical_replicates_labels),not applicable)
override technical_replicates_labels:=
override technical_replicates:=
endif
endif

# backwards compatibility
ifdef technical.replicates
technical_replicates?=$(technical.replicates)
endif

# ex. technical.replicates=SE1,SE2;PE1,PE3,PE4;SE3
# means that there are two groups of tech. replicates (separated by;), group 1 composed by SE1 and SE2  and group2 composed by PE1,PE3 and PE4, SE3 does not have tech. replicates
ifndef technical_replicates
technical_replicates=
 $(info *	technical_replicates=NONE)
else
 $(info *	technical_replicates=$(technical_replicates))
 # validate the libraries names
 comma=,
 quote="
 $(foreach  l,$(subst $(quote), ,$(subst $(comma), ,$(subst ;, ,$(technical_replicates)))),$(call check_param_ok,$(strip $(l))))
#"
## check - should be improved...
num_tr=$(words $(sort $(subst $(quote), ,$(subst $(comma), ,$(subst ;, ,$(technical_replicates))))))
num_libs=$(words $(sort $(se) $(pe)))

ifdef tr_validation
ifneq ($(num_tr),$(num_libs))
$(error technical_replicates should include all libs in se and pe: found $(num_tr) labels but expected $(num_libs))
endif
## technical replicate names 
ifndef technical_replicates_labels
$(error technical_replicates_labels should be defined when technical_replicates is defined)
endif
endif

##
## number of technical replicates labels needs to match the number of groups defined in technical.replicates
# technical.replicates.labels=G1;G2;G3
ng1=$(words $(sort $(subst ;, ,$(technical_replicates))))
ng2=$(words $(sort $(subst ;, ,$(subst $(space),_,$(technical_replicates_labels)))))
ifdef tr_validation
ifneq ($(ng1),$(ng2))
$(error technical_replicates inconsistent with technical_replicates_labels: found $(ng2) labels but expected $(ng1) labels in technical_replicates_labels)
endif
endif
$(info *	number of technical replicate groups/labels=$(ng1))
endif



#####################
# Other Optional parameters
#####################
$(info * )
$(info * Optional Parameters:)

## RNASeq type
SUPPORTED_RNASEQ_TYPES:=blk sc
ifeq (,$(filter $(rnaseq_type),$(SUPPORTED_RNASEQ_TYPES)))
$(call p_info,[ERROR] Invalid rnaseq_type - valid values are $(SUPPORTED_RNASEQ_TYPES))
$(error Invalid rnaseq_type)
endif
$(info *	rnaseq_type=$(rnaseq_type))

ifeq ($(rnaseq_type),sc)
## single cell protocol
SUPPORTED_SCP:=none smart-seq2 smart-seq smart drop-seq 10xV1 10xV2 10xV1a

# drop-seq 10x
ifeq (,$(filter $(sc_protocol),$(SUPPORTED_SCP)))
$(call p_info,[ERROR] Invalid sc_protocol - valid values are $(SUPPORTED_SCP))
$(error Invalid sc_protocol)
endif
$(info *	sc_protocol=$(sc_protocol))

endif

################################################################
## single cell specific parameters

CELL_TAG?=CR

## Simple cell filter
## Cells with less or equal than $(sc_non_zero_rows) are discarded
sc_non_zero_rows=1
## increasing this filter can make the intermidiate files smaller but there is no track of the cells lost

## parameters passed to bam_umi_count (umi_count quant. option)
bam_umi_count_params?=--min_reads 1 --multi_mapped --min_umis 1

#

## may be used by different methods
## currently only bam_umi_count takes advantage of the following options (to reduce memory usage)
## max. number of cell barcodes
sc_max_cells=800000
## max. number of features quantified (if protein coding only genes are considered then this value can be reduced)
sc_max_features=80000
## average number of features expected to be expressed per cell
sc_feat_cell=5000

## --uniq_mapped --multi_mapped
umis_params?=--cb_cutoff 2

## use sample barcodes... if available (y|n)
sc_use_sample_barcode?=n

########################
## single cell filtering
cell_filt_min_features?=15   # minimum number of features expressed as a percentage of the total number of features

cell_filt_max_ERCC?=30  # maximum percentage of expression that may be atributed to ERCC spike-ins

# pre-blacklisted cells
cell_filt_controls?=   # file with the a known list of cells that should not be used in downstream analysis

cell_filt_outliers?=y   # filter outliers based on the total number of counts/expr (y|n)
## Exclude outliers based on the 5*median absolute difference (like scater).
cell_outliers_mad=5

## minimum expression per feature
cell_filt_min_expression?=1

cell_filt_min_tot_expr?=1000 # minimum number of counts per cell

sc_quant_viz?=tsne

## tsne plot parameters
# filter genes based on the number of cells where they are expressed
tsne_min_cells?=1
# filter cells based on the number of genes expressed
tsne_min_genes?=1
# max perplexity to consider
tsne_max_perplexity?=35

########################
## Clustering
## Only used in sc mode
min_clusters?=2
max_clusters?=2
clustering_method:=sc3


# Optional two column tsv file with cell annotations 
# cells_samples_annotation=
#

# not used yet
# valid_clustering_methods=sc3 none

ifeq ($(sc_use_sample_barcode),y)

ifeq ($(quant_method),umis)
$(error umis does not use sample barcodes. sc_use_sample_barcode is set to $(sc_use_sample_barcode))
endif
else
bam_umi_count_params+= --ignore_sample
endif

# dge files will be in mtx (MM) format
ifeq ($(quant_method),umi_count)
expr_format?=mtx
expr_ext?=mtx.gz
else
ifeq ($(quant_method),umis)
expr_format?=mtx
expr_ext?=mtx
else
expr_format?=tsv
expr_ext?=tsv
endif
endif

#********
# Threads
#********
ifndef max_threads
 max_threads=$(def_max_threads)
endif
$(info *	max_threads=$(max_threads))

#********************
# Temporary directory
#********************
# 
ifndef tmp_dir
 tmp_dir=$(data_dir)/tmp
endif

#$(shell mkdir -p tmp_dir)
$(info *	tmp_dir=$(tmp_dir) (temporary directory))

#******************
# Quality filtering
#******************
# on|off|report=none(deprecated)
ifdef qc
qual_filtering:=$(qc)
endif

ifndef qual_filtering
 qual_filtering=$(def_qual_filtering)
endif

ifndef cont_mapper
cont_mapper=$(def_cont_mapper)
endif

# alias
ifndef qc
qc=$(qual_filtering)
endif

# on - perform filtering of the reads
# report - just collect QC stats
# off - skip filtering and collection of QC stats
qc_modes=on off report
ifeq ($(qc_modes),$(filter-out $(qc),$(qc_modes)))
$(call p_error,Invalid qc/qual_filtering value)
endif

$(info *	qc/filter=$(qc))


##############
# Fine tune QC

# enable disable first QC stage (poly-AT and quality based trimming and filtering)
# qc_stage1=y|n
qc_stage1?=y

# Trim poly-A/T? y|n
trim_poly_at?=$(def_trim_poly_at)
# minimum poly-at length
# by default, if a read has at least 10 consecutive A or T in the edges then it will be trimmed. This option is only used if trim_poly_at is set to y
trim_poly_at_len?=$(def_trim_poly_at_len)

ifndef trim_reads
 trim_reads=$(def_trim_reads)
endif

#*************
# Min. quality
#*************
ifdef min_read_quality
 ifeq ($(strip $(min_read_quality)),)
  undefine min_read_quality
 endif
endif

ifndef min_read_quality
 min_read_quality=$(def_min_read_quality)
endif


# Maximum (percentage) of uncalled bases acceptable in a read
# max_n=100 disables the filtering
max_n?=0

#*******************
# Contamination file

# disable contamination index
# cont_index=no
ifndef cont_index
 cont_index=$(def_cont_index)
endif

#**************
# Mapper to use
#**************
ifndef mapper
 mapper=$(def_mapper)
endif
$(info *	mapper=$(mapper))

# gems
SUPPORTED_MAPPERS=tophat1 tophat2 smalt gsnap soapsplice bwa1 bwa2 bowtie1 bowtie2 gem star osa mapsplice hisat2
ifeq ($(rnaseq_type),sc)
SUPPORTED_MAPPERS+= kallisto
ifeq ($(mapper),kallisto)
# aligns to the transcriptome
mapper_splicing=no
endif
endif

ifeq (,$(filter $(mapper),none $(SUPPORTED_MAPPERS)))
$(call p_info,[ERROR] Invalid mapper)
$(error $(mapper) not implemented)
endif

mapper:=$(strip $(mapper))
#***************************
# Quantification/Transcr. assembly program
#***************************

ifndef quant_norm_method
quant_norm_method=$(def_quant_norm_method)
endif

ifndef quant_norm_tool
quant_norm_tool=$(def_quant_norm_tool)
endif

ifndef quant_norm_count_method
quant_norm_count_method=$(def_quant_norm_count_method)
endif

ifndef quant_norm_count_tool
quant_norm_count_tool=$(def_quant_norm_count_tool)
endif

###########################
# Method to count reads mapped to features (genes, exons, ...)
###########################
ifndef quant_method
 quant_method:=$(def_quant_method)
endif

SUPPORTED_QUANT_METHODS=htseq1 htseq2 cufflinks1 cufflinks2 cufflinks1_nd cufflinks2_nd  flux_cap nurd stringtie stringtie_nd rsem kallisto salmon umi_count umis featurecounts

# methods that produce transcript level quantification by default
TRANS_QUANT_METHODS=flux_cap cufflinks1 cufflinks2 cufflinks1_nd cufflinks2_nd nurd stringtie stringtie_nd rsem kallisto salmon umi_count bitseq 
#umi_count kallisto_umi 
# umis?

#rsem isoem sailfish bitseq
ifeq (,$(filter $(quant_method),$(SUPPORTED_QUANT_METHODS) none))
$(call p_info,[ERROR] quant_method)
$(error $(quant_method) not supported)
endif

quant_method:=$(strip $(quant_method))
$(info *	quant_method=$(quant_method))

ifndef exon_quant
exon_quant=$(def_exon_quant)
else
exon_quant:=$(strip $(exon_quant))
endif

ifndef exon_quant_method
exon_quant_method=dexseq
endif


## use a method that supports transcript quantification
ifndef transcript_quant
transcript_quant=$(def_transcript_quant)
ifneq (,$(filter $(quant_method),$(TRANS_QUANT_METHODS)))
$(info * Enabling transcript_quant since  $(quant_method) supports it by default)
override transcript_quant:=y
endif
endif

ifneq ($(transcript_de_method),none)
ifeq ($(transcript_quant),n)
$(info * Enabling transcript_quant since DE transcript is expression is set to y)
override transcript_quant:=y
endif
endif
## a transcript quantification method may be used but we do not care
## about the transcript expression - just flag this to skip a few
## steps
transcript_expr?=n
ifeq ($(transcript_quant),y)
override transcript_expr=y
endif

$(info *	exon_quant=$(exon_quant))
ifeq ($(exon_quant),y)
$(info *	exon_quant_method=$(exon_quant_method))
endif
$(info *	transcript_quant=$(transcript_quant))
$(info *	transcript_expr=$(transcript_expr))

## Needed for transcript quantification
mapTrans2gene=$(auxdata_toplevel_folder)/$(gtf_file_basename).mapTrans2Gene.tsv

ifeq ($(transcript_quant),y)
ifneq ($(dt_fc),n)
override gen_riu:=n
endif
endif

# disable dt and riu since the code still does not support mtx files efficiently
ifeq ($(expr_ext),mtx)
override gen_riu:=n
override dt_fc:=n
endif

#####################
# Exon quantification
SUPPORTED_EXON_QUANT_METHODS=stringtie dexseq htseq1 htseq2 featurecounts


ifeq ($(strip $(exon_quant)),y)
ifneq (,$(filter $(quant_method),cufflinks1 cufflinks2 cufflinks2_nd cufflinks1_nd nurd stringtie stringtie_nd rsem kallisto salmon))
$(error Sorry... currently $(quant_method) cannot be used together with exon_quantification)
endif

# dexseq can only be used with methods that do not support exon quant
ifeq ($(exon_quant_method),dexseq) 
ifneq (,$(filter $(quant_method),cufflinks1 cufflinks2 cufflinks2_nd cufflinks1_nd nurd stringtie stringtie_nd rsem kallisto salmon))
$(error Sorry... currently DEXSeq cannot be used together with $(quant_method))
endif
else 
ifneq (,$(filter $(exon_quant_method),cufflinks1 cufflinks2 cufflinks2_nd cufflinks1_nd nurd rsem kallisto salmon))
$(error $(exon_quant_method) does not provide exon level quantification)
endif
#$(info Exon quantification requires that all exons in the gtf file have an exon_id attribute)
endif
#exon_quant_method=$(quant_method)
#$(info *	exon_quant_method=$(exon_quant_method))
endif

ifndef gtf_file_wexonid
gtf_file_wexonid=$(gtf_file_abspath)
#$(gtf_file_wexonid)=$(gtf_file_abspath).exon_id.gtf
endif



##################################
# method to normalize the counts
##################################

SUPPORTED_NORM_TOOLS=cufflinks1 cufflinks2 cufflinks1_nd cufflinks2_nd flux_cap  nurd stringtie irap kallisto scran
ifeq (,$(filter $(quant_norm_tool),none $(SUPPORTED_NORM_TOOLS)))
$(call p_error,quant_norm_tool '$(quant_norm_tool)' invalid)
endif

ifeq ($(quant_norm_method),rpkm)
$(call p_error,quant_norm_method Please use fpkm instead of RPKM)
endif

# Normalization methods
SUPPORTED_NORM_METHODS=fpkm uq-fpkm fpkm-uq deseq_nlib tpm scran_gene scran_spike

ifeq (,$(filter $(quant_norm_method),none $(SUPPORTED_NORM_METHODS)))
$(call p_error,quant_norm_method '$(quant_norm_method)' invalid)
endif

ifneq (,$(filter $(quant_method),cufflinks1 cufflinks2 cufflinks1_nd cufflinks2_nd nurd stringtie stringtie_nd flux_cap))
ifndef quant_norm_tool
quant_norm_tool=$(quant_method)
override quant_norm_method=fpkm
$(info * Enabling generation of table of FPKMs since $(quant_method) produces F/RPKMs by default)
endif
ifeq ($(quant_norm_tool),none)
quant_norm_tool=$(quant_method)
override quant_norm_method=fpkm
$(info * Enabling generation of table of FPKMs since $(quant_method) produces F/RPKMs by default)
endif
endif

$(info *	quant_norm_tool=$(quant_norm_tool))
$(info *	quant_norm_method=$(quant_norm_method))

ifeq ($(quant_norm_tool),irap)
## Sequence of biotypes, separated by comma, that define the set of genes that are used to compute the total number of reads (rpkm/fpkm/tpm)
## default: all biotypes
## e.g.
## quant_norm_mass_biotypes=protein_coding,lncRNA,pseudogenes
quant_norm_mass_biotypes?=
$(info *	quant_norm_mass_biotypes=$(quant_norm_mass_biotypes))
endif

irap_raw2metric_params?=
ifneq ($(irap_raw2metric_params),)
$(info *	irap_raw2metric_params=$(irap_raw2metric_params))
endif

# Method to normalize the counts to normalised counts 

SUPPORTED_NORM_COUNT_METHODS=deseq_nlib scran_gene scran_spike
ifeq (,$(filter $(quant_norm_count_method),none $(SUPPORTED_NORM_COUNT_METHODS)))
$(call p_error,quant_norm_count_method '$(quant_norm_count_method)' invalid)
endif

#******************************
# Use unspliced mapping: yes/no 
#******************************
# this option is going to be removed in the future
ifndef mapper_splicing
 mapper_splicing=$(def_mapper_splicing)
endif
$(info *	mapper_splicing=$(mapper_splicing))

#********************
# Gene level DE
#********************
ifndef de_method
ifdef compare
de_method=$(compare)
endif
endif

ifndef de_method
 de_method=$(def_de_method)
endif

SUPPORTED_DE_METHODS=cuffdiff1 cuffdiff2 cuffdiff1_nd cuffdiff2_nd deseq edger voom deseq2 ebseq
ifeq (,$(filter $(de_method),none $(SUPPORTED_DE_METHODS)))
$(call p_info,[ERROR] de_method)
$(error $(de_method) not supported)
endif

de_method:=$(strip $(de_method))
$(info *	de_method=$(de_method))

ifndef de_pvalue_cutoff=
 de_pvalue_cutoff=$(def_de_pvalue_cutoff)
endif
$(info *	de_pvalue_cutoff=$(de_pvalue_cutoff))

ifndef de_num_genes_per_table
de_num_genes_per_table=$(def_de_num_genes_per_table)
endif

# TSV file
# format (fields may be missing with NA)- (gene) ID is mandatory and should match the one given in the gtf file
# "ID","Name","locus","source","lname","GO","GOterm","KEGG"
ifndef annot_tsv
annot_tsv=$(def_annot_tsv)
#else
# TODO: check if file format is ok and file !=auto
endif

ifeq (auto,$(annot_tsv)) 
override annot_tsv:=$(auxdata_toplevel_folder)/$(reference_basename).gene.annot.tsv
endif


ifeq (off,$(annot_tsv)) 
override annot_tsv:=$(auxdata_toplevel_folder)/empty.gene.annot.tsv
endif


ifndef de_annot_genes_only
de_annot_genes_only=$(def_de_annot_genes_only)
endif

#$(info annot_tsv=$(annot_tsv))

# used by de_seq, edger, voom
gene_de_min_count?=0
ifdef de_min_count
## backward compatibility
gene_de_min_count?=$(de_min_count)
endif

transcript_de_min_count?=0
exon_de_min_count?=0

# 1 - feature
get_gtf_for_feature=$(if $(subst exon,,$(1)),$(gtf_file_wexonid),$(gtf_file_abspath))

#********************
# Transcript level DE
#********************
# by default transcript DE is disabled

# reuse some gene level quantification methods...not ideal
SUPPORTED_TRANSCRIPT_DE_METHODS=edger voom deseq2 cuffdiff1 cuffdiff2 ebseq

ifeq (,$(filter $(transcript_de_method),none $(SUPPORTED_TRANSCRIPT_DE_METHODS)))
$(call p_info,[ERROR] transcript_de_method)
$(error $(transcript_de_method) not supported)
endif

$(info *	transcript_de_method=$(transcript_de_method))

ifndef transcript_de_pvalue_cutoff=
 transcript_de_pvalue_cutoff=$(def_de_pvalue_cutoff)
endif
$(info *	transcript_de_pvalue_cutoff=$(transcript_de_pvalue_cutoff))

ifndef de_num_transcripts_per_table
de_num_transcripts_per_table=$(def_de_num_transcripts_per_table)
endif

ifndef de_annot_transcripts_only
de_annot_transcripts_only=$(def_de_annot_transcripts_only)
endif

transcript_de_min_count?=10


#********************************
# by default exon DE is disabled
SUPPORTED_EXON_DE_METHODS=dexseq

ifeq (,$(filter $(exon_de_method),none $(SUPPORTED_EXON_DE_METHODS)))
$(call p_info,[ERROR] exon_de_method)
$(error $(exon_de_method) not supported)
endif

$(info *	exon_de_method=$(exon_de_method))

ifndef exon_de_pvalue_cutoff=
exon_de_pvalue_cutoff=$(def_de_pvalue_cutoff)
endif
$(info *	exon_de_pvalue_cutoff=$(exon_de_pvalue_cutoff))

ifndef de_num_exons_per_table
de_num_exons_per_table=$(def_de_num_exons_per_table)
endif

ifndef de_annot_exons_only
de_annot_exons_only=$(def_de_annot_exons_only)
endif

exon_de_min_count?=10

###############################################
# isl enabled -> stage3 targets=stage4
ifeq ($(isl_mode),y)
$(info *	isl_mode enabled)
endif

#************
# Constraints
#************
ifdef compare_mappers
 ifdef conditions
  $(error Incompatible options: compare_mappers and conditions)
 endif
endif



ifeq ($(mapper_splicing),no)
 files_indexed=$(trans_abspath)
else
# word 1 = reference 
 files_indexed=$(reference_prefix) $(gtf_file_abspath)
endif

ifeq (y,$(de_annot_genes_only))
 ifeq ($(de_method),cuffdiff1)
  $(error Incompatible options: de_annot_genes_only & $(quant_method))
 endif
 ifeq ($(de_method),cuffdiff2)
  $(error Incompatible options: de_annot_genes_only & $(quant_method))
 endif
endif

##########
# other options

# alias
ifdef html
gen_html_report=$(html)
endif

# maximum number of mappings allowed for each read
ifndef max_hits
 max_hits=$(def_max_hits)
endif

# fix NH flag in BAM files (y,n)
ifndef fix_NH
 fix_NH=$(def_fix_NH)
endif

#
# Htseq - produce a sam file with annotations?y/n
ifndef htseq_sam_output
htseq_sam_output=$(htseq_sam_output_def)
endif

ifeq ($(htseq_sam_output),y)
htseq_sam_output=
endif
#
# mem in MB
ifndef max_mem
max_mem=$(def_max_mem)
endif

# max memory in GB
max_mem_gb:=$(shell expr $(max_mem) \/ 1000)

# set memory (in bytes) to 65% of the max. memory available
#ifndef SAMTOOLS_SORT_MEM
# SAMTOOLS_SORT_MEM:=`expr $(max_mem) \* 1000000 \* 65 \/ 100`
#endif


# samtools 1.x
ifndef SAMTOOLS_SORT_MEM
 SAMTOOLS_SORT_MEM:=$(shell bash -c "expr  $(max_mem_gb) \* 90 \/ 100 \/ $(max_threads) ")G
ifeq ($(SAMTOOLS_SORT_MEM),0G)
override SAMTOOLS_SORT_MEM:=1G
endif
endif


# convert Gb to Mb - smatools sort controls better the memory usage with M
# include an extra 1MB
MIN_MEM?=1
ifndef SAMTOOLS_SORT_MEM_MT
 SAMTOOLS_SORT_MEM_MT:=$(shell echo "($(max_mem_gb)-$(MIN_MEM))*0.95/$(max_threads)*1000000000+1200000" | bc )
#shell bash -c "expr \( $(max_mem_gb) - $(MIN_MEM) \) \* 50 \/ 100 \/ $(max_thread
endif


ifndef feat_length
feat_length=$(auxdata_toplevel_folder)/$(gtf_file_basename).lengths.Rdata
endif

ifndef exon_length
exon_length=$(auxdata_toplevel_folder)/$(gtf_file_basename).lengths.Rdata
endif


##
## barcodes (umi, cell, sample)
barcode_post_process_bam?=n
$(info *	barcode_post_process_bam=$(barcode_post_process_bam))


# 1 - in bam
# 2 - out bam
ifeq ($(barcode_post_process_bam),y)
## alignments to the genome
define do_post_process_bam_cmd=
bam_add_tags  --inbam $(1) --outbam - | bam_annotate.sh -b - -e $(auxdata_toplevel_folder)/$(reference_basename).exons.bed   -i $(auxdata_toplevel_folder)/$(reference_basename).introns.bed  -g $(auxdata_toplevel_folder)/$(reference_basename).genes.bed6 -t $(auxdata_toplevel_folder)/$(reference_basename).transcripts.bed6 > $(2).tmp && mv $(2).tmp $(2)
endef
## Alignments to the transcriptome - assumes that the chr/seq contains the transcript ID
define do_post_process_trans_bam_cmd=
bam_add_tags --tx_2_gx $(mapTrans2gene) --tx  --inbam $(1) --outbam - > $(2).tmp && mv $(2).tmp $(2)
endef

else
define do_post_process_bam_cmd=
true 
endef
define do_post_process_trans_bam_cmd=
true 
endef
endif

###############################################################################
# When the number of libraries is greater than BIG_LIM then pass the arguments
# to some scripts from stdin
ifndef BIG_LIM
BIG_LIM:=400
endif

ifeq ($(shell expr $(words $(se) $(pe)) \<  $(BIG_LIM)),0)
# too many libs to be able to pass them as an argument
$(call p_info, Big number of libraries mode (>$(BIG_LIM)))
BIG_EXP=400
# cmd=$1
# out_file=$2
# args=$3
ifeq ($(dry_run),n)
define pass_args_stdin=
cat $(2).in | $(1) -stdin && rm -f $(2).in
endef
else
define pass_args_stdin=
$(call args2file,$(2).in,$(3)) cat $(2).in | $(1) -stdin && rm -f $(2).in
endef
endif

# 1- infiles
# 2- outfile
define stdin_cat=
$(call args2file,$(2).in,$(1)) xargs -a $(2).in cat > $(2)
endef

else
BIG_EXP=0
#$(call p_info, Small number of libraries mode)
# cmd=$1
# out_file=$2
# args=$3
define pass_args_stdin=
$(1) $(3) 
endef

# $(1) input files
# 2- outfile
define stdin_cat=
cat $(1) > $(2)
endef

endif


# *****************************************************
# Paths
qc_toplevel_folder1:=$(subst /.,,$(name)/$(qc_label))
qc_toplevel_folder:=$(qc_toplevel_folder1)/$(qc_method)

mapper_toplevel_folder1:=$(subst /.,,$(qc_toplevel_folder)/$(mapper_label))
mapper_toplevel_folder:=$(mapper_toplevel_folder1)/$(mapper)

quant_toplevel_folder1:=$(subst /.,,$(mapper_toplevel_folder)/$(quant_label))
quant_toplevel_folder:=$(quant_toplevel_folder1)/$(quant_method)

de_toplevel_folder1:=$(subst /.,,$(quant_toplevel_folder)/$(de_label))
de_toplevel_folder:=$(de_toplevel_folder1)/$(de_method)

tde_toplevel_folder1:=$(de_toplevel_folder1)
tde_toplevel_folder:=$(tde_toplevel_folder1)/$(transcript_de_method)

ede_toplevel_folder1:=$(de_toplevel_folder1)
ede_toplevel_folder:=$(ede_toplevel_folder1)/$(exon_de_method)

################################################################################
# Make stuff
phony_targets=
silent_targets= 
precious_targets=

################################################################################
# AUXILIARY FUNCTIONS
################################################################################

# FALSE=empty

# return a string in the form of group1,group2,...
define groups2str=
$(call spaces2commas,$(call groupsnames))
endef

define remove_spaces=
$(subst $(space),,$(1))
endef

define spaces2commas=
$(shell echo $(1)|tr " " ",")
endef

#$(sort $(strip $(foreach l,$(groups),$(foreach g,$($(l)),$(g)))))
define groupsnames=
$(sort $(strip $(foreach l,$(groups),$(l) )))
endef

define factorsnames=
$(sort $(strip $(foreach l,$(factors),$(l) )))
endef

define factors2str=
$(call spaces2commas,$(call factorsnames))
endef

define factorsdef2str=
$(call remove_spaces,$(foreach g,$(call factorsnames),$(call spaces2commas,$(strip $($(g))));))
endef

# return groups definition
define groupsdef2str=
$(call remove_spaces,$(foreach g,$(call groupsnames),$(call spaces2commas,$(strip $($(g))));))
endef

define libname2ofastq=
$($(1))
endef

# 1 - lib
# 2 - pe
define get_fastq_prefix=
$(if $(filter $(2),pe),$(subst _2,,$(subst _1,,$(basename $(word 1,$($(1)))))),$(basename $(word 1,$($(1)))))
endef

# fix the libname by excluding _1 and _2 from PE files
define fix_libname=
$(if $(call valid_libname,$(1)),$(1),$(call get_lib_name,$(1)))
endef

define valid_libname=
$(if $(filter $(strip $(1)),$(pe) $(se)),y)
endef

# 1 - libname
# return y if is PE '' otherwise
define is_pe_lib=
$(if $(filter $(strip $(1)),$(pe)),y)
endef

# 1-libname
# path to the filtered fastq files 
define libname2fastq=
$(if $(call is_pe_lib,$(1)),$(call lib2filt_folder,$(1))$(1)_1.f.fastq.gz $(call lib2filt_folder,$(1))$(1)_2.f.fastq.gz,$(call lib2filt_folder,$(1))$(1).f.fastq.gz)
endef

# 1-libname
define bam_file_for_lib=
$(if $(call is_pe_lib,$(1)),$(1).pe.hits.bam,$(1).se.hits.bam)
endef


define get_lib_name=
$(if $(filter $(1),$(se)),$(1),$(if $(findstring $(patsubst %_1,%,$(1)),$(pe)),$(patsubst %_1,%,$(1)),$(if $(findstring $(patsubst %_2,%,$(1)),$(pe)),$(patsubst %_2,%,$(1)),$(1))))
endef

# 1 -lib
define path2lib_bam=
$(call lib2bam_folder,$(1))/$(call bam_file_for_lib,$(1))
endef

# 1-group
define  get_group_bam_files=
$(foreach l,$($(1)), $(call path2lib_bam,$(l)))
endef

# save all variables defined in Make (filename)
define vars2file=
$(file > $(1))
$(foreach ivar,$(.VARIABLES),$(if $(filter $(origin $(ivar)),'file' 'command line',$(file >> $(1),$(ivar)=$($(ivar))))))
endef
#$(call vars2file,/tmp/lixo)

# the bam file names of each group are concatened using ,
define  get_contrast_bam_files=
$(foreach c,$($(1)), $(shell echo $(call get_group_bam_files,$(c)) | sed "s/ /,/g"))
endef

# 
define get_de_annot=
$(if $(annot_tsv),--annotation $(annot_tsv),)
endef

# DE options
define get_de_annot_genes_only=
$(if $(subst y,,$(de_annot_genes_only)),,--annot-genes-only)
endef

#1
define filename2libname=
$(foreach l,$(se) $(pe),$(if $(filter $(l).,$(1)),$(l)))
endef

# 1 tsv filename
#  exon, gene or CDS?
# DE filename to Analysis level
define DEfilename2AL=
$(if $(findstring .genes_de,$(1)),gene,$(if $(findstring .exons_de,$(1)),exon,CDS))
endef

define DEfilepath2demethod=
$(shell basename `dirname $(1)`)
endef

# sam/bam cat
# avoid samtools cat if there is only a file
# $(1)=files
define samcat=
$(if $(strip $(word 2,$(1))),samtools cat -o - $(1),cat $(1))
endef

# two line empty file (ensure that size !=0 and # lines >1)
# 1 - filename
define empty_file=
echo Empty	Empty > $(1) && echo Empty	Empty  >> $(1)
endef

# 1 - filename
# return y if it is an empty file or '' otherwise
define is_empty_file=
$(if  $(realpath $(1)),$(if $(shell head -n 1 $(1) | grep   Empty),y))
endef

define zero_lines_file=
$(if $(realpath $(1)),$(if  $(shell wc -l $(1)|cut -f 1 -d\ |sed "s/^0$//",,no)),)
endef
#1-file
#2-args
define args2file=
 $(file > $(1),$(2))
 $(shell echo >> $(1))
endef
# add the new line in the end

wave_file?=
#1 - args
define wave_echo=
$(if $(wave_file),$(file > $(wave_file),$(1)),echo $(1))
endef
####################
###################
# 
# $(1) var name
cached_vars=

cached_vars_file=$(name)/cached_vars.mk
ifndef use_cached_vars
use_cached_vars=n
endif

define cached_var=
$(if $(call is_defined,$(1)), $(call p_debug,cache hit $(1)), $(call set_$(1)) $(call p_debug,cache miss $(1)) $(eval cached_vars+= $(1))) $($(1))
endef

ifeq (use_cached_vars,y)
$(call file_exists,$(cached_vars_file))	
include $(cached_vars_file)
$(info [INFO] Loaded cached variables.)
endif

#
define list_cached_vars=
	$(foreach var,$(cached_vars) cached_vars, echo $(var)=$(strip $(call cached_var,$(var)));)
endef

################################################################################
# Minimum set of targets for each stage
STAGE0_TARGETS=
STAGE1_TARGETS=
STAGE2_TARGETS=
STAGE3_TARGETS=
STAGE4_TARGETS=
STAGE5_TARGETS=

################################################################################
# Files produced at each stage
CLEAN_UP_TARGETS=
BOOTSTRAP_TARGETS=
# stage0
SETUP_DATA_FILES=
STAGE1_OUT_FILES=
STAGE2_OUT_FILES=
STAGE3_OUT_FILES=
STAGE4_OUT_FILES=
STAGE5_OUT_FILES=

################################################################################
# wave* are used with a job scheduler in a HPC (the current code -
# irap_lsf - is for LSF)

WAVEB_TARGETS?=
WAVE0_TARGETS?=
WAVE1_TARGETS?=
# BAM+basic stats
WAVE2_TARGETS?=
# $(STAGE3_TARGETS)
WAVE3_TARGETS?=
WAVE3_s_TARGETS?=
WAVE4_TARGETS?=
WAVE5_TARGETS?=
WAVE6_TARGETS?=

################################################################################
# STAGE3 library level targets
STAGE3_S_TARGETS?=
# STAGE3 library level output files
STAGE3_S_OFILES?=
STAGE1_S_TARGETS?=

##
#STAGE1_OUT_FILES+=$(foreach p,$(se),$(call lib2filt_folder,$(p))$(p).f.fastq.gz) $(foreach p,$(pe),$(call lib2filt_folder,$(p))$(p)_1.f.fastq.gz)

ifneq ($(mapper),none)

bam_files:=$(foreach p,$(pe), $(call lib2bam_folder,$(p))$(p).pe.hits.byname.bam) $(foreach s,$(se), $(call lib2bam_folder,$(s))$(s).se.hits.bam)
STAGE2BYNAME_OUT_FILES:=$(subst .hits.bam,.hits.byname.bam,$(bam_files))

else

bam_files=
STAGE2_OUT_FILES=
STAGE2BYNAME_OUT_FILES=

endif


ifeq (umi_count,$(quant_method))
WAVE3_p_TARGETS+=$(subst .hits.bam,.hits.bytag_$(CELL_TAG).bam,$(bam_files))
endif



################################################################################
# Load extra code
# shared code
# QC stats
include $(irap_path)/../aux/mk/irap_qc_stats.mk
# Mapping
include $(irap_path)/../aux/mk/irap_map.mk
# JBrowse
include $(irap_path)/../aux/mk/irap_jbrowse.mk
# Software used (versions and citations)
include $(irap_path)/../aux/mk/irap_citations.mk
# Quantification
include $(irap_path)/../aux/mk/irap_quant.mk
include $(irap_path)/../aux/mk/irap_sc_quant.mk
# Normalization
include $(irap_path)/../aux/mk/irap_norm.mk
# DE
include $(irap_path)/../aux/mk/irap_de.mk
# Gene set enrichment analysis
include $(irap_path)/../aux/mk/irap_gse.mk


# Fusion
include $(irap_path)/../aux/mk/irap_fusion.mk

# QC - single cell
include $(irap_path)/../aux/mk/irap_sc_qc.mk
include $(irap_path)/../aux/mk/irap_clustering.mk
include $(irap_path)/../aux/mk/irap_sc_vis.mk

# Junctions
include $(irap_path)/../aux/mk/irap_junction.mk


# HTML Reporting
include $(irap_path)/../aux/mk/irap_report.mk
#
include $(irap_path)/../aux/mk/irap_status.mk

ifdef irap_devel
$(call p_info,Loading code under development)
# include under development features
include $(irap_path)/../aux/mk/irap_snp_indel_calling.mk
endif

# Atlas (atlas specific stuff)
include $(irap_path)/../aux/mk/irap_atlas.mk

# Check if the options provided are valid
ifeq (invalid,$(shell irap_paths $(mapper) $(quant_method) $(quant_norm_tool) $(quant_norm_method) $(de_method) $(transcript_de_method) $(exon_de_method) $(gse_tool) $(has_stranded_data) $(rnaseq_type) $(sc_protocol)))
  $(error invalid combination mapper:$(mapper) -> quant_method:$(quant_method) -> quant_norm method:$(quant_norm_method) quant_norm_tool:$(quant_norm_tool) -> de_method:$(de_method) transcriptDE:$(transcript_de_method) exonDE:$(exon_de_method) rnaseq_type:$(rnaseq_type) sc_protocol:$(sc_protocol) for the given data)
endif

$(info *========================================================)
pprint_path=$(subst /./,/,$(1))
$(info * Files/folders: )
$(info ** Top level folder: $(name))
$(info ** Aux data: $(auxdata_toplevel_folder))
$(info ** QC for input files: $(call pprint_path,$(qc_toplevel_folder)))
ifneq ($(mapper),none)
$(info ** Alignment files: $(call pprint_path,$(mapper_toplevel_folder)))
endif
ifneq ($(quant_method),none)
$(info ** Quantification matrices: $(call pprint_path,$(quant_toplevel_folder)))
endif
ifneq ($(de_method),none)
$(info ** Gene level DE: $(call pprint_path,$(de_toplevel_folder)))
endif
ifneq ($(transcript_de_method),none)
$(info ** Transcript level DE: $(call pprint_path,$(tde_toplevel_folder)))
endif
ifneq ($(exon_de_method),none)
$(info ** Exon level DE: $(call pprint_path,$(ede_toplevel_folder)))
endif
$(info --------------------------------------------------------)
#################################################################################
ifneq ($(mapper),none)
index_files:=$(call $(mapper)_index_filenames,$(word 1,$(files_indexed)),$(word 1,$(files_indexed)))
else
index_files=
endif

#*********************
# print all variables
ifdef debug
ifeq ($(debug),1) 
$(info * DEBUG)
VARS2PRINT=reference_prefix gtf_file_abspath index_files files_indexed feat_mapping_files
$(foreach v,$(VARS2PRINT),$(info $v=$($v)))
endif
endif
$(call p_info,[DONE] Initialization)


###################
# Quality Filtering
###################
read_qual_filter_common_params=tmp_dir=$(tmp_dir)  threads=$(max_threads)  qual_filtering=$(qual_filtering)  min_qual=$(min_read_quality) trim=$(trim_reads) cont_index=$(cont_index) mapper=$(cont_mapper) max_n=$(max_n) max_mem=$(max_mem) poly_at_len=$(trim_poly_at_len) trim_poly_at=$(trim_poly_at) qc_stage1=$(qc_stage1) trim_n_trailing_bases=$(trim_n_trailing_bases)

# get a param value pair iff the value passed is not empty and not undef
# 1 - param name
# 2 - value
get_param_value_pair=$(if $(call not_empty_not_undef,$(2)),$(1)$(2),)
not_empty_not_undef=$(if $(subst undef,,$(1)),1,)

# 1 - lib
define get_opt_barcode_params=
barcode_min_qual=$(barcode_min_qual) \
$(call get_param_value_pair,umi_read=,$(subst read,,$($(1)_umi_read))) $(call get_param_value_pair,umi_size=,$($(1)_umi_size))  $(call get_param_value_pair,umi_offset=,$($(1)_umi_offset)) \
$(call get_param_value_pair,cell_read=,$(subst read,,$($(1)_cell_read))) $(call get_param_value_pair,cell_size=,$($(1)_cell_size))  $(call get_param_value_pair,cell_offset=,$($(1)_cell_offset)) \
$(call get_param_value_pair,sample_read=,$(subst read,,$($(1)_sample_read))) $(call get_param_value_pair,sample_size=,$($(1)_sample_size))  $(call get_param_value_pair,sample_offset=,$($(1)_sample_offset)) \
 $(call get_param_value_pair,known_umi_file=,$($(1)_known_umi_file)) \
 $(call get_param_value_pair,known_cells_file=,$($(1)_known_cells_file)) \
 $(call get_param_value_pair,index1=,$(notdir $($(1)_index1))) \
 $(call get_param_value_pair,index2=,$(notdir $($(1)_index2))) \
 $(call get_param_value_pair,index3=,$(notdir $($(1)_index3))) \
 $(call get_param_value_pair,read1_offset=,$($(1)_read1_offset))  $(call get_param_value_pair,read2_offset=,$($(1)_read2_offset)) \
 $(call get_param_value_pair,read1_size=,$($(1)_read1_size)) $(call get_param_value_pair,read2_size=,$($(1)_read2_size))
endef


define not_empty=
	$(if $(call file_exists,$(1),),$(call p_error,File $(1) not found))
	$(if $(shell wc -l $(1)|cut -f 1 -d\ ),,$(call p_error,Empty file $(1)))
endef

################################################################################
#
################################################################################
# Generic file extension rules

# generate the fastq from BAM
# $(data_dir)/raw_data/$(species)/%.fastq: $(data_dir)/raw_data/$(species)/%.bam
# 	bamToFastq -i $^ -fq $@.tmp  && mv $@.tmp $@
# # pair-end data
# # bam files need to be sorted by name
# $(data_dir)/raw_data/$(species)/%_1.fastq $(data_dir)/raw_data/$(species)/%_2.fastq: $(data_dir)/raw_data/$(species)/%.bam
# 	$(call p_info,"Converting bam to fastq...note that the bam file needs to be sorted by name.")
# 	bamToFastq -i $^ -fq $*_1.fastq.tmp -fq2 $*_2.fastq.tmp && mv $*_1.fastq.tmp $*_1.fastq && mv $*_2.fastq.tmp $*_2.fastq

# uncompress rules
%.fa: %.fa.gz
	gunzip -c $< > $@.tmp && mv $@.tmp $@

%.fasta: %.fasta.gz
	gunzip -c $< > $@.tmp && mv $@.tmp $@

%.gtf: %.gtf.gz
	gunzip -c $< > $@.tmp && mv $@.tmp $@

$(reference_prefix).fa: $(reference_prefix)
	ln -s `basename $<`  $@

ifeq ($(user_trans),auto)
$(user_trans_abspath): $(gtf_file_abspath) $(reference_abspath)
	irap_gtf_to_fasta -g $(gtf_file_abspath) --genome $(reference_abspath) --types '$(user_trans_biotypes)' --out $@.tmp && mv $@.tmp $@
endif

ifdef spikein_data
$(trans_abspath): $(subst .gz,,$(user_trans_abspath)) $(subst .gz,,$(spikein_fasta_abspath))
	cat $^ > $@.tmp && mv $@.tmp $@


$(reference_abspath): $(subst .gz,,$(user_reference_abspath)) $(subst .gz,,$(spikein_fasta_abspath))
	cat $^ > $@.tmp && mv $@.tmp $@

$(spikein_gtf_file): $(spikein_fasta_abspath) 
	spikein_fasta2gtf.pl $< > $@.tmp && mv $@.tmp $@

$(gtf_file_abspath): $(user_gtf_abspath) $(spikein_gtf_file)
	cat $^ > $@.tmp && mv $@.tmp $@
endif


# give an error (should never happen)
#%.gtf:
#	$(call p_error, Missing gtf file $@)
# nf
%.gtf.irap.stats: %.gtf
%.gtf.irap.stats: %.gtf.gz



%.gtf.checked: %.gtf
	irap_check_gtf $(gtf_file_abspath) &&  touch $@

%.gtf.bed: %.gtf
	gtf2bed.pl $< > $@

%.bam.bedGraph: %.bam $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt
	bedtools  genomecov -ibam $< -bg -g $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt > $@.tmp && mv $@.tmp $@

# gtf file with exon_id attribute
%.gtf.exon_id.gtf: %.gtf %.gtf.checked
	gtf_add_exon_id.pl $< > $@.tmp && mv $@.tmp $@


# sort a bam file by name and index
%.byname.bam: %.bam
	rm -f $*.byname.tmp.{0,1,2,3,4,5,6,7,8,9}*.bam && samtools sort  -n -m $(SAMTOOLS_SORT_MEM) -T $*.byname.tmp -o $*.byname.tmp.bam $< && mv $*.byname.tmp.bam $@ 

# 1 - tag
define make-sort-bam-by-tag-rule=
# sort a bam file by a tag
%.bytag_$(1).bam: %.bam
	mkdir -p $$* && rm -f $$*/* && \
	samtools sort --thread $$(max_threads) -T $$*/ -m $$(SAMTOOLS_SORT_MEM_MT) -t $(1) $$< -o $$@.tmp && mv $$@.tmp $$@
endef
# rules to sort a bam file by the CB or CR tags
$(foreach tag,CR CB,$(eval $(call make-sort-bam-by-tag-rule,$(tag))))


# index a bam file
%.bam.bai: %.bam
	samtools index -@ $(max_threads) $<
# samtools 1.1 does not support the second argument :(
#	samtools index $< $@.tmp && mv $@.tmp $@

%.cram: %.bam
	samtools view --threads $(max_threads) -C -T $(reference_abspath) $< > $@.tmp && mv $@.tmp $@

#
# bigWig from bed 
# bed file needs to be sorted and converted to bedGraph
%.bw: %.bed $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt
	set -o pipefail && tail -n +2 $< | cut -f 1,2,3,5 | sort -k1,1 -k2,2n  > $<.sorted.bed && \
	bedGraphToBigWig $<.sorted.bed $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt $@.tmp && mv $@.tmp $@

# bigwig from bedgraph
%.bw: %.bedGraph $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt
	bedGraphToBigWig $< $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt $@.tmp && mv $@.tmp $@

define quant_levelFromFilename=
$(if $(findstring exons,$(1)),exon,$(if $(findstring genes,$(1)),gene,mRNA))
endef

# TSV (with feature value) is converted to bedGraph (http://genome.ucsc.edu/goldenPath/help/bedgraph.html)
# bedGraph generated is sorted
%.bedGraph: %.$(expr_ext) $(gff3_file_abspath).csv $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt
	set -o pipefail && tsv2bed.R $<  $(call quant_levelFromFilename,$*) $(gff3_file_abspath).csv $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt | \
	sort -k1,1 -k2,2n | \
	bedtools merge -scores mean -i - > $@.tmp &&\
	mv $@.tmp $@

# fail if file needs to be generated
%.cdna.fa:
	$(call p_error,Missing cdna file $@)

%.mapping_trans.tsv: %.gtf
	irap_gtf2mapping.pl --gtf $< --feature transcript > $@.tmp && mv $@.tmp $@

%.mapping_exons.tsv: %.gtf
	irap_gtf2mapping.pl --gtf $< --feature exon > $@.tmp && mv $@.tmp $@


# rename does not work on some distros :(
#	rename $*.mapping.tmp $*.mapping $*.mapping.*

################################################################################
# Does nothing...for now
quickcheck:

################################################################################
# stage 0 - setup/initialization
setup: setup_dirs setup_files

phony_targets+= setup setup_files

################################################################################
# Setup initial files
# file with the length of the features (gene, isoform, exon)
BOOTSTRAP_TARGETS+= setup_dirs $(trans_abspath) $(gtf_file_abspath).checked  $(gff3_file_abspath).filt.gff3   $(name)/version

SETUP_DATA_FILES+= setup_data_files2 $(auxdata_toplevel_folder)/$(gtf_file_basename).gene_class.txt $(index_files)   $(gtf_file_abspath).exon_id.gtf $(juncs_file_abspath)   $(annot_tsv)  $(auxdata_toplevel_folder)/$(reference_basename).introns.bed  $(auxdata_toplevel_folder)/$(reference_basename).genes.bed6 $(auxdata_toplevel_folder)/$(reference_basename).transcripts.bed6 $(auxdata_toplevel_folder)/$(reference_basename).genes.bed $(auxdata_toplevel_folder)/$(reference_basename).exons.bed $(feat_mapping_files) 

WAVE0_TARGETS+= $(auxdata_toplevel_folder)/$(gtf_file_basename).gene_class.txt $(word 1,$(index_files))   $(gtf_file_abspath).exon_id.gtf $(juncs_file_abspath)   $(annot_tsv)  $(auxdata_toplevel_folder)/$(reference_basename).introns.bed  $(auxdata_toplevel_folder)/$(reference_basename).genes.bed6 $(auxdata_toplevel_folder)/$(reference_basename).transcripts.bed6  $(feat_mapping_files)

setup2_files:=$(auxdata_toplevel_folder)/$(reference_basename).$(gtf_file_basename).data_info.tsv 
WAVE1_TARGETS+= $(setup2_files)  $(word 2,$(index_files))

setup_data_files2: $(setup2_files) 
#phony_targets+=setup_data_files2

# TODO: move $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.txt to $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt
# No need to include these files since they are generated by some rule
# $(auxdata_toplevel_folder)/introns.bed -> $(auxdata_toplevel_folder)/genes.bed $(auxdata_toplevel_folder)/exons.bed 
# $(gff3_file_abspath).csv -> $(gff3_file_abspath) 

setup_files: $(SETUP_DATA_FILES)  $(BOOTSTRAP_TARGETS)

print_stage0_files: setup_dirs
	echo $(SETUP_DATA_FILES)

print_stage1_files:
	echo $(STAGE1_OUT_FILES)

print_stage1_s_files:
	echo $(STAGE1_S_TARGETS)

###############################################################

print_stage1_folder:
	echo $(qc_toplevel_folder)

print_stage2_folder:
	echo $(mapper_toplevel_folder)

print_stage3_folder:
	echo $(quant_toplevel_folder)

###############################################################

ifneq ($(mapper),none)
# index the reference
$(word 1,$(index_files)): $(files_indexed)
	$(call p_info,INDEXING)
	$(call run_$(mapper)_index,$(word 1,$(files_indexed)))
ifneq ($(word 2,$(index_files)),)
# index the annotation/transcriptome
$(word 2,$(index_files)): $(files_indexed)
	$(call p_info,INDEXING 2)
	$(call run_$(mapper)_index_annot,$(word 1,$(files_indexed)))
endif
endif

# Create the file with gene classification by GeneId
$(refgeneclass_file): $(gtf_file_abspath)
	gtf2geneclass.sh $< > $@.tmp && mv $@.tmp  $@
#	$(error Missing required file $@)

# use tophat to generate the juncs file
# .gtf.juncs: .gtf
# note: the file may be empty
# in this case gtf2juncs returns an error that is ignored
$(juncs_file_abspath): $(gtf_file_abspath)
	tophat2_gtf_juncs $< > $@.tmp 
	mv $@.tmp $@

#
# .gff3.csv: .gff3
$(gff3_file_abspath): $(gtf_file_abspath)
	set -o pipefail && gtf2gff3.pl $< | sort -k1,1 -k4,4n > $@.tmp && mv $@.tmp $@

# Filter the gff3 file to contain only the chr that are in the fasta file
# use the faidx file to ensure that the ordering is the same
$(gff3_file_abspath).filt.gff3: $(gff3_file_abspath) $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.bed 
	bedtools intersect -wa -a $(gff3_file_abspath) -b $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.bed  > $@.tmp &&  bedtools sort -faidx $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.bed  -i $@.tmp > $@.tmp2 && mv $@.tmp2 $@ && rm -f $@.tmp


# 
$(gff3_file_abspath).csv: $(gff3_file_abspath)
	gff32csv.R $<   >   $@.tmp && mv $@.tmp $@

#$(auxdata_toplevel_folder)/$(reference_basename).gene_class.txt: $(refgeneclass_file)
$(auxdata_toplevel_folder)/$(gtf_file_basename).gene_class.txt: $(refgeneclass_file)
	cp $< $@.tmp && mv $@.tmp $@

# file with chr\tchr_size
$(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt: $(reference_abspath).fai
	cut -f 1,2 $< > $@.tmp && mv $@.tmp $@

#######################################################
# Collect some stats about the reference and annotation
data_stats: $(auxdata_toplevel_folder)/$(reference_basename).$(gtf_file_basename).data_info.tsv

print_data_stats_file:
	@echo $(auxdata_toplevel_folder)/$(reference_basename).$(gtf_file_basename).data_info.tsv

$(auxdata_toplevel_folder)/$(reference_basename).$(gtf_file_basename).data_info.tsv: $(auxdata_toplevel_folder)/$(reference_basename).irap_stats.tsv $(auxdata_toplevel_folder)/$(gtf_file_basename).irap_stats.tsv
	echo "Species|$(species)" | tr "|" "\n" > $@.tmp2
	echo "Genome|$(reference_basename)" | tr "|" "\n" > $@.tmp3
	paste $@.tmp2 $@.tmp3 $^   > $@.tmp && mv $@.tmp $@ && rm -f $@.tmp2 $@.tmp3

$(auxdata_toplevel_folder)/$(reference_basename).irap_stats.tsv: $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt
	irap_genome_stats -i $< -o $@.tmp && mv $@.tmp $@

$(auxdata_toplevel_folder)/$(gtf_file_basename).irap_stats.tsv: $(auxdata_toplevel_folder)/$(gtf_file_basename).lengths.Rdata
	irap_annot_stats -i $< -o $@.tmp && mv $@.tmp $@

#########################################################

# temporary rule:
# remove $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.txt
# and rename to  $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt
$(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.txt: $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.txt
	sort -k 1,1 $<  > $@.tmp  && mv $@.tmp $@

$(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.bed: $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.txt
	set -o pipefail && cat $< | awk 'BEGIN {OFS="\t";} {print  $$1,"1",$$2;}' > $@.tmp && mv $@.tmp $@

$(reference_abspath).fai: $(reference_abspath)
	samtools faidx $< && sleep 1

#ifeq ($(realpath ,)
# automatically generated annot file
# avoid generating file if it is in the reference directory
# Create a Rdata file to speedup loading the matrix 
$(refgeneannot_file) $(refgeneannot_file).Rdata: $(gtf_file_abspath)
	irap_gtf2annot --gtf $< -s $(species) --cores $(max_threads) --rdata -o $(refgeneannot_file).tmp && mv $(refgeneannot_file).tmp $(refgeneannot_file) && mv $(refgeneannot_file).tmp.Rdata $(refgeneannot_file).Rdata

$(auxdata_toplevel_folder)/$(reference_basename).gene.annot.tsv: $(refgeneannot_file) 
	cp $< $@.tmp && mv $@.tmp $@

$(auxdata_toplevel_folder)/$(reference_basename).gene.annot.tsv.Rdata: $(refgeneannot_file).Rdata
	cp $< $@.tmp && mv $@.tmp $@
#endif

$(auxdata_toplevel_folder)/empty.gene.annot.tsv: 
	echo '"ID","Name","locus","source","lname","GO","GOterm","KEGG","biotype"' | tr "," "\t" > $@

# collect genes, transcripts and exons lengths
# generate only once
$(gtf_file_abspath).lengths.Rdata: $(gtf_file_abspath)
	irap_gtf2featlength --gtf $< -o $@.tmp --cores $(max_threads) && mv $@.tmp.Rdata $@ && mv $@.tmp.gene_length.tsv $@.gene_length.tsv && mv $@.tmp.exon_length.tsv $@.exon_length.tsv && mv $@.tmp.trans_length.tsv $@.trans_length.tsv


$(auxdata_toplevel_folder)/$(gtf_file_basename).lengths.Rdata: $(gtf_file_abspath).lengths.Rdata
	sleep 2 && cp $< $@.tmp && mv $@.tmp $@ 

precious_targets+=$(auxdata_toplevel_folder)/$(gtf_file_basename).lengths.Rdata


$(DEXSEQ_GFF).lengths.Rdata: $(DEXSEQ_GFF)
	irap_DexSeqExonLen --gff $< -o $@.tmp --cores $(max_threads) && mv $@.tmp.Rdata $@

# sleep to ensure that the file will not have the same timestamp
$(auxdata_toplevel_folder)/$(gtf_file_basename).dexseq.lengths.Rdata: $(gtf_file_abspath).DEXSeq.gff.lengths.Rdata
	sleep 2 && cp $< $@.tmp && mv $@.tmp $@



# $(trans_file): 
# 	$(call p_error,Missing trans_file = $(trans_file))

################################################################################
# Setup directory structure
phony_targets+= setup_dirs

setup_dirs: $(tmp_dir)  $(auxdata_toplevel_folder)/ \
	$(qc_toplevel_folder)/ \
	$(if $(mapper),$(mapper_toplevel_folder)/) \
	$(if $(quant_method),$(quant_toplevel_folder)/)  \
	$(if $(de_method),$(de_toplevel_folder)/)
	$(call p_info,[DONE] Directory structure created)

bootstrap: $(BOOTSTRAP_TARGETS)

$(tmp_dir):
	@mkdir -p $@

%/:
	@mkdir -p $@




################################################################################
# Initial fastq QC
################################################################################
# alias: stage1=qc=quality_filtering_and_report
qc: quality_filtering_and_report 

phony_targets+= qc quality_filtering_and_report
# 
# TODO:	signature=filtering parameters (if they are the same then avoid recomputation between experiments)  minlen=$(min_read_len) min_qual=$(min_qual) qual_perc=$(min_qu
quality_filtering_and_report: setup	$(STAGE1_OUT_FILES)
	$(call p_info,[DONE] Quality filtering)

phony_targets+= do_qc
do_qc: $(STAGE1_OUT_FILES)

################################################################################
# Mapping
################################################################################
# All mapping files have the same name (but are placed in different folders)
# The mappings are placed in $(mapper_toplevel_folder) with the suffix (se/pe)hits.bam

# ** TODO **
# 1. Use singleton reads if available  (f.sing.fastq file)
mapping: $(STAGE1_S_TARGETS) $(mapper_toplevel_folder)/ $(mapper)_mapping
	$(call p_info,[DONE] Mapping)

print_stage2_files:
	echo $(STAGE2_OUT_FILES)

print_bam_filenames:
	echo $(bam_files)

#*************
# Generic rule
#*************
#####################
phony_targets+= $(mapper)_mapping mapping stage2_tracks print_bam_filenames

outbams=
ifeq ($(mapper),none)
$(mapper)_mapping:
else
outbams:=$(foreach p,$(pe), $(call lib2bam_folder,$(p))$(p).pe.hits.bam) $(foreach s,$(se), $(call lib2bam_folder,$(s))$(s).se.hits.bam)
STAGE2_OUT_FILES+=$(outbams)
STAGE2_TARGETS+=$(outbams)

WAVE2_TARGETS+=$(outbams)
$(mapper)_mapping: $(index_files) $(outbams)

endif

# stage2_tracks_targets is empty if report generation is disabled
stage2_tracks_targets:=$(call rep_browse,$(subst .bam,.bam.tracks,$(outbams)))	

STAGE3_OUT_FILES+=$(stage2_tracks_targets)
STAGE3_TARGETS+=$(stage2_tracks_targets)
WAVE3_TARGETS+=$(stage2_tracks_targets)

stage2_tracks: $(stage2_tracks_targets)
	$(call p_info,[DONE] Generated stage 2 tracks)

stage2_upload_tracks: $(subst .tracks,.tracks.uploaded,$(stage2_tracks_targets))
	$(call p_info,[DONE] Uploaded stage 2 tracks)

# 
mappingbyname: $(mapper)_mappingbyname

$(mapper)_mappingbyname: $(index_files) $(STAGE2BYNAME_OUT_FILES)

# interleaved fastq file
# special rule for GEM with PE
# 1-mapper - broken - discontinued support for GEM
define mapper_ifiles=
$(if $(subst gem,,$(1)),$(3),$(name)$(qc_method)/$(2)_int.f.fastq)
endef

##############################################
# Create the rules on the fly for each library
# $1 = libname


define make-se-bam-rule=
ifeq ($(mapper),star)
$(call lib2bam_folder,$(1))$(1).se.hits.bam.trans.bam: $(call lib2bam_folder,$(1))$(1).se.hits.bam
endif
$(call lib2bam_folder,$(1))$(1).se.hits.bam: $(call lib2filt_folder,$(1))$(1).f.fastq.gz $(index_files) $(gtf_file_abspath) $(reference_prefix)
	mkdir -p $(call lib2bam_folder,$(1)) && \
	$(call run_$(mapper)_map,$(1),$$<,$$@)
endef



define make-pe-bam-rule=
ifeq ($(mapper),star)
$(call lib2bam_folder,$(1))$(1).pe.hits.bam.trans.bam: $(call lib2bam_folder,$(1))$(1).pe.hits.bam
endif
$(call lib2bam_folder,$(1))$(1).pe.hits.bam: $(call mapper_ifiles,$(mapper),$(1),$(call lib2filt_folder,$(1))$(1)_1.f.fastq.gz $(call lib2filt_folder,$(1))$(1)_2.f.fastq.gz)   $(index_files)
	mkdir -p $(call lib2bam_folder,$(1)) && \
	$(call run_$(mapper)_map,$(1),$(call mapper_ifiles,$(mapper),$(1),$(call lib2filt_folder,$(1))$(1)_1.f.fastq.gz $(call lib2filt_folder,$(1))$(1)_2.f.fastq.gz),$$@)
endef


# create the output directories
ifneq ($(mapper),none)
#$(foreach l,$(se) $(pe),$(eval $(shell mkdir -p $(call lib2bam_folder,$(l)))))

# rules for SE libraries
$(foreach l,$(se),$(eval $(call make-se-bam-rule,$(l))))
# rules for PE libraries
$(foreach l,$(pe),$(eval $(call make-pe-bam-rule,$(l))))
endif


# interleaved fastq file
%_int.f.fastq: %_1.f.fastq %_2.f.fastq
	fastq2interleaved.pl $^ $@.tmp && mv $@.tmp $@

%.$(mapper).index: %.fa


$(mapper_toplevel_folder)/alignments.bam: $(foreach s,$(se), $(call lib2bam_folder,$(s))$(s).se.hits.bam)
	set -o pipefail && $(call samcat,$^) | samtools sort -m $(SAMTOOLS_SORT_MEM) -T $@.sorted -o $@.sorted.bam -  && mv  $@.sorted.bam $@  && samtools index $@


$(mapper_toplevel_folder)/alignments.bam.paired.bam: $(foreach s,$(pe), $(call lib2bam_folder,$(s))$(s).pe.hits.bam)
	set -o pipefail && $(call samcat,$^) | samtools sort  -m $(SAMTOOLS_SORT_MEM) -T $@.tmp -o $@.tmp.bam -  && mv $@.tmp.bam $@ && samtools index $@


$(reference_abspath)_files: $(reference_abspath)
	irap_fasta_split.pl $< $@



################################################################################
# STAGE 3
################################################################################


STAGE3_TSV_FILES+=$(call exons_quant_files)\
	       $(call transcripts_quant_files)\
	       $(call genes_quant_files)

print_stage3_files:
	echo $(sort $(STAGE3_OUT_FILES))

print_stage3_s_files:
	echo $(sort $(STAGE3_S_OFILES))

print_stage3_targets:
	echo $(sort $(STAGE3_TARGETS))

print_stage3_s_targets: 
	echo $(STAGE3_S_TARGETS)

################################################################################
# Stage 4
################################################################################

phony_targets+= print_stage4_files print_stage4_targets

print_stage4_files:
	echo $(STAGE4_OUT_FILES)
print_stage4_targets:
	echo $(STAGE4_TARGETS)

################################################################################
# Stage 5
################################################################################

phony_targets+= print_stage5_files  print_stage5_targets
print_stage5_files:
	echo $(STAGE5_OUT_FILES)

print_stage5_targets:
	echo $(STAGE5_TARGETS)

#############################################################
# GSE
###########################################
# IRAP targets
# todo: mv this to the irap_gse file

#############################################################
## Cleanup
phony_targets+= clean full_clean clean_data_files


clean: $(CLEANUP_TARGETS)

full_clean: clean_data_files
	rm -fr $(name)/

clean_data_files:
	rm -rf $(SETUP_DATA_FILES) $(index_files) $(reference_abspath).fa


# delete all files related to the libraries up to quantification
lib_full_clean: 
	rm -rf $(STAGE1_OUT_FILES)
	$(foreach l,$(se) $(pe), $(if $(subst .,,$($(l)_dir)), rm -rf $(qc_toplevel_folder)/$($(l)_dir)))
	rm -rf $(STAGE2_OUT_FILES)
	##$(foreach l,$(se) $(pe), $(if $(subst .,,$($(l)_dir)), rm -rf $(mapper_toplevel_folder))/$($(l)_dir)))
	rm -rf $(STAGE3_OUT_FILES)
	##	$(foreach l,$(se) $(pe), $(if  $(subst .,,$($(l)_dir)), rm -rf $(quant_toplevel_folder)/$($(l)_dir)))

# TODO: archive (delete everything except the "main" output files for each stage
#
#############################################################
################################################################################
phony_targets+= stage0 stage1 stage2 stage3 stage3a stage3as stage3b stage4 stage5
silent_targets+= 

# alias - make it easier for the user
stage0: setup
stage1: setup quality_filtering_and_report
stage2: setup mapping
stage3: setup $(STAGE3_OUT_FILES) $(STAGE3_S_OFILES) $(STAGE3_S_OFILES)
stage3a: setup $(quant_method)_quant 
stage3as: setup quantification_s
stage3s: $(STAGE3_S_OFILES)
#stage3b: setup stage3a
# deprecated
stage3b: setup $(shell rm -f $(quant_toplevel_folder)/rawcounts.$(quant_method).$(expr_ext)) stage3a
stage4: setup $(STAGE4_OUT_FILES)
stage5: setup $(STAGE5_OUT_FILES)

# deprecated: to be removed in the future
assemble: quantification


# *************
# Transcriptome (WIP)
# *************
%.trans.fa: $(reference_abspath) %.ref.gtf
	extract_transcripts $^ > $@

%.cdna_ncrna.fa: %.cdna.all.fa %.ncrna.fa
	cat $^ > $@

%.ref.gtf: %.gtf %.ncrna.fa
	filterGTF.rb  $*.ncrna.fa $< > $@

%.ref.gtf: %.ref.gff
	ensembl_gtf_to_gff.pl $< > $@

##############################################################################
# Reporting

phony_targets+= report
report: report_all_targets

###################################################
# Execution status
# Creates a status.tsv file with the completion status of each stage of the pipeline based on the files created
phony_targets+= status_html status stage2_tracks stage3_tracks stage4_tracks
status_html: $(name)/report/status.html
	$(call p_info,Created $<)

status: $(name)/$(name).status.tsv
	$(call p_info,Created $<)
	@cat $<

### Tracks
define exclude_empty=
$(foreach f,$(1), $(if $(call is_empty_file,$(f)),,$(f)))
endef

define stage3_tracks_targets=
$(subst .$(expr_ext),.$(expr_ext).tracks,$(call exclude_empty,$(STAGE3_TSV_FILES)))
endef

stage3_tracks: $(call stage3_tracks_targets)
	$(call p_info,[DONE] Generated stage 3 tracks)

stage3_upload_tracks: $(subst .tracks,.tracks.uploaded,$(stage3_tracks_targets))
	$(call p_info,[DONE] Uploaded stage 3 tracks)

stage4_tracks_targets=$(subst .$(expr_ext),.$(expr_ext).tracks,$(STAGE4_OUT_FILES))

stage4_tracks: $(stage4_tracks_targets)
	$(call p_info,[DONE] Generated stage 4 tracks)

stage4_upload_tracks: $(subst .tracks,.tracks.uploaded,$(stage4_tracks_targets))
	$(call p_info,[DONE] Uploaded stage 4 tracks)


# target to generate some tracks
get_tracks: stage2_tracks stage3_tracks stage4_tracks 


###################################################
# 
#$(name)/report/all_options.txt:
#	( $(foreach v, $(call interesting_vars), $(v) = $($(v)))) ) > /dev/stdout

#define pprint_var
#
#echo "$(1) = $($(1))"
#
#endef
#define interesting_vars=
#$(filter-out %_sd %_rs $(se) $(pe),$(foreach v,$(.VARIABLES),$(if $(subst file,,$(#origin $(v))),,$(v) )))
#endef

###################################################
# Alternative targets (lightweight)
do_index: setup_dirs $(word 1,$(index_files))

# TODO: rename the BAM file and add the headers if necessary
do_mapping: $(mapper)_mapping

###################################################
phony_targets+= save_cache 
save_cache: 
	($(call list_cached_vars)) > $(cached_vars_file).tmp && mv $(cached_vars_file).tmp $(cached_vars_file)


######################################
## job scheduler
## targets are placed into bins/waves

phony_targets+= run_wave_b run_wave_0 run_wave_1 run_wave_2 run_wave_3 run_wave_4 run_wave_5 run_wave3_s run_wave3_p run_wave_6 run_wave_7

# Note
# targets variables should contain the name of targets that may be run in parallel 
# Up to 
WAVEB_TARGETS+=$(BOOTSTRAP_TARGETS)
WAVE0_TARGETS+=
WAVE1_TARGETS+=$(STAGE1_TARGETS)
# BAM+basic stats
WAVE2_TARGETS+=
# $(STAGE3_TARGETS)
WAVE3_TARGETS+=
# targets that should be complete before running stage3
# e.g., post alignment (sorting bams)
WAVE3_p_TARGETS+=
WAVE3_s_TARGETS+=
WAVE4_TARGETS+=
WAVE5_TARGETS+=
WAVE6_TARGETS+=
WAVE7_TARGETS+=

run_wave_b: $(WAVEB_TARGETS)
run_wave_0: $(WAVE0_TARGETS)
run_wave_1: $(WAVE1_TARGETS)
run_wave_2: $(WAVE2_TARGETS)
run_wave_3: $(WAVE3_TARGETS)
run_wave_3_s: $(WAVE3_s_TARGETS)
run_wave_3_p: $(WAVE3_p_TARGETS)
run_wave_4: $(WAVE4_TARGETS)
run_wave_5: $(WAVE5_TARGETS)
run_wave_6: $(WAVE6_TARGETS)
run_wave_7: $(WAVE7_TARGETS)

print_wave_b_targets:
	$(call wave_echo,$(sort $(WAVEB_TARGETS)))

print_wave_0_targets:
	$(call wave_echo,$(sort $(WAVE0_TARGETS)))

print_wave_1_targets:
	$(call wave_echo,$(sort $(WAVE1_TARGETS)))

print_wave_2_targets:
	$(call wave_echo,$(sort $(WAVE2_TARGETS)))

print_wave_3_targets:
	$(call wave_echo,$(sort $(WAVE3_TARGETS)))

print_wave_3_s_targets:
	$(call wave_echo,$(sort $(WAVE3_s_TARGETS)))

print_wave_3_p_targets:
	$(call wave_echo,$(sort $(WAVE3_p_TARGETS)))

print_wave_4_targets:
	$(call wave_echo,$(sort $(WAVE4_TARGETS)))

print_wave_5_targets:
	$(call wave_echo,$(sort $(WAVE5_TARGETS)))

print_wave_6_targets:
	$(call wave_echo,$(sort $(WAVE6_TARGETS)))

print_wave_7_targets:
	$(call wave_echo,$(sort $(WAVE7_TARGETS)))

###################################################
# libraries
print_libs:
	@echo -n $(file > /dev/stdout,$(pe) $(se))

print_se_libs:
	@echo -n $(file > /dev/stdout,$(se))

print_pe_libs:
	@echo -n $(file > /dev/stdout,$(pe)) 

phony_targets+= print_libs lib_isl
###################################################
lib_isl: $(STAGE1_S_TARGETS) stage2 stage3as $(WAVE3_s_TARGETS) $(STAGE3_S_OFILES)

lib_isl_stage1: $(STAGE1_S_TARGETS)

###################################################
# Keep the versions used in the top level folder
$(name)/version: $(IRAP_DIR)/version
	touch $@ && cat $@ $(IRAP_DIR)/version |sort -u > $@.tmp && mv $@.tmp $@
###################################################
# FORCE the program to run even if files haven't changed
FORCE: 

#PHONY: performance improvement. Tell MAKE that those targets don't generate any files. 
.PHONY:  $(phony_targets)
.SILENT: $(silent_targets)
# disable deletion of temporary files
.SECONDARY: $(precious_targets)
.PRECIOUS: $(precious_targets)
###################################################
