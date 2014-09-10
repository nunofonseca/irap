#; -*- mode: Makefile;-*-
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
#    $Id: 0.1.3 Nuno Fonseca Wed Dec 26 16:16:19 2012$
# =========================================================
# 1: bam_file 2: fixed_bam_file
define bam_fix_nh=
$(if $(subst y,,$(fix_NH)),mv $(1) $(2),bam_fix_NH $(1) $(2))
endef

ifndef SAMTOOLS_SORT_MEM
SAMTOOLS_SORT_MEM=5000000000
endif

run_picard=java -jar $(abspath  $(irap_path)/../bin/picard-tools/$(1).jar)

########################################################################
# Quality encoding parameters
#
# Bowtie1, Bowtie2, Tophat
define qual_option1=
$(shell if [ "$(1)" == "33" ]; then echo "--solexa-quals"; else echo "--phred64-quals"; fi)
endef

define qual_option33_64=
$(shell if [ "$(1)" == "33" ]; then echo "--phred33-quals"; else echo "--phred64-quals"; fi)
endef

# BWA
define qual_bwa=
$(if $(findstring 33,$($(1)_qual)),,-I)
endef

# Soapsplice
define qual_soap=
$(shell if [ "$(1)" == "33" ]; then echo "-q 1"; else echo "-q 0"; fi)
endef

# mapsplice
define qual_mapsplice=
$(shell if [ "$(1)" == "33" ]; then echo "--qual-scale phred33"; else echo "--qual-scale solexa64"; fi)
endef

# passion
define qual_passion=
$(shell if [ "$(1)" == "33" ]; then echo "--phred33"; else echo "--phred64"; fi)
endef

# default values
# options that are passed directly to the mapper
tophat1_map_options?=
tophat2_map_options?=
bowtie1_map_options?=
bowtie2_map_options?=
bwa1_map_options?=
bwa2_map_options?=
gsnap_map_options?=
smalt_map_options?=
gem_map_options?=
gem_index_options?=
gems_map_options?=
gems_index_options?=
soapsplice_map_options?=
mapsplice_map_options?=
star_map_options?=
star_index_options?=
osa_index_options?=

tophat1_aln_options?=
tophat2_aln_options?=
bowtie1_aln_options?=
bowtie2_aln_options?=
bwa1_aln_options?=
bwa2_aln_options?=
gsnap_aln_options?=
smalt_index_options?=
soapsplice_aln_options?=
mapsplice_aln_options?=
star_aln_options?=
osa_aln_options?=

#####################################################
# Bowtie1
# ifeq ($(mapper),bowtie1) 
#  mapper_splicing=no
# endif

define run_bowtie1_index=
	irap_map.sh bowtie1  bowtie-build --offrate 3 $(1) $(1)
endef

#
# same arguments used for *_index
define bowtie1_index_filename=
$(2).1.ebwt
endef

# -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities
#  -q                 query input files are FASTQ .fq/.fastq (default)
# --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 64)
# -a/--all           report all alignments per read (much slower than low -k)
#--best             hits guaranteed best stratum; ties broken by quality
#  --strata           hits in sub-optimal strata aren't reported (requires --best)
# -m <int>           suppress all alignments if > <int> exist (def: no limit)
# -k <int>           report up to <int> good alignments per read (default: 1)
#bowtie1_map_params+=  --fullref --sam -q -a --best --strata   --chunkmbs 256 -m 100 -p $(max_threads) $(bowtie1_map_options)
bowtie1_map_params+=  --fullref --sam -q --best --strata -k $(max_hits)   --chunkmbs 256 -m `expr $(max_hits) + 1` -p $(max_threads) $(bowtie1_map_options)

# Set mean insert size and standard deviation for bowtie1 and bowtie2
define bowtie_ins_sd_params=
    $(if $(findstring $(1),$(pe)), -I 0 -X `expr $($(1)_ins) + $($(1)_sd)`)
endef

define bowtie1_file_params=
	$(if $(findstring $(1),$(pe)), $(call bowtie_ins_sd_params,$(1)) -1 $(word 1,$(2)) -2 $(word 2,$(2)), $(call tophat_qual_option,$($(1)_qual)) $(2))
endef

define run_bowtie1_map=
	irap_map.sh  bowtie1 bowtie  $(bowtie1_map_params) $(subst .1.ebwt,,$(index_files)) $(call bowtie1_file_params,$(1),$(2))  | \
	samtools view -T $(reference_abspath) -F 0xC -bS - > $(3).tmp.bam && \
	$(call bam_fix_nh,$(3).tmp.bam,-) | \
	samtools sort -m $(SAMTOOLS_SORT_MEM) - $(3).tmp && \
	rm -f $(3).tmp2.bam && \
	mv $(3).tmp.bam $(3)
endef

#####################################################
# Bowtie2

# ifeq ($(mapper),bowtie2) 
#  mapper_splicing=no
# endif

define run_bowtie2_index=
	irap_map.sh bowtie2  bowtie2-build --offrate 3 $(1) $(1)
endef

# same arguments used for *_index
define bowtie2_index_filename=
$(2).1.bt2
endef

define bowtie2_file_params=
	$(if $(findstring $(1),$(pe)), $(call bowtie_ins_sd_params,$(1)) -1 $(word 1,$(2)) -2 $(word 2,$(2)), $(call tophat_qual_option,$($(1)_qual)) -U $(2))
endef

bowtie2_map_params+= --end-to-end -k $(max_hits) -p $(max_threads)  $(bowtie2_map_options)
# 0xC filter correctly paired alignemnts
define run_bowtie2_map=
	irap_map.sh  bowtie2 bowtie2  $(bowtie2_map_params) -x $(subst .1.bt2,,$(index_files)) $(call bowtie2_file_params,$(1),$(2))  | \
	samtools view -T $(reference_abspath) -F 0xC -bS - > $(3).tmp.bam && \
	$(call bam_fix_nh,$(3).tmp.bam,-) | \
	samtools sort -m $(SAMTOOLS_SORT_MEM) - $(3).tmp && \
	rm -f $(3).tmp2.bam && \
	mv $(3).tmp.bam $(3)
endef

#####################################################
# Tophat 1 & 2

#--solexa-quals
define tophat_qual_option=
	$(shell if [ "$(1)" == "33" ]; then echo ""; else echo "--solexa-quals"; fi)
endef

define tophat_seglength_option=
	$(shell if [ $(1) \< 60 ]; then echo "--segment-length 20"; else echo ""; fi)
endef


# override: the reads need to be trimmed when tophat is used 
ifeq ($(mapper),tophat1)
	trim_reads=y
endif
ifeq ($(mapper),tophat2)
	trim_reads=y
endif

min_intron_len=6
tophat1_map_params= --min-intron-length $(min_intron_len) $(tophat1_map_options)
tophat2_map_params= --no-coverage-search --min-intron-length $(min_intron_len) $(tophat2_map_options)
tophat_no_splicing= --no-gtf-juncs --no-novel-juncs --transcriptome-only --no-novel-indels

ifeq ($(mapper_splicing),no)
	tophat1_map_params+= $(tophat_no_splicing)
	tophat2_map_params+= $(tophat_no_splicing)
endif


# strand option
# irap
# first  -> --library-type=fr-firststrand
# second -> --library-type=fr-secondstrand 
irap_strand2tophatoption=$(if $(findstring $(1),first),fr-firststrand,$(if $(findstring $(1),second),fr-secondstrand,fr-unstranded))

# 1 - libname
define tophat_strand_params=
	$(if $(call lib_strand_info,$(1)),--library-type=$(call irap_strand2tophatoption,$($(1)_strand)),)
endef


define tophat_ins_sd_params=
	$(if $(findstring $(1),$(pe)), --mate-inner-dist $($(1)_ins) --mate-std-dev $($(1)_sd))
endef

define run_tophat1_index=
        $(call run_bowtie1_index,$(1),$(1))
endef

# generate the transcriptome once (v.2.0.10 or above)\
#	
define run_tophat2_index=
        $(call run_bowtie2_index,$(1),$(1)) && \
	mkdir -p $(call tophat2_trans_index_filename,$(1),$(1)) &&  \
	irap_map.sh tophat2 tophat -G $(gtf_file_abspath) --transcriptome-index $(call tophat2_trans_index_filename,$(1),$(1))/ $(tophat_reference_prefix)
endef

# same arguments used for *_index
define tophat1_index_filename=
	$(call bowtie1_index_filename,$(1),$(1))
endef
define tophat2_index_filename=
	$(call bowtie2_index_filename,$(1),$(1))
endef
define tophat2_trans_index_filename=
	$(subst .fa,,$(1))_th2_trans
endef


# Warning: tophat does not like reads with different sizes in the same file
# splice mismatches -m0 (0-2)"
# --transcriptome-index
define tophat_setup_dirs=
	$(shell if [ ! -h $(reference_prefix).fa ] ; then  ln -s $(reference_prefix) $(reference_prefix).fa; fi)
	$(shell if [ ! -e $(call lib2bam_folder,$(1))$(1)/tmp ] ; then mkdir -p $(call lib2bam_folder,$(1))$(1)/tmp; fi)
endef

tophat_reference_prefix=$(reference_prefix)

# cuffdiff complains about the order..
# Error: sort order of reads in BAMs must be the same
# TODO: test and do the same to tophat2
define run_tophat1_map=
        $(call tophat_setup_dirs,$(1))
	irap_map.sh tophat1 tophat  -p $(max_threads) $(call tophat_seglength_option,$($(1)_rs),$(1)) $(call tophat_qual_option,$($(1)_qual)) $(call tophat_strand_params,$(1)) $(tophat1_map_params) $(call tophat_ins_sd_params,$(1)) -G $(gtf_file_abspath) --tmp-dir $(call lib2bam_folder,$(1))$(1)/tmp -o $(call lib2bam_folder,$(1))$(1)	 $(tophat_reference_prefix) $(2) &&\
	samtools sort -m $(SAMTOOLS_SORT_MEM) $(call lib2bam_folder,$(1))$(1)/accepted_hits.bam  $(3).tmp && \
	mv $(3).tmp.bam $(3)
endef



#picard is necessary to fix some issues with PE aligned reads
#  ("Malformed SAM line: MRNM == '*' although flag bit &0x0008 cleared", 'line 67')
#$(if $(findstring $(1),$(pe)), $(call run_picard,FixMateInformation) INPUT=$(call lib2bam_folder,$(1))$(1)/accepted_hits.bam ASSUME_SORTED=false VALIDATION_STRINGENCY=LENIENT TMP_DIR=$(call lib2bam_folder,$(1))$(1) && ) 
# bam_tophat2_pe_fix fix unmapped reads flags
define run_tophat2_map=
        $(call tophat_setup_dirs,$(1))
	irap_map.sh tophat2 tophat2  -p $(max_threads) $(call tophat_seglength_option,$($(1)_rs),$(1)) $(call tophat_qual_option,$($(1)_qual)) $(call tophat_strand_params,$(1)) $(tophat2_map_params) $(call tophat_ins_sd_params,$(1)) -G $(gtf_file_abspath) --tmp-dir $(call lib2bam_folder,$(1))$(1)/tmp -o $(call lib2bam_folder,$(1))$(1) --transcriptome-index $(call tophat2_trans_index_filename,$(file_indexed),$(file_indexed))/  --rg-id $(1) --rg-sample $(1)  $(tophat_reference_prefix) $(2) && \
	samtools merge - $(call lib2bam_folder,$(1))$(1)/accepted_hits.bam $(call lib2bam_folder,$(1))$(1)/unmapped.bam | \
	bam_tophat2_fix - - | \
	samtools sort -m $(SAMTOOLS_SORT_MEM) - $(call lib2bam_folder,$(1))$(1)/$(1) &&\
	mv $(call lib2bam_folder,$(1))$(1)/$(1).bam $(3)
endef


####################################################
# SMALT
smalt_map_params= -n $(max_threads)  -f samsoft  $(smalt_map_options)
smalt_index_params= -s 3 $(smalt_index_options)
# Note: -p option broken in 0.6.4 when used in conjunction with PE data
# only use -p with SE

# reference, index.prefix
define run_smalt_index=
	irap_map.sh smalt smalt_x86_64 index $(smalt_index_params) $(call smalt_index_filename,$(1)) $(1) && touch  $(call smalt_index_filename,$(1))
endef

# PE
# maximum insert size -i
# minimum insert size -k
# -l type of read pair library (pe|mp|pp)
# enable -x?
define smalt_ins_param=
$(if $(findstring $(1),$(pe)), -j 0 -i `expr $($(1)_ins) + $($(1)_sd)`,-p)
endef

# same arguments used for *_index
define smalt_index_filename=
$(1).smalt.index
endef

# pe
# insert size is the distance between the 5'ends of the mapped reads
# 5'-3'<---....---> 3'-5' 
define run_smalt_map=
	irap_map.sh smalt smalt_x86_64 check $(2) && \
	irap_map.sh smalt smalt_x86_64 map $(smalt_map_params) $(call smalt_ins_param,$(1)) -o $(3).unsrt.sam $(index_files) $(2) && \
	samtools view -T $(reference_abspath) -bS $(3).unsrt.sam > $(3).tmp.bam  && \
	$(call bam_fix_nh,$(3).tmp.bam,$(3).unsrt.bam) && \
	samtools sort -m $(SAMTOOLS_SORT_MEM)  $(3).unsrt.bam $(3).tmp2 && \
	rm -f $(3).tmp.bam $(3).unsrt.* && \
	mv $(3).tmp2.bam $(3)
endef

######################################################
# GSNAP

# -N 1 : novel splicing
# -m, --max-mismatches=FLOAT
# TODO: tune parameters!!!
gsnap_map_params=  -A sam  -t $(max_threads)  $(gsnap_map_options)
gsnap_index_options= 

ifeq ($(mapper_splicing),no)
 gsnap_map_params+= -N 0
else
 gsnap_map_params+= -N 1
endif

define gsnap_ins_param=
$(if $(findstring $(1),$(pe)), --pairexpect=$($(1)_ins)  --pairdev=$($(1)_sd))
endef

define run_gsnap_index=
	irap_map.sh gsnap gmap_build $(gsnap_index_options)  -D $(reference_dir) -d `basename $(call gsnap_index_dirname,$(1))` -T  $(reference_dir) $(1)  && touch $(call gsnap_index_filename,$(1))
endef

define gsnap_index_filename=
$(1).gsnap.index
endef

define gsnap_index_dirname=
$(1).gsnap
endef


define run_gsnap_map=
	irap_map.sh gsnap gsnap $(gsnap_map_params) $(call gsnap_ins_param,$(1))  -D $(reference_dir) -d `basename $(index_files)|sed "s/.index//"` $(2)  | samtools view -T $(reference_abspath) -bS - | \
	samtools sort -m $(SAMTOOLS_SORT_MEM) - $(3).tmp && \
	mv $(3).tmp.bam $(3) 
endef

######################################################
# SOAPSPLICE

# -S strand (1,2,3)
# -q quality (1/33,2/64)
# -p max threads
# -f 2  (Sam output)
# -m <int>    Maximum mismatch for one-segment alignment, <= 5, 3 (default)
#  -g <int>    Maximum indel for one-segment alignment, <= 2, 2 (default)
#  -i <int>    Length of tail that can be ignored in one-segment alignment, 7 (default)
# -l <int>    The minimum distance between paired-end reads, 50 (default)
#  -I <int>    The insert length of paired-end reads
soapsplice_map_params=  -S 3  -p $(max_threads)  -f 2  -g 1000 $(soapsplice_map_options)
soapsplice_index_options= 

ifeq ($(mapper_splicing),no)
 soapsplice_map_params+=
$(error mapper_splicing option must be on when using soapsplice)
else
 soapsplice_map_params+=
endif

define run_soapsplice_index=
	sed 's/ .*//g' $(1) > $(1).soapsplice.fa  && irap_map.sh soap_splice  2bwt-builder $(1).soapsplice.fa $(call soapsplice_index_filename,$(1)) && touch $(call soapsplice_index_filename,$(1)) 
endef

define soapsplice_index_filename=
$(1).soapsplice.index
endef

# $(1)=libname
define soapsplice_ins_param=
$(if $(findstring $(1),$(pe)),-I $($(1)_ins)) 
endef

# $(1)=libname $(2)=filename(s)
define soapsplice_file_param=
$(if $(findstring $(1),$(pe)),-1 $(word 1,$(2)) -2 $(word 2,$(2)),-1 $(2))
endef

# WARNING: cufflinks requires the spliced alignments with the XS record defined (soapsplice does not include the XS field in the SAM file
define run_soapsplice_map=
	irap_map.sh soap_splice soapsplice $(soapsplice_map_params)  -d $(index_files).index -o $(3).unsrt $(call soapsplice_file_param,$(1),$(2)) $(call soapsplice_ins_param,$(1))  -q $($(1)_qual) &&\
	samtools view -T $(reference_abspath) -bS $(3).unsrt.sam > $(3).unsrt.bam  && \
	$(call bam_fix_nh,$(3).unsrt.bam,$(3).unsrt2.bam) && rm -f $(3).unsrt.* && \
	echo "@HD	VN:1.0	SO:coordinate" > $(3).tmp.H && \
	samtools view -H  $(3).unsrt2.bam >> $(3).tmp.H && \
	samtools reheader $(3).tmp.H $(3).unsrt2.bam | samtools sort -m $(SAMTOOLS_SORT_MEM) - $(3).tmp && \
	rm -f $(3).tmp.H  $(3).unsrt2.bam && \
	mv $(3).tmp.bam $(3)
endef

########################################################################
# BWA 1 
bwa1_index_options= -a bwtsw
bwa1_aln_params:= -t $(max_threads) $(bwa1_aln_options)
bwa1_map_params:= $(bwa1_map_options)

define bwa1_index_filename=
$(1).bwa1.amb
endef

define run_bwa1_index=
	sed 's/ .*//g' $(1) > $(1).bwa1.fa  && irap_map.sh bwa bwa index  $(bwa1_index_options) -p $(1).bwa1 $(1).bwa1.fa 
endef


define bwa1_map_se=
	irap_map.sh bwa bwa aln $(bwa1_aln_params)  $(call qual_bwa,$(1)) $(subst .amb,,$(index_files)) $(2) | \
	irap_map.sh bwa bwa samse $(bwa1_map_params)  $(subst .amb,,$(index_files)) /dev/fd/1 $(2) > $(3).sam
endef

define bwa1_map_pe=
	irap_map.sh bwa bwa aln $(bwa1_aln_params) -f $(word 1,$(2)).bwa1 $(call qual_bwa,$(1)) $(subst .amb,,$(index_files)) $(word 1,$(2))  && \
	irap_map.sh bwa bwa aln $(bwa1_aln_params) -f $(word 2,$(2)).bwa1 $(call qual_bwa,$(1)) $(subst .amb,,$(index_files)) $(word 2,$(2))  && \
	irap_map.sh bwa bwa sampe -a `expr $($(1)_ins) + $($(1)_sd)` $(call qual_bwa,$(1)) $(bwa1_map_params) $(subst .amb,,$(index_files))  $(word 1,$(2)).bwa1 $(word 2,$(2)).bwa1  $(2) > $(3).sam
endef

#-i 0  -e 10000
define run_bwa1_map=
	$(if $(findstring $(1),$(pe)),$(call bwa1_map_pe,$(1),$(2),$(3)),$(call bwa1_map_se,$(1),$(2),$(3))) &&\
	samtools view -T $(reference_abspath) -bS $(3).sam > $(3).tmp.bam &&\
	$(call bam_fix_nh,$(3).tmp.bam,-) | \
	samtools sort -m $(SAMTOOLS_SORT_MEM) - $(3).tmp && \
	rm -f $(3).tmp2.bam && \
	mv $(3).tmp.bam $(3)
endef

# ########################################################################
# BWA 2 
bwa2_index_options= $(bwa1_index_options)
bwa2_aln_params:= -t $(max_threads) -q 3 $(bwa2_aln_options)
bwa2_map_params:= $(bwa2_map_options)

# reuse bwa2 index file if possible/necessary
define bwa2_index_filename=
$(1).bwa2.amb
endef

define run_bwa2_index=
	sed 's/ .*//g' $(1) > $(1).bwa2.fa  && irap_map.sh bwa bwa index  $(bwa2_index_options) -p $(1).bwa2 $(1).bwa2.fa 
endef

define run_bwa2_map=
	irap_map.sh bwa bwa bwasw $(bwa2_map_params)  $(subst .amb,,$(index_files)) $(2) | \
	samtools view -T $(reference_abspath) -bS - > $(3).tmp.bam &&\
	$(call bam_fix_nh,$(3).tmp.bam,-)  | \
	samtools sort -m $(SAMTOOLS_SORT_MEM) - $(3).tmp && \
	rm -f $(3).tmp2.bam && \
	mv $(3).tmp.bam $(3)
endef

########################################################################
# GEM (general)
# [Single-end alignment]
#     --mismatch-alphabet <symbols>            (default="ACGT")
#     -m <max_mismatches>|<%_mismatches>       (default=0.04)
#     -e <max_edit_distance>|<%_differences>   (default=0.00)
#     --min-matched-bases <number>|<%>         (default=0.80)
#     --max-big-indel-length <number>          (default=15)
#     -s|--strata-after-best <number>          (default=0)
#     --fast-mapping <number>|'adaptive'       (default=false)
#     --unique-mapping                         (default=false)
#     --allow-incomplete-strata <number>|<%>   (default=0.00)
#   [Selecting alignments for output (single-end mode) or pairing (paired-end mode)]
#     -d|--max-decoded-matches <number>|'all'  (default=20)
#     -D|--min-decoded-strata <number>         (default=1)
#   [Paired-end alignment]
#     -p|--paired-end-alignment                (default=false)
#     -b|--map-both-ends                       (default=false)
#     --min-insert-size <number>               (default=0)
#     --max-insert-size <number>               (default=1000)
#     -E <max_edit_distance>|<%_differences>   (default=0.30)
#     --max-extendable-matches <number>|'all'  (default=20)
#     --max-extensions-per-match <number>      (default=1)
#     --unique-pairing                         (default=false)
# ADD
# --max-decoded-matches
gem_map_params= --threads $(max_threads) $(gem_map_options)
gem_index_params= --threads $(max_threads) --max-memory unlimited $(gem_index_options)
#gem_index_params= -t $(max_threads)--max-memory unlimited $(gem_index_options) --for-rna-mapper

define gem_qual_option=
	$(shell if [ "$(1)" != "33" ]; then echo "-q 'offset-64'"; else echo "-q offset-33"; fi)
endef

define gem_index_filename=
$(1).gem.index
endef

# TODO: --max-extendable-matches=0
# !!! PE maybe flagged has invalid (due to distance) 
##     --min-insert-size <number>               (default=0)
#     --max-insert-size <number>               (default=100000) ...max intron size
define gem_pairing_param=
$(if $(findstring $(1),$(pe)), --paired-end-alignment --min-insert-size `expr $($(1)_ins) - $($(1)_sd)`  --max-insert-size `expr $($(1)_ins) + $($(1)_sd)`)
endef

define gem2sam_pairing_param=
$(if $(findstring $(1),$(pe)),   --expect-paired-end-reads,  --expect-single-end-reads)
endef

# (1) ref, index file
define run_gem_index=
	irap_map.sh GEM gem-indexer -i $(1) -o $(call gem_index_filename,$(1)) $(gem_index_params) && touch $(call gem_index_filename,$(1))
endef

# TODO: SAM/BAM file reference is wrong
# SE reads have paired flags enabled
# sed code is to replaces the names of the sequences from 1 dna, 2 dna, ... to 1,2,...
# lib,fastq,out
# PE mode requires interlaced fastq file
#	 sed -i  -e 's/\([a-zA-Z0-9]\+\) dna:/\1:/g' $(3).gem.map &&
define run_gem_map=
	 irap_map.sh GEM gem-mapper $(gem_map_params) $(call gem_pairing_param,$(1)) $(call gem_qual_option,$($(1)_qual)) -I $(index_files).gem -i $(2) -o $(3).gem && \
	 irap_map.sh GEM  gem-2-sam -i $(3).gem.map -I $(index_files) --emit-correct-flags -T $(max_threads) $(call gem2sam_pairing_param,$(1)) $(call gem_qual_option,$($(1)_qual)) -o /dev/fd/1 | sed  -e 's/\([a-zA-Z0-9]\+\) dna/\1/g'   | \
	 samtools  view -T $(reference_abspath) -bS - > $(3).tmp2.bam &&\
	 bam_fix_se_flag $(3).tmp2.bam - | \
	 samtools sort -m $(SAMTOOLS_SORT_MEM) - $(3).tmp   && \
	 mv $(3).tmp.bam $(3) && rm -rf $(3).tmp2* $(3).tmp.*
endef

# old
#	 sed -i -e 's/\([a-zA-Z0-9]\+\) dna:/\1:/g' $(3).gem.map && \
# 	 $(call bam_fix_nh,$(3).tmp.bam,$(3).tmp21.bam) &&   : no need, NH flag is now fixed
#-o $(3).tmp.sam && \
	 # samtools  view -T $(reference_abspath) -bS $(3).tmp.sam > $(3).tmp.bam &&\

# TODO: bam_fix -NH -se?

############
# GEM split (not working)
#--min-split-size <min_split_size> ==tophat
# -m mismatches
# TODO: not fully working (junction format unknown, gtf2junction mentioned in the code is not distributed!?)
gems_map_params= --min-split-size 20  $(gems_map_options)
gems_index_params= -t $(max_threads) --max-memory unlimited --for-rna-mapper $(gems_index_options) 

ifeq ($(mapper),gems)
ifeq ($(mapper_splicing),no)
 $(error mapper_splicing option must be enabled when using gems)
endif
endif

define gems_index_filename=
$(1).gems.index
endef

define run_gems_index=
	irap_map.sh GEM gem-indexer -i $(1) -o $(call gems_index_filename,$(1)) $(gems_index_params) && touch $(call gems_index_filename,$(1))
endef

#  -J $(juncs_file_abspath)
define run_gems_map=
	 irap_map.sh GEM gem-rna-mapper $(gems_map_params) $(call gem_qual_option,$($(1)_qual))  -I $(index_files).gem -i $(2) -o $(3).gems && \
	 irap_map.sh GEM  gem-2-sam -i $(3).gems.map -I $(index_files) --emit-correct-flags -T $(max_threads) $(call gem2sam_pairing_param,$(1)) $(call gem_qual_option,$($(1)_qual)) -o /dev/fd/1 | sed  -e 's/\([a-zA-Z0-9]\+\) dna/\1/g'   | \
	 samtools  view -T $(reference_abspath) -bS - > $(3).tmp2.bam &&\
	 bam_fix_se_flag $(3).tmp2.bam - | \
	 samtools sort -m $(SAMTOOLS_SORT_MEM) - $(3).tmp   && \
	 mv $(3).tmp.bam $(3) && rm -rf $(3).tmp2* $(3).tmp.*

endef

###########################################################################
# Razers3 TODO
razers3_map_options+= --indels --max-hits $(max_hits)
#razers3_map_options= -of 4 --global-alignment -tc $(THREADS) -gn 1 --percent-identity 90
#razers3_index_options=  

# Missing: PE
#-of 2 (eland)
#--max-hits
#--unique
#-id (allow indels,default is off)
#-i %id
# razers3 does not like reads' names with spaces -> problem in sam
# (1) ref, index file, touch file
define razers3_index=
	$(call run_mapper,razers3,"razers3_index_$(1)","sed 's/ .*//g' $(1) > $(2).razers3.fa") && touch $(3)

endef
# --seed-length 16 (default)
# define run_razers3=
# 	 $(call run_mapper,razers3,"razers3_map_$(2)","sed 's/ /_/g' $(2) > $(2).razers3.fa &&\
# 	razers3 $(1).razers3.fa $(2).razers3.fa $(razers3_map_options) -o $(3).sam &&\
# 	samtools view -T $(1).razers3.fa -bS $(3).sam > $(3).tmp" ) && mv $(3).tmp $(3) 
# endef

#$(ref)/
#%.razers3.bam: %.fastq $(ref).razers3.index 
#	$(call run_razers3,$(ref),$<,$@)

# index
%.fa.razers3.index: %.fa
	$(call razers3_index,$<,$<,$@)

########################################################################
# BFAST
# (1) ref, index file, touch file
# -A 0/1
#

# WARNING!!! add / to the name of the temporary directory otherwise bfast index fails
bfast_A:=0
define bfast_index=
	rm -f $(1).*.1.1.bif
	rm -f $(1).tmp/*
	mkdir -p $(1).tmp
	$(call run_mapper,bfast,"bfast_index_$(1)","bfast fasta2brg $(bfast_index_options1) -f $(1) && bfast index $(bfast_index_options2) -T $(1).tmp/ -f $(1)") && touch $(3)
endef
bfast_index_options2= -w 8 -m '111101111011101111' -A $(bfast_A)
bfast_index_options1=  -A $(bfast_A)
bfast_map_options= -A $(bfast_A)

# refindex, fastq files, outfile.bam
# bfast easyalign will run bfast match, bfast localalign, and bfast postprocess
#define run_bfast=
#	$(call run_mapper,bfast,"bfast_map_$(2)","bfast easyalign $(bfast_map_options) -T $(1).tmp/ -f $(1) -r $(2) > $(3).sam && samtools  view -T $(ref) -bS $(3).sam > $(3).tmp ")  && mv $(3).tmp $(3)
#endef


#$(ref).bfast.index
#%.bfast.bam: %.fastq  %.fa.bfast.index
#	$(call run_bfast,$(ref).bfast.fa,$<,$@)

# index
%.fa.bfast.index: %.fa.bfast.fa
	$(call bfast_index,$<,$<,$@)

%.fa.bfast.fa: %.fa
	ln -s $< $@

#######################################################################
# Star

# Notes: 
#  1) The mapping quality MAPQ (column 5) is 255 for uniquely mapping reads, and int(-log10(1-1/Nmap)) for
#     multi-mapping reads. This scheme is same as the one used by Tophat and is compatible with Cufflinks.
#  2)  --outSAMstrandField intronMotif option, which will generate the XS strand on  all alignments that contain splice junctions. 
 
# outFilterMultimapNmax 10
# int: read alignments will be output only if the read maps fewer than thisvalue, otherwise no alignments will be output

# outFilterMismatchNoverLmax 0.3
#int: alignment will be output only if its ratio of mismatches to mapped length is less than this value

# splice junctions
# --alignIntronMax is the max intron size for a splice that occurs inside each mate, i.e. the one recorded as xxxN in the .sam CIGAR.

# --genomeLoad NoSharedMemory
# --genomeLoad LoadAndRemove
# note: --outSAMattributes Standard  --outSAMstrandField intronMotif (order is important?)
#   --outSAMstrandField intronMotif - unstranded data only
# stranded RNA-seq data: do not need to use any specific STAR options
star_map_params=  --genomeLoad NoSharedMemory --runThreadN $(max_threads) --outSAMunmapped Within --outFilterMultimapNmax $(max_hits) --outSAMattributes Standard  --outSAMstrandField intronMotif   $(star_map_options) 
star_index_params= --runThreadN $(max_threads) 


define star_qual_option=
	$(shell if [ "$(1)" != "33" ]; then echo "-q 'offset-64'"; else echo "-q offset-33"; fi)
endef

define star_index_filename=
$(1).star.index
endef

define star_index_dirname=
$(1).star
endef

# STAR requires one directory with a file per chr
# no tabs and spaces allowed in chr's names
# (1) ref, index file
define run_star_index=
	irap_fasta_split.pl $(1) $(call star_index_dirname,$(1)) && \
	sed -i 's/ .*//' $(call star_index_dirname,$(1))/*.fa && \
	irap_map.sh star star --runMode genomeGenerate $(star_index_params) --genomeDir $(call star_index_dirname,$(1)) --genomeFastaFiles $(call star_index_dirname,$(1))/*.fa && \
	touch $(call star_index_filename,$(1))
endef

# TODO: if splicing use 
ifeq ($(mapper_splicing),no)
 star_map_params+= --sjdbOverhang  0
else
# TODO: replace 20 by readlength-1 
 star_map_params+= --sjdbOverhang  20
endif

# TODO: fix
# Seg. fault 
# irap conf=test3.conf mapper=star quant_method=basic   exon_quant=y  transcript_quant=y de_method=none  stage2
define run_star_map=
	 irap_map.sh star star $(star_map_params) --genomeDir $(call star_index_dirname,$(file_indexed)) \
	--readFilesIn $(2) --outFileNamePrefix $(3) --sjdbFileChrStartEnd  $(juncs_file_abspath) &&\
	samtools view -bS $(3)Aligned.out.sam  > $(3).tmp.bam && \
	$(call bam_fix_nh,$(3).tmp.bam,$(3).tmp21.bam) &&  \
	bam_fix_se_flag $(3).tmp21.bam - | \
	samtools sort -m $(SAMTOOLS_SORT_MEM) - $(3).tmp  && \
	mv $(3).tmp.bam $(3) && rm -f $(3)Aligned.out.sam $(3).tmp21.bam
endef


#
# --sjdbFileChrStartEnd  -
# string: path to the file with genomic coordinates (chr <tab> start
# <tab> end) for the introns
# --sjdbOverhang  0
# int>=0: length of the donor/acceptor sequence on each side of the
# junctions, ideally = (mate_length - 1)
# if =0, splice junction database is not used


#######################################################################
# OSA

# Notes: 

osa_map_params=  --genomeLoad NoSharedMemory --runThreadN $(max_threads) --outSAMunmapped Within --outFilterMultimapNmax $(max_hits) --outSAMattributes Standard  --outSAMstrandField intronMotif --SearchNovelExonJunction True   $(osa_map_options) 
osa_index_params= --runThreadN $(max_threads)  $(osa_index_options)


define osa_index_filename=
$(1).osa.index
endef

define osa_index_dirname=
$(1).osa
endef

define osa_gene_model_name=
$(subst .gtf,,$(basename $(gtf_file)))
endef

define osa_ref_lib_name=
$(subst .fasta,,$(subst .fa,,$(shell basename $(1))))
endef

# when an error occurs the exit status is 0!? :(
#
define run_osa_index=
	irap_map.sh osa osa.exe --buildref `dirname $(1)` $(1) $(call osa_ref_lib_name,$(1)) &&\
	irap_map.sh osa osa.exe --buildgm $(call osa_index_dirname,$(1)) $(gtf_file_abspath) $(call osa_ref_lib_name,$(1)) $(call osa_gene_model_name) &&\
	ls $(shell dirname $(1))/ReferenceLibrary/$(call osa_ref_lib_name,$(1)).gindex1 &&\
	touch $(call osa_index_filename,$(1))
endef

# TODO: if splicing use 
ifeq ($(mapper_splicing),no)
 osa_map_params+= --sjdbOverhang  0
else
# TODO: replace 20 by readlength-1 
 osa_map_params+= --sjdbOverhang  20
endif

#
# generate a  configuration file and return the file name
#irap_gen_osa_conf.sh file1 threads expression outfile [file2 ins sd] 
# 1 - libname
# 2 - fastq_files
# 3 - outfile
define osa_conf_file=
	irap_gen_osa_conf.sh $(word 1,$(2)) $(max_threads) none $(3) $(if $(findstring $(1),$(pe)), $(word 2,$(2)) $($(1)_ins)  $($(1)_sd))
endef
# // Possible values: True, False. Default value=False"
# // possible values: None, TPM, RPKM, TPM_Transcript, RPKM_Transcript. Default value =None

# outputs BAM
# BAM includes unmapped

# osa requires the options to be passed in a conf. file...not user friendly :(
# gene model must be absolute path?! (01/2012)
# osa ignores the OutputName for the bam file :(
# - output name: remove the _1/2 and .fastq from the input filename
define run_osa_map=
	$(call osa_conf_file,$(1),$(2),`dirname $(3)`/$(1)_tmp) > $(3).conf && \
	 irap_map.sh osa osa.exe -alignrna `dirname $(file_indexed)` $(call osa_ref_lib_name,$(file_indexed)) \
	 $(call osa_index_dirname,$(file_indexed))/ReferenceLibrary/$(call osa_ref_lib_name,$(file_indexed))_GeneModels/$(call osa_gene_model_name) $(3).conf && \
	samtools sort -m $(SAMTOOLS_SORT_MEM) `dirname $(3)`/$(if $(findstring $(1),$(pe)),$(1)_f.bam,$(1).f.bam) $(3).tmp &&\
	mv $(3).tmp.bam $(3) && rm -f `dirname $(3)`/$(1)_f.bam  && rm -f $(3).conf
endef

#	 irap_map.sh osa osa.exe -alignrna `dirname $(file_indexed)` $(call osa_ref_lib_name,$(file_indexed))  `dirname $(file_indexed)`/$(call osa_ref_lib_name,$(file_indexed))_GeneModels/$(call osa_gene_model_name) $(3).conf && \ $(call osa_index_dirname,$(1)) $(gtf_file_abspath) $(call osa_ref_lib_name,$(1)) $(call osa_gene_model_name) &&\

######################################################
# Mapsplice

# global alignement
# --min-map-len 0
# max number of alignments
# -k 
# -m / --splice-mis <int> 	Maximum number of mismatches that are allowed in the first/last segment crossing a splice junction in the range of [0, 2]. Default is 1.
#(Maximum number of mismatches that are allowed in the middle segment crossing a splice junction is always fixed at 2.)
# --bam
# -o outputdir
# --threads thread
# -qual-scale phred33,phred64,solexa64
#-I / --max-intron <int>
#-i / --min-intron <int>
# Fusion
#--gene-gtf <string> 
#--fusion | --fusion-non-canonical  	--fusion: Search for canonical and semi-canonical fusion junctions.
#--fusion-non-canonical: Search for canonical, semi-canonical, and non-canonical fusion junctions. 
# mapsplice fails with smaller values
mapsplice_min_intron=11
mapsplice_map_params=   --min-intron $(mapsplice_min_intron)  $(mapsplice_map_options) --max-intron 10000 --max-hits $(max_hits) --del $(min_intron_len) --ins $(min_intron_len)
mapsplice_index_params= --offrate 3 

define mapsplice_index_filename=
$(1).mapsplice.index
endef

define mapsplice_index_prefix=
$(1).mapsplice
endef

define mapsplice_file_params=
	$(if $(findstring $(1),$(pe)), -1 $(word 1,$(2)) -2 $(word 2,$(2)) $(call qual_mapsplice,$($(1)_qual)), -1 $(2) $(call qual_mapsplice,$($(1)_qual)))
endef

#
define run_mapsplice_index=
	sed 's/ .*//g' $(1) > $(call mapsplice_index_prefix,$(1)).fa && \
	irap_fasta_split.pl $(call mapsplice_index_prefix,$(1)).fa $(call mapsplice_index_prefix,$(1)) && \
	irap_map.sh mapsplice  bowtie-build $(mapsplice_index_params) $(call mapsplice_index_prefix,$(1)).fa $(call mapsplice_index_prefix,$(1)) && \
	touch $(call mapsplice_index_filename,$(1))
endef


#-c The directory containing the sequence files of reference genome. All sequence files are required to:
# BAM file does not contain the NH flag and the mate information is not ok (htseq fails)
define run_mapsplice_map=
	 irap_map.sh mapsplice python $(IRAP_DIR)/bin/mapsplice/mapsplice.py  $(mapsplice_map_params) --threads	 $(max_threads) --bam -o  $(call lib2bam_folder,$(1))$(1) -c $(call mapsplice_index_prefix,$(file_indexed)) -x $(call mapsplice_index_prefix,$(file_indexed)) $(call mapsplice_file_params,$(1),$(2)) &&\
	samtools fixmate  $(call lib2bam_folder,$(1))$(1)/alignments.bam $(3).fix && \
	$(call bam_fix_nh,$(3).fix,-) | \
	samtools sort -m $(SAMTOOLS_SORT_MEM) - $(3).tmp && \
	mv $(3).tmp.bam $(3) && rm -f $(3).fix
endef

########################################################################
ifneq ($(mapper),tophat2)
define $(mapper)_index_filenames=
	$(call $(mapper)_index_filename,$(1),$(1)) 
endef
else
define tophat2_index_filenames=
	$(call tophat2_index_filename,$(1),$(1)) 	$(call tophat2_trans_index_filename,$(1),$(1))
endef
endif
