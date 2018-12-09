#; -*- mode: Makefile;-*-
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

# QC statistics
# Generation of a single file with the stats collected from different stages of the analysis


######################
# do_quality_filtering_and_report(1,2,3,4)
# 1 - libname
# 2 - input folder
# 3 - output folder
# 4 - extra options
#s
define do_quality_filtering_and_report=
	irap_fastq_qc -j 1 $(read_qual_filter_common_params) input_dir=$(2) read_size=$($(1)_rs)   qual=$($(1)_qual) f="$(notdir $($(1)))" out_prefix=$(1) is_pe=$(call is_pe_lib,$(1)) out_dir=$(3)  $(4) $(call get_opt_barcode_params,$(1)) 
endef
#	irap_fastq_qc $(read_qual_filter_common_params) data_dir=$(raw_data_dir)$($(1)_dir) read_size=$($(1)_rs)   qual=$($(1)_qual) f="$($(1))" out_prefix=$(1) is_pe=$(#call is_pe_lib,$(1)) out_dir=$(auxdata_toplevel_folder)/$($(1)_dir)  $(2)


################################################################################
# Filtering
################################################################################
# reuse the data filtered between experiments?


ifdef min_read_len
	read_qual_filter_common_params+=min_len=$(min_read_len)
endif


####################################################################
# Preprocessing of the input files (.fastq/.bam), QC, filtering
define rep_is_on=
$(if $(filter off,$(qc)),,$(1))
endef
#$(info $(call rep_is_on))
intermediate_targets=
define make-pe-qc-rule=
$(call lib2filt_folder,$(1))$(1)_1.f.fastq%gz $(call lib2filt_folder,$(1))$(1)_2.f.fastq%gz $(call lib2filt_folder,$(1))$(1)_1.f.fastqc%tsv $(call lib2filt_folder,$(1))$(1)_2.f.fastqc%tsv $(call lib2filt_folder,$(1))$(1)_1.f%csv $(call lib2filt_folder,$(1))$(1)_2.f%csv $(call rep_is_on,$(call lib2filt_folder,$(1))$(1)_1.f.fastqc%zip) $(call rep_is_on,$(call lib2filt_folder,$(1))$(1)_2.f.fastqc%zip): $(raw_data_dir)$($(1)_dir)/$(notdir $(word 1,$($(1))))  $(raw_data_dir)$($(1)_dir)/$(notdir $(word 2,$($(1))))
	$$(call p_info,Filtering $(call fix_libname,$(1)))
	$$(call do_quality_filtering_and_report,$(call fix_libname,$(1)),$(raw_data_dir)$($(1)_dir),$$(@D),) || (rm -f $(call lib2filt_folder,$(1))$(1)_1.f.fastq.gz && exit 1)
endef

define make-se-qc-rule=

$(call lib2filt_folder,$(1))$(1).f.fastq%gz $(call lib2filt_folder,$(1))$(1).f.fastqc%tsv $(call lib2filt_folder,$(1))$(1).f%csv $(call rep_is_on,$(call lib2filt_folder,$(1))$(1).f.fastqc%zip): $(raw_data_dir)$($(1)_dir)/$(notdir $($(1)))
	$$(call p_info,Filtering $(call fix_libname,$(1)))
	$$(call do_quality_filtering_and_report,$(1),$(raw_data_dir)$($(1)_dir),$$(@D),) || (rm -f $(call lib2filt_folder,$(1))$(1).f.fastq.gz && exit 1)

endef

ifneq ($(deps_check),nocheck)
# rules for SE libraries
$(foreach l,$(se),$(eval $(call make-se-qc-rule,$(l))))
# rules for PE libraries
$(foreach l,$(pe),$(eval $(call make-pe-qc-rule,$(l))))


qc_files:=$(foreach p,$(se),$(call lib2filt_folder,$(p))$(p).f.fastq.gz) $(foreach p,$(pe),$(call lib2filt_folder,$(p))$(p)_1.f.fastq.gz)
STAGE1_OUT_FILES+=$(qc_files)
STAGE1_TARGETS+=$(qc_files)
STAGE1_S_TARGETS+=$(qc_files)

endif

CLEANUP_TARGETS+= clean_quality_filtering_and_report
phony_targets+= clean_quality_filtering_and_report
# Cleanup
clean_quality_filtering_and_report:
	$(foreach p,$(se) $(pe),$(call do_quality_filtering_and_report,$(p),$(raw_data_dir)$($(p)_dir)/,$(dir $(call lib2filt_folder,$(p))),clean))


#########################################################################
# a single file with the mapping stats
MAPPING_REPORT_PRE_STATS=
ifneq ($(mapper),none)
ifneq ($(deps_check),nocheck)
ifneq ($(transcriptome_alignment_only),n)
# only get stats (gene, exon, ...) when the alignments were made against the genome
MAPPING_REPORT_PRE_STATS:=$(foreach s,$(se),$(call lib2bam_folder,$(s))$(s).se.hits.bam.stats.csv) $(foreach s,$(pe),$(call lib2bam_folder,$(s))$(s).pe.hits.bam.stats.csv  )
endif
endif

LIBS_QC_REQS=$(mapper_toplevel_folder)/stats_raw.tsv $(mapper_toplevel_folder)/stats_perc.tsv  

ifeq ($(transcriptome_alignment_only),n)
LIBS_QC_REQS+=$(mapper_toplevel_folder)/featstats_raw.tsv  $(mapper_toplevel_folder)/featstats_perc.tsv   $(mapper_toplevel_folder)/genestats_raw.tsv
$(mapper_toplevel_folder)/libs_qc.tsv: $(LIBS_QC_REQS)
	irap_append2tsv --in "$(mapper_toplevel_folder)/stats_raw.tsv $(mapper_toplevel_folder)/featstats_raw.tsv $(mapper_toplevel_folder)/genestats_raw.tsv" --exclude_aggr  --cols_not_sorted --out $@.1.tmp &&\
	irap_append2tsv --in "$(mapper_toplevel_folder)/stats_perc.tsv $(mapper_toplevel_folder)/featstats_perc.tsv $(mapper_toplevel_folder)/genestats_perc.tsv" --exclude_aggr --add_row_suffix "_perc" --cols_not_sorted --out $@.2.tmp && \
	irap_append2tsv --in "$@.1.tmp $@.2.tmp" --exclude_aggr --transpose --out $@.tmp && mv $@.tmp $@ &&\
	rm -f $@.1.tmp $@.2.tmp
else
## alignment to transcripts
$(mapper_toplevel_folder)/libs_qc.tsv: $(LIBS_QC_REQS)
	irap_append2tsv --in "$(mapper_toplevel_folder)/stats_raw.tsv " --exclude_aggr  --cols_not_sorted --out $@.1.tmp &&\
	irap_append2tsv --in "$(mapper_toplevel_folder)/stats_perc.tsv  " --exclude_aggr --add_row_suffix "_perc" --cols_not_sorted --out $@.2.tmp && \
	irap_append2tsv --in "$@.1.tmp $@.2.tmp" --exclude_aggr --transpose --out $@.tmp && mv $@.tmp $@ &&\
	rm -f $@.1.tmp $@.2.tmp

endif
## endif($(transcriptome_alignment_only),n)

else
# empty file
$(mapper_toplevel_folder)/libs_qc.tsv:
	touch $@

## mapper
endif

## the reference in the BAMs generated with kallisto is not the mapper
ifneq ($(mapper),kallisto)
ifneq ($(mapper),none)
ifneq ($(deps_check),nocheck)
MAPPING_REPORT_PRE_STATS+=$(foreach s,$(se), $(call lib2bam_folder,$(s))$(s).se.hits.bam.gene.stats  $(call lib2bam_folder,$(s))$(s).se.hits.bam.stats) $(foreach s,$(pe), $(call lib2bam_folder,$(s))$(s).pe.hits.bam.gene.stats $(call lib2bam_folder,$(s))$(s).pe.hits.bam.stats)
endif
endif	
endif
WAVE3_s_TARGETS+=$(MAPPING_REPORT_PRE_STATS)
ifneq ($(mapper),none)
WAVE3_TARGETS+=$(mapper_toplevel_folder)/libs_qc.tsv
endif
WAVE4_TARGETS+=

# merge into a single file the statistics collected from the BAMs 
$(mapper_toplevel_folder)/stats_raw%tsv $(mapper_toplevel_folder)/stats_perc%tsv: $(foreach s,$(se),$(call lib2bam_folder,$(s))$(s).se.hits.bam) $(foreach p,$(pe),$(call lib2bam_folder,$(p))$(p).pe.hits.bam)
	$(call pass_args_stdin,irap_bams2tsv,$(mapper_toplevel_folder)/stats_raw.tsv, --pe "$(call remove_spaces,$(foreach p,$(pe),;$(call lib2bam_folder,$(p))$(p).pe.hits.bam))" --se "$(call remove_spaces,$(foreach s,$(se),;$(call lib2bam_folder,$(s))$(s).se.hits.bam))"  --pe_labels "$(call remove_spaces,$(foreach p,$(pe),;$(p)))" --se_labels "$(call remove_spaces,$(foreach s,$(se),;$(s)))" --out $(mapper_toplevel_folder)/$(mapper)) && mv $(mapper_toplevel_folder)/$(mapper)_mapping_stats_raw.tsv $(mapper_toplevel_folder)/stats_raw.tsv && mv $(mapper_toplevel_folder)/$(mapper)_mapping_stats_perc.tsv $(mapper_toplevel_folder)/stats_perc.tsv

#
$(mapper_toplevel_folder)/featstats_raw%tsv $(mapper_toplevel_folder)/featstats_perc%tsv:  $(foreach s,$(se),$(call lib2bam_folder,$(s))$(s).se.hits.bam.stats) $(foreach p,$(pe),$(call lib2bam_folder,$(p))$(p).pe.hits.bam.stats)
	$(call pass_args_stdin,merge_featstats,$(mapper_toplevel_folder)/featstats_raw.tsv, --header --stats "$(call remove_spaces,$(foreach p,$(pe),;$(call lib2bam_folder,$(p))$(p).pe.hits.bam.stats))$(call remove_spaces, $(foreach s,$(se),;$(call lib2bam_folder,$(s))$(s).se.hits.bam.stats))"  --labels "$(call remove_spaces,$(foreach p,$(pe) $(se),;$(p)))"  --out $(mapper_toplevel_folder)/$(mapper).ftmp) && mv $(mapper_toplevel_folder)/$(mapper).ftmp_featstats_raw.tsv $(mapper_toplevel_folder)/featstats_raw.tsv && mv $(mapper_toplevel_folder)/$(mapper).ftmp_featstats_perc.tsv $(mapper_toplevel_folder)/featstats_perc.tsv

$(mapper_toplevel_folder)/genestats_raw%tsv $(mapper_toplevel_folder)/genestats_perc%tsv:   $(foreach s,$(se),$(call lib2bam_folder,$(s))$(s).se.hits.bam.gene.stats) $(foreach p,$(pe),$(call lib2bam_folder,$(p))$(p).pe.hits.bam.gene.stats)
	$(call pass_args_stdin,merge_featstats,$(mapper_toplevel_folder)/genestats_raw.tsv, --stats "$(call remove_spaces,$(foreach p,$(pe),;$(call lib2bam_folder,$(p))$(p).pe.hits.bam.gene.stats))$(call remove_spaces,$(foreach p,$(se),;$(call lib2bam_folder,$(p))$(p).se.hits.bam.gene.stats))"  --labels "$(call remove_spaces,$(foreach p,$(pe) $(se),;$(p)))"  --out $(mapper_toplevel_folder)/$(mapper).gtmp) && mv $(mapper_toplevel_folder)/$(mapper).gtmp_featstats_perc.tsv $(mapper_toplevel_folder)/genestats_perc.tsv && mv $(mapper_toplevel_folder)/$(mapper).gtmp_featstats_raw.tsv $(mapper_toplevel_folder)/genestats_raw.tsv 

################
#
# statistics per bam file
#
ifeq ($(transcriptome_alignment_only),n)
# 
%.bam.gff3: %.bam $(gff3_file_abspath).filt.gff3 $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.txt
	bedtools coverage -counts -sorted  -a $(gff3_file_abspath).filt.gff3 -b $< -g  $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.txt > $@.tmp  &&\
	mv $@.tmp $@

%.bam.stats: %.bam.gff3 
	mapping_feature_stats --in $< --out $@.tmp -c "`basename $*`" && mv $@.tmp $@

%.bam.stats.csv: %.bam 
	irapBAM2stats bam=$< || ( rm -f $@ && exit 1)

%.bam.gene.stats: %.bam $(auxdata_toplevel_folder)/$(reference_basename).exons.bed $(auxdata_toplevel_folder)/$(reference_basename).introns.bed $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.txt
	set -o pipefail && echo -n "Exons	" > $@.tmp &&\
	bedtools intersect -sorted -g $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.txt -abam $<  -b $(auxdata_toplevel_folder)/$(reference_basename).exons.bed |samtools view -c - >> $@.tmp && echo >> $@ &&\
	echo -n "Introns	" >> $@.tmp &&\
	bedtools intersect -sorted -g $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.txt -abam $<  -b $(auxdata_toplevel_folder)/$(reference_basename).introns.bed |samtools view -c - >> $@.tmp && echo >> $@ && \
	expr `wc -l $@.tmp | cut -f 1 -d\ ` == 2 && \
	mv $@.tmp $@

else

%.bam.stats: %.bam
	echo WIP
	touch $@

endif

# bed files required to get some extra stats
# exons.bed
$(auxdata_toplevel_folder)/$(reference_basename).exons.bed: $(gff3_file_abspath).filt.gff3 $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.bed
	set -o pipefail && cat $< | awk 'BEGIN{OFS="\t";} $$3=="exon" {print $$1,$$4,$$5,$$6,$$6,$$7}' | sort -u| bedtools sort -faidx $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.bed -i /dev/stdin > $@.tmp.bed && \
	bedtools sort -faidx $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.bed -i $@.tmp.bed > $@.tmp && \
	mv $@.tmp $@ && rm -f $@.tmp.bed

# genes.bed
$(auxdata_toplevel_folder)/$(reference_basename).genes.bed: $(gff3_file_abspath).filt.gff3 $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.bed
	set -o pipefail &&  cat $< | awk 'BEGIN{OFS="\t";} $$3=="gene" {print $$1,$$4,$$5}' | sort -k1,1 -k2,2n -u  > $@.tmp.bed &&\
	bedtools merge -i $@.tmp.bed  | bedtools sort -faidx $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.bed -i /dev/stdin  > $@.tmp && \
	mv $@.tmp $@ && rm -f $@.tmp.bed

# 
$(auxdata_toplevel_folder)/$(reference_basename).genes.bed6: $(gtf_file_abspath)
	set -o pipefail &&  sed -E 's/[^\t]*gene_id "([^;]+)".*$$/\1/' $< | awk 'BEGIN{OFS="\t";} $$3=="exon" {print $$1,$$4,$$5,$$9,$$6,$$7}' |  sort -u| bedtools sort -faidx $(auxdata_toplevel_folder)/$(reference_basename).chr_sizes.sorted.bed -i /dev/stdin > $@.tmp.bed &&\
	mv $@.tmp.bed $@ && rm -f $@.tmp.bed


# introns
$(auxdata_toplevel_folder)/$(reference_basename).introns.bed: $(auxdata_toplevel_folder)/$(reference_basename).genes.bed $(auxdata_toplevel_folder)/$(reference_basename).exons.bed
	 sort -k1,1 -k2,2n -u $(auxdata_toplevel_folder)/$(reference_basename).exons.bed > $@.tmp.bed && bedtools subtract -sorted -a $< -b $@.tmp.bed > $@.tmp && if [ `wc -l $@.tmp |cut -f 1 -d\ ` == 0 ]; then echo -e 'dummy_entry\t1\t1' > $@.tmp; fi && mv $@.tmp $@


## transcripts
$(auxdata_toplevel_folder)/$(reference_basename).transcripts.bed6:  $(gtf_file_abspath)
	set -o pipefail && grep -E "(exon)" $< | sed -E 's/[^\t]*transcript_id "([^;]+)".*$$/\1/'|awk 'BEGIN{OFS="\t";} {print $$1,$$4,$$5,$$9,$$6,$$7}' |  sort -u| bedtools sort -i /dev/stdin > $@.tmp.bed &&\
	mv $@.tmp.bed $@ && rm -f $@.tmp.bed


#
ifneq ($(mapper),none)
STAGE4_OUT_FILES+=$(mapper_toplevel_folder)/libs_qc.tsv
STAGE4_TARGETS+=$(mapper_toplevel_folder))/libs_qc.tsv
endif
######################################################################


######################################################################
# qc=none|on|off
STAGE1_OUT_FILES+=$(qc_toplevel_folder)/fastq_qc_report.tsv
STAGE1_TARGETS+=$(qc_toplevel_folder)/fastq_qc_report.tsv

STAGE1_OUT_FILES+=$(qc_toplevel_folder)/qc.tsv
STAGE1_TARGETS+=$(qc_toplevel_folder)/qc.tsv

ifeq ($(qc),off)
$(qc_toplevel_folder)/qc.html $(qc_toplevel_folder)/qc.tsv: 
	touch $@

# empty file
FASTQC_REPORT_FILES=
$(qc_toplevel_folder)/fastq_qc_report.tsv:
	touch $@

## nothing
print_qc_dirs_files:

## end qc=off
else

ifeq ($(HUGE_NUM_LIBS),0)
## qc=on |report
$(qc_toplevel_folder)/qc%html $(qc_toplevel_folder)/qc%tsv: $(conf) $(call must_exist,$(auxdata_toplevel_folder)/)  $(qc_toplevel_folder)/fastq_qc_report.tsv $(qc_toplevel_folder)/qc_report.csv $(lib_info)
	irap_report_qc $(IRAP_REPORT_MAIN_OPTIONS) --conf $(conf) --out_dir $(qc_toplevel_folder) --qc_dir $(qc_toplevel_folder) --css $(CSS_FILE) $(call get_lib_info_option) || ( rm -f $(qc_toplevel_folder)/qc.tsv && exit 1)

else
# too many libs - this code needs to be optimized to support more than 20k libs
$(qc_toplevel_folder)/qc%html $(qc_toplevel_folder)/qc%tsv:
	touch $(qc_toplevel_folder)/qc.html $(qc_toplevel_folder)/qc.tsv
endif

ifneq ($(qc),off)
ifneq ($(deps_check),nocheck)
FASTQC_REPORT_FILES=$(foreach p,$(se),$(call lib2filt_folder,$(p))$(p).f.fastqc.tsv) $(foreach p,$(pe),$(call lib2filt_folder,$(p))$(p)_1.f.fastqc.tsv $(call lib2filt_folder,$(p))$(p)_2.f.fastqc.tsv)

QC_CSV_FILES=$(foreach p,$(se),$(call lib2filt_folder,$(p))$(p).f.csv) $(foreach p,$(pe),$(call lib2filt_folder,$(p))$(p)_1.f.csv $(call lib2filt_folder,$(p))$(p)_2.f.csv)

ZIP_FILES=$(foreach p,$(se),$(qc_toplevel_folder)/$($(p)_dir)$(p).f.fastqc.zip) $(foreach p,$(pe),$(qc_toplevel_folder)/$($(p)_dir)$(p)_1.f.fastqc.zip $(qc_toplevel_folder)/$($(p)_dir)$(p)_2.f.fastqc.zip)
else
ZIP_FILES=
QC_CSV_FILES=
FASTQC_REPORT_FILES=
endif
else
ZIP_FILES=
QC_CSV_FILES=
FASTQC_REPORT_FILES=
endif

STAGE1_S_TARGETS+=$(QC_CSV_FILES) $(FASTQC_REPORT_FILES) $(ZIP_FILES)

## work around the limitation on the number of arguments
print_qc_dirs_files:  
	@echo -n $(file > /dev/stdout, $(QC_CSV_FILES) $(FASTQC_REPORT_FILES) $(ZIP_FILES))


$(qc_toplevel_folder)/fastq_qc_report.tsv:  $(FASTQC_REPORT_FILES)
	$(call pass_args_stdin,irap_merge2tsv,$@.tmp, --in "$(subst $(space),;,$^)"  --out $@.tmp) && mv $@.tmp $@

$(qc_toplevel_folder)/qc_report.csv:  $(QC_CSV_FILES)
	$(call stdin_cat,$^,$@.tmp) && mv $@.tmp $@

## qc=on |report
endif

ifeq ($(deps_check),nocheck)
$(call p_info,irap_qc_stats.mk loaded)
endif





