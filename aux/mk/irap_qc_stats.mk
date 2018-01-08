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
	irap_fastq_qc $(read_qual_filter_common_params) input_dir=$(2) read_size=$($(1)_rs)   qual=$($(1)_qual) f="$(notdir $($(1)))" out_prefix=$(1) is_pe=$(call is_pe_lib,$(1)) out_dir=$(3)  $(4) $(call get_opt_barcode_params,$(1))
endef
#	irap_fastq_qc $(read_qual_filter_common_params) data_dir=$(raw_data_dir)$($(1)_dir) read_size=$($(1)_rs)   qual=$($(1)_qual) f="$($(1))" out_prefix=$(1) is_pe=$(#call is_pe_lib,$(1)) out_dir=$(name)/data/$($(1)_dir)  $(2)


################################################################################
# Filtering
################################################################################
# reuse the data filtered between experiments?


ifdef min_read_len
	read_qual_filter_common_params+=min_len=$(min_read_len)
endif


####################################################################
# Preprocessing of the input files (.fastq/.bam), QC, filtering

define make-pe-qc-rule=
$(call lib2filt_folder,$(1))$(1)_1.f.fastq.gz $(call lib2filt_folder,$(1))$(1)_2.f.fastq.gz $(call lib2filt_folder,$(1))$(1)_1.f.fastqc.tsv: $(raw_data_dir)$($(1)_dir)/$(notdir $(word 1,$($(1))))  $(raw_data_dir)$($(1)_dir)/$(notdir $(word 2,$($(1))))
	$$(call p_info,Filtering $(call fix_libname,$(1)))
	$$(call do_quality_filtering_and_report,$(call fix_libname,$(1)),$(raw_data_dir)$($(1)_dir),$$(@D),) || (rm -f $(call lib2filt_folder,$(1))$(1)_1.f.fastq.gz && exit 1)
endef

define make-se-qc-rule=
$(call lib2filt_folder,$(1))$(1).f.fastq.gz $(call lib2filt_folder,$(1))$(1).f.fastqc.tsv: $(raw_data_dir)$($(1)_dir)/$(notdir $($(1)))
	$$(call p_info,Filtering $(call fix_libname,$(1)))
	$$(call do_quality_filtering_and_report,$(1),$(raw_data_dir)$($(1)_dir),$$(@D),) || (rm -f $(call lib2filt_folder,$(1))$(1).f.fastq.gz && exit 1)

endef
# rules for SE libraries
$(foreach l,$(se),$(eval $(call make-se-qc-rule,$(l))))
# rules for PE libraries
$(foreach l,$(pe),$(eval $(call make-pe-qc-rule,$(l))))


qc_files:=$(foreach p,$(se),$(call lib2filt_folder,$(p))$(p).f.fastq.gz) $(foreach p,$(pe),$(call lib2filt_folder,$(p))$(p)_1.f.fastq.gz)
STAGE1_OUT_FILES+=$(qc_files)
STAGE1_TARGETS+=$(qc_files)


CLEANUP_TARGETS+= clean_quality_filtering_and_report
phony_targets+= clean_quality_filtering_and_report
# Cleanup
clean_quality_filtering_and_report:
	$(foreach p,$(se) $(pe),$(call do_quality_filtering_and_report,$(p),$(raw_data_dir)$($(p)_dir)/,$(dir $(call lib2filt_folder,$(p))),clean))


# TODO: deprecated/rm
#print_qc_dirs_files:
#	echo	$(foreach l,$(se) $(pe),$(name)/report/riq/$($(l)_dir) )


#########################################################################
# a single file with the mapping stats
$(name)/$(mapper)/libs_qc.tsv: $(name)/$(mapper)/stats_raw.tsv $(name)/$(mapper)/stats_perc.tsv  $(name)/$(mapper)/featstats_raw.tsv   $(name)/$(mapper)/genestats_raw.tsv 
	irap_append2tsv --in "$(name)/$(mapper)/stats_raw.tsv $(name)/$(mapper)/featstats_raw.tsv $(name)/$(mapper)/genestats_raw.tsv" --exclude_aggr  --cols_not_sorted --out $@.1.tmp &&\
	irap_append2tsv --in "$(name)/$(mapper)/stats_perc.tsv $(name)/$(mapper)/featstats_perc.tsv $(name)/$(mapper)/genestats_perc.tsv" --exclude_aggr --add_row_suffix "_perc" --cols_not_sorted --out $@.2.tmp && \
	irap_append2tsv --in "$@.1.tmp $@.2.tmp" --exclude_aggr --transpose --out $@.tmp && mv $@.tmp $@ &&\
	rm -f $@.1.tmp $@.2.tmp


MAPPING_REPORT_PRE_STATS:=$(foreach s,$(se),$(call lib2bam_folder,$(s))$(s).se.hits.bam.stats.csv $(call lib2bam_folder,$(s))$(s).se.hits.bam.stats $(call lib2bam_folder,$(s))$(s).se.hits.bam.gene.stats) $(foreach s,$(pe),$(call lib2bam_folder,$(s))$(s).pe.hits.bam.stats.csv $(call lib2bam_folder,$(s))$(s).pe.hits.bam.stats $(call lib2bam_folder,$(s))$(s).pe.hits.bam.gene.stats)


WAVE3_s_TARGETS+=$(MAPPING_REPORT_PRE_STATS)
WAVE3_TARGETS+=$(name)/$(mapper)/libs_qc.tsv
WAVE4_TARGETS+=

# merge into a single file the statistics collected from the BAMs 
$(name)/$(mapper)/stats_raw.tsv $(name)/$(mapper)/stats_perc.tsv: $(foreach s,$(se),$(call lib2bam_folder,$(s))$(s).se.hits.bam) $(foreach p,$(pe),$(call lib2bam_folder,$(p))$(p).pe.hits.bam)
	$(call pass_args_stdin,irap_bams2tsv,$(name)/$(mapper)/stats_raw.tsv, --pe "$(call remove_spaces,$(foreach p,$(pe),;$(call lib2bam_folder,$(p))$(p).pe.hits.bam))" --se "$(call remove_spaces,$(foreach s,$(se),;$(call lib2bam_folder,$(s))$(s).se.hits.bam))"  --pe_labels "$(call remove_spaces,$(foreach p,$(pe),;$(p)))" --se_labels "$(call remove_spaces,$(foreach s,$(se),;$(s)))" --out $(name)/$(mapper)/$(mapper)) && mv $(name)/$(mapper)/$(mapper)_mapping_stats_raw.tsv $(name)/$(mapper)/stats_raw.tsv && mv $(name)/$(mapper)/$(mapper)_mapping_stats_perc.tsv $(name)/$(mapper)/stats_perc.tsv

#
$(name)/$(mapper)/featstats_raw.tsv $(name)/$(mapper)/featstats_perc.tsv:  $(foreach s,$(se),$(call lib2bam_folder,$(s))$(s).se.hits.bam.stats) $(foreach p,$(pe),$(call lib2bam_folder,$(p))$(p).pe.hits.bam.stats)
	$(call pass_args_stdin,merge_featstats,$(name)/$(mapper)/featstats_raw.tsv, --header --stats "$(call remove_spaces,$(foreach p,$(pe),;$(call lib2bam_folder,$(p))$(p).pe.hits.bam.stats))$(call remove_spaces, $(foreach s,$(se),;$(call lib2bam_folder,$(s))$(s).se.hits.bam.stats))"  --labels "$(call remove_spaces,$(foreach p,$(pe) $(se),;$(p)))"  --out $(name)/(mapper)/$(mapper).ftmp) && mv $(name)/(mapper)/$(mapper).ftmp_featstats_raw.tsv $@ && mv $(name)/(mapper)/$(mapper).ftmp_featstats_perc.tsv $(name)/(mapper)/$(mapper)/featstats_perc.tsv

$(name)/$(mapper)/genestats_raw.tsv $(name)/$(mapper)/genestats_perc.tsv:   $(foreach s,$(se),$(call lib2bam_folder,$(s))$(s).se.hits.bam.gene.stats) $(foreach p,$(pe),$(call lib2bam_folder,$(p))$(p).pe.hits.bam.gene.stats)
	$(call pass_args_stdin,merge_featstats,$(name)/$(mapper)/genestats_raw.tsv, --stats "$(call remove_spaces,$(foreach p,$(pe),;$(call lib2bam_folder,$(p))$(p).pe.hits.bam.gene.stats))$(call remove_spaces,$(foreach p,$(se),;$(call lib2bam_folder,$(p))$(p).se.hits.bam.gene.stats))"  --labels "$(call remove_spaces,$(foreach p,$(pe) $(se),;$(p)))"  --out $(name)/$(mapper)/$(mapper).gtmp) && mv $(name)/$(mapper)/$(mapper).gtmp_featstats_raw.tsv $(name)/$(mapper)/genestats_raw.tsv && mv $(name)/$(mapper)/$(mapper).gtmp_featstats_perc.tsv $(name)/$(mapper)/genestats_perc.tsv

################

# statistics per bam file
%.bam.gff3: %.bam $(gff3_file_abspath).filt.gff3 $(name)/data/$(reference_basename).chr_sizes.sorted.txt
	bedtools coverage -counts -sorted  -a $(gff3_file_abspath).filt.gff3 -b $< -g  $(name)/data/$(reference_basename).chr_sizes.sorted.txt > $@.tmp  &&\
	mv $@.tmp $@

%.bam.stats: %.bam.gff3 
	mapping_feature_stats --in $< --out $@.tmp -c "`basename $*`" && mv $@.tmp $@

%.bam.stats.csv: %.bam 
	irapBAM2stats bam=$< || ( rm -f $@ && exit 1)

%.bam.gene.stats: %.bam $(name)/data/$(reference_basename).exons.bed $(name)/data/$(reference_basename).introns.bed $(name)/data/$(reference_basename).chr_sizes.sorted.txt
	echo -n "Exons	" > $@.tmp &&\
	bedtools intersect -sorted -g $(name)/data/$(reference_basename).chr_sizes.sorted.txt -abam $<  -b $(name)/data/$(reference_basename).exons.bed |samtools view -c - >> $@.tmp && echo >> $@ &&\
	echo -n "Introns	" >> $@.tmp &&\
	bedtools intersect -sorted -g $(name)/data/$(reference_basename).chr_sizes.sorted.txt -abam $<  -b $(name)/data/$(reference_basename).introns.bed |samtools view -c - >> $@.tmp && echo >> $@ && \
	expr `wc -l $@.tmp | cut -f 1 -d\ ` == 2 && \
	mv $@.tmp $@

# bed files required to get some extra stats
# exons.bed
$(name)/data/$(reference_basename).exons.bed: $(gff3_file_abspath).filt.gff3 $(name)/data/$(reference_basename).chr_sizes.sorted.bed
	cat $< | awk 'BEGIN{OFS="\t";} $$3=="exon" {print $$1,$$4,$$5,$$6,$$6,$$7}' | sort -u| bedtools sort -i /dev/stdin > $@.tmp.bed && \
	bedtools merge -i $@.tmp.bed | bedtools sort -faidx $(name)/data/$(reference_basename).chr_sizes.sorted.bed -i /dev/stdin > $@.tmp && \
	mv $@.tmp $@ && rm -f $@.tmp.bed

# genes.bed
$(name)/data/$(reference_basename).genes.bed: $(gff3_file_abspath).filt.gff3 $(name)/data/$(reference_basename).chr_sizes.sorted.bed
	cat $< | awk 'BEGIN{OFS="\t";} $$3=="gene" {print $$1,$$4,$$5}' |  sort -u| bedtools sort -i /dev/stdin > $@.tmp.bed &&\
	bedtools merge -i $@.tmp.bed  | bedtools sort -faidx $(name)/data/$(reference_basename).chr_sizes.sorted.bed -i /dev/stdin  > $@.tmp && \
	mv $@.tmp $@ && rm -f $@.tmp.bed

# 
$(name)/data/$(reference_basename).genes.bed6: $(gtf_file_abspath)
	sed -E 's/[^\t]*gene_id "([^;]+)".*$$/\1/' $< | awk 'BEGIN{OFS="\t";} $$3=="exon" {print $$1,$$4,$$5,$$9,$$6,$$7}' |  sort -u| bedtools sort -i /dev/stdin > $@.tmp.bed &&\
	mv $@.tmp.bed $@ && rm -f $@.tmp.bed


# introns
$(name)/data/$(reference_basename).introns.bed: $(name)/data/$(reference_basename).genes.bed $(name)/data/$(reference_basename).exons.bed
	bedtools subtract -sorted -a $< -b $(name)/data/$(reference_basename).exons.bed > $@.tmp && if [ `wc -l $@.tmp |cut -f 1 -d\ ` == 0 ]; then echo -e 'dummy_entry\t1\t1' > $@.tmp; fi && mv $@.tmp $@


## transcripts
$(name)/data/$(reference_basename).transcripts.bed6:  $(gtf_file_abspath)
	grep -E "(exon)" $< | sed -E 's/[^\t]*transcript_id "([^;]+)".*$$/\1/'|awk 'BEGIN{OFS="\t";} {print $$1,$$4,$$5,$$9,$$6,$$7}' |  sort -u| bedtools sort -i /dev/stdin > $@.tmp.bed &&\
	mv $@.tmp.bed $@ && rm -f $@.tmp.bed


#
STAGE4_OUT_FILES+=$(name)/$(mapper)/libs_qc.tsv
STAGE4_TARGETS+=$(name)/$(mapper)/libs_qc.tsv

######################################################################


######################################################################
# qc=none|on|off
ifeq ($(qc),none)
$(name)/report/qc.html $(name)/report/qc.tsv: 

else
$(name)/report/qc.html $(name)/report/qc.tsv: $(conf) $(call must_exist,$(name)/data/)  $(name)/report/fastq_qc_report.tsv
	irap_report_qc $(IRAP_REPORT_MAIN_OPTIONS) --conf $(conf) --rep_dir $(name)/report || ( rm -f $@ && exit 1)
endif

STAGE2_OUT_FILES+=$(name)/report/qc.tsv
STAGE2_TARGETS+=$(name)/report/qc.tsv

ifeq ($(qc),none)
# empty file
FASTQC_REPORT_FILES=
$(name)/report/fastq_qc_report.tsv:
	touch $@

else

FASTQC_REPORT_FILES:=$(foreach p,$(se),$(call lib2filt_folder,$(p))$(p).f.fastqc.tsv) $(foreach p,$(pe),$(call lib2filt_folder,$(p))$(p)_1.f.fastqc.tsv)


STAGE2_OUT_FILES+=$(name)/report/fastq_qc_report.tsv
STAGE2_TARGETS+=$(name)/report/fastq_qc_report.tsv

ifeq  ($(qc),off)
#FASTQC_REPORT_FILES=$(foreach p,$(pe),$(name)/report/riq/$($(p)_dir)/$(call get_fastq_prefix,$(p),pe)_1.fastqc.tsv $(name)/report/riq/$($(p)_dir)$(call get_fastq_prefix,$(p),pe)_2.fastqc.tsv) $(foreach p,$(se),$(name)/report/riq/$($(p)_dir)$(call get_fastq_prefix,$(p),se).fastqc.tsv)

$(name)/report/fastq_qc_report.tsv:  $(FASTQC_REPORT_FILES)
	$(call pass_args_stdin,irap_merge2tsv,$@.tmp, --in='$(subst $(space),;,$^)'  --out $@.tmp) && mv $@.tmp $@

%.fastqc.tsv: %.fastqc.zip
	irap_fastqc2tsv $< > $@.tmp && mv $@.tmp $@

else


$(name)/report/fastq_qc_report.tsv:  $(FASTQC_REPORT_FILES)
	$(call pass_args_stdin,irap_merge2tsv,$@.tmp, --in='$(subst $(space),;,$^)'  --out $@.tmp) && mv $@.tmp $@

# SE
%.f.fastqc.tsv: %.f.fastqc.zip
	irap_fastqc2tsv $< | sed "1s/.f$$//" > $@.tmp && mv $@.tmp $@

# PE
%_1.f.fastqc.tsv: %_1.f.fastqc.zip 
	irap_fastqc2tsv $< | sed "1s/.f$$//" > $@.tmp && mv $@.tmp $@

%_2.f.fastqc.tsv: %_2.f.fastqc.zip 
	irap_fastqc2tsv $< | sed "1s/.f$$//" > $@.tmp && mv $@.tmp $@

endif
endif





