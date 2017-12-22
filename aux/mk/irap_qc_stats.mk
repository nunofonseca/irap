#; -*- mode: Makefile;-*-
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
$(call lib2filt_folder,$(1))$(1)_1.f.fastq.gz $(call lib2filt_folder,$(1))$(1)_2.f.fastq.gz: $(raw_data_dir)$($(1)_dir)/$(notdir $(word 1,$($(1))))  $(raw_data_dir)$($(1)_dir)/$(notdir $(word 2,$($(1))))
	$$(call p_info,Filtering $(call fix_libname,$(1)))
	$$(call do_quality_filtering_and_report,$(call fix_libname,$(1)),$(raw_data_dir)$($(1)_dir),$$(@D),) || (rm -f $$@ && exit 1)
endef

define make-se-qc-rule=
$(call lib2filt_folder,$(1))$(1).f.fastq.gz: $(raw_data_dir)$($(1)_dir)/$(notdir $($(1)))
	$$(call p_info,Filtering $(call fix_libname,$(1)))
	$$(call do_quality_filtering_and_report,$(1),$(raw_data_dir)$($(1)_dir),$$(@D),) || (rm -f $$@ && exit 1)

endef
# rules for SE libraries
$(foreach l,$(se),$(eval $(call make-se-qc-rule,$(l))))
# rules for PE libraries
$(foreach l,$(pe),$(eval $(call make-pe-qc-rule,$(l))))


CLEANUP_TARGETS+= clean_quality_filtering_and_report
phony_targets+= clean_quality_filtering_and_report
# Cleanup
clean_quality_filtering_and_report:
	$(foreach p,$(se) $(pe),$(call do_quality_filtering_and_report,$(p),$(raw_data_dir)$($(p)_dir)/,$(dir $(call lib2filt_folder,$(p))),clean))


# TODO: deprecated
print_qc_dirs_files:
	echo	$(foreach l,$(se) $(pe),$(name)/report/riq/$($(l)_dir) )



####
#############################
# a single file with the mapping stats
$(name)/report/libs_qc.tsv: $(name)/$(mapper)/stats_raw.tsv $(name)/$(mapper)/stats_perc.tsv  $(name)/$(mapper)/featstats_raw.tsv $(name)/$(mapper)/featstats_perc.tsv  $(name)/$(mapper)/genestats_raw.tsv 
	irap_append2tsv --in "$(name)/$(mapper)/stats_raw.tsv $(name)/$(mapper)/featstats_raw.tsv $(name)/$(mapper)/genestats_raw.tsv" --exclude_aggr  --cols_not_sorted --out $@.1.tmp &&\
	irap_append2tsv --in "$(name)/$(mapper)/stats_perc.tsv $(name)/$(mapper)/featstats_perc.tsv $(name)/$(mapper)/genestats_perc.tsv" --exclude_aggr --add_row_suffix "_perc" --cols_not_sorted --out $@.2.tmp && \
	irap_append2tsv --in "$@.1.tmp $@.2.tmp" --exclude_aggr --transpose --out $@.tmp && mv $@.tmp $@ &&\
	rm -f $@.1.tmp $@.2.tmp

# qc=none|on|off
ifeq ($(qc),none)
$(name)/report/qc.html $(name)/report/qc.tsv: 

else
$(name)/report/qc.html $(name)/report/qc.tsv: $(conf) $(call must_exist,$(name)/data/)  $(name)/report/fastq_qc_report.tsv
	irap_report_qc $(IRAP_REPORT_MAIN_OPTIONS) --conf $(conf) --rep_dir $(name)/report || ( rm -f $@ && exit 1)
endif

ifeq ($(qc),none)
# empty file
FASTQC_REPORT_FILES=
$(name)/report/fastq_qc_report.tsv:
	touch $@

else
FASTQC_REPORT_FILES=$(foreach p,$(pe),$(name)/report/riq/$($(p)_dir)raw_data/$(p)_1.f.fastqc.tsv $(name)/report/riq/$($(p)_dir)raw_data/$(p)_2.f.fastqc.tsv) $(foreach p,$(se),$(name)/report/riq/$($(p)_dir)raw_data/$(p).f.fastqc.tsv)

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





