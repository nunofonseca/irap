#; -*- mode: Makefile;-*-
# =========================================================
# Copyright 2012-2015,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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

#*******************
# Indel/snp calling
#*******************

#-u tells it to output into an uncompressed bcf file (rather than compressed)
#-D tells it to keep read depth for each sample
#-r tells it which chromosome region to call SNPs for (you can omit this if you wan#t to do the whole genome, but in the interest of speed, we picked a 50kb region)
#-f tells it that the next argument is going to be the reference genome file
#-r 2L:100,000-150,000

#[9:49:13 AM] Nuno Fonseca:  samtools mpileup -r 1:10000-100000 -g -u -f ~nf/storage3/ref.fa ../../georgepanos/tophat2/SM.X1549.120928.8.se.hits.bam | bcftools view -cgbv - > test.bcf


ifdef test_snp 
indel_snp_calling_method=samtools

endif

# Currently only one tool is supported
indel_snp_calling_method_valid_tools=samtools none
def_indel_snp_calling_method=none

# initialize
mpileup_params?=
bcf_call_params?=
vcf_params?=

mpileup_params+= -v -u
bcf_call_params+=  -c -O b
vcf_params+= 

# default values
ifndef indel_snp_calling_method
indel_snp_calling_method=$(def_indel_snp_calling_method)
endif

################
# print the parameters
$(info *	indel_snp_calling_method=$(indel_snp_calling_method))
ifneq ($(indel_snp_calling_method),none) 
$(info *	mpileup_params=$(mpileup_params))
endif

indel_snp_calling_method_valid_methods=samtools none
ifeq (,$(filter $(indel_snp_calling_method),$(indel_snp_calling_method_valid_tools)))
$(call p_info,[ERROR] indel_snp_calling_method)
$(error $(indel_snp_calling_method) not supported)
endif

################################
# samtools
snp_dir=$(name)/$(mapper)/snp

lib2snp_folder=$(snp_dir)/$($(1)_dir)

ifeq ($(indel_snp_calling_method),samtools)


indel_snp_calling_setup: $(snp_dir) $(reference_abspath).fai

SETUP_DATA_FILES+=indel_snp_calling_setup

bcftools_cmd=bcftools

$(snp_dir):
	mkdir -p $@

%.vcf.gz: %.bcf
	$(bcftools_cmd) view $< | gzip -c > $@.tmp && mv $@.tmp $@
endif

phony_targets+= snp_indel_calling_stage

ifeq ($(indel_snp_calling_method),none)
indel_snp_calling_setup:
snp_indel_calling_stage: 

else
BCF_FILES=$(subst /$(mapper)/,/$(mapper)/snp/,$(subst .hits.bam,.$(indel_snp_calling_method).bcf,$(STAGE2_OUT_FILES)))
VCF_FILES=$(subst .bcf,.vcf.gz,$(BCF_FILES))

snp_indel_calling_stage: indel_snp_calling_setup $(VCF_FILES)


define make-snp-rule=
$(call lib2snp_folder,$(1))$(2).$(indel_snp_calling_method).bcf: $(call lib2bam_folder,$(1))$(2).hits.bam 
	mkdir -p $$(@D) && samtools mpileup $$(mpileup_params) -f $$(reference_abspath) $$<  | gzip -c | $(bcftools_cmd) call $$(bcf_call_params) /dev/stdin | > $$@.bcf.tmp && mv $$@.bcf.tmp $$@
endef

$(foreach l,$(se),$(eval $(call make-snp-rule,$(l),$(l).se)))
$(foreach l,$(pe),$(eval $(call make-snp-rule,$(l),$(l).pe)))

endif

SNP_INDEL_CALLING_FILES:
	echo $(BCF_FILES)

#########
# reports

