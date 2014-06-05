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

mpileup_params+= $(mpileup_params) -gu
bcf_params+= $(bcf_params) -cgvb
vcf_params+= $(vcf_params)

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

ifeq ($(indel_snp_calling_method),samtools)


indel_snp_calling_setup: $(snp_dir) $(reference_abspath).fai

$(snp_dir):
	mkdir -p $@

$(snp_dir)/%.$(indel_snp_calling_method).bcf: $(name)/$(mapper)/%.bam 
	samtools mpileup $(mpileup_params) -f $(reference_abspath) $< | bcftools view $(bcf_params) - > $@.bcf.tmp && mv $@.bcf.tmp $@

%.vcf: %.bcf
	bcftools view $< > $@.tmp && mv $@.tmp $@
endif

phony_targets+= snp_indel_calling_stage

ifeq ($(indel_snp_calling_method),none)
indel_snp_calling_setup:
snp_indel_calling_stage: 

else
BCF_FILES=$(subst /$(mapper)/,/$(mapper)/snp/,$(subst .bam,.$(indel_snp_calling_method).bcf,$(STAGE2_OUT_FILES)))
VCF_FILES=$(subst .bcf,.vcf,$(BCF_FILES))

snp_indel_calling_stage: indel_snp_calling_setup $(VCF_FILES)

endif

SNP_INDEL_CALLING_FILES:
	echo $(BCF_FILES)
#########
# reports
$(name)/report/$(mapper)/$(quant_method)/$(de_method)/%.gse.$(gse_tool).$(gse_method).go.html: $(name)/$(mapper)/$(quant_method)/$(de_method)/%.gse.$(gse_tool).$(gse_method).go.tsv
	$(call run_gse_report,$<,$@,,"$(mapper)x$(quant_method)x$(de_method)",$*)


