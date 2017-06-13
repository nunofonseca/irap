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

# normalized expression values can be computed
# by the quantification tool (e.g., cufflinks produces r/fpkms)
# or by irap code iff there are raw counts

# filenames
# <expr. level>.<norm.method>.<quant.method>.<norm.tool>.tsv
# expr.level=genes|transcripts|exons
# norm.method=deseq_nlib|rpkm|...
# norm.tool=irap or the quantification method


# add feature length to stage0 
ifeq ($(quant_norm_method),fpkm)
SETUP_DATA_FILES+=$(feat_length)
endif
ifeq ($(quant_norm_method),tpm)
SETUP_DATA_FILES+=$(feat_length)
endif

######
# irap
# FPKMs irap computes the FPKMs based on the raw counts
$(name)/$(mapper)/$(quant_method)/genes.fpkm.$(quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).tsv $(feat_length)
	irap_raw2metric --tsv $<  --lengths $(feat_length) --feature gene --metric fpkm --out $@.tmp && mv $@.tmp $@	

ifeq ($(exon_quant),y)
$(name)/$(mapper)/$(quant_method)/exons.fpkm.$(exon_quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/exons.raw.$(exon_quant_method).tsv $(feat_length)
	irap_raw2metric --tsv $<  --lengths $(exon_length) --feature exon --metric fpkm --out $@.tmp && mv $@.tmp $@	
endif

$(name)/$(mapper)/$(quant_method)/transcripts.fpkm.$(quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).tsv $(feat_length)
	irap_raw2metric --tsv $<  --lengths $(feat_length) --feature transcript --metric fpkm --out $@.tmp && mv $@.tmp $@	


# per library
$(name)/$(mapper)/$(quant_method)/%.genes.fpkm.$(quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/%.genes.raw.$(quant_method).tsv
	irap_raw2metric --tsv $<  --lengths $(feat_length) --feature gene --metric fpkm --out $@.tmp && mv $@.tmp $@	

$(name)/$(mapper)/$(quant_method)/%.transcripts.fpkm.$(quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/%.transcripts.raw.$(quant_method).tsv
	irap_raw2metric --tsv $<  --lengths $(feat_length) --feature transcript --metric fpkm --out $@.tmp && mv $@.tmp $@	


$(name)/$(mapper)/$(quant_method)/%.exons.fpkm.$(exon_quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/%.exons.raw.$(exon_quant_method).tsv
	irap_raw2metric --tsv $<  --lengths $(exon_length) --feature exon --metric fpkm --out $@.tmp && mv $@.tmp $@	

# 
ifdef atlas_run
norm_files:=$(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.genes.fpkm.$(quant_method).irap.tsv $(call lib2quant_folder,$(p))$(p).pe.genes.tpm.$(quant_method).irap.tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.fpkm.$(quant_method).irap.tsv $(call lib2quant_folder,$(s))$(s).se.genes.tpm.$(quant_method).irap.tsv)
STAGE3_S_TARGETS+=$(norm_files)
STAGE3_S_OFILES+=$(norm_files)

ifeq ($(exon_quant),y)
norm_files:=$(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.exons.fpkm.$(exon_quant_method).irap.tsv $(call lib2quant_folder,$(p))$(p).pe.exons.tpm.$(exon_quant_method).irap.tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.exons.fpkm.$(exon_quant_method).irap.tsv $(call lib2quant_folder,$(s))$(s).se.exons.tpm.$(exon_quant_method).irap.tsv)
STAGE3_S_TARGETS+=$(norm_files)
STAGE3_S_OFILES+=$(norm_files)
endif
endif

#####################
# 
# TPM
$(name)/$(mapper)/$(quant_method)/genes.tpm.$(quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).tsv $(feat_length)
	irap_raw2metric --tsv $<  --lengths $(feat_length) --feature gene --metric tpm --out $@.tmp && mv $@.tmp $@	

ifeq ($(exon_quant),y)
$(name)/$(mapper)/$(quant_method)/exons.tpm.$(exon_quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/exons.raw.$(exon_quant_method).tsv $(feat_length)
	irap_raw2metric --tsv $<  --lengths $(exon_length) --feature exon --metric tpm --out $@.tmp && mv $@.tmp $@	
endif

$(name)/$(mapper)/$(quant_method)/transcripts.tpm.$(quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).tsv $(feat_length)
	irap_raw2metric --tsv $<  --lengths $(feat_length) --feature transcript --metric tpm --out $@.tmp && mv $@.tmp $@	


# per library
$(name)/$(mapper)/$(quant_method)/%.genes.tpm.$(quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/%.genes.raw.$(quant_method).tsv
	irap_raw2metric --tsv $<  --lengths $(feat_length) --feature gene --metric tpm --out $@.tmp && mv $@.tmp $@	

$(name)/$(mapper)/$(quant_method)/%.transcripts.tpm.$(quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/%.transcripts.raw.$(quant_method).tsv
	irap_raw2metric --tsv $<  --lengths $(feat_length) --feature transcript --metric tpm --out $@.tmp && mv $@.tmp $@	


$(name)/$(mapper)/$(quant_method)/%.exons.tpm.$(exon_quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/%.exons.raw.$(exon_quant_method).tsv
	irap_raw2metric --tsv $<  --lengths $(exon_length) --feature exon --metric tpm --out $@.tmp && mv $@.tmp $@	

#########
# DEseq - normalize by library size
# 
# deseq_min_reads:=5
# Note: counts of technical replicates  have to be summed up into a single column
$(name)/$(mapper)/$(quant_method)/genes.deseq_nlib.$(quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).tsv
	irap_deseq_norm $<  > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/transcripts.deseq_nlib.$(quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).tsv
	irap_deseq_norm $<  > $@.tmp && mv $@.tmp $@

ifeq ($(exon_quant),y)
$(name)/$(mapper)/$(quant_method)/exons.deseq_nlib.$(exon_quant_method).irap.tsv: $(name)/$(mapper)/$(quant_method)/exons.raw.$(exon_quant_method).tsv
	irap_deseq_norm $<  > $@.tmp && mv $@.tmp $@
endif

##################################
# tpm
#$(name)/$(mapper)/$(quant_method)/exons.tpm.$(exon_quant_method).irap.tsv:
#	$(call p_error, Under implementation)

#$(name)/$(mapper)/$(quant_method)/genes.tpm.$(quant_method).irap.tsv:
#	$(call p_error, Under implementation)

#$(name)/$(mapper)/$(quant_method)/transcripts.tpm.$(quant_method).irap.tsv:
#	$(call p_error, Under implementation)

#################################################################
# Disabled: just copy the file with the raw quantification values
# $(name)/$(mapper)/$(quant_method)/%.raw.$(quant_method).tsv 
$(name)/$(mapper)/$(quant_method)/%.nlib.none.tsv: 
	$(call p_info,Quantification normalization disabled. Please set the quant_norm_method parameter if you do not want this behavior.)
	@$(call empty_file,$@)


################################################################################
# Normalization
################################################################################
ifneq ($(quant_method),none)
ifneq ($(quant_norm_tool),none)
ifneq ($(quant_norm_method),none)

nquant_files=$(foreach m,$(quant_norm_method),$(name)/$(mapper)/$(quant_method)/genes.$(m).$(quant_method).$(quant_norm_tool).tsv)


ofiles:=$(foreach m,$(quant_norm_method),$(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.genes.$(m).$(quant_method).$(quant_norm_tool).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.$(m).$(quant_method).$(quant_norm_tool).tsv))
STAGE3_S_OFILES+= $(ofiles)
nofiles=$(ofiles)

ifeq ($(transcript_quant),y)
nquant_files+=$(foreach m,$(quant_norm_method),$(name)/$(mapper)/$(quant_method)/transcripts.$(m).$(quant_method).$(quant_norm_tool).tsv)
ofiles:=$(foreach m,$(quant_norm_method),$(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.$(m).$(quant_method).$(quant_norm_tool).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.$(m).$(quant_method).$(quant_norm_tool).tsv))
STAGE3_S_OFILES+= $(ofiles)
nofiles+=$(ofiles)
endif

ifeq ($(exon_quant),y)
nquant_files+=$(foreach m,$(quant_norm_method),$(name)/$(mapper)/$(quant_method)/exons.$(m).$(exon_quant_method).$(quant_norm_tool).tsv)
ofiles=$(foreach m,$(quant_norm_method),$(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.exons.$(m).$(exon_quant_method).$(quant_norm_tool).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.exons.$(m).$(exon_quant_method).$(quant_norm_tool).tsv))
STAGE3_S_OFILES+= $(ofiles)
nofiles+=$(ofiles)


endif


################################################
## collect QC stats
## $(1) = metric
## $(2) = quant_norm_tool
define make-qc-norm-rules=
$(name)/$(mapper)/$(quant_method)/%.genes.$(1).$(quant_method).$(2).quant_qc.tsv: $(name)/$(mapper)/$(quant_method)/%.genes.$(1).$(quant_method).$(2).tsv
	irap_quant_qc --tsv $$< --feature gene --metric $(1) --gtf $$(gtf_file_abspath) --out $$@.tmp && mv $$@.tmp $$@

$(name)/$(mapper)/$(quant_method)/%.exons.$(1).$(exon_quant_method).$(2).quant_qc.tsv: $(name)/$(mapper)/$(quant_method)/%.exons.$(1).$(exon_quant_method).$(2).tsv
	irap_quant_qc --tsv $$< --feature exon --metric $(1) --gtf $$(gtf_file_abspath) --out $$@.tmp && mv $$@.tmp $$@

$(name)/$(mapper)/$(quant_method)/%.transcripts.$(1).$(quant_method).$(2).quant_qc.tsv: $(name)/$(mapper)/$(quant_method)/%.transcripts.$(1).$(quant_method).$(2).tsv
	irap_quant_qc --tsv $$< --feature transcript --metric $(1) --gtf $$(gtf_file_abspath) --out $$@.tmp && mv $$@.tmp $$@

endef

## generate the rules for each norm. method
$(foreach nm,$(SUPPORTED_NORM_METHODS),$(eval $(call make-qc-norm-rules,$(nm),$(quant_norm_tool))))


####################################
##
ifeq ($(isl_mode),y)
STAGE3_S_TARGETS+= $(nofiles:.tsv=.quant_qc.tsv)
STAGE3_S_OFILES+=$(nofiles:.tsv=.quant_qc.tsv)
endif


STAGE4_OUT_FILES+= $(nquant_files) 

phony_targets+= norm_quant
norm_quant: $(quant_method)_quant $(nquant_files)
endif
else
norm_quant: 
endif
else
norm_quant: 
endif
