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


ifneq ($(no_deps_check),nocheck)

# 
ifdef atlas_run
norm_files1=$(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.genes.fpkm.$(quant_method).irap.$(expr_ext) $(call lib2quant_folder,$(p))$(p).pe.genes.tpm.$(quant_method).irap.$(expr_ext)) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.fpkm.$(quant_method).irap.$(expr_ext) $(call lib2quant_folder,$(s))$(s).se.genes.tpm.$(quant_method).irap.$(expr_ext))
STAGE3_S_TARGETS+=$(norm_files1)
STAGE3_S_OFILES+=$(norm_files1)

ifeq ($(exon_quant),y)
norm_files2=$(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.exons.fpkm.$(exon_quant_method).irap.$(expr_ext) $(call lib2quant_folder,$(p))$(p).pe.exons.tpm.$(exon_quant_method).irap.$(expr_ext)) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.exons.fpkm.$(exon_quant_method).irap.$(expr_ext) $(call lib2quant_folder,$(s))$(s).se.exons.tpm.$(exon_quant_method).irap.$(expr_ext))
STAGE3_S_TARGETS+=$(norm_files2)
STAGE3_S_OFILES+=$(norm_files2)
endif
endif

endif

#########
# DEseq - normalize by library size
# 
# deseq_min_reads:=5
# Note: counts of technical replicates  have to be summed up into a single column
$(quant_toplevel_folder)/genes.deseq_nlib.$(quant_method).irap.tsv: $(quant_toplevel_folder)/genes.raw.$(quant_method).tsv
	irap_deseq_norm $<  > $@.tmp && mv $@.tmp $@

$(quant_toplevel_folder)/transcripts.deseq_nlib.$(quant_method).irap.tsv: $(quant_toplevel_folder)/transcripts.raw.$(quant_method).tsv
	irap_deseq_norm $<  > $@.tmp && mv $@.tmp $@

ifeq ($(exon_quant),y)
$(quant_toplevel_folder)/exons.deseq_nlib.$(exon_quant_method).irap.tsv: $(quant_toplevel_folder)/exons.raw.$(exon_quant_method).tsv
	irap_deseq_norm $<  > $@.tmp && mv $@.tmp $@
endif

# scran- single-cell aware normalisation

scran_params=
scran_dep=
ifdef spikein_fasta
scran_params:=-s $(spikein_gtf_file)
scran_dep:=$(spikein_gtf_file)
endif

$(quant_toplevel_folder)/genes.scran_gene.$(quant_method).$(expr_ext): $(quant_toplevel_folder)/genes.raw.filtered.$(quant_method).$(expr_ext) $(scran_dep)
	irap_scran_normalise -i $< -o $@.tmp --tsv -n genes $(scran_params) && mv $@.tmp $@

$(quant_toplevel_folder)/transcripts.scran_gene.$(quant_method).$(expr_ext): $(quant_toplevel_folder)/transcripts.raw.filtered.$(quant_method).$(expr_ext) $(scran_dep)
	irap_scran_normalise -i $< -o $@.tmp --tsv -n genes $(scran_params) && mv $@.tmp $@

# Return normalised count matrices, possibly in addition to another normalisation

normalised_count_matrices=$(patsubst %,$(quant_toplevel_folder)/genes.$(quant_norm_libsize_method).$(quant_method).%,$(expr_ext) mtx.gz mtx_cols.gz mtx_rows.gz)
ifeq ($(transcript_expr),y)
normalised_count_matrices+=$(patsubst %,$(quant_toplevel_folder)/transcripts.$(quant_norm_libsize_method).$(quant_method).%,$(expr_ext) mtx.gz mtx_cols.gz mtx_rows.gz)
endif


#################################################################
# Disabled: just copy the file with the raw quantification values
# $(quant_toplevel_folder)/%.raw.$(quant_method).tsv 
$(quant_toplevel_folder)/%.nlib.none.tsv: 
	$(call p_info,Quantification normalization disabled. Please set the quant_norm_method parameter if you do not want this behavior.)
	@$(call empty_file,$@)


################################################################################
# Normalization
################################################################################
ifneq ($(quant_method),none)
ifneq ($(quant_norm_tool),none)
ifneq ($(quant_norm_method),none)

nquant_files=$(foreach m,$(quant_norm_method),$(quant_toplevel_folder)/genes.$(m).$(quant_method).$(quant_norm_tool).tsv)


ofiles1=$(foreach m,$(quant_norm_method),$(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.genes.$(m).$(quant_method).$(quant_norm_tool).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.$(m).$(quant_method).$(quant_norm_tool).tsv))
STAGE3_S_OFILES+= $(ofiles1)
nofiles=$(ofiles1)

ifeq ($(transcript_expr),y)
nquant_files+=$(foreach m,$(quant_norm_method),$(quant_toplevel_folder)/transcripts.$(m).$(quant_method).$(quant_norm_tool).tsv)
ofiles2=$(foreach m,$(quant_norm_method),$(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.$(m).$(quant_method).$(quant_norm_tool).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.$(m).$(quant_method).$(quant_norm_tool).tsv))
STAGE3_S_OFILES+= $(ofiles2)
nofiles+=$(ofiles2)
endif

ifeq ($(exon_quant),y)
nquant_files+=$(foreach m,$(quant_norm_method),$(quant_toplevel_folder)/exons.$(m).$(exon_quant_method).$(quant_norm_tool).tsv)
ofiles3=$(foreach m,$(quant_norm_method),$(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.exons.$(m).$(exon_quant_method).$(quant_norm_tool).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.exons.$(m).$(exon_quant_method).$(quant_norm_tool).tsv))
STAGE3_S_OFILES+= $(ofiles3)
nofiles+=$(ofiles3)

endif

# 1 - feature 
irap_raw2metric_mass_params=$(if $(quant_norm_mass_biotypes), --mass_biotypes $(quant_norm_mass_biotypes) --gtf $(call get_gtf_for_feature,$(1)))

################################################
## qc stats + normalization
## $(1) = metric
## $(2) = quant_norm_tool
## Note: this is a bit inconsistent now since quant_norm_tools will generate tsv files
define make-norm-rules=
$(quant_toplevel_folder)/%.genes.$(1).$(quant_method).$(2).quant_qc.tsv: $(quant_toplevel_folder)/%.genes.$(1).$(quant_method).$(2).$(expr_ext)
	irap_quant_qc --ifile $$< --feature gene --metric $(1)  --$(expr_format) --gtf $$(gtf_file_abspath) --out $$@ || (rm -f $$@ && exit 1)

$(quant_toplevel_folder)/%.exons.$(1).$(exon_quant_method).$(2).quant_qc.tsv: $(quant_toplevel_folder)/%.exons.$(1).$(exon_quant_method).$(2).$(expr_ext)
	irap_quant_qc --ifile $$< --feature exon --metric $(1)  --$(expr_format) --gtf $$(gtf_file_abspath) --out $$@  || (rm -f $$@ && exit 1)

$(quant_toplevel_folder)/%.transcripts.$(1).$(quant_method).$(2).quant_qc.tsv: $(quant_toplevel_folder)/%.transcripts.$(1).$(quant_method).$(2).$(expr_ext)
	irap_quant_qc --ifile $$< --feature transcript --metric $(1)  --$(expr_format) --gtf $$(gtf_file_abspath) --out $$@ || (rm -f $$@ && exit 1)


$(quant_toplevel_folder)/genes.$(1).$(quant_method).irap.$(expr_ext): $(quant_toplevel_folder)/genes.raw.$(quant_method).$(expr_ext) $(feat_length)
	irap_raw2metric $(irap_raw2metric_params) --ifile $$<  --lengths $(feat_length)   --$(expr_format) --feature gene --metric $(1) $(call irap_raw2metric_mass_params,gene) --out $$@  || (rm -f $$@ && exit 1)

ifeq ($(exon_quant),y)
$(quant_toplevel_folder)/exons.$(1).$(exon_quant_method).irap.$(expr_ext): $(quant_toplevel_folder)/exons.raw.$(exon_quant_method).$(expr_ext) $(feat_length)
	irap_raw2metric $(irap_raw2metric_params) --ifile $$<  --lengths $(exon_length) --feature exon  --$(expr_format) --metric $(1) $(call irap_raw2metric_mass_params,exon) --out $$@  || (rm -f $$@ && exit 1)
endif

$(quant_toplevel_folder)/transcripts.$(1).$(quant_method).irap.$(expr_ext): $(quant_toplevel_folder)/transcripts.raw.$(quant_method).$(expr_ext) $(feat_length)
	irap_raw2metric $(irap_raw2metric_params) --ifile $$<  --lengths $(feat_length)  --$(expr_format) --feature transcript --metric $(1) $(call irap_raw2metric_mass_params,transcript) --out $$@  || (rm -f $$@ && exit 1)


# per library
$(quant_toplevel_folder)/%.genes.$(1).$(quant_method).irap.$(expr_ext): $(quant_toplevel_folder)/%.genes.raw.$(quant_method).$(expr_ext) $(feat_length)
	irap_raw2metric $(irap_raw2metric_params) --ifile $$<   --$(expr_format) --lengths $(feat_length) --feature gene --metric $(1) $(call irap_raw2metric_mass_params,gene) --out $$@ || (rm -f $$@ && exit 1)

$(quant_toplevel_folder)/%.transcripts.$(1).$(quant_method).irap.$(expr_ext): $(quant_toplevel_folder)/%.transcripts.raw.$(quant_method).$(expr_ext)
	irap_raw2metric $(irap_raw2metric_params) --ifile $$<  --lengths $(feat_length)  --$(expr_format) --feature transcript --metric $(1) $(call irap_raw2metric_mass_params,transcript) --out $$@ || (rm -f $$@ && exit 1)


$(quant_toplevel_folder)/%.exons.$(1).$(exon_quant_method).irap.$(expr_ext): $(quant_toplevel_folder)/%.exons.raw.$(exon_quant_method).$(expr_ext)
	irap_raw2metric $(irap_raw2metric_params) --ifile $$<  --lengths $(exon_length) --$(expr_format) --feature exon --metric $(1) $(call irap_raw2metric_mass_params,exon) --out $$@ || (rm -f $$@ && exit 1)

endef

## generate the rules for each norm. method
$(foreach nm,fpkm tpm fpkm-uq uq-fpkm,$(eval $(call make-norm-rules,$(nm),$(quant_norm_tool))))


####################################
##
ifeq ($(isl_mode),y)
STAGE3_S_TARGETS+= $(nofiles:.tsv=.quant_qc.tsv)
STAGE3_S_OFILES+=$(nofiles:.tsv=.quant_qc.tsv)
endif
ifeq ($(isl_mode),y)
STAGE3_S_TARGETS+= $(nofiles:.tsv=.quant_qc.tsv)
STAGE3_S_OFILES+=$(nofiles:.tsv=.quant_qc.tsv)
endif


STAGE4_OUT_FILES+= $(nquant_files) 
WAVE5_TARGETS+=$(nquant_files)

phony_targets+= norm_quant
norm_quant: $(quant_method)_quant $(nquant_files)
endif
else
norm_quant: 
endif
else
norm_quant: 
endif
