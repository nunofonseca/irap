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

#**************
# Gene level DE
#**************
de_targets=

ifneq (none,$(de_method))
de_targets+=$(foreach cont,$(contrasts),$(de_toplevel_folder)/$(cont).genes_de.tsv)
STAGE4_OUT_FILES+=$(de_targets)
WAVE5_TARGETS+=$(de_targets)
endif

phony_targets+= de_files DE

de_files:
	echo $(de_targets)

DE: stage3 $(de_targets) 
	$(call p_info,[DONE] Differential analysis)

#*********************
# Transcript level DE
#*********************
trans_de_targets=

ifneq (none,$(transcript_de_method))
trans_de_targets+=$(foreach cont,$(contrasts),$(tde_toplevel_folder)/$(cont).transcripts_de.tsv)
STAGE4_OUT_FILES+=$(trans_de_targets)
endif

phony_targets+= transcript_de_files trans_DE

transcript_de_files:
	echo $(trans_de_targets)

trans_DE: stage3 $(trans_de_targets)
	$(call p_info,[DONE] Differential (transcript) analysis)

#*********************
# Exon level DE
#*********************
exon_de_targets=

ifneq (none,$(exon_de_method))
exon_de_targets+=$(foreach cont,$(contrasts),$(ede_toplevel_folder)/$(cont).exons_de.tsv)
STAGE4_OUT_FILES+=$(exon_de_targets)
endif
$(call p_debug,exon_de_targets=$(exon_de_targets))

exon_de_files:
	echo $(exon_de_targets)

exon_DE: stage3 $(exon_de_targets)
	$(call p_info,[DONE] Differential (exon) analysis)

phony_targets+= exon_de_files exon_DE


########################################

# % == cond
#1=contrast
define get_contrast_labels=
$(shell echo $($(1))|sed "s/ /,/g")
endef
define get_contrast_def=
$(shell echo $(foreach c,$($(1)),$(call get_contrast_labels,$(c)))|sed "s/ /;/g")
endef

#************
# DESeq
#************
# 1=counts file
# 2=contrast
# 3=de tsv file
# 4=feature type
define run_deseq=
irap_DE_deseq --tsv $(1)  --contrasts "$(call get_contrast_def,$(2))" --labels "$(call get_contrast_labels,$(2))" --out $(3) $(call get_de_annot) $(call get_de_annot_genes_only) $(if $(technical.replicates),--tech-rep "$(technical.replicates)") --feature $(4) --min $($(4)_de_min_count)
endef

#************
# DESeq2
#************
# 1=counts file
# 2=contrast
# 3=de tsv file
# 4=feature type
# TODO: include support for technical replicates
ifndef deseq2_params
deseq2_params=--independent-filtering
endif
#

define run_deseq2=
irap_DE_deseq2 --tsv $(1) --contrasts "$(call get_contrast_def,$(2))" --labels "$(call get_contrast_labels,$(2))" $(deseq2_params) --out $(3) $(call get_de_annot) $(call get_de_annot_genes_only) $(if $(technical.replicates),--tech-rep "$(technical.replicates)")  --feature $(4) --min $($(4)_de_min_count)
endef

#************
# edgeR
#************

define run_edger=
irap_DE_edgeR --tsv $(1) --contrasts "$(call get_contrast_def,$(2))" --labels "$(call get_contrast_labels,$(2))" --out $(3) $(call get_de_annot) $(call get_de_annot_genes_only) $(if $(technical.replicates),--tech-rep "$(technical.replicates)")  --feature $(4) --min $($(4)_de_min_count)
endef

#************
# VOOM
#************

define run_voom=
irap_DE_voom --tsv $(1)  --contrasts "$(call get_contrast_def,$(2))" --labels "$(call get_contrast_labels,$(2))" --out $(3) $(call get_de_annot) $(call get_de_annot_genes_only) $(if $(technical.replicates),--tech-rep "$(technical.replicates)")  --feature $(4) --min $($(4)_de_min_count)
endef

#************
# EBSeq
#************

define run_ebseq=
irap_DE_ebseq --tsv $(1)  --contrasts "$(call get_contrast_def,$(2))" --labels "$(call get_contrast_labels,$(2))" --out $(3) $(call get_de_annot) $(call get_de_annot_genes_only) $(if $(technical.replicates),--tech-rep "$(technical.replicates)") --feature $(4) --g2t $(mapTrans2gene) --min $($(4)_de_min_count)
endef

#************
# DEXSeq
#************
# 4 - always exon
define run_dexseq=
irap_DE_dexseq --tsv $(1) --min $(exon_de_min_count) --contrasts "$(call get_contrast_def,$(2))" --labels "$(call get_contrast_labels,$(2))" --out $(3) $(call get_de_annot) $(call get_de_annot_genes_only) $(if $(technical.replicates),--tech-rep "$(technical.replicates)") --feature $(4) 
endef

################################################################################
# Differential Analysis
################################################################################
# Cuffdiff
$(de_toplevel_folder1)/cuffdiff1/%.genes_de.tsv: $(de_toplevel_folder1)/$(de_method)/%.$(de_method)/gene_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

$(de_toplevel_folder1)/cuffdiff2/%.genes_de.tsv: $(de_toplevel_folder1)/$(de_method)/%.$(de_method)/gene_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

$(de_toplevel_folder1)/cuffdiff1_nd/%.genes_de.tsv: $(de_toplevel_folder1)/$(de_method)/%.$(de_method)/gene_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

$(de_toplevel_folder1)/cuffdiff2_nd/%.genes_de.tsv: $(de_toplevel_folder1)/$(de_method)/%.$(de_method)/gene_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

#
$(tde_toplevel_folder1)/cuffdiff1/%.transcripts_de.tsv: $(de_toplevel_folder1)/$(de_method)/%.$(de_method)/isoform_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

$(tde_toplevel_folder1)/cuffdiff2/%_transcripts_de.tsv: $(de_toplevel_folder1)/$(de_method)/%.$(de_method)/isoform_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

$(tde_toplevel_folder1)/cuffdiff1_nd/%_transcripts_de.tsv: $(de_toplevel_folder1)/$(de_method)/%.$(de_method)/isoform_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

$(tde_toplevel_folder1)/cuffdiff2_nd/%_transcripts_de.tsv: $(de_toplevel_folder1)/$(de_method)/%.$(de_method)/isoform_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

# b) Differential analysis with gene and transcript discovery
$(de_toplevel_folder1)/cuffdiff1/%.cuffdiff1 $(de_toplevel_folder1)/cuffdiff1/%.cuffdiff1/gene_exp.diff $(de_toplevel_folder1)/cuffdiff1/%.cuffdiff1/isorform_exp.diff: $(reference_abspath)   $(de_toplevel_folder1)/$(name).merged.gtf
	$(call run_cuffdiff,cufflinks1,cuffdiff --no-update-check -u -o $(de_toplevel_folder1)/cuffdiff1/$*.$(de_method).tmp -b $(reference_abspath) -p $(max_threads)  $(de_toplevel_folder1)/$(name).merged.gtf  $(call get_contrast_bam_files,$*)) && mv $(de_toplevel_folder1)/cuffdiff1/$*.$(de_method).tmp $(de_toplevel_folder1)/cuffdiff1/$*.$(de_method)


$(de_toplevel_folder1)/cuffdiff2/%.cuffdiff2 $(de_toplevel_folder1)/cuffdiff2/%.cuffdiff2/gene_exp.diff $(de_toplevel_folder1)/cuffdiff2/%.cuffdiff2/isoform_exp.diff: $(reference_abspath)   $(de_toplevel_folder1)/$(name).merged.gtf
	$(call run_cuffdiff,cufflinks2,cuffdiff --no-update-check -u -o $(de_toplevel_folder1)/cuffdiff2/$*.$(de_method).tmp -b $(reference_abspath) -p $(max_threads)  $(de_toplevel_folder1)/$(name).merged.gtf  $(call get_contrast_bam_files,$*)) && mv $(de_toplevel_folder1)/cuffdiff2/$*.$(de_method).tmp $(de_toplevel_folder1)/cuffdiff2/$*.$(de_method)

# a) Differential analysis without gene and transcript discovery  (use gtf from annotation)
$(de_toplevel_folder1)/cuffdiff1_nd/%.cuffdiff1_nd $(de_toplevel_folder1)/cuffdiff1_nd/%.cuffdiff1_nd/gene_exp.diff $(de_toplevel_folder1)/cuffdiff1_nd/%.cuffdiff1_nd/isoform_exp.diff: $(reference_abspath)  $(gtf_file_abspath)
	$(call run_cuffdiff,cufflinks1,cuffdiff --no-update-check -u -o $(de_toplevel_folder1)/cuffdiff1_nd/$*.$(de_method).tmp -b $(reference_abspath) -p $(max_threads)  $(gtf_file_abspath)  $(call get_contrast_bam_files,$*)) && mv $(de_toplevel_folder1)/cuffdiff1_nd/$*.$(de_method).tmp $(de_toplevel_folder1)/cuffdiff1_nd/$*.$(de_method)

# $(de_toplevel_folder1)/$(name).merged.gtf
$(de_toplevel_folder1)/cuffdiff2_nd/%.cuffdiff2_nd $(de_toplevel_folder1)/cuffdiff2_nd/%.cuffdiff2_nd/gene_exp.diff $(de_toplevel_folder1)/cuffdiff2_nd/%.cuffdiff2_nd/isoform_exp.diff: $(reference_abspath)  $(gtf_file_abspath)
	$(call run_cuffdiff,cufflinks2,cuffdiff --no-update-check -u -o $(de_toplevel_folder1)/cuffdiff2_nd/$*.$(de_method).tmp -b $(reference_abspath) -p $(max_threads)  $(gtf_file_abspath)  $(call get_contrast_bam_files,$*)) && mv $(de_toplevel_folder1)/cuffdiff2_nd/$*.$(de_method).tmp $(de_toplevel_folder1)/cuffdiff2_nd/$*.$(de_method)

###############
# DE

# DESEQ
$(de_toplevel_folder1)/deseq/%.genes_de.tsv $(de_toplevel_folder1)/deseq/%.genes_de.Rdata: $(quant_toplevel_folder)/genes.raw.$(quant_method).tsv $(annot_tsv) 
	$(call run_deseq,$<,$*,$(@D)/$*.deseq,gene) && mv $(@D)/$*.deseq/de.tsv $(@D)/$*.genes_de.tsv && mv $(@D)/$*.deseq/de.Rdata $(@D)/$*.genes_de.Rdata

$(tde_toplevel_folder1)/deseq/%.transcripts_de.tsv $(tde_toplevel_folder1)/deseq/%.transcripts_de.Rdata: $(quant_toplevel_folder)/transcripts.raw.$(quant_method).tsv $(annot_tsv) 
	$(call run_deseq,$<,$*,$(@D)/$*.deseq,transcript) && mv $(@D)/$*.deseq/de.tsv $(@D)/$*.transcripts_de.tsv && mv $(@D)/$*.deseq/de.Rdata $(@D)/$*.transcripts_de.Rdata

#ex. irap_DE_deseq matrix_file.tsv  "colnameA,colnameB,colnameC;colnameD,colnameF" "grouplabel1,grouplabel2" oprefix [annot.file.tsv] [tech.replicates.def]

# DESEQ2
$(de_toplevel_folder1)/deseq2/%.genes_de.tsv $(de_toplevel_folder1)/deseq2/%.genes_de.Rdata: $(quant_toplevel_folder)/genes.raw.$(quant_method).tsv $(annot_tsv) 
	$(call run_deseq2,$<,$*,$(@D)/$*.deseq2,gene) && mv $(@D)/$*.deseq2/de.tsv $(@D)/$*.genes_de.tsv && mv $(@D)/$*.deseq2/de.Rdata $(@D)/$*.genes_de.Rdata

$(tde_toplevel_folder1)/deseq2/%.transcripts_de.tsv $(tde_toplevel_folder1)/deseq2/%.transcripts_de.Rdata: $(quant_toplevel_folder)/transcripts.raw.$(quant_method).tsv $(annot_tsv) 
	$(call run_deseq2,$<,$*,$(@D)/$*.deseq2,transcript) && mv $(@D)/$*.deseq2/de.tsv $(@D)/$*.transcripts_de.tsv && mv $(@D)/$*.deseq2/de.Rdata $(@D)/$*.transcripts_de.Rdata


# EDGER
$(de_toplevel_folder1)/edger/%.genes_de.tsv $(de_toplevel_folder1)/edger/%.genes_de.Rdata: $(quant_toplevel_folder)/genes.raw.$(quant_method).tsv $(annot_tsv) 
	$(call run_edger,$<,$*,$(@D)/$*.edger,gene) && mv $(@D)/$*.edger/de.tsv $(@D)/$*.genes_de.tsv && mv $(@D)/$*.edger/de.Rdata $(@D)/$*.genes_de.Rdata

$(tde_toplevel_folder1)/edger/%.transcripts_de.tsv $(tde_toplevel_folder1)/edger/%.transcripts_de.Rdata: $(quant_toplevel_folder)/transcripts.raw.$(quant_method).tsv $(annot_tsv)  $(mapTrans2gene)
	$(call run_edger,$<,$*,$(@D)/$*.edger,transcript) && mv $(@D)/$*.edger/de.tsv $(@D)/$*.transcripts_de.tsv && mv $(@D)/$*.edger/de.Rdata $(@D)/$*.transcripts_de.Rdata

# VOOM
$(de_toplevel_folder1)/voom/%.genes_de.tsv $(de_toplevel_folder1)/voom/%.genes_de.Rdata: $(quant_toplevel_folder)/genes.raw.$(quant_method).tsv $(annot_tsv) 
	$(call run_voom,$<,$*,$(@D)/$*.voom,gene) && mv $(@D)/$*.voom/de.tsv $(@D)/$*.genes_de.tsv && mv $(@D)/$*.voom/de.Rdata $(@D)/$*.genes_de.Rdata

$(tde_toplevel_folder1)/voom/%.transcripts_de.tsv $(tde_toplevel_folder1)/voom/%.transcripts_de.Rdata: $(quant_toplevel_folder)/transcripts.raw.$(quant_method).tsv $(annot_tsv) 
	$(call run_voom,$<,$*,$(@D)/$*.voom,transcript) && mv $(@D)/$*.voom/de.tsv $(@D)/$*.transcripts_de.tsv && mv $(@D)/$*.voom/de.Rdata $(@D)/$*.transcripts_de.Rdata

# ebseq
$(de_toplevel_folder1)/ebseq/%.genes_de.tsv $(de_toplevel_folder1)/ebseq/%.genes_de.Rdata: $(quant_toplevel_folder)/genes.raw.$(quant_method).tsv $(annot_tsv)  $(mapTrans2gene)
	$(call run_ebseq,$<,$*,$(@D)/$*.ebseq,gene) && mv $(@D)/$*.ebseq/de.tsv $(@D)/$*.genes_de.tsv && mv $(@D)/$*.ebseq/de.Rdata $(@D)/$*.genes_de.Rdata

$(tde_toplevel_folder1)/ebseq/%.transcripts_de.tsv $(tde_toplevel_folder1)/ebseq/%.transcripts_de.Rdata: $(quant_toplevel_folder)/transcripts.raw.$(quant_method).tsv $(annot_tsv)  $(mapTrans2gene)
	$(call run_ebseq,$<,$*,$(@D)/$*.ebseq,transcript) && mv $(@D)/$*.ebseq/de.tsv $(@D)/$*.transcripts_de.tsv && mv $(@D)/$*.ebseq/de.Rdata $(@D)/$*.transcripts_de.Rdata


# dexseq
$(ede_toplevel_folder1)/dexseq/%.exons_de.tsv $(ede_toplevel_folder1)/dexseq/%.exons_de.Rdata: $(quant_toplevel_folder)/exons.raw.$(exon_quant_method).tsv $(DEXSEQ_GFF)
	$(call run_dexseq,$<,$*,$(@D)/$*.dexseq,exon) && mv $(@D)/$*.dexseq/de.tsv $(@D)/$*.exons_de.tsv && mv $(@D)/$*.dexseq/de.Rdata $(@D)/$*.exons_de.Rdata


# BAYSEQ
$(de_toplevel_folder1)/bayseq/%.genes_de.tsv: 
	echo TODO:RUN_Bayseq_FOR_DE
