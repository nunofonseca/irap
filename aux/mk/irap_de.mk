#; -*- mode: Makefile;-*-
# =========================================================
# Copyright 2012,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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

#************
# DE
#************

# used by de_seq, edger, voom
ifndef de_min_count
de_min_count=0
endif

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
# TODO: include support for technical replicates
define run_deseq=
irap_DE_deseq --tsv $(1) --min $(de_min_count) --contrasts "$(call get_contrast_def,$(2))" --labels "$(call get_contrast_labels,$(2))" --out $(3) $(call get_de_annot) $(call get_de_annot_genes_only)
endef

#************
# edgeR
#************

define run_edger=
irap_DE_edgeR --tsv $(1) --min $(de_min_count) --contrasts "$(call get_contrast_def,$(2))" --labels "$(call get_contrast_labels,$(2))" --out $(3) $(call get_de_annot) $(call get_de_annot_genes_only)
endef

#************
# VOOM
#************

define run_voom=
irap_DE_voom --tsv $(1) --min $(de_min_count) --contrasts "$(call get_contrast_def,$(2))" --labels "$(call get_contrast_labels,$(2))" --out $(3) $(call get_de_annot) $(call get_de_annot_genes_only)
endef

################################################################################
# Differential Analysis
################################################################################

# Cuffdiff
$(name)/$(mapper)/$(quant_method)/cuffdiff1/%.genes_de.tsv: $(name)/$(mapper)/$(quant_method)/$(de_method)/%.$(de_method)/gene_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/cuffdiff2/%.genes_de.tsv: $(name)/$(mapper)/$(quant_method)/$(de_method)/%.$(de_method)/gene_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/cuffdiff1_nd/%.genes_de.tsv: $(name)/$(mapper)/$(quant_method)/$(de_method)/%.$(de_method)/gene_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/cuffdiff2_nd/%.genes_de.tsv: $(name)/$(mapper)/$(quant_method)/$(de_method)/%.$(de_method)/gene_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

#
$(name)/$(mapper)/$(quant_method)/cuffdiff1/%.transcripts_de.tsv: $(name)/$(mapper)/$(quant_method)/$(de_method)/%.$(de_method)/isoform_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/cuffdiff2/%_transcripts.de.tsv: $(name)/$(mapper)/$(quant_method)/$(de_method)/%.$(de_method)/isoform_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/cuffdiff1_nd/%_transcripts.de.tsv: $(name)/$(mapper)/$(quant_method)/$(de_method)/%.$(de_method)/isoform_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/cuffdiff2_nd/%_transcripts.de.tsv: $(name)/$(mapper)/$(quant_method)/$(de_method)/%.$(de_method)/isoform_exp.diff
	cp $< $@.tmp && mv $@.tmp $@

# b) Differential analysis with gene and transcript discovery
$(name)/$(mapper)/$(quant_method)/cuffdiff1/%.cuffdiff1 $(name)/$(mapper)/$(quant_method)/cuffdiff1/%.cuffdiff1/gene_exp.diff $(name)/$(mapper)/$(quant_method)/cuffdiff1/%.cuffdiff1/isorform_exp.diff: $(reference_abspath)   $(name)/$(mapper)/$(quant_method)/$(name).merged.gtf
	$(call run_cuffdiff,cufflinks1,cuffdiff --no-update-check -u -o $(name)/$(mapper)/$(quant_method)/cuffdiff1/$*.$(de_method).tmp -b $(reference_abspath) -p $(max_threads)  $(name)/$(mapper)/$(quant_method)/$(name).merged.gtf  $(call get_contrast_bam_files,$*)) && mv $(name)/$(mapper)/$(quant_method)/cuffdiff1/$*.$(de_method).tmp $(name)/$(mapper)/$(quant_method)/cuffdiff1/$*.$(de_method)


$(name)/$(mapper)/$(quant_method)/cuffdiff2/%.cuffdiff2 $(name)/$(mapper)/$(quant_method)/cuffdiff2/%.cuffdiff2/gene_exp.diff $(name)/$(mapper)/$(quant_method)/cuffdiff2/%.cuffdiff2/isoform_exp.diff: $(reference_abspath)   $(name)/$(mapper)/$(quant_method)/$(name).merged.gtf
	$(call run_cuffdiff,cufflinks2,cuffdiff --no-update-check -u -o $(name)/$(mapper)/$(quant_method)/cuffdiff2/$*.$(de_method).tmp -b $(reference_abspath) -p $(max_threads)  $(name)/$(mapper)/$(quant_method)/$(name).merged.gtf  $(call get_contrast_bam_files,$*)) && mv $(name)/$(mapper)/$(quant_method)/cuffdiff2/$*.$(de_method).tmp $(name)/$(mapper)/$(quant_method)/cuffdiff2/$*.$(de_method)

# a) Differential analysis without gene and transcript discovery  (use gtf from annotation)
$(name)/$(mapper)/$(quant_method)/cuffdiff1_nd/%.cuffdiff1_nd $(name)/$(mapper)/$(quant_method)/cuffdiff1_nd/%.cuffdiff1_nd/gene_exp.diff $(name)/$(mapper)/$(quant_method)/cuffdiff1_nd/%.cuffdiff1_nd/isoform_exp.diff: $(reference_abspath)  $(gtf_file_abspath)
	$(call run_cuffdiff,cufflinks1,cuffdiff --no-update-check -u -o $(name)/$(mapper)/$(quant_method)/cuffdiff1_nd/$*.$(de_method).tmp -b $(reference_abspath) -p $(max_threads)  $(gtf_file_abspath)  $(call get_contrast_bam_files,$*)) && mv $(name)/$(mapper)/$(quant_method)/cuffdiff1_nd/$*.$(de_method).tmp $(name)/$(mapper)/$(quant_method)/cuffdiff1_nd/$*.$(de_method)

# $(name)/$(mapper)/$(quant_method)/$(name).merged.gtf
$(name)/$(mapper)/$(quant_method)/cuffdiff2_nd/%.cuffdiff2_nd $(name)/$(mapper)/$(quant_method)/cuffdiff2_nd/%.cuffdiff2_nd/gene_exp.diff $(name)/$(mapper)/$(quant_method)/cuffdiff2_nd/%.cuffdiff2_nd/isoform_exp.diff: $(reference_abspath)  $(gtf_file_abspath)
	$(call run_cuffdiff,cufflinks2,cuffdiff --no-update-check -u -o $(name)/$(mapper)/$(quant_method)/cuffdiff2_nd/$*.$(de_method).tmp -b $(reference_abspath) -p $(max_threads)  $(gtf_file_abspath)  $(call get_contrast_bam_files,$*)) && mv $(name)/$(mapper)/$(quant_method)/cuffdiff2_nd/$*.$(de_method).tmp $(name)/$(mapper)/$(quant_method)/cuffdiff2_nd/$*.$(de_method)

###############
# gene level DE

# DESEQ
$(name)/$(mapper)/$(quant_method)/deseq/%.genes_de.tsv $(name)/$(mapper)/$(quant_method)/deseq/%.genes_de.Rdata: $(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).tsv $(annot_tsv) 
	$(call run_deseq,$<,$*,$(@D)/$*.deseq) && mv $(@D)/$*.deseq/de.tsv $(@D)/$*.genes_de.tsv && mv $(@D)/$*.deseq/de.Rdata $(@D)/$*.genes_de.Rdata

#ex. irap_DE_deseq matrix_file.tsv  "colnameA,colnameB,colnameC;colnameD,colnameF" "grouplabel1,grouplabel2" oprefix [annot.file.tsv] [tech.replicates.def]

# EDGER
$(name)/$(mapper)/$(quant_method)/edger/%.genes_de.tsv $(name)/$(mapper)/$(quant_method)/edger/%.genes_de.Rdata: $(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).tsv $(annot_tsv) 
	$(call run_edger,$<,$*,$(@D)/$*.edger) && mv $(@D)/$*.edger/de.tsv $(@D)/$*.genes_de.tsv && mv $(@D)/$*.edger/de.Rdata $(@D)/$*.genes_de.Rdata

# VOOM
$(name)/$(mapper)/$(quant_method)/voom/%.genes_de.tsv $(name)/$(mapper)/$(quant_method)/voom/%.genes_de.Rdata: $(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).tsv $(annot_tsv) 
	$(call run_voom,$<,$*,$(@D)/$*.voom) && mv $(@D)/$*.voom/de.tsv $(@D)/$*.genes_de.tsv && mv $(@D)/$*.voom/de.Rdata $(@D)/$*.genes_de.Rdata

# BAYSEQ
$(name)/$(mapper)/$(quant_method)/bayseq/%.genes_de.tsv: 
	echo TODO:RUN_Bayseq_FOR_DE
