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

#************
# GSE
#************

# Currently only one tool is supported
def_gse_method=fisher
def_gse_pvalue=0.05
def_gse_minsize=3

# default values
ifndef gse_tool
gse_tool=$(def_gse_tool)
endif

ifndef gse_method
gse_method=$(def_gse_method)
endif

ifndef gse_pvalue
gse_pvalue=$(def_gse_pvalue)
endif

# minimum number of genes in a set
ifndef gse_minsize
gse_minsize=$(def_gse_minsize)
endif

################
# plots
ifndef gse_minedge
gse_minedge=3
endif

ifndef gse_top
gse_top=Inf
endif

################

# print the parameters
$(info *	gse_tool=$(gse_tool))
ifneq ($(gse_tool),none) 
$(info *	gse_method=$(gse_method))
$(info *	gse_pvalue=$(gse_pvalue))
$(info *	gse_minsize=$(gse_minsize)  (minimum number of genes))
$(info *	gse_minedge=$(gse_minedge)  (used in the GSE related plots))
$(info *	gse_top=$(gse_top)  (used in the GSE related plots))
endif

#gene2go_mapping
#gene2pathway_mapping


gse_valid_tools=piano none
ifeq (,$(filter $(gse_tool),$(gse_valid_tools)))
$(call p_info,[ERROR] gse_tool)
$(error $(gse_tool) not supported)
endif

################################
# Piano

ifeq ($(gse_tool),piano)

# foldchange column corresponds to logfoldchange

# if the user defines a file with the mapping between genes and go terms then
# use it instead of the annot. file (gene ids are in the first column, go term id or description in the second
ifdef gene2go_mapping
$(call file_exists,$(gene2go_mapping))
$(info *	Using gene to go terms mapping file $(gene2go_mapping))
gse_go_mapping=$(gene2go_mapping)
gse_go_annot_col=2
gse_go_map_file_option=--go
else
# use gene annot file
gse_go_mapping=$(annot_tsv)
gse_go_annot_col=GOterm
gse_go_map_file_option=--annotation
endif

# findstring - empty if not found
foldchange_col=$(if $(findstring cuffdiff,$(de_method)),10,6)


define run_piano_goterm=
irap_GSE_piano --tsv $1 --out $(subst .tsv,,$2) --foldchange-col $(foldchange_col) --annotation_col $(gse_go_annot_col) $(gse_go_map_file_option) $(gse_go_mapping) --pvalue $(gse_pvalue) --minsize $(gse_minsize) --method $(gse_method) --minedge $(gse_minedge) --top $(gse_top)
endef

#########
# Pathway
ifdef gene2pathway_mapping
$(call file_exists,$(gene2pathway_mapping))
$(info *	Using gene to go terms mapping file $(gene2pathway_mapping))
gse_pathway_mapping=$(gene2pathway_mapping)
gse_pathway_annot_col=2
gse_pathway_map_file_option=--go
else
gse_pathway_mapping=$(annot_tsv)
gse_pathway_annot_col=KEGG
gse_pathway_map_file_option=--annotation
endif

define run_piano_kegg=
irap_GSE_piano --tsv $1 --out $(subst .tsv,,$2) --foldchange-col $(foldchange_col) --annotation_col $(gse_pathway_annot_col) $(gse_pathway_map_file_option) $(gse_pathway_mapping) --pvalue $(gse_pvalue) --minsize $(gse_minsize) --method $(gse_method) --minedge $(gse_minedge) --top $(gse_top)
endef
################################
# validate the options
gse_piano_valid_methods=mean median sum fisher fisher-exact stouffer tailStrength wilcoxon reporter page
ifeq (,$(filter $(gse_method),$(gse_piano_valid_methods)))
$(call p_info,[ERROR] gse_method)
$(error $(gse_method) not supported)
endif


$(name)/$(mapper)/$(quant_method)/$(de_method)/%.gse.$(gse_tool).$(gse_method).go.tsv: $(name)/$(mapper)/$(quant_method)/$(de_method)/%_de.tsv $(gse_go_mapping)
	$(call run_piano_goterm,$<,$@)

$(name)/$(mapper)/$(quant_method)/$(de_method)/%.gse.$(gse_tool).$(gse_method).kegg.tsv: $(name)/$(mapper)/$(quant_method)/$(de_method)/%_de.tsv $(gse_pathway_mapping)
	$(call run_piano_kegg,$<,$@)

else
# nothing to do
gse_stage:  


endif


phony_targets+= gse_stage


#########
# reports
$(name)/report/$(mapper)/$(quant_method)/$(de_method)/%.gse.$(gse_tool).$(gse_method).go.html: $(name)/$(mapper)/$(quant_method)/$(de_method)/%.gse.$(gse_tool).$(gse_method).go.tsv
	$(call run_gse_report,$<,$@,,"$(mapper)x$(quant_method)x$(de_method)",$*)

$(name)/report/$(mapper)/$(quant_method)/$(de_method)/%.gse.$(gse_tool).$(gse_method).kegg.html: $(name)/$(mapper)/$(quant_method)/$(de_method)/%.gse.$(gse_tool).$(gse_method).kegg.tsv
	$(call run_gse_report,$<,$@,--pathway,"$(mapper)x$(quant_method)x$(de_method)",$*)
