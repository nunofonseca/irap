#; -*- mode: Makefile;-*-
# =========================================================
# Copyright 2012-2013,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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

# default values
ifndef gse_tool
gse_tool=none
endif

ifndef gse_method
gse_method=fisher
endif

ifndef gse_pvalue
gse_pvalue=0.05
endif

# minimum number of genes in a set
ifndef gse_minsize
gse_minsize=3
endif

# print the parameters
$(info *	gse_tool=$(gse_tool))
ifneq ($(gse_tool),none) 
$(info *	gse_method=$(gse_method))
$(info *	gse_pvalue=$(gse_pvalue))
$(info *	gse_minsize=$(gse_minsize)  (minimum number of genes))
endif


gse_valid_tools=piano none
ifeq (,$(filter $(gse_tool),$(gse_valid_tools)))
$(call p_info,[ERROR] gse_tool)
$(error $(gse_tool) not supported)
endif

################################
# Piano

ifeq ($(gse_tool),piano)

define run_piano_goterm=
irap_GSE_piano --tsv $1 --out $(subst .tsv,,$2) --annotation_col GOterm --annotation $(annot_tsv) --pvalue $(gse_pvalue) --minsize $(gse_minsize) --method $(gse_method) 
endef

define run_piano_kegg=
irap_GSE_piano --tsv $1 --out $(subst .tsv,,$2) --annotation_col KEGG --annotation $(annot_tsv) --pvalue $(gse_pvalue) --minsize $(gse_minsize) --method $(gse_method) 
endef

################################
# validate the options
gse_piano_valid_methods=mean median sum fisher stouffer tailStrength wilcoxon reporter page
ifeq (,$(filter $(gse_method),$(gse_piano_valid_methods)))
$(call p_info,[ERROR] gse_method)
$(error $(gse_method) not supported)
endif


$(name)/$(mapper)/$(quant_method)/$(de_method)/%.$(gse_tool).$(gse_method).go.tsv: $(name)/$(mapper)/$(quant_method)/$(de_method)/%_de.tsv $(annot_tsv)
	$(call run_piano_goterm,$<,$@)

$(name)/$(mapper)/$(quant_method)/$(de_method)/%.$(gse_tool).$(gse_method).kegg.tsv: $(name)/$(mapper)/$(quant_method)/$(de_method)/%_de.tsv $(annot_tsv)
	$(call run_piano_kegg,$<,$@)

else
# nothing to do
gse_stage:  

endif

phony_targets+= gse_stage
