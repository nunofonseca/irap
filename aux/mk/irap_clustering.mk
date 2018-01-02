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

## sample/cell clustering

##
ifeq ($(rnaseq_type),sc)

clustering_files=
all_clustering_files=
ifneq ($(clustering_method),none)
# always based on gene expression
clustering_files=$(name)/$(mapper)/$(quant_method)/$(clustering_method)/genes.raw.filtered.$(quant_method).clusters.tsv
ifeq ($(transcript_expr),y)
# if transcript quantification is available then also generate clusters based on transcript quantification
clustering_files+=$(name)/$(mapper)/$(quant_method)/$(clustering_method)/transcripts.raw.filtered.$(quant_method).clusters.tsv
endif
endif




# generate the QC files
clustering: $(clustering_files)


## SC3 rules
# also generates _marker_genes.tsv for each k
# and multiple files with the coordinates of TSNe based on different perplexity values - tsne_perp_PERP_VAL.tsv
$(name)/$(mapper)/$(quant_method)/sc3/genes.raw.filtered.$(quant_method).clusters.tsv: $(name)/$(mapper)/$(quant_method)/genes.raw.filtered.$(quant_method).$(expr_ext) 
	mkdir -p $(@D) && irap_sc3 -i $< --$(expr_format) --out $(name)/$(mapper)/$(quant_method)/sc3/genes.raw.filtered.$(quant_method).irap --min_clusters $(min_clusters) --max_clusters $(max_clusters) --max_threads $(max_threads) && mv $(name)/$(mapper)/$(quant_method)/sc3/genes.raw.filtered.$(quant_method).irap_clusters.tsv $@ || ( rm -f $@* && exit 1)

$(name)/$(mapper)/$(quant_method)/sc3/transcripts.raw.filtered.$(quant_method).clusters.tsv: $(name)/$(mapper)/$(quant_method)/transcripts.raw.filtered.$(quant_method).$(expr_ext) 
	mkdir -p $(@D) &&  irap_sc3 -i $< --$(expr_format) --out $(name)/$(mapper)/$(quant_method)/sc3/transcripts.raw.$(quant_method).irap --min_clusters $(min_clusters) --max_clusters $(max_clusters) --max_threads $(max_threads)  || ( rm -f $@* && exit 1)


# STAGE4_OFILES+=
# STAGE4_TARGETS+=


endif
