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

## visualization related files
ifeq ($(rnaseq_type),sc)


## 2017-11-01: TSNE
sc_visualization_files=

ifneq ($(sc_quant_viz),none)
# always based on gene expression
sc_visualization_files+=$(quant_toplevel_folder)/genes.raw.filtered.$(quant_method).tsne.tsv
ifeq ($(transcript_expr),y)
# if transcript quantification is available then also generate clusters based on transcript quantification
sc_visualization_files+=$(quant_toplevel_folder)/transcripts.raw.filtered.$(quant_method).tsne.tsv
endif
endif


# generate the QC files
sc_visualization: $(sc_visualization_files)

## TSNE
# multiple files with the coordinates of TSNe based on different perplexity values - tsne_perp_PERP_VAL.tsv
%.tsne.tsv: %.$(expr_ext) 
	irap_tsne -i $< --$(expr_format) --out $@ --max_threads $(max_threads) -C $(tsne_min_cells) -G $(tsne_min_genes) && cp $@_tsne_perp_10.tsv $@ || ( rm -f $@* && exit 1)

# requires filtered+normalized expression
STAGE5_OUTFILES+=$(sc_visualization_files)
STAGE5_TARGETS+=$(sc_visualization_files)

WAVE7_TARGETS+=$(sc_visualization_files)
endif
