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

## QC for single cell data

##
ifeq ($(rnaseq_type),sc)

cell_qc_files=$(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).qc.tsv
ifeq ($(transcript_expr),y)
cell_qc_files+=$(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).qc.tsv
endif

# generate the QC files
cell_qc: $(cell_qc_files)

irap_cell_qc_params=$(if $(cell_filt_controls),--controls "$(cell_filt_controls)") --min_features $(cell_filt_min_features) $(if $(subst y,,$(cell_filt_outliers)),,--outliers) --max_ERCC $(cell_filt_max_ERCC) --min_counts $(cell_filt_min_cell_expr) --min_expression $(cell_filt_min_expression) 

$(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).qc.tsv: $(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).$(expr_ext) $(cell_filt_controls)
	irap_cell_qc --$(expr_format) -i $< --out $@.tmp $(irap_cell_qc_params) && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).qc.tsv: $(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).$(expr_ext) $(cell_filt_controls)
	irap_cell_qc --$(expr_format) -i $< --out $@.tmp $(irap_cell_qc_params) && mv $@.tmp $@


## out files
filtered_expr_matrices=$(name)/$(mapper)/$(quant_method)/genes.raw.filtered.$(quant_method).$(expr_ext)
ifeq ($(transcript_expr),y)
filtered_expr_matrices+=$(name)/$(mapper)/$(quant_method)/transcripts.raw.filtered.$(quant_method).$(expr_ext)
endif

# apply the QC filtering
filter_cells: $(filtered_expr_matrices)

## filter the expression matrix based on the QC outcome
$(name)/$(mapper)/$(quant_method)/genes.raw.filtered.$(quant_method).$(expr_ext): $(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).$(expr_ext) $(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).qc.tsv
	irap_filter_cols  -i $< -o $@ --$(expr_format) --qc $(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).qc.tsv || (rm -f $@ && exit 1)

$(name)/$(mapper)/$(quant_method)/transcripts.raw.filtered.$(quant_method).$(expr_ext): $(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).$(expr_ext) $(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).qc.tsv
	irap_filter_cols  -i $< -o $@ --$(expr_format) --qc $(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).qc.tsv  || (rm -f $@ && exit 1)


# STAGE4_OFILES+=
# STAGE4_TARGETS+=


endif
