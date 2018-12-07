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
clustering_files=
##
ifeq ($(rnaseq_type),sc)
all_clustering_files=
ifneq ($(clustering_method),none)
# always based on gene expression
clustering_files=$(quant_toplevel_folder)/$(clustering_method)/genes.$(quant_norm_libsize_method).$(quant_method).clusters.tsv
ifeq ($(transcript_expr),y)
# if transcript quantification is available then also generate clusters based on transcript quantification
clustering_files+=$(quant_toplevel_folder)/$(clustering_method)/transcripts.$(quant_norm_libsize_method).$(quant_method).clusters.tsv
endif
endif


SUPPORTED_CLUSTERING_METHODS=sc3 none

# generate the QC files
clustering: $(clustering_files)

cluster_spike_param=
cluster_spike_dep=
ifdef spikein_fasta
cluster_spike_param:=-s $(spikein_gtf_file)
cluster_spike_dep:=$(spikein_gtf_file)
endif

## SC3 rules
# also generates _marker_genes.tsv for each k
# and multiple files with the coordinates of TSNe based on different perplexity values - tsne_perp_PERP_VAL.tsv
$(quant_toplevel_folder)/sc3/genes.$(quant_norm_libsize_method).$(quant_method).clusters.tsv: $(quant_toplevel_folder)/genes.raw.filtered.$(quant_method).$(expr_ext) $(quant_toplevel_folder)/genes.$(quant_norm_libsize_method).$(quant_method).$(expr_ext) $(cluster_spike_dep) 
	mkdir -p $(@D) && irap_sc3 -i $< -j $(word 2,$^) --$(expr_format) --out $(quant_toplevel_folder)/sc3/genes.$(quant_norm_libsize_method).$(quant_method).irap --min_clusters $(min_clusters) --max_clusters $(max_clusters) --max_threads $(max_threads) $(cluster_spike_param) && mv $(quant_toplevel_folder)/sc3/genes.$(quant_norm_libsize_method).$(quant_method).irap_clusters.tsv $@ || ( rm -f $@* && exit 1)

$(quant_toplevel_folder)/sc3/transcripts.$(quant_norm_libsize_method).$(quant_method).clusters.tsv: $(quant_toplevel_folder)/transcripts.raw.filtered.$(quant_method).$(expr_ext) $(quant_toplevel_folder)/transcripts.$(quant_norm_libsize_method).$(quant_method).$(expr_ext) $(cluster_spike_dep)
	mkdir -p $(@D) && irap_sc3 -i $< -j $(word 2,$^) --$(expr_format) --out $(quant_toplevel_folder)/sc3/transcripts.$(quant_norm_libsize_method).$(quant_method).irap --min_clusters $(min_clusters) --max_clusters $(max_clusters) --max_threads $(max_threads) $(cluster_spike_param) || ( rm -f $@* && exit 1)


STAGE5_OUTFILES+=$(clustering_files)
STAGE5_TARGETS+=$(clustering_files)

WAVE6_TARGETS+=$(clustering_files)


endif
