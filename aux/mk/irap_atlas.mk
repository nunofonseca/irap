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
# deprecated
ATLAS_SC_FILES?=
EXPR_FILES?=
ifdef atlas_run
$(info * Atlas mode enabled)

SETUP_DATA_FILES+=$(feat_mapping_files)
endif

ifeq ($(rnaseq_type),sc)
##ALL_TARGETS+=atlas_sc_wrap_up
## expression matrix
## QC
ATLAS_SC_FILES+=$(qc_toplevel_folder)/qc.tsv   $(cell_qc_files) $(filtered_expr_matrices)  $(clustering_files) $(report_toplevel_folder)/software.tsv   $(report_toplevel_folder)/versions.tsv
ifneq ($(sc_quant_viz),none)
ATLAS_SC_FILES+=$(sc_visualization_files)
endif
ifneq ($(mapper),none)
ATLAS_SC_FILES+=$(mapper_toplevel_folder)/libs_qc.tsv
endif
# mtx - always
EXPR_FILES=$(quant_toplevel_folder)/genes.raw.$(quant_method).mtx.gz $(quant_toplevel_folder)/genes.uraw.$(quant_method).mtx.gz $(quant_toplevel_folder)/genes.raw.$(quant_method).mtx_cols.gz $(quant_toplevel_folder)/genes.raw.$(quant_method).mtx_cols.gz $(quant_toplevel_folder)/genes.$(quant_norm_method).$(quant_method).$(quant_norm_tool).$(expr_ext) $(filtered_expr_matrices)
ifeq ($(transcript_quant),y)
EXPR_FILES+=$(quant_toplevel_folder)/transcripts.raw.$(quant_method).$(expr_ext) $(quant_toplevel_folder)/transcripts.uraw.$(quant_method).$(expr_ext) $(quant_toplevel_folder)/transcripts.$(quant_norm_method).$(quant_method).$(quant_norm_tool).$(expr_ext)
endif
# non-10x and no drop-seq
ifneq ($(expr_format),mtx)
EXPR_FILES+=$(quant_toplevel_folder)/genes.$(quant_norm_method).$(quant_method).$(quant_norm_tool).mtx.gz $(quant_toplevel_folder)/genes.$(quant_norm_method).$(quant_method).$(quant_norm_tool).mtx_cols.gz $(quant_toplevel_folder)/genes.$(quant_norm_method).$(quant_method).$(quant_norm_tool).mtx_rows.gz $(subst .tsv,.mtx.gz,$(filtered_expr_matrices)) $(subst .tsv,.mtx_rows.gz,$(filtered_expr_matrices)) $(subst .tsv,.mtx_cols.gz,$(filtered_expr_matrices))
endif
endif


phony_targets+= atlas_sc_wrap_up atlas_bundle
atlas_sc_wrap_up: atlas_bundle

sc_bundle_dir=$(name)/sc_bundle
## single cell bundle
atlas_bundle: $(sc_bundle_dir)/ $(ATLAS_SC_FILES) $(EXPR_FILES)
	cp -ar $(ATLAS_SC_FILES) $(sc_bundle_dir) 
	cp -ar $(sort $(EXPR_FILES)) $(sc_bundle_dir)
ifneq ($(sc_quant_viz),none)
	cp -ar $(quant_toplevel_folder)/genes.raw.filtered.$(quant_method).*tsne_perp*.tsv $(sc_bundle_dir)
	cp -ar $(subst .clusters.tsv,*_marker_genes.tsv,$(clustering_files)) $(sc_bundle_dir)

$(sc_bundle_dir)/: 
	mkdir -p $@


endif 

%.mtx.gz %.mtx_cols.gz %.mtx_rows.gz: %.tsv
	irap_tsv2mtx --tsv $< --out $*.tmp.mtx && mv $*.tmp.mtx.gz $*.mtx.gz  && mv $*.tmp.mtx_cols.gz  $*.mtx_cols.gz  && mv $*.tmp.mtx_rows.gz  $*.mtx_rows.gz 

# Reduce the resolution of some images
ATLAS_IMAGES2CONVERT=$(shell ls --color=never -1 $(report_toplevel_folder)/read_filtering_plot.png $(if $(call GEN_REPORT_QC_ONLY),,$(report_toplevel_folder)/mapping/$(mapper)*.png) 2>/dev/null | grep -v orig.png | grep -v scaled.png )
ATLAS_SCALED_IMAGES=$(subst .png,.scaled.png,$(ATLAS_IMAGES2CONVERT))

atlas_wrap_up: $(name)/atlas_html.tar.gz

atlas_wrap_up_clean:
	rm -f $(name)/atlas_html.tar.gz $(ATLAS_SCALED_IMAGES)

#######################
# Get a tarball with the plots and HTML files without irap's menu, HTML head and CSS.


# TODO: add the plots for mapping...check if all files are there!
ATLAS_REPORT_FILES=$(report_toplevel_folder)/qc.html  $(report_toplevel_folder)/qc.tsv $(report_toplevel_folder)/read_filtering_plot.png.eps $(report_toplevel_folder)/software.tsv $(report_toplevel_folder)/info.html $(report_toplevel_folder)/irap.conf $(report_toplevel_folder)/qc.html $(shell find $(report_toplevel_folder)/riq/ -name "fastqc_report.html" -print) $(shell find $(report_toplevel_folder)/riq/ -type d -name "Images"  -print) $(shell find  $(report_toplevel_folder)/riq/ -type d -name "Icons" -print) $(ATLAS_IMAGES2CONVERT) $(ATLAS_SCALED_IMAGES)  $(if $(call GEN_REPORT_QC_ONLY),,$(report_toplevel_folder)/mapping/$(mapper).html  $(shell ls --color=never $(report_toplevel_folder)/mapping/$(mapper)*.png.eps))



$(name)/atlas_html.tar.gz: $(report_toplevel_folder) $(report_toplevel_folder)/software.tsv $(report_toplevel_folder)/irap.conf $(ATLAS_SCALED_IMAGES)
	$(file >>.atlas_wrap_up_files) $(foreach O,$(ATLAS_REPORT_FILES),$(file >>.atlas_wrap_up_files,$(O))) tar czvf $(name)/tmp.atlas.tgz   --files-from  .atlas_wrap_up_files 
	mkdir -p $(name)/atlas && rm -rf $(name)/atlas/*  && \
	cd $(name)/atlas && tar xzvf ../tmp.atlas.tgz && rm ../tmp.atlas.tgz && mv $(report_toplevel_folder)/* . && rm -rf $(name) .atlas_wrap_up_files  && \
	find . -name "*.scaled.png"  -exec rename .scaled.png .png {} \;  && \
	for html_file in `find . -name "*.htm*"`; do atlas_clean_html.sh $$html_file; done && \
	cd .. && tar czvf $(@F).tmp  atlas/ && mv $(@F).tmp $(@F) && rm -rf atlas/ &&\
	echo Created $@
# remove header, footer and menu
# rename the scaled images

$(report_toplevel_folder)/irap.conf: $(conf)
	cp $< $@.tmp && mv $@.tmp $@

atlas_test: $(ATLAS_IMAGES2CONVERT)
	echo $^

# reduce the size of the images
%.orig.png: %.png
	cp $< $@

%.scaled.png: %.orig.png
	convert -scale 75% $< $*.scaled.png

##########################
# Useful temporary targets
atlas_clean_report:
	$(call p_info, Cleaning some HTML files)
	@rm -f $(report_toplevel_folder)/qc.html $(report_toplevel_folder)/mapping/$(mapper).html
	$(call p_info, Files deleted)
	$(call p_info, Please rerun 'irap conf=... report' to update the files)


