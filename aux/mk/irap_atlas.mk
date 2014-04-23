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
#    $Id$
# =========================================================
# Atlas specific 
# all variable, functions, and rules should contain a atlas suffix or prefix

atlas_wrap_up: $(name)/atlas_html.tar.gz


#######################
# Get a tarball with the plots and HTML files without irap's menu, HTML head and CSS.

# TODO: add the plots for mapping...check if all files are there!
ATLAS_REPORT_FILES=$(name)/report/qc.html $(name)/report/read_filtering_plot.png $(name)/report/qc.tsv $(name)/report/read_filtering_plot.png.eps $(name)/report/mapping/$(mapper).html $(name)/report/mapping/$(mapper)*.png $(name)/report/mapping/$(mapper)*.png.eps $(name)/report/mapping/$(mapper)*.tsv $(name)/report/qc.html $(name)/report/riq/raw_data/*/fastqc_report.html $(name)/report/riq/raw_data/*/Images $(name)/report/riq/raw_data/*/Icons/ $(name)/report/software.tsv $(name)/report/info.html $(name)/report/irap.conf 

$(name)/atlas_html.tar.gz: $(name)/report $(name)/report/software.tsv $(name)/report/irap.conf
	mkdir -p $(name)/atlas && rm -rf $(name)/atlas/* && \
	tar czvf $(name)/tmp.atlas.tgz  $(ATLAS_REPORT_FILES) &&\
	cd $(name)/atlas && tar xzvf ../tmp.atlas.tgz && rm ../tmp.atlas.tgz && mv $(name)/report/* . && rm -rf $(name) && \
	for html_file in `find . -name "*.htm*"`; do atlas_clean_html.sh $$html_file; done && \
	cd .. && tar czvf $(@F).tmp  atlas/ && mv $(@F).tmp $(@F) &&\
	echo Created $@
# remove header, footer and menu

$(name)/report/irap.conf: $(conf)
	cp $< $@.tmp && mv $@.tmp $@


##########################
# Useful temporary targets
atlas_clean_report:
	$(call p_info, Cleaning some HTML files)
	@rm -f $(name)/report/qc.html $(name)/report/mapping/$(mapper).html
	$(call p_info, Files deleted)
	$(call p_info, Please rerun 'irap conf=... report' to update the files)
