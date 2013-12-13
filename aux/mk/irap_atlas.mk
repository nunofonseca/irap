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

atlas_wrap_up: get_tracks $(name)/atlas_html.tar.gz


#######################
# Get a tarball with the plots and HTML files without irap's menu, HTML head and CSS.
ATLAS_REPORT_FILES=$(name)/report/qc.html $(name)/report/mapping/$(mapper).html $(name)/report/mapping/$(mapper)/*.png $(name)/report/qc.html $(name)/report/read_filtering_plot.png $(name)/report/riq/raw_data/*/fastqc_report.html $(name)/report/riq/raw_data/*/Images $(name)/report/riq/raw_data/*/Icons/

$(name)/atlas_html.tar.gz: $(name)/report
	mkdir -p $(name)/atlas && rm -rf $(name)/atlas/*
	tar czvf $(name)/tmp.atlas.tgz  $(ATLAS_REPORT_FILES)
	cd $(name)/atlas && tar xzvf ../tmp.atlas.tgz && rm ../tmp.atlas.tgz && mv $(name)/report/* . && rm -rf $(name)
# remove header, footerm and menu
	for html_file in `find . -name "*.htm*"`; do atlas_clean_html.sh $$html_file; done
	cd .. && tar czvf $(@F).tmp  atlas/ && mv $(@F).tmp $(@F)



