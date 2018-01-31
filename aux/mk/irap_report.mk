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
# Rules for creating the reports/HTML pages
#
#
BROWSER_DIR=jbrowse

# default CSS file
CSS_FILE?=irap.css
PATH2CSS_FILE?=$(IRAP_DIR)/aux/css/$(CSS_FILE)

IRAP_REPORT_MAIN_OPTIONS?=
ifdef reuse_menu
IRAP_REPORT_MAIN_OPTIONS += --reuse-menu
endif

#
report_qc_only?=n

GEN_REPORT_QC_ONLY:=$(if $(filter $(strip $(report_qc_only)),y),y,)

# 1- pat
define quiet_ls=
$(shell ls --color=never -1 $(1) 2>/dev/null)
endef
# 1 - pat
# return only the the first file (most recent file)
define quiet_ls1=
$(shell ls --color=never -1 -t $(1) 2>/dev/null| head -n 1)
endef


define file_not_empty=
$(if $(call is_empty_file,$(1)),,$(1))
endef

#define mapping_dirs=
#$(strip $(call cached_var,MAPPING_DIRS))
#endef

define quant_dirs=
$(strip $(call cached_var,QUANT_DIRS))
endef

#$(strip $(call cached_var,DE_DIRS))
define de_dirs=

endef
####################################################
#ifdef report_find_files @deprecated

# all targets
REPORT_TARGETS=

phony_targets+= report_all_targets

ifeq ($(gen_html_report),n)

#$(info $(REPORT_TARGETS))
report_all_targets:

else
# report enabled

info_targets=$(report_toplevel_folder)/info.html $(report_toplevel_folder)/versions.html

REPORT_TARGETS+=report_setup $(info_targets)  $(if $(call GEN_REPORT_QC_ONLY), qc_report, qc_report mapping_report quant_report de_report gse_report )   end_report $(report_toplevel_folder)/about.html

####################################################
#
de_html_files:= $(foreach c,$(contrasts),$(patsubst %.tsv,%.html,$(call quiet_ls,$(de_toplevel_folder)/$(c)*_de.tsv)))
gse_html_files:=$(foreach c,$(contrasts),$(patsubst %.tsv,%.html,$(call quiet_ls,$(de_toplevel_folder)/$(c)*.gse.*.tsv)))

ge_files:=$(shell find $(quant_toplevel_folder) -regextype egrep -type f -regex '.*/genes.*.(tsv|tsv.gz|mtx|mtx.gz)$$'  2>/dev/null )
te_files:=$(filter-out .riu. .dt.,$(shell find $(quant_toplevel_folder) -type f -regextype egrep -regex '.*/transcripts.*.(tsv|tsv.gz|mtx|mtx.gz)$$'  2>/dev/null ))
ee_files:=$(shell find $(quant_toplevel_folder) -regextype egrep -type f -regex '*./exons.*.(tsv|tsv.gz|mtx|mtx.gz)$$'  2>/dev/null )
all_quant_files:=$(ge_files) $(te_files) $(ee_files)

quant_html_files:=$(patsubst %.tsv,%.html,$(patsubst %.mtx.gz,%.html,$(patsubst %.tsv.gz,%.html,$(all_quant_files))))

mapping_dirs:=$(mapper_toplevel_folder)

$(call p_debug, cached_vars=$(cached_vars))

quant2report_folder:=$(shell realpath --relative-to=$(quant_toplevel_folder) $(report_toplevel_folder))
mapper2report_folder:=$(shell realpath --relative-to=$(mapper_toplevel_folder) $(report_toplevel_folder))
qc2report_folder:=$(shell realpath --relative-to=$(qc_toplevel_folder) $(report_toplevel_folder))
de2report_folder:=$(shell realpath --relative-to=$(de_toplevel_folder) $(report_toplevel_folder))


########################################################################
## DE
# 1 metric
# 2 TSV file
# 3 out dir
# 4 out file
# 5 title
# --anotation ....
#
define DE_tsv2html=
	tsvDE2html --flavour $(1) --tsv $(2) --out $(3)/$(4) --cut-off $(de_pvalue_cutoff) --species $(species) --feature $(call DEfilename2AL,$(2)) --browser  $(de2report_folder)/$(BROWSER_DIR)/ --css $(de2report_folder)/$(CSS_FILE) --title "$(5)" -a $(annot_tsv) -m $(de_num_genes_per_table)
endef

## GSE
# input,output,options,pipeline,contrast
define run_gse_report=
irap_report_gse --tsv $1 --out $2 $3  --gse_method "$(gse_tool):$(gse_method)" --pipeline $4 --contrast $5 --pvalue $(gse_pvalue)  --css ../../../$(CSS_FILE) 
endef


# GE
# 1 metric
# 2 TSV file
# 3 out dir
# 4 out file
# 5 title
# 6 feature
# --anotation ....
define GE_tsv2html=
	tsvGE2html -m $(1) --ifile $(2) --out $(3)/$(4) --species $(species)  --browser $(quant2report_folder)/$(BROWSER_DIR)/ --css $(quant2report_folder)/$(CSS_FILE) --title "$(5)" -a $(annot_tsv)  --gdef "$(call groupsdef2str)" --gnames "$(call groups2str)" -f $(6) --feat_mapping $(word 1,$(feat_mapping_files))
endef

#-x min value
#-r replicates
#-f feature (gene,exon,CDS)

#1 DEST FILE
#2 OUTDIR
#3 TSV FILE
define  gen_htseq_report=
	$(if $(3),irap_htseq_report.R $(2) $(3) $(de_min_count) && touch $(1),)
endef


# 
clean_report: $(name)/
	@find $(report_toplevel_folder)/ $(quant_toplevel_folder) $(mapper_toplevel_folder)/ $(de_toplevel_folder) -maxdepth 1 -name "*.html*"   -type f -exec rm -f {} \; ; 
	$(call p_info,Report data partially cleaned up)

##############################################################################
# Produce a HTML report
#report: $(report_toplevel_folder)/index.html mapping_report quant_report de_report
phony_targets+=report_setup clean_report


report_setup: $(call must_exist,$(report_toplevel_folder)) $(call must_exist,$(report_toplevel_folder)/mapping/) $(call must_exist,$(report_toplevel_folder)/de/) $(call must_exist,$(report_toplevel_folder)/quant/) $(call rep_browse,report_browser_setup) $(call must_exist,$(report_toplevel_folder)/irap.css) $(call must_exist,$(report_toplevel_folder)/menu.css) $(feat_mapping_files)

SETUP_DATA_FILES+=report_setup


$(report_toplevel_folder): $(report_toplevel_folder)/

$(report_toplevel_folder)/irap.css: $(PATH2CSS_FILE)
	cp -f $< $@

$(report_toplevel_folder)/menu.css: $(IRAP_DIR)/aux/css/menu.css
	cp -f $< $@

#############################
# QC
phony_targets+=qc_report

qc_html_files=$(qc_toplevel_folder)/qc.html

qc_report: $(qc_html_files)


#############################
# TODO: info.html
phony_targets+=info_report


info_report: report_setup $(info_targets)
status_report: report_setup $(report_toplevel_folder)/status.html

$(report_toplevel_folder)/info.html: $(report_toplevel_folder)/$(call notdir,$(conf))
	irap_report_expinfo --conf $<  --css $(CSS_FILE) --out $@.tmp && mv $@.tmp $@

#
$(report_toplevel_folder)/versions.html: $(report_toplevel_folder)/software.tsv $(conf) 
	tsvSoftware2html -i $< -o $@.tmp && mv $@.tmp.html $@


$(report_toplevel_folder)/$(call notdir,$(conf)): $(conf)
	cp $< $@.tmp && mv $@.tmp $@

#############################
phony_targets+=mapping_report quant_report mapping_report_req

mapping_report_targets=$(foreach m,$(call mapping_dirs), $(mapper_toplevel_folder)/$(shell basename $(m)).html) 

#$(report_toplevel_folder)/mapping/comparison.html 


mapping_report_files:
	echo $(call mapping_report_targets)
	echo $(call mapping_dirs)

print_mapping_dirs:
	echo $(MAPPING_DIRS)

mapping_report: report_setup $(mapper_toplevel_folder)/$(mapper).html


#$(foreach m,$(mapping_dirs),$(mapper_toplevel_folder)/$(m).html_req)
print_mapping_report_req: $(mapper_toplevel_folder)/$(mapper).html_req

$(mapper_toplevel_folder)/%.html_req:
	echo $(MAPPING_REPORT_PRE_STATS)

$(mapper_toplevel_folder)/%.html_doreq: $(MAPPING_REPORT_PRE_STATS)
	@echo done $@

# files required to produce the mapping report
mapping_report_req: $(MAPPING_REPORT_PRE_STATS)
#	@echo "done"

# Mapping report for a specific mapper
%/$(mapper).html:  $(conf)  %/stats_raw.tsv %/stats_perc.tsv  %/featstats_raw.tsv %/featstats_perc.tsv  %/genestats_raw.tsv 
	$(call pass_args_stdin,irap_report_mapping,$@, --out $(subst .html,,$@).1.html --mapper $(mapper) --bam_stats $(@D)/stats_raw.tsv --bam_statsp $(@D)/stats_perc.tsv --bam_fstats $(@D)/featstats_raw.tsv --bam_fstatsp $(@D)/featstats_perc.tsv --bam_gstats $(@D)/genestats_raw.tsv --css $(mapper2report_folder)/$(CSS_FILE) --cores $(max_threads) ) && mv $(subst .html,,$@).1.html  $@

################

# M
# 
define only_existing_files=
$(foreach f,$(1),$(if $(realpath $(f)),$(f) ,))
endef

define mappersFromReportPath=
$(subst /align_overall_comparison.png.tsv,,$(subst $(report_toplevel_folder)/mapping/,,$(1)))
endef

# find dirs with stats_raw.tsv 
# only perform the comparison on the existing TSV files
mappers_folders:=$(shell find $(name) -name stats_raw.tsv -exec dirname {} \;  2> /dev/null)

$(report_toplevel_folder)/mapping_comparison.html: $(foreach f,$(mappers_folders), $f/align_overall_comparison.png.tsv)
	mappers_comp_sum.R --tsv "$^" --labels "$(foreach f,$(mappers_folders), $(basename $f))" --out $(@D)/comparison --css  ../../$(CSS_FILE) && touch $@


phony_targets+= 

########################
phony_targets+=quant_report quant_report_files
silent_targets+=quant_report quant_report_files

#
quant_report: report_setup $(quant_html_files)

quant_report_files: 
	echo $(quant_html_files)

#######################################
# Quant. at gene level

define ge_html2level=
$(subst transcript,CDS,$(patsubst %s,%,$(word 1,$(subst ., ,$(notdir $*)))))
endef

define ge_html2metric=
$(word 2,$(subst ., ,$(notdir $1)))
endef

# 
#$(_toplevel_folder)/%.html: 
#	$(call GE_tsv2html,$(call ge_html2metric,$*),$(call quiet_ls1,$(name)/$(subst _x_,/,$*).tsv),$(@D),$(notdir $*).t,$(subst _x_, x ,$(subst /,,$(dir $*))),$(call ge_html2level,$*)) && \
#	cp $(subst .html,,$@).t.html $@

$(quant_toplevel_folder)/%.html: $(quant_toplevel_folder)/%.tsv $(annot_tsv) $(feat_mapping_files)
	$(call GE_tsv2html,$(call ge_html2metric,$*),$<,$(@D),$(notdir $*).t,$(subst _x_, x ,$(subst /,,$(dir $*))),$(call ge_html2level,$*)) && \
	cp $(subst .html,,$@).t.html $@

$(quant_toplevel_folder)/%.html: $(quant_toplevel_folder)/%.tsv.gz  $(annot_tsv) $(feat_mapping_files)
	$(call GE_tsv2html,$(call ge_html2metric,$*),$<,$(@D),$(notdir $*).t,$(subst _x_, x ,$(subst /,,$(dir $*))),$(call ge_html2level,$*)) && \
	cp $(subst .html,,$@).t.html $@


############################
# DE
phony_targets+=de_report de_report_files
silent_targets+=de_report_files

# TODO
de_report: report_setup $(de_html_files)

# just print the name of the files that will be produced
de_report_files:
	echo $(call de_html_files,$(name))


%.genes_de.html: %.genes_de.tsv $(annot_tsv)
	mkdir -p $(@D)
	$(call DE_tsv2html,$(subst _nd,,$(call DEfilepath2demethod,$@)),$<,$(@D),$(subst .html,,$(shell basename $@)),$(subst /, x ,$*))

############################
# GSE
phony_targets+=gse_report gse_report_files
silent_targets+=gse_report_files

#$(call gse_html_files,$(name))
# only generates the html iff the respective GSE tsv file exist
gse_report: report_setup 

# just print the name of the files that will be produced
gse_report_files:
	echo $(call gse_html_files,$(name))




$(report_toplevel_folder)/about.html: 
	cat $(IRAP_DIR)/aux/html/page.header.html $(IRAP_DIR)/aux/html/about.html  $(IRAP_DIR)/aux/html/page.footer.html >  $@


#########################
#mapping_report de_report
# TODO: remove/fix this in the future (currently necessary to update the menu)
phony_targets+=end_report

end_report: $(report_toplevel_folder)/index.html $(call must_exist,$(report_toplevel_folder)/irap.css)


# TODO: replace versions.html by info_report
# TODO $(call must_exist,$(report_toplevel_folder)/status.html)a
ifeq ($(report_qc_only),y)
$(report_toplevel_folder)/index.html: $(conf) $(info_targets) $(qc_html_files) $(call rep_browse,$(report_toplevel_folder)/jbrowse/index.html)  $(report_toplevel_folder)/about.html $(call must_exist,$(report_toplevel_folder)/irap.css) $(call must_exist,$(report_toplevel_folder)/menu.css)
	cp  $(report_toplevel_folder)/info.html $@ &&
	irap_report_main $(IRAP_REPORT_MAIN_OPTIONS) --conf $(conf) --rep_dir $(report_toplevel_folder) -m "" -q "" -d "" &&
	sleep 2 &&
	touch $@
else
$(report_toplevel_folder)/index.html: $(conf) $(info_targets)  $(call quant_html_files) $(qc_html_files) $(call mapping_report_targets) $(call de_html_files,$(name)) $(call gse_html_files,$(name))  $(call rep_browse,$(report_toplevel_folder)/jbrowse/index.html)  $(report_toplevel_folder)/about.html $(call must_exist,$(report_toplevel_folder)/irap.css) $(call must_exist,$(report_toplevel_folder)/menu.css)
	cp  $(report_toplevel_folder)/info.html $@ &&
	irap_report_main $(IRAP_REPORT_MAIN_OPTIONS) --conf $(conf) --rep_dir $(report_toplevel_folder) -m "$(call mapping_dirs)" -q "$(call quant_dirs,$(name))" -d "$(call de_dirs,$(name))" &&
	sleep 2 &&
	touch $@
endif

endif
