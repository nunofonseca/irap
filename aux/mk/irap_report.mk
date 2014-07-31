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
#    $Id: 0.1.3 Nuno Fonseca Fri Dec 21 11:56:56 2012$
# =========================================================
# Rules for creating the reports

BROWSER_DIR=jbrowse

ifndef CSS_FILE
CSS_FILE=irap.css
endif
PATH2CSS_FILE=$(IRAP_DIR)/aux/css/irap.css

# useful functions
define mapping_dirs=
$(shell find  $(1) -name "*.hits.bam" 2>/dev/null  | sed "s|/[^/]*.bam||" | sed -E "s|($(1)/[^/]+)(/.*)|\1|" | sort -u | grep -E "$(shell echo $(SUPPORTED_MAPPERS)| sed 's/ /|/g')"   2>/dev/null )
endef

define bam_files=
$(shell find  $(1) -name "*.hits.bam"  2>/dev/null | sort -u)
endef

#define de_dirs_orig=
#$(shell ls -d -1 $(1)/*/*/* 2>/dev/null| grep -E "($(shell echo $(SUPPORTED_DE_MET#HODS)| sed 's/ /|/g'))$$")
#endef

define de_dirs=
$(shell ls -d -1 $(1)/{$(shell echo $(SUPPORTED_MAPPERS)| sed 's/ /,/g')}/{$(shell echo $(SUPPORTED_QUANT_METHODS)| sed 's/ /,/g')}/{$(shell echo $(SUPPORTED_DE_METHODS)| sed 's/ /,/g')}/ 2>/dev/null)
endef


#define quant_dirs_orig=
#$(shell ls -d -1 $(1)/*/* 2>/dev/null| grep -E "($(shell echo $(SUPPORTED_QUANT_ME#THODS)| sed 's/ /|/g'))$$")
#endef

define quant_dirs=
$(shell ls -d -1 $(1)/{$(shell echo $(SUPPORTED_MAPPERS)| sed 's/ /,/g')}/{$(shell echo $(SUPPORTED_QUANT_METHODS)| sed 's/ /,/g')}/ 2>/dev/null)
endef

# 1 - exp name
define de_html_files=
$(subst $(name)/,$(name)/report/,$(foreach d,$(call de_dirs,$(1)),$(subst .tsv,.html,$(call quiet_ls,$(d)/*_de.tsv))))
endef

# 1 - exp name
define gse_html_files=
$(subst $(name)/,$(name)/report/,$(foreach d,$(call de_dirs,$(1)),$(subst .tsv,.html,$(call quiet_ls,$(d)/*.gse.*.tsv))))
endef

# mapper quant raw|nlib|rpkm gene|exon|trans
define quant_target=
$(if $(call quiet_ls1,$(name)/$(1)/$(2)/$(4)s.$(3).*.tsv), $(name)/report/quant/$(1)_x_$(2)/$(4).$(3).html, )
endef


# 1 metric
# 2 TSV file
# 3 out dir
# 4 out file
# 5 title
# --anotation ....
#
define DE_tsv2html=
	tsvDE2html --flavour $(1) --tsv $(2) --out $(3)/$(4) --cut-off $(de_pvalue_cutoff) --species $(species) --feature $(call DEfilename2AL,$(2)) --browser ../../../$(BROWSER_DIR)/ --css ../../../$(CSS_FILE) --title "$(5)" -a $(annot_tsv) -m $(de_num_genes_per_table)
endef


# input,output,options,pipeline,contrast
define run_gse_report=
irap_report_gse --tsv $1 --out $2 $3  --gse_method "$(gse_tool):$(gse_method)" --pipeline $4 --contrast $5 --pvalue $(gse_pvalue)  --css ../../../$(CSS_FILE)
endef


# 1 metric
# 2 TSV file
# 3 out dir
# 4 out file
# 5 title
# 6 feature
# --anotation ....
define GE_tsv2html=
	tsvGE2html -m $(1) --tsv $(2) --out $(3)/$(4) --species $(species)  --browser ../../../../$(BROWSER_DIR)/ --css ../../../../$(CSS_FILE) --title "$(5)" -a $(annot_tsv)  --gdef "$(call groupsdef2str)" --gnames "$(call groups2str)" -f $(6) --feat_mapping $(feat_mapping_file)
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

# 1- pat
define quiet_ls=
$(shell ls -1 $(1) 2>/dev/null)
endef

# 1 - pat
# return only the the first file (most recent file)
define quiet_ls1=
$(shell ls -1 -t $(1) 2>/dev/null| head -n 1)
endef

ifndef IRAP_REPORT_MAIN_OPTIONS
IRAP_REPORT_MAIN_OPTIONS=
endif
# 
ifdef reuse_menu
IRAP_REPORT_MAIN_OPTIONS += --reuse-menu
endif

must_exist=$(if  $(realpath $(1)),,$(1))


clean_report: 
	@find $(name)/report/mapping/ $(name)/report/quant/ $(name)/report/de/  -maxdepth 1 -type f -exec rm -f {} \; 
	$(call p_info,Report folder partially cleaned up)

##############################################################################
# Produce a HTML report
#report: $(name)/report/index.html mapping_report quant_report de_report
phony_targets+=report_setup clean_report


report_setup: $(call must_exist,$(name)/report) $(call must_exist,$(name)/report/mapping/) $(call must_exist,$(name)/report/de/) $(call must_exist,$(name)/report/quant/) $(call rep_browse,report_browser_setup) $(call must_exist,$(name)/report/irap.css) $(call must_exist,$(name)/report/menu.css)

$(name)/report/:
	mkdir -p $@

$(name)/report/mapping/:
	mkdir -p $@

$(name)/report/de/:
	mkdir -p $@

$(name)/report/quant/:
	mkdir -p $@

$(name)/report/irap.css: $(PATH2CSS_FILE)
	cp -f $< $@

$(name)/report/menu.css: $(IRAP_DIR)/aux/css/menu.css
	cp -f $< $@

#############################
# QC
phony_targets+=qc_report

qc_html_files=$(name)/report/qc.html

qc_report: $(qc_html_files)

$(name)/report/qc.html: $(conf) $(call must_exist,$(name)/data/)
	irap_report_qc $(IRAP_REPORT_MAIN_OPTIONS) --conf $(conf) --rep_dir $(name)/report 

#############################
# TODO: info.html
phony_targets+=info_report
info_targets=$(name)/report/info.html $(name)/report/versions.html

info_report: report_setup $(info_targets)

$(name)/report/info.html: $(name)/report/$(call notdir,$(conf))
	irap_report_expinfo --conf $<  --css $(CSS_FILE) --out $@.tmp && mv $@.tmp $@

#
$(name)/report/versions.html: $(name)/report/software.tsv $(conf) 
	tsvSoftware2html -i $< -o $@.tmp && mv $@.tmp.html $@

$(name)/report/status.html:

$(name)/report/$(call notdir,$(conf)): $(conf)
	cp $< $@.tmp && mv $@.tmp $@

#############################
phony_targets+=mapping_report quant_report

define mapping_report_targets=
$(foreach m,$(call mapping_dirs,$(name)), $(name)/report/mapping/$(shell basename $(m)).html) $(name)/report/mapping/comparison.html 
endef

#mapping_report_targets=$(foreach m,$(call mapping_dirs,$(name)), $(name)/report/mapping/$(shell basename $(m)).html) $(name)/report/mapping/comparison.html 

mapping_report_files:
	echo $(call mapping_report_targets)
	echo $(call mapping_dirs,$(name))

mapping_report: report_setup $(call mapping_report_targets)


$(name)/report/mapping/%.html: $(name)/%/   $(foreach p,$(pe),$(name)/%/$($(p)_dir)$(p).pe.hits.bam) $(foreach s,$(se),$(name)/%/$($(s)_dir)$(s).se.hits.bam) $(foreach p,$(pe),$(name)/%/$($(p)_dir)$(p).pe.hits.bam.stats) $(foreach s,$(se),$(name)/%/$($(s)_dir)$(s).se.hits.bam.stats) $(foreach p,$(pe),$(name)/%/$($(p)_dir)$(p).pe.hits.bam.gene.stats) $(foreach s,$(se),$(name)/%/$($(s)_dir)$(s).se.hits.bam.gene.stats) $(conf) $(call must_exist,$(name)/report/mapping/)
	$(call pass_args_stdin,irap_report_mapping,$@, --out $(subst .html,,$@).1.html --mapper $* --pe $(subst $(space),,$(foreach p,$(pe),;$(name)/$*/$($(p)_dir)$(p).pe.hits.bam)) --se $(subst $(space),,$(foreach s,$(se),;$(name)/$*/$($(s)_dir)$(s).se.hits.bam))  --pe_labels $(subst  $(space),,$(foreach p,$(pe),;$(p))) --se_labels $(subst $(space),,$(foreach s,$(se),;$(s)))" --css ../$(CSS_FILE) $@ ) && mv $(subst .html,,$@).1.html  $@

# statistics per bam file
%.bam.stats: %.bam $(gff3_file_abspath)
	bedtools coverage -abam $< -counts -b $(gff3_file_abspath) > $@.tmp && \
	mv $@.tmp $@

#$(name)/report/mapping/%.stats: $(call must_exist,$(name)/report/mapping/) $(call must_exist,$(name)/%/)  $(foreach p,$(pe),$(name)/%/$($(p)_dir)$(p).pe.hits.bam.stats) $(foreach s,$(se),$(name)/%/$($(s)_dir)$(s).se.hits.bam.stats) $(foreach p,$(pe),$(name)/%/$($(p)_dir)$(p).pe.hits.bam.gene.stats) $(foreach s,$(se),$(name)/%/$($(s)_dir)$(s).se.hits.bam.gene.stats) 

%.bam.gene.stats: %.bam $(name)/data/exons.bed $(name)/data/introns.bed
	echo -n "Exons	" > $@.tmp &&\
	bedtools intersect -abam $<  -b $(name)/data/exons.bed |samtools view -c - >> $@.tmp && echo >> $@ &&\
	echo -n "Introns	" >> $@.tmp &&\
	bedtools intersect -abam $<  -b $(name)/data/introns.bed |samtools view -c - >> $@.tmp && echo >> $@ && \
	mv $@.tmp $@

# bed files required to get some extra stats
# exons.bed
$(name)/data/exons.bed: $(gff3_file_abspath) 
	cat $< | awk 'BEGIN{OFS="\t";} $$3=="exon" {print $$1,$$4,$$5}' | bedtools sort -i /dev/stdin | bedtools merge -i /dev/stdin > $@.tmp && \
	mv $@.tmp $@

# genes.bed
$(name)/data/genes.bed: $(gff3_file_abspath)
	cat $< | awk 'BEGIN{OFS="\t";} $$3=="gene" {print $$1,$$4,$$5}' |  bedtools sort -i /dev/stdin | bedtools merge -i /dev/stdin > $@.tmp && \
	mv $@.tmp $@

# introns
$(name)/data/introns.bed: $(name)/data/genes.bed $(name)/data/exons.bed
	bedtools subtract -a $< -b $(name)/data/exons.bed > $@.tmp && mv $@.tmp $@


# M
# 
define only_existing_files=
$(foreach f,$(1),$(if $(realpath $(f)),$(f) ,))
endef

define mappersFromReportPath=
$(subst /align_overall_comparison.png.tsv,,$(subst $(name)/report/mapping/,,$(1)))
endef

# only perform the comparison on the existing TSV files
$(name)/report/mapping/comparison.html: $(call only_existing_files,$(foreach m,$(call mapping_dirs,$(name)), $(name)/report/mapping/$(shell basename $(m))/align_overall_comparison.png.tsv))
	mappers_comp_sum.R --tsv "$^" --labels "$(foreach f,$^, $(call mappersFromReportPath,$(f)))" --out $(@D)/comparison --css  ../../$(CSS_FILE) && touch $@
#	mappers_comp_sum.R --tsv "$^" --labels "$(foreach m,$(call mapping_dirs,$(name)), $(shell basename $(m)))" --out $(@D)/comparison --css  ../irap.css && touch $@


# detailed info per bam
$(name)/report/mapping/$(mapper)/%/index.html: FORCE
	bam_report.R --bam $(name)/$(mapper)/$(call bam_file_for_lib,$*) -d $(@D) --fastq "$(call libname2fastq,$*)" --cores $(max_threads)

########################
phony_targets+=quant_report quant_report_files
silent_targets+=quant_report quant_report_files

quant_html_files=$(foreach q,$(SUPPORTED_QUANT_METHODS),$(foreach m,$(SUPPORTED_MAPPERS),$(foreach f,gene exon transcript,$(foreach metric,raw nlib rpkm,$(call quant_target,$(m),$(q),$(metric),$(f)) ))))


quant_report: report_setup $(quant_html_files)

quant_report_files: 
	echo $(quant_html_files)

#######################################
# Quant. at gene level
$(name)/report/quant/%/gene.raw.html: 
	$(call GE_tsv2html,"raw",$(call quiet_ls1,$(name)/$(subst _x_,/,$*)/genes.raw.$(shell echo $*|sed "s/.*_x_//").tsv),$(@D),gene.raw.t,$(subst _x_, x ,$*),"gene") && \
	cp $(subst .html,,$@).t.html $@

$(name)/report/quant/%/gene.rpkm.html: 
	$(call GE_tsv2html,"rpkm",$(call quiet_ls1,$(name)/$(subst _x_,/,$*)/genes.rpkm.*.tsv),$(@D),gene.rpkm.t,$(subst _x_, x ,$*),"gene") && \
	cp $(subst .html,,$@).t.html $@

$(name)/report/quant/%/gene.nlib.html: 
	$(call GE_tsv2html,"nlib",$(call quiet_ls1,$(name)/$(subst _x_,/,$*)/genes.nlib.*.tsv),$(@D),gene.nlib.t,$(subst _x_, x ,$*),"gene") && \
	cp $(subst .html,,$@).t.html $@


# Transcript
$(name)/report/quant/%/transcript.raw.html: 
	$(call GE_tsv2html,"raw",$(call quiet_ls1,$(name)/$(subst _x_,/,$*)/transcripts.raw.$(shell echo $*|sed "s/.*_x_//").tsv),$(@D),transcript.raw.t,$(subst _x_, x ,$*),"CDS") && \
	cp $(subst .html,,$@).t.html $@

$(name)/report/quant/%/transcript.rpkm.html: 
	$(call GE_tsv2html,"rpkm",$(call quiet_ls1,$(name)/$(subst _x_,/,$*)/transcripts.rpkm.*.tsv),$(@D),transcript.rpkm.t,$(subst _x_, x ,$*),"CDS") && \
	cp $(subst .html,,$@).t.html $@

$(name)/report/quant/%/transcript.nlib.html: 
	$(call GE_tsv2html,"nlib",$(call quiet_ls1,$(name)/$(subst _x_,/,$*)/transcripts.nlib.*.tsv),$(@D),transcript.nlib.t,$(subst _x_, x ,$*),"CDS") && \
	cp $(subst .html,,$@).t.html $@

# exon
$(name)/report/quant/%/exon.raw.html: 
	$(call GE_tsv2html,"raw",$(call quiet_ls1,$(name)/$(subst _x_,/,$*)/exons.raw.$(shell echo $*|sed "s/.*_x_//").tsv),$(@D),exon.raw.t,$(subst _x_, x ,$*),"exon") && \
	cp $(subst .html,,$@).t.html $@

$(name)/report/quant/%/exon.rpkm.html: 
	$(call GE_tsv2html,"rpkm",$(call quiet_ls1,$(name)/$(subst _x_,/,$*)/exons.rpkm.*.tsv),$(@D),exon.rpkm.t,$(subst _x_, x ,$*),"exon") && \
	cp $(subst .html,,$@).t.html $@

$(name)/report/quant/%/exon.nlib.html: 
	$(call GE_tsv2html,"nlib",$(call quiet_ls1,$(name)/$(subst _x_,/,$*)/exons.nlib.*.tsv),$(@D),exon.nlib.t,$(subst _x_, x ,$*),"exon") && \
	cp $(subst .html,,$@).t.html $@

##########################
# one rule by quant option
$(name)/report/quant/%_x_htseq1.html: $(name)/report/quant/%_x_htseq1/gene.raw.html $(name)/report/quant/%_x_htseq1/gene.rpkm.html $(name)/report/quant/%_x_htseq1/gene.nlib.html
	touch $@

$(name)/report/quant/%_x_htseq2.html: $(name)/report/quant/%_x_htseq2/gene.raw.html $(name)/report/quant/%_x_htseq2/gene.rpkm.html  $(name)/report/quant/%_x_htseq2/gene.nlib.html
	touch $@


$(name)/report/quant/%_x_flux_cap.html: $(name)/report/quant/%_x_flux_cap/gene.raw.html $(name)/report/quant/%_x_flux_cap/gene.rpkm.html $(name)/report/quant/%_x_flux_cap/gene.nlib.html
	touch $@

$(name)/report/quant/%_x_basic.html: $(name)/report/quant/%_x_basic/gene.raw.html $(name)/report/quant/%_x_basic/gene.rpkm.html $(name)/report/quant/%_x_basic/gene.nlib.html
	touch $@

$(name)/report/quant/%_x_cufflinks1_nd.html: $(name)/report/quant/%_x_cufflinks1_nd/gene.raw.html $(name)/report/quant/%_x_cufflinks1_nd/gene.rpkm.html $(name)/report/quant/%_x_cufflinks1_nd/gene.nlib.html
	touch $@

$(name)/report/quant/%_x_cufflinks1.html: $(name)/report/quant/%_x_cufflinks1/gene.raw.html $(name)/report/quant/%_x_cufflinks1/gene.rpkm.html $(name)/report/quant/%_x_cufflinks1/gene.nlib.html
	touch $@

$(name)/report/quant/%_x_cufflinks2.html: $(name)/report/quant/%_x_cufflinks2/gene.raw.html $(name)/report/quant/%_x_cufflinks2/gene.rpkm.html $(name)/report/quant/%_x_cufflinks2/gene.nlib.html
	touch $@

$(name)/report/quant/%_x_cufflinks2_nd.html: $(name)/report/quant/%_x_cufflinks2/gene.raw.html $(name)/report/quant/%_x_cufflinks2/gene.rpkm.html $(name)/report/quant/%_x_cufflinks2/gene.nlib.html
	touch $@

$(name)/report/quant/%_x_scripture.html: 
	$(info TODO: complete $@)
	touch $@

$(name)/report/quant/%_x_nurd.html: $(name)/report/quant/%_x_nurd/gene.raw.html $(name)/report/quant/%_x_nurd/gene.rpkm.html $(name)/report/quant/%_x_nurd/gene.nlib.html
	touch $@

## quant report: transcripts
$(name)/report/quant/%_x_htseq1.transcript.html: $(name)/report/quant/%_x_htseq1/transcript.raw.html $(name)/report/quant/%_x_htseq1/transcript.rpkm.html  $(name)/report/quant/%_x_htseq1/transcript.nlib.html
	touch $@

$(name)/report/quant/%_x_htseq2.transcript.html: $(name)/report/quant/%_x_htseq2/transcript.raw.html $(name)/report/quant/%_x_htseq2/transcript.rpkm.html  $(name)/report/quant/%_x_htseq2/transcript.nlib.html
	touch $@

$(name)/report/quant/%_x_flux_cap.transcript.html: $(name)/report/quant/%_x_flux_cap/transcript.raw.html $(name)/report/quant/%_x_flux_cap/transcript.rpkm.html  $(name)/report/quant/%_x_flux_cap/transcript.nlib.html
	touch $@

$(name)/report/quant/%_x_cufflinks1.transcript.html: $(name)/report/quant/%_x_cufflinks1/transcript.raw.html $(name)/report/quant/%_x_cufflinks1/transcript.rpkm.html  $(name)/report/quant/%_x_cufflinks1/transcript.nlib.html
	touch $@

$(name)/report/quant/%_x_cufflinks1_nd.transcript.html: $(name)/report/quant/%_x_cufflinks1_nd/transcript.raw.html $(name)/report/quant/%_x_cufflinks1_nd/transcript.rpkm.html  $(name)/report/quant/%_x_cufflinks1_nd/transcript.nlib.html
	touch $@

$(name)/report/quant/%_x_cufflinks2.transcript.html: $(name)/report/quant/%_x_cufflinks2/transcript.raw.html $(name)/report/quant/%_x_cufflinks2/transcript.rpkm.html  $(name)/report/quant/%_x_cufflinks2/transcript.nlib.html
	touch $@

$(name)/report/quant/%_x_cufflinks2_nd.transcript.html: $(name)/report/quant/%_x_cufflinks2_nd/transcript.raw.html $(name)/report/quant/%_x_cufflinks2_nd/transcript.rpkm.html  $(name)/report/quant/%_x_cufflinks2_nd/transcript.nlib.html
	touch $@

$(name)/report/quant/%_x_nurd.transcript.html: $(name)/report/quant/%_x_nurd/transcript.raw.html $(name)/report/quant/%_x_nurd/transcript.rpkm.html  $(name)/report/quant/%_x_nurd/transcript.nlib.html
	touch $@

## quant report: exons
$(name)/report/quant/%_x_htseq1.exon.html: $(name)/report/quant/%_x_htseq1/exon.raw.html $(name)/report/quant/%_x_htseq1/exon.rpkm.html  $(name)/report/quant/%_x_htseq1/exon.nlib.html
	touch $@
$(name)/report/quant/%_x_htseq2.exon.html: $(name)/report/quant/%_x_htseq2/exon.raw.html $(name)/report/quant/%_x_htseq2/exon.rpkm.html  $(name)/report/quant/%_x_htseq2/exon.nlib.html
	touch $@
$(name)/report/quant/%_x_flux_cap.exon.html: $(name)/report/quant/%_x_flux_cap/exon.raw.html $(name)/report/quant/%_x_flux_cap/exon.rpkm.html  $(name)/report/quant/%_x_flux_cap/exon.nlib.html
	touch $@

$(name)/report/quant/%_x_cufflinks1.exon.html: $(name)/report/quant/%_x_cufflinks1/exon.raw.html $(name)/report/quant/%_x_cufflinks1/exon.rpkm.html  $(name)/report/quant/%_x_cufflinks1/exon.nlib.html
	touch $@

$(name)/report/quant/%_x_cufflinks1_nd.exon.html: $(name)/report/quant/%_x_cufflinks1_nd/exon.raw.html $(name)/report/quant/%_x_cufflinks1_nd/exon.rpkm.html  $(name)/report/quant/%_x_cufflinks1_nd/exon.nlib.html
	touch $@

$(name)/report/quant/%_x_cufflinks2.exon.html: $(name)/report/quant/%_x_cufflinks2/exon.raw.html $(name)/report/quant/%_x_cufflinks2/exon.rpkm.html  $(name)/report/quant/%_x_cufflinks2/exon.nlib.html
	touch $@

$(name)/report/quant/%_x_cufflinks2_nd.exon.html: $(name)/report/quant/%_x_cufflinks2_nd/exon.raw.html $(name)/report/quant/%_x_cufflinks2_nd/exon.rpkm.html  $(name)/report/quant/%_x_cufflinks2_nd/exon.nlib.html
	touch $@

$(name)/report/quant/%_x_nurd.exon.html: $(name)/report/quant/%_x_nurd/exon.raw.html $(name)/report/quant/%_x_nurd/exon.rpkm.html  $(name)/report/quant/%_x_nurd/exon.nlib.html
	touch $@

############################
# DE
phony_targets+=de_report de_report_files
silent_targets+=de_report_files

de_report: report_setup $(call de_html_files,$(name))

# just print the name of the files that will be produced
de_report_files:
	echo $(call de_html_files,$(name))


$(name)/report/%.genes_de.html: $(name)/%.genes_de.tsv $(annot_tsv)
	mkdir -p $(@D)
	$(call DE_tsv2html,$(subst _nd,,$(call DEfilepath2demethod,$@)),$<,$(@D),$(subst .html,,$(shell basename $@)),$(subst /, x ,$*))

############################
# GSE
phony_targets+=gse_report gse_report_files
silent_targets+=gse_report_files

# only generates the html iff the respective GSE tsv file exist
gse_report: report_setup $(call gse_html_files,$(name))
# just print the name of the files that will be produced
gse_report_files:
	echo $(call gse_html_files,$(name))


############################
# all targets
phony_targets+= report_all_targets
report_all_targets:  report_setup $(info_targets) qc_report mapping_report quant_report de_report gse_report  end_report $(name)/report/about.html



$(name)/report/about.html: 
	cat $(IRAP_DIR)/aux/html/page.header.html $(IRAP_DIR)/aux/html/about.html  $(IRAP_DIR)/aux/html/page.footer.html >  $@


#########################
#mapping_report de_report
# TODO: remove/fix this in the future (currently necessary to update the menu)
phony_targets+=end_report

end_report: $(name)/report/index.html $(call must_exist,$(name)/report/irap.css)


# TODO: replace versions.html by info_report
# TODO $(call must_exist,$(name)/report/status.html)a
$(name)/report/index.html: $(conf) $(info_targets)  $(quant_html_files) $(qc_html_files) $(call mapping_report_targets) $(call de_html_files,$(name)) $(call gse_html_files,$(name))  $(call rep_browse,$(name)/report/jbrowse/index.html)  $(name)/report/about.html $(call must_exist,$(name)/report/irap.css) $(call must_exist,$(name)/report/menu.css)
	cp  $(name)/report/info.html $@ &&
	irap_report_main $(IRAP_REPORT_MAIN_OPTIONS) --conf $(conf) --rep_dir $(name)/report -m "$(call mapping_dirs,$(name))" -q "$(call quant_dirs,$(name))" -d "$(call de_dirs,$(name))" &&
	sleep 2 &&
	touch $@
