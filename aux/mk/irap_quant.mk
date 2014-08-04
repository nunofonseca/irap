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
#    $Id: 0.1.3 Nuno Fonseca Wed Dec 26 16:16:19 2012$
# =========================================================

#*****************
# Cufflinks 1 & 2
#*****************
ifndef cufflinks1_params
	cufflinks1_params=--max-bundle-frags 90000000
endif
ifndef cufflinks2_params
	cufflinks2_params=--max-bundle-frags 90000000
endif

#--FDR $(de_pvalue_cuffoff)
#--multi-read-correct more accurate
cufflinks1_params+= --no-update-check --min-isoform-fraction 0.05 --multi-read-correct
cufflinks2_params+= $(cufflinks1_params)

# cuffdiff params
ifndef cuffdiff1_params=
cuffdiff1_params=
endif

ifndef cuffdiff2_params=
cuffdiff2_params=
endif
# params
# 1- bam file
# 2- library name 
# 3- .cuff.gtf filename
# 4- .isoforms.fpkm_tracking
# 5- .genes.fpkm_tracking
# Use irap_mapper.sh for now to set up the environment
define cufflinks_g_option=
$(if $(findstring _nd,$(quant_method)),-G,-g)
endef

define run_cufflinks1=
	irap_wrapper.sh cufflinks1 cufflinks -p $(max_threads) -o $(call lib2quant_folder,$(2))$(2) $(call cufflinks_g_option) $(gtf_file_abspath) $(cufflinks1_params) $(1) && mv $(call lib2quant_folder,$(2))$(2)/transcripts.gtf $(3) && mv $(call lib2quant_folder,$(2))$(2)/isoforms.fpkm_tracking $(4) && mv $(call lib2quant_folder,$(2))$(2)/genes.fpkm_tracking $(5)
endef

define run_cufflinks2=
	irap_wrapper.sh cufflinks2 cufflinks -p $(max_threads) -o $(call lib2quant_folder,$(2))$(2)  $(call cufflinks_g_option) $(gtf_file_abspath) $(cufflinks2_params) $(1) && mv $(call lib2quant_folder,$(2))$(2)/transcripts.gtf $(3) && mv $(call lib2quant_folder,$(2))$(2)/isoforms.fpkm_tracking $(4) && mv $(call lib2quant_folder,$(2))$(2)/genes.fpkm_tracking $(5)
endef

# run_cuffmerge cufflinks1|cufflinks2 
define run_cuffmerge=
        irap_wrapper.sh $(subst _nd,,$(1)) cuffmerge -o $(name)/$(mapper)/$(quant_method)/ -g $(gtf_file_abspath) -s $(reference_abspath) -p $(max_threads) $(name)/$(mapper)/$(quant_method)/$(name).assemblies.txt  && mv $(name)/$(mapper)/$(quant_method)/merged.gtf $(2)
endef

# run_cuffdiff cufflinks1|cufflinks2  cmd
define run_cuffdiff=
        irap_wrapper.sh $(1) $(2)
endef


###########################################################
#
ifndef ireckon_params=
	ireckon_params=
endif

#--multi-read-correct more accurate
ireckon_params+= -m bwa  -n $(max_threads)

# 1-BAM
# 2-ref
# 3-bed  preprocessed with FormatTools
define run_ireckon=
	irap_wrapper.sh ireckon ireckon $(1)  $(2) -o $(name)/$(mapper)/ireckon/$(2) $(ireckon_params) &&\
	mv $(name)/$(mapper)/ireckon/$(2)/result.gtf $(5)
#output_dir/
# $(call ireckon_lib_option) $(gtf_file_abspath) $(cufflinks1_params)
endef


#******************
# Flux capacitator
#******************
ifndef fluxcap_mem
fluxcap_mem=$(max_mem_gb)G
endif

# make it simple
define get_read_descriptor=
$(if $(filter $(1),PAIRED),$(1),SIMPLE)
endef
#--annotation-mapping $(3) --read-descriptor $(call get_read_descriptor,$(3)) 
# --annotation-mapping (PAIRED|SINGLE)
#1 bam file
#2 gtf
#3 SINGLE|PAIRED
#4 output
# bam file the needs to have the bam file indexed (sorted by reference ID and then the leftmost coordinate before indexing)
# flux_cap requires a file with some parameters :(
define run_flux_cap=
	rm -f $(1).par &&\
	echo COVERAGE_STATS true > $(1).par &&\
	echo COVERAGE_FILE `basename $(1).cov` >> $(1).par &&\
	echo "transcript	locus	gene	reads	length	RPKM"> $(4).tmp2 &&\
	TMP_DIR=$(tmp_dir) FLUX_MEM=$(fluxcap_mem) flux-capacitor --force --annotation $(2)  -i $(1) --threads $(max_threads) -p $(1).par -o $(4).gtf && \
	sed  "s/.*transcript_id/transcript_id/" $(4).gtf | sed "s/\"//g" | sed "s/;/\t/g" | sed -E "s/([a-zA-Z]*_id|reads|length|RPKM| )//g" >> $(4).tmp2 && \
	mv $(4).tmp2 $(4) 
endef


#************
# HTSeq-count
#************

# 
ifndef htseq_params
htseq_params= -q  --stranded=no
endif


define htseq_sam_output_param=
$(if $(filter y,$(htseq_sam_output)),--samout=$(subst .tsv,.anot.sam,$(1)),)
endef

define htseq_id_param=
$(if $(filter $(1),gene),--idattr gene_id,$(if $(findstring $(1),exon),--idattr exon_id,--idattr transcript_id))
endef 

define htseq_mode_param=
$(if $(filter $(1),htseq1),-m union,-m intersection-nonempty)
endef
# BAM file needs to be sorted by read name
# run_htseq, bam file, gtf, output file
# mode:union
#
#1- BAM file
#2- GTF/GFF file
#3- output file 
#4- mode
#5- exon|gene|transctript 
define run_htseq=
	samtools view $(1) | \
	htseq-count $(call htseq_sam_output_param,$(3))  $(call htseq_mode_param,$(4))  $(call htseq_id_param,$(5))  $(htseq_params) - $(2)  > $(3).tmp && \
	tail -n -5 $(3).tmp > $(3).$(4).stats &&\
	head -n -5 $(3).tmp > $(3).tmp2 &&\
	mv $(3).tmp2 $(3) && rm -f $(3).tmp
endef

# 
#1- BAM file
#2- GTF/GFF file
#3- output file 
#4- mode
#5- exon|gene|transtript 
define run_htseq1=
	$(call run_htseq,$(1),$(2),$(3),htseq1,$(4))
endef

#mode: intersection-nonempty
define run_htseq2=
	$(call run_htseq,$(1),$(2),$(3),htseq2,$(4))
endef

#******************
# NURD
#******************

ifndef nurd_params
nurd_params= 
endif

# 1 bam file
# 2 gtf
# 3 nurd outgenes.reads.out
# bam file needs to be sorted by name
define run_nurd=
        samtools view $(1) > $(1).sam && \
	NURD $(nurd_params) -O $(3).dir -G $(2) -S $(1).sam && \
	rm -f $(1).sam && \
	mv $(3).dir/`basename $(1).sam`.nurd.all_expr $(3) 
endef
# ideally nurd should support BAM
# or accept SAM from stdin...

#******************
# IsoEM
#******************

ifndef isoem_params
isoem_params= 
endif

#sample/genome_aligned_reads.iso_estimates     		- estimated FPKMs (Fragments Per Kilobase per Million reads) for all the isoforms in the GTF file
#sample/genome_aligned_reads.gene_estimates    		- estimated FPKMs for genes
#sample/genome_aligned_reads_iso_read_coverage.bed	- isoforms coverage by reads

#1 bam file
#2 gtf
#3 SINGLE|PAIRED
#4 output
# bam file needs to be sorted by name
define run_isoem=
endef


#************
# 1=bam 2=gtf 3=out gene_file 4=out exon_file
define run_naive_count=
	irap_naive_count.pl  $(2) $(1) $(4).tmp $(3).tmp && mv $(3).tmp  $(3) && mv $(4).tmp  $(4)
endef


#*****************
# Cufflinks 1 & 2
#*****************
# Assemble and quantify the transcripts
# 1-prefix
# 2-odir
# 3-tracking file
# 4-cufflinks_flavour
define cufflinks2counts=
	irap_cufflinksFPKMs2counts -c $(3) --label "$(1)" -r $($(subst .pe,,$(subst .se,,$(1)))_rs) -o $(2)/$(1).genes.raw.$(4).tsv -i $(2)/$(1).transcripts.raw.$(4).tsv
endef

# call we avoid two functions?
define cufflinks2counts_g=
	irap_cufflinksFPKMs2counts -c $(3) --label "$(1)" -r $($(subst .pe,,$(subst .se,,$(1)))_rs) -o $(2)/$(1).genes.raw.$(4).tsv -i /dev/null
endef

define cufflinks2counts_t=
	irap_cufflinksFPKMs2counts -c $(3) --label "$(1)" -r $($(subst .pe,,$(subst .se,,$(1)))_rs) -o /dev/null -i $(2)/$(1).transcripts.raw.$(4).tsv
endef


#######
# The SAM file supplied to Cufflinks must be sorted by reference position. 
# -M/--mask-file <mask.(gtf/gff)> 
# novel genes=-g gtf
# only known=-G gtf
# $1 - lib
# $2 - cufflinks1|cufflinks1_nd|cuffflinks2|...
# $3 - bam file prefix (includes .se|.pe)
define make-cufflinks-quant-rule=
$(call lib2quant_folder,$(1))$(3).cuff.gtf $(call lib2quant_folder,$(1))$(3).cuff.genes.fpkm_tracking $(call lib2quant_folder,$(1))$(3).cuff.isoforms.fpkm_tracking: $(call lib2bam_folder,$(1))$(3).hits.bam
	$(call run_$(subst _nd,,$(2)),$$<,$(1),$$@,$(call lib2quant_folder,$(1))$(3).cuff.isoforms.fpkm_tracking,$(call lib2quant_folder,$(1))$(3).cuff.genes.fpkm_tracking)


$(call lib2quant_folder,$(1))$(3).transcripts.raw.$(quant_method).tsv: $(call lib2quant_folder,$(1))$(3).cuff.isoforms.fpkm_tracking
	$$(call p_info,Warning! Cufflinks produces FPKMS and does not provide counts. Generating pseudo counts $$@.)
	$$(call cufflinks2counts_t,$(3),$$(@D),$$^,$(quant_method))

$(call lib2quant_folder,$(1))$(3).genes.raw.$(quant_method).tsv: $(call lib2quant_folder,$(1))$(3).cuff.isoforms.fpkm_tracking
	$$(call p_info,Warning! Cufflinks produces FPKMS and does not provide counts. Generating pseudo counts $$@.)
	$$(call cufflinks2counts_g,$(3),$$(@D),$$^,$(quant_method))

endef

# generate the rules for htseq
ifeq ($(patsubst cufflinks%,,$(quant_method)),)
$(foreach l,$(se),$(eval $(call make-cufflinks-quant-rule,$(l),$(quant_method),$(l).se,gene)))	
$(foreach l,$(pe),$(eval $(call make-cufflinks-quant-rule,$(l),$(quant_method),$(l).pe,gene)))

# cufflinks* specific rules

$(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) )  > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/%xons.raw.$(quant_method).tsv: 
	$(call p_info,Warning! Cufflinks produces FPKMS and does not provide counts. Generating empty file $@.)
	@$(call empty_file,$@)

$(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) ) > $@.tmp && mv $@.tmp $@

# Run Cuffmerge on all your assemblies to create a single merged transcriptome annotation
$(name)/$(mapper)/$(quant_method)/$(name).merged.gtf:  $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.cuff.gtf) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.cuff.gtf) $(name)/$(mapper)/$(quant_method)/$(name).assemblies.txt
	$(call run_cuffmerge,$(quant_method),$@)

# cuffmerge' assemblies.txt - lists the assembly file for each sample. 
$(name)/$(mapper)/$(quant_method)/$(name).assemblies.txt: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.cuff.gtf) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.cuff.gtf)
	rm -f $@
	touch $@.tmp && for f in $^; do echo $$f >> $@.tmp ; done && mv $@.tmp $@

endif

#*****************
# HTSeq
#*****************

# * HTSeq1
# counts per gene
$(name)/$(mapper)/htseq1/genes.raw.htseq1.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@


# counts per transcript
$(name)/$(mapper)/htseq1/transcripts.raw.htseq1.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).tsv)	
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

# counts per exon
$(name)/$(mapper)/htseq1/exons.raw.htseq1.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.exons.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.exons.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

## htseq bam file needs to be sorted by name
# $1 - lib
# $2 - htseq1|htseq2
# $3 - bam file prefix (includes .se|.pe)
# $4 - quantification of exon|gene|trans
# $5 - gtf file
define make-htseq-quant-rule=
$(call lib2quant_folder,$(1))$(3).$(4)s.raw.$(2).tsv: $(call lib2bam_folder,$(1))$(3).hits.byname.bam $(5)
	mkdir -p $$(@D) && $$(call run_$(2),$$<,$(5),$$@,$(4))
endef

# generate the rules for htseq
ifeq ($(patsubst htseq%,,$(quant_method)),)
$(foreach l,$(se),$(eval $(call make-htseq-quant-rule,$(l),$(quant_method),$(l).se,gene,$(gtf_file_abspath))))	
$(foreach l,$(pe),$(eval $(call make-htseq-quant-rule,$(l),$(quant_method),$(l).pe,gene,$(gtf_file_abspath))))
ifeq ($(exon_quant),y)
$(foreach l,$(se),$(eval $(call make-htseq-quant-rule,$(l),$(quant_method),$(l).se,exon,$(gtf_file_abspath).exon_id.gtf)))	
$(foreach l,$(pe),$(eval $(call make-htseq-quant-rule,$(l),$(quant_method),$(l).pe,exon,$(gtf_file_abspath).exon_id.gtf)))
endif
ifeq ($(transcript_quant),y)
$(foreach l,$(se),$(eval $(call make-htseq-quant-rule,$(l),$(quant_method),$(l).se,transcript,$(gtf_file_abspath))))	
$(foreach l,$(pe),$(eval $(call make-htseq-quant-rule,$(l),$(quant_method),$(l).pe,transcript,$(gtf_file_abspath))))
endif
#$(foreach l,$(se) $(pe),$(info $(call make-htseq-quant-rule,$(l),htseq1,$(l).pe,gene)))
endif

# * HTSeq2

# counts per gene
$(name)/$(mapper)/htseq2/genes.raw.htseq2.tsv: $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.genes.raw.htseq2.tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@
# counts per transcript
$(name)/$(mapper)/htseq2/transcripts.raw.htseq2.tsv: $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

# counts per exon
$(name)/$(mapper)/htseq2/exons.raw.htseq2.tsv: $(foreach p, $(pe),$(call lib2quant_folder,$(p))$(p).pe.exons.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.exons.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@


#*****************
# Simple count
#*****************

# TODO: from exon2gene
$(name)/$(mapper)/basic/genes.raw.basic.tsv: $(foreach p,$(pe),$(name)/$(mapper)/$(quant_method)/$(p).genes.raw.$(quant_method).tsv) $(foreach s,$(se), $(name)/$(mapper)/$(quant_method)/$(s).genes.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/basic/transcripts.raw.basic.tsv: 
	$(call p_info, Warning! Basic method does not produce quantification at transcript level (in TODO). Generating empty file $@.)
	@$(call empty_file,$@)

$(name)/$(mapper)/basic/exons.raw.basic.tsv: 
	$(call p_info, Warning! Basic method does not produce quantification at exon level (in TODO). Generating empty file $@.)
	@$(call empty_file,$@)

#$(name)/$(mapper)/basic/exons.rpkm.basic.tsv: $(name)/$(mapper)/basic/exons.raw.basic.tsv $(feat_length)
#	irap_raw2metric --tsv $<  --lengths $(feat_length) --feature exon --metric rpkm --out $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/basic/%.genes.raw.basic.tsv $(name)/$(mapper)/basic/%.exons.raw.basic.tsv: $(name)/$(mapper)/%.se.hits.byname.bam $(gtf_file_abspath).exon_id.gtf
	$(call run_naive_count,$<,$(gtf_file_abspath).exon_id.gtf,$(@D)/$*.genes.raw.basic.tsv,$(@D)/$*.exons.raw.basic.tsv)

$(name)/$(mapper)/basic/%.genes.raw.basic.tsv $(name)/$(mapper)/basic/%.exons.raw.basic.tsv: $(name)/$(mapper)/%.pe.hits.byname.bam $(gtf_file_abspath).exon_id.gtf
	$(call run_naive_count,$<,$(gtf_file_abspath).exon_id.gtf,$(@D)/$*.genes.raw.basic.tsv,$(@D)/$*.exons.raw.basic.tsv)

$(name)/$(mapper)/basic/%.transcripts.raw.basic.tsv: $(name)/$(mapper)/%.pe.hits.byname.bam
	$(call p_info, Warning! Basic method does not produce quantification at transcript level (in TODO). Generating empty file $@.)
	@$(call empty_file,$@)

$(name)/$(mapper)/basic/%.transcripts.raw.basic.tsv: $(name)/$(mapper)/%.se.hits.byname.bam
	$(call p_info, Warning! Basic method does not produce quantification at transcript level (in TODO). Generating empty file $@.)
	@$(call empty_file,$@)

#*****************
# Flux-capacitor
#*****************

$(name)/$(mapper)/flux_cap/genes.raw.flux_cap.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/flux_cap/exons.raw.flux_cap.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.exons.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.exons.raw.$(quant_method).tsv)
	@$(call empty_file,$@)

$(name)/$(mapper)/flux_cap/transcripts.raw.flux_cap.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

#
# Run flux
#
# $1 - lib
# $2 - SINGLE|PE
# $3 - bam file prefix (includes .se|.pe)
define make-fluxcap-quant-rule=
$(call lib2quant_folder,$(1))$(3).flux_cap.tsv: $(call lib2bam_folder,$(1))$(3).hits.byname.bam $(gtf_file_abspath) $(call lib2bam_folder,$(1))$(3).hits.byname.bam $(call lib2bam_folder,$(1))$(3).hits.bam.bai 
	mkdir -p $$(@D) && $$(call run_flux_cap,$$<,$$(gtf_file_abspath),$(2),$$@)
endef

ifeq (flux_cap,$(quant_method))
$(foreach l,$(se),$(eval $(call make-fluxcap-quant-rule,$(l),SINGLE,$(l).se)))	
$(foreach l,$(pe),$(eval $(call make-fluxcap-quant-rule,$(l),PAIRED,$(l).pe)))
endif

# Process flux output
%.genes.raw.flux_cap.tsv: %.flux_cap.tsv 
	tsv.aggrbygene $< counts $@.tmp && mv $@.tmp $@ 

%.genes.rpkm.flux_cap.tsv: %.flux_cap.tsv
	tsv.aggrbygene $< RPKM $@.tmp && mv  $@.tmp $@

%.transcripts.raw.flux_cap.tsv: %.flux_cap.tsv
	cut -f 1,4 $< |tail -n +2 > $@.tmp && mv  $@.tmp $@

%.transcripts.rpkm.flux_cap.tsv: %.flux_cap.tsv
	cut -f 1,6 $< | tail -n +2  > $@.tmp && mv  $@.tmp $@

%.exons.raw.flux_cap.tsv: 
	$(call p_info, Warning! Flux_capacitor does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)

%.exons.rpkm.flux_cap.tsv: 
	$(call p_info, Warning! Flux_capacitor does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)


#*****************
# NURD
#*****************

$(name)/$(mapper)/nurd/genes.raw.nurd.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/nurd/exons.raw.nurd.tsv: 
	$(call p_info, Warning! NURD does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)

$(name)/$(mapper)/nurd/transcripts.raw.nurd.tsv: 
#	$(call p_info, Warning! NURD does not produce raw quantification at transcript level. Generating empty file $@.)
	@$(call empty_file,$@)

#
# Run flux
#
# $1 - lib
# $2 - bam file prefix (includes .se|.pe)
define make-nurd-quant-rule=
$(call lib2quant_folder,$(1))$(2).nurd.tsv: $(call lib2bam_folder,$(1))$(2).hits.bam $(gtf_file_abspath) 
	$$(call run_nurd,$$<,$$(gtf_file_abspath),$$@)
# Process nurd output
$(call lib2quant_folder,$(1))$(2).genes.raw.nurd.tsv: $(call lib2quant_folder,$(1))$(2).nurd.tsv 
	cut -f 1,3 $$< > $$@.tmp && mv $$@.tmp $$@

$(call lib2quant_folder,$(1))$(2).genes.rpkm.nurd.tsv: $(call lib2quant_folder,$(1))$(2).nurd.tsv
	cut -f 1,6 $$< > $$@.tmp && mv $$@.tmp $$@


# transcripts/isoforms
$(call lib2quant_folder,$(1))$(2).transcripts.raw.nurd.tsv: $(call lib2quant_folder,$(1))$(2).nurd.tsv
	irap_nurd2tsv --tsv $$< --out $$@.tmp && mv $$@.tmp $$@

$(call lib2quant_folder,$(1))$(2).transcripts.rpkm.nurd.tsv: $(call lib2quant_folder,$(1))$(2).nurd.tsv
	irap_nurd2tsv --tsv $$< --out $$@.tmp && mv $$@.tmp $$@

endef

ifeq (nurd,$(quant_method))
$(foreach l,$(se),$(eval $(call make-nurd-quant-rule,$(l),$(l).se)))	
$(foreach l,$(pe),$(eval $(call make-nurd-quant-rule,$(l),$(l).pe)))
endif


# exons
%.exons.raw.nurd.tsv: 
	$(call p_info, Warning! NURD does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)

%.exons.rpkm.nurd.tsv: 
	$(call p_info, Warning! NURD does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)


