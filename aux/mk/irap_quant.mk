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

## Notes:
## 1 - Always produce gene expression quantification

## Needed for transcript quantification
mapTrans2gene=$(name)/data/$(gtf_file_basename).mapTrans2Gene.tsv

$(mapTrans2gene): $(gtf_file_abspath)
	genMapTrans2Gene -i $< -o $@.tmp -c $(max_threads) && mv $@.tmp $@

## Output files produced
STAGE3_S_OFILES+= $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).$(expr_ext)) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).$(expr_ext)) 

## Quantification statistics
a_quant_qc_stats:=
quant_qc_stats:=$(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).quant_qc.tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).quant_qc.tsv)

a_quant_qc_stats+= $(quant_qc_stats)

ifeq ($(transcript_expr),y)
quant_qc_statst:= $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).transcripts.raw.$(quant_method).quant_qc.tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).quant_qc.tsv)
a_quant_qc_stats+= $(quant_qc_statst)

STAGE3_S_OFILES+= $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).$(expr_ext)) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).$(expr_ext)) 
endif

ifeq ($(exon_quant),y)
quant_qc_statse:= $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).exons.raw.$(exon_quant_method).quant_qc.tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.exons.raw.$(exon_quant_method).quant_qc.tsv)
a_quant_qc_stats+= $(quant_qc_statse)

STAGE3_S_OFILES+= $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.exons.raw.$(exon_quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.exons.raw.$(exon_quant_method).tsv) 
endif


#*****************
# Cufflinks 1 & 2
#*****************
ifndef cufflinks1_params
	cufflinks1_params=--max-bundle-frags 190000000 --min-isoform-fraction 0.05 --multi-read-correct
endif
ifndef cufflinks2_params
	cufflinks2_params=--max-bundle-frags 190000000 --min-isoform-fraction 0.05 --multi-read-correct
endif

#--FDR $(de_pvalue_cuffoff)
#--multi-read-correct more accurate
cufflinks1_params+= --no-update-check 
cufflinks2_params+= --no-update-check 

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
	irap_wrapper.sh cufflinks1 cufflinks -p $(max_threads) -o $(call lib2quant_folder,$(2))$(2) $(call cufflinks_g_option) $(gtf_file_abspath) $(call tophat_strand_params,$(2)) $(cufflinks1_params) $(1) && mv $(call lib2quant_folder,$(2))$(2)/transcripts.gtf $(3) && mv $(call lib2quant_folder,$(2))$(2)/isoforms.fpkm_tracking $(4) && mv $(call lib2quant_folder,$(2))$(2)/genes.fpkm_tracking $(5)
endef

define run_cufflinks2=
	irap_wrapper.sh cufflinks2 cufflinks -p $(max_threads) -o $(call lib2quant_folder,$(2))$(2)  $(call cufflinks_g_option) $(gtf_file_abspath) $(call tophat_strand_params,$(2)) $(cufflinks2_params) $(1) && mv $(call lib2quant_folder,$(2))$(2)/transcripts.gtf $(3) && mv $(call lib2quant_folder,$(2))$(2)/isoforms.fpkm_tracking $(4) && mv $(call lib2quant_folder,$(2))$(2)/genes.fpkm_tracking $(5)
endef

# run_cuffmerge cufflinks1|cufflinks2 
define run_cuffmerge=
        irap_wrapper.sh $(subst _nd,,$(1)) cuffmerge -o $(name)/$(mapper)/$(quant_method)/ -g $(gtf_file_abspath) -s $(reference_abspath) -p $(max_threads) $(name)/$(mapper)/$(quant_method)/$(name).assemblies.txt  && mv $(name)/$(mapper)/$(quant_method)/merged.gtf $(2)
endef

# run_cuffdiff cufflinks1|cufflinks2  cmd
define run_cuffdiff=
        irap_wrapper.sh $(1) $(2)
endef


#*****************
# Stringtie
#*****************
stringtie_params?=

stringtie_params+= 

## library specific parameters
##--rf	Assumes a stranded library fr-firststrand.
##--fr	Assumes a stranded library fr-secondstrand.
# 1 - libname
stringtie_lib_params=$(if $(findstring $($(1)_strand),first),--rf,$(if $(findstring $($(1)_strand),second),--fr,))


# params
# 1- bam file
# 2- library name 
# 3- quant parameter -e for quantification only
# 4- out: gtf_filename
# Use irap_mapper.sh for now to set up the environment
# output:
# $(2).t_data.ctab - fpkm per transcript|gene
# $(2).gtf - fpkm per transcript/exon
# use -l CUFF to reuse the cufflinks code for novel transcripts
define run_stringtie=
	mkdir -p $(call lib2quant_folder,$(2))$(2) && irap_wrapper.sh stringtie stringtie  -p $(max_threads) $(3)  -o $(call lib2quant_folder,$(2))$(2)/$(2).tmp.gtf -B   -G $(gtf_file_abspath) $(call stringtie_lib_params,$(2))  $(stringtie_params) $(1)  && rename .tmp. . $(call lib2quant_folder,$(2))$(2)/$(2)*.tmp* && mv $(call lib2quant_folder,$(2))$(2)/$(2).gtf $(4)
endef
#-C $(call lib2quant_folder,$(2))$(2)/$(2).cov.tmp.gtf

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
	echo "transcript	locus	gene	reads	length	FPKM"> $(4).tmp2 &&\
	TMP_DIR=$(tmp_dir) FLUX_MEM=$(fluxcap_mem) flux-capacitor --force --annotation $(2)  -i $(1) --threads $(max_threads) -p $(1).par -o $(4).gtf && \
	sed  "s/.*transcript_id/transcript_id/" $(4).gtf | sed "s/\"//g" | sed "s/;/\t/g" | sed -E "s/([a-zA-Z]*_id|reads|length|FPKM| )//g" >> $(4).tmp2 && \
	mv $(4).tmp2 $(4) 
endef


#************************************
# bulkRNA-seq 
#************************************

#************
# HTSeq-count
#************

# 
htseq_params?=
htseq_params+= -q  


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
#6- lib
define run_htseq=
	samtools view -F 4 $(1) | \
	htseq-count  $(call htseq_sam_output_param,$(3))  $(call htseq_mode_param,$(4)) $(call htseq_strand_params,$(6)) $(call htseq_id_param,$(5))  $(htseq_params) - $(2)  > $(3).tmp && \
	tail -n -5 $(3).tmp > $(3).$(4).stats &&\
	head -n -5 $(3).tmp > $(3).tmp2 &&\
	mv $(3).tmp2 $(3) && rm -f $(3).tmp
endef

# 
#1- BAM file
#2- GTF/GFF file
#3- output file 
#4- exon|gene|transtript 
#5- libname
define run_htseq1=
	$(call run_htseq,$(1),$(2),$(3),htseq1,$(4),$(5))
endef

#mode: intersection-nonempty
define run_htseq2=
	$(call run_htseq,$(1),$(2),$(3),htseq2,$(4),$(5))
endef

# strand option
# irap
# first  -> --library-type=fr-firststrand 
# second -> --library-type=fr-secondstrand riu
irap_strand2htseqoption=$(if $(findstring $(1),first),reverse,$(if $(findstring $(1),second),yes,no))


# 1 - libname
define htseq_strand_params=
	$(if $(call lib_strand_info,$(1)),--stranded=$(call irap_strand2htseqoption,$($(1)_strand)),)
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
# 1.0.9 or earlier
# define run_nurd=
#         samtools view $(1) > $(1).sam && \
# 	NURD $(nurd_params) -O $(3).dir -G $(2) -S $(1).sam && \
# 	rm -f $(1).sam && \
# 	mv $(3).dir/`basename $(1).sam`.nurd.all_expr $(3) 
# endef
define run_nurd=
        samtools view $(1) > $(1).sam && \
	NURD $(nurd_params) -O $(3).dir -G $(2) -S $(1).sam && \
	rm -f $(1).sam && \
	mv $(3).dir/`basename $(1).sam`.nurd.rpkm $(3).nurd.tsv  && \
	mv $(3).dir/`basename $(1).sam`.nurd.rdcnt $(3).raw.nurd.tsv 
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

# call to avoid two functions?
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
	$(call run_$(subst _nd,,$(2)),$$<,$(1),$(call lib2quant_folder,$(1))$(3).cuff.gtf,$(call lib2quant_folder,$(1))$(3).cuff.isoforms.fpkm_tracking,$(call lib2quant_folder,$(1))$(3).cuff.genes.fpkm_tracking)


$(call lib2quant_folder,$(1))$(3).transcripts.raw.$(quant_method).tsv: $(call lib2quant_folder,$(1))$(3).cuff.isoforms.fpkm_tracking
	$$(call p_info,Warning! Cufflinks produces FPKMS and does not provide counts. Generating pseudo counts $$@.)
	$$(call cufflinks2counts_t,$(3),$$(@D),$$^,$(quant_method))

$(call lib2quant_folder,$(1))$(3).genes.raw.$(quant_method).tsv: $(call lib2quant_folder,$(1))$(3).cuff.isoforms.fpkm_tracking
	$$(call p_info,Warning! Cufflinks produces FPKMS and does not provide counts. Generating pseudo counts $$@.)
	$$(call cufflinks2counts_g,$(3),$$(@D),$$^,$(quant_method))

endef

# generate the rules for cufflinks
ifeq ($(patsubst cufflinks%,,$(quant_method)),)
$(foreach l,$(se),$(eval $(call make-cufflinks-quant-rule,$(l),$(quant_method),$(l).se,gene)))	
$(foreach l,$(pe),$(eval $(call make-cufflinks-quant-rule,$(l),$(quant_method),$(l).pe,gene)))

# cufflinks* specific rules
$(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@

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

############
# cufflinks
ifeq ($(patsubst cufflinks%,,$(quant_norm_tool)),)

$(name)/$(mapper)/$(quant_method)/%.genes.fpkm.$(quant_method).$(quant_norm_tool).tsv: $(name)/$(mapper)/$(quant_method)/%.cuff.genes.fpkm_tracking
	cut -f 1,10 $<  | tail -n +2  |grep -v CUFF > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/%.exons.fpkm.$(quant_method).$(quant_norm_tool).tsv: $(name)/$(mapper)/$(quant_method)/%.cuff.gtf
	grep -w exon $<  | grep -v CUFF | sed -e 's/.*gene_id .\([^\"]*\).; transcript_id .\([^\"]*\).; exon_number .\([^\"]*\).; FPKM .\([^\"]*\)..*/\1.\2.\3\t\4/' > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/%.transcripts.fpkm.$(quant_method).$(quant_norm_tool).tsv: $(name)/$(mapper)/$(quant_method)/%.cuff.isoforms.fpkm_tracking
	cut -f 1,10 $<  | tail -n +2  |grep -v CUFF  > $@.tmp && mv $@.tmp $@


$(name)/$(mapper)/$(quant_method)/genes.fpkm.cufflinks%.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.fpkm.$(quant_method).$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.genes.fpkm.$(quant_method).$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) ) > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/transcripts.fpkm.cufflinks%.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.fpkm.$(quant_method).$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.transcripts.fpkm.$(quant_method).$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) ) > $@.tmp && mv $@.tmp $@

ifeq ($(exon_quant),y)
$(name)/$(mapper)/$(quant_method)/exons.fpkm.$(exon_quant_method).cufflinks%.tsv:
	$(call p_info, Warning! Cufflinks does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)
endif
endif

endif

#*****************
# HTSeq
#*****************

# * HTSeq1
# counts per gene
$(name)/$(mapper)/htseq1/genes.raw.htseq1.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

# * HTSeq2
# counts per gene
$(name)/$(mapper)/htseq2/genes.raw.htseq2.tsv: $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.genes.raw.htseq2.tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@


## htseq bam file needs to be sorted by name
# $1 - lib
# $2 - htseq1|htseq2
# $3 - bam file prefix (includes .se|.pe)
# $4 - quantification of exon|gene|trans
# $5 - gtf file
define make-htseq-quant-rule=
$(call lib2quant_folder,$(1))$(3).$(4)s.raw.$(2).tsv: $(call lib2bam_folder,$(1))$(3).hits.byname.bam $(5)
	mkdir -p $$(@D) && $$(call run_$(2),$$<,$(5),$$@,$(4),$(1))
endef


# generate the rules for htseq
ifeq ($(patsubst htseq%,,$(quant_method)),)

# gene level quantification
$(foreach l,$(se),$(eval $(call make-htseq-quant-rule,$(l),$(quant_method),$(l).se,gene,$(gtf_file_abspath))))	
$(foreach l,$(pe),$(eval $(call make-htseq-quant-rule,$(l),$(quant_method),$(l).pe,gene,$(gtf_file_abspath))))

ifeq ($(exon_quant),y)
# exon level quantification
ifeq ($(patsubst htseq%,,$(exon_quant_method)),)
$(foreach l,$(se),$(eval $(call make-htseq-quant-rule,$(l),$(exon_quant_method),$(l).se,exon,$(gtf_file_wexonid))))	
$(foreach l,$(pe),$(eval $(call make-htseq-quant-rule,$(l),$(exon_quant_method),$(l).pe,exon,$(gtf_file_wexonid))))

# counts per exon
$(name)/$(mapper)/$(quant_method)/exons.raw.$(exon_quant_method).tsv: $(foreach p, $(pe),$(call lib2quant_folder,$(p))$(p).pe.exons.raw.$(exon_quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.exons.raw.$(exon_quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

endif

endif

# transcript level quantification
ifeq ($(transcript_quant),y)
$(foreach l,$(se),$(eval $(call make-htseq-quant-rule,$(l),$(quant_method),$(l).se,transcript,$(gtf_file_abspath))))	
$(foreach l,$(pe),$(eval $(call make-htseq-quant-rule,$(l),$(quant_method),$(l).pe,transcript,$(gtf_file_abspath))))

# counts per transcript
$(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).tsv: $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@


endif

endif


#*****************
# Simple count
#*****************

ifeq (basic,$(quant_method))
# TODO: from exon2gene
$(name)/$(mapper)/basic/genes.raw.basic.tsv: $(foreach p,$(pe),$(name)/$(mapper)/$(quant_method)/$(p).genes.raw.$(quant_method).tsv) $(foreach s,$(se), $(name)/$(mapper)/$(quant_method)/$(s).genes.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/basic/transcripts.raw.basic.tsv: 
	$(call p_info, Warning! Basic method does not produce quantification at transcript level (in TODO). Generating empty file $@.)
	@$(call empty_file,$@)

$(name)/$(mapper)/basic/exons.raw.basic.tsv: 
	$(call p_info, Warning! Basic method does not produce quantification at exon level (in TODO). Generating empty file $@.)
	@$(call empty_file,$@)

#$(name)/$(mapper)/basic/exons.fpkm.basic.tsv: $(name)/$(mapper)/basic/exons.raw.basic.tsv $(feat_length)
#	irap_raw2metric --tsv $<  --lengths $(feat_length) --feature exon --metric fpkm --out $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/basic/%.genes.raw.basic.tsv $(name)/$(mapper)/basic/%.exons.raw.basic.tsv: $(name)/$(mapper)/%.se.hits.byname.bam $(gtf_file_wexonid)
	$(call run_naive_count,$<,$(gtf_file_wexonid),$(@D)/$*.genes.raw.basic.tsv,$(@D)/$*.exons.raw.basic.tsv)

$(name)/$(mapper)/basic/%.genes.raw.basic.tsv $(name)/$(mapper)/basic/%.exons.raw.basic.tsv: $(name)/$(mapper)/%.pe.hits.byname.bam $(gtf_file_wexonid)
	$(call run_naive_count,$<,$(gtf_file_wexonid),$(@D)/$*.genes.raw.basic.tsv,$(@D)/$*.exons.raw.basic.tsv)

$(name)/$(mapper)/basic/%.transcripts.raw.basic.tsv: $(name)/$(mapper)/%.pe.hits.byname.bam
	$(call p_info, Warning! Basic method does not produce quantification at transcript level (in TODO). Generating empty file $@.)
	@$(call empty_file,$@)

$(name)/$(mapper)/basic/%.transcripts.raw.basic.tsv: $(name)/$(mapper)/%.se.hits.byname.bam
	$(call p_info, Warning! Basic method does not produce quantification at transcript level (in TODO). Generating empty file $@.)
	@$(call empty_file,$@)

endif

#*****************
# Flux-capacitor
#*****************

#
# Run flux
#
# $1 - lib
# $2 - SINGLE|PE
# $3 - bam file prefix (includes .se|.pe)
define make-fluxcap-quant-rule=
$(call lib2quant_folder,$(1))$(3).flux_cap.tsv: $(call lib2bam_folder,$(1))$(3).hits.byname.bam $(gtf_file_abspath) $(call lib2bam_folder,$(1))$(3).hits.byname.bam $(call lib2bam_folder,$(1))$(3).hits.bam.bai 
	mkdir -p $$(@D) && $$(call run_flux_cap,$$<,$$(gtf_file_abspath),$(2),$$@)

$(call lib2quant_folder,$(1))$(3).genes.fpkm.flux_cap.flux_cap.tsv: $(call lib2quant_folder,$(1))$(3).flux_cap.tsv
	tsv.aggrbygene $$< FPKM $$@.tmp && mv  $$@.tmp $$@
endef

ifeq (flux_cap,$(quant_method))
$(foreach l,$(se),$(eval $(call make-fluxcap-quant-rule,$(l),SINGLE,$(l).se)))	
$(foreach l,$(pe),$(eval $(call make-fluxcap-quant-rule,$(l),PAIRED,$(l).pe)))


# Process flux output
%.genes.raw.flux_cap.tsv: %.flux_cap.tsv 
	tsv.aggrbygene $< counts $@.tmp && mv $@.tmp $@



%.transcripts.raw.flux_cap.tsv: %.flux_cap.tsv
	cut -f 1,4 $< |tail -n +2 > $@.tmp && mv  $@.tmp $@

%.transcripts.fpkm.flux_cap.flux_cap.tsv: %.flux_cap.tsv
	cut -f 1,6 $< | tail -n +2  > $@.tmp && mv  $@.tmp $@

%.exons.raw.flux_cap.tsv: 
	$(call p_info, Warning! Flux_capacitor does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)

%.exons.fpkm.flux_cap.flux_cap.tsv: 
	$(call p_info, Warning! Flux_capacitor does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)



$(name)/$(mapper)/flux_cap/genes.raw.flux_cap.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).genes.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).genes.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/flux_cap/exons.raw.flux_cap.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).exons.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).exons.raw.$(quant_method).tsv)
	@$(call empty_file,$@)

$(name)/$(mapper)/flux_cap/transcripts.raw.flux_cap.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).transcripts.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).transcripts.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/genes.fpkm.flux_cap.flux_cap.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).genes.fpkm.$(quant_method).$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).genes.fpkm.$(quant_method).$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/transcripts.fpkm.flux%.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).transcripts.fpkm.$(quant_method).$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).transcripts.fpkm.$(quant_method).$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@

ifeq ($(exon_quant),y)
$(name)/$(mapper)/$(quant_method)/exons.fpkm.flux_cap.flux_cap.tsv: 
	$(call p_info, Warning! Flux-capacitor does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)
endif

endif
#*****************
# NURD
#*****************

ifeq (nurd,$(quant_method))
$(name)/$(mapper)/nurd/genes.raw.nurd.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@


$(name)/$(mapper)/nurd/transcripts.raw.nurd.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

#$(name)/$(mapper)/nurd/transcripts.raw.nurd.tsv: 
#	$(call p_info, Warning! NURD does not produce raw quantification at transcript level. Generating empty file $@.)
#	@$(call empty_file,$@)

#
# Run flux
#
# $1 - lib
# $2 - bam file prefix (includes .se|.pe)
define make-nurd-quant-rule=
$(call lib2quant_folder,$(1))$(2).nurd.tsv $(call lib2quant_folder,$(1))$(2).raw.nurd.tsv : $(call lib2bam_folder,$(1))$(2).hits.bam $(gtf_file_abspath) 
	$$(call run_nurd,$$<,$$(gtf_file_abspath),$(call lib2quant_folder,$(1))$(2))
# Process nurd output
$(call lib2quant_folder,$(1))$(2).genes.raw.nurd.tsv: $(call lib2quant_folder,$(1))$(2).nurd.tsv 
	cut -f 1,3 $$< > $$@.tmp && mv $$@.tmp $$@

$(call lib2quant_folder,$(1))$(2).genes.fpkm.nurd.nurd.tsv: $(call lib2quant_folder,$(1))$(2).nurd.tsv
	cut -f 1,6 $$< > $$@.tmp && mv $$@.tmp $$@


# transcripts/isoforms
$(call lib2quant_folder,$(1))$(2).transcripts.raw.nurd.tsv: $(call lib2quant_folder,$(1))$(2).raw.nurd.tsv
	irap_nurd2tsv --tsv $$< --out $$@.tmp && mv $$@.tmp $$@

$(call lib2quant_folder,$(1))$(2).transcripts.fpkm.nurd.nurd.tsv: $(call lib2quant_folder,$(1))$(2).nurd.tsv
	irap_nurd2tsv --tsv $$< --out $$@.tmp && mv $$@.tmp $$@

endef


$(foreach l,$(se),$(eval $(call make-nurd-quant-rule,$(l),$(l).se)))	
$(foreach l,$(pe),$(eval $(call make-nurd-quant-rule,$(l),$(l).pe)))


# exons
%.exons.raw.nurd.tsv: 
	$(call p_info, Warning! NURD does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)

%.exons.fpkm.nurd.nurd.tsv: 
	$(call p_info, Warning! NURD does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)


# nurd
$(name)/$(mapper)/nurd/genes.fpkm.nurd.nurd.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.fpkm.$(quant_method).nurd.tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.fpkm.$(quant_method).nurd.tsv)
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) ) > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/nurd/transcripts.fpkm.nurd.nurd.tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.fpkm.$(quant_method).nurd.tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.transcripts.fpkm.$(quant_method).nurd.tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) ) > $@.tmp && mv $@.tmp $@

ifeq ($(exon_quant),y)
$(name)/$(mapper)/nurd/exons.fpkm.$(exon_quant_method).nurd.tsv: 
	$(call p_info, Warning! NURD does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)

$(name)/$(mapper)/nurd/exons.raw.$(exon_quant_method).tsv: 
	$(call p_info, Warning! NURD does not produce quantification at exon level. Generating empty file $@.)
	@$(call empty_file,$@)
endif

$(name)/$(mapper)/$(quant_method)/%.nlib.nurd.tsv:
	$(call p_info, Warning! Unable to generate  nlib file $@ with $(quant_method).)
	@$(call empty_file,$@)
$(name)/$(mapper)/$(quant_method)/%.nlib.nurd.tsv:
	$(call p_info, Warning! Unable to generate nlib file $@ with $(quant_method).)
	@$(call empty_file,$@)

endif

#################################################
# Stringtie

# generate the rules for stringtie
ifeq ($(patsubst stringtie%,,$(quant_method)),)
# 
define stringtie2counts_g=
	irap_stringtie2counts -c $(3) --label "$(1)" -r $($(subst .pe,,$(subst .se,,$(1)))_rs) -p $(if $(call is_pe_lib,$($(subst .pe,,$(subst .se,,$(1))))),y,n) -o $(5).tmp -i /dev/null && mv $(5).tmp $(5)
endef

define stringtie2counts_t=
	irap_stringtie2counts -c $(3) --label "$(1)" -r $($(subst .pe,,$(subst .se,,$(1)))_rs) -p $(if $(call is_pe_lib,$(subst .pe,,$(subst .se,,$(1)))),y,n) -o $(5).tmp -i $(2)/$(1).transcripts.raw.$(4).tsv && mv $(5).tmp $(5)
endef


define make-stringtie-quant-rule=

$(call lib2quant_folder,$(1))$(3).$(quant_method).gtf: $(call lib2bam_folder,$(1))$(3).hits.bam
	$(call run_stringtie,$(call lib2bam_folder,$(1))$(3).hits.bam,$(1),$(if $(filter $(quant_method),stringtie_nd),-e,-l CUFF),$$@)

$(call lib2quant_folder,$(1))$(3).transcripts.fpkm.$(quant_method).$(quant_norm_tool).tsv: $(call lib2quant_folder,$(1))$(3).$(quant_method).gtf $(call lib2quant_folder,$(1))$(1)/t_data.ctab 
	cut -f 6,12 $$(@D)/$(1)/t_data.ctab | tail -n +2 > $$@.tmp && mv $$@.tmp $$@

$(call lib2quant_folder,$(1))$(3).genes.fpkm.$(quant_method).$(quant_norm_tool).tsv: $(call lib2quant_folder,$(1))$(3).$(quant_method).gtf $(call lib2quant_folder,$(1))$(1)/t_data.ctab 
	tsv.aggrbygene.stringtie $(call lib2quant_folder,$(1))$(1)/t_data.ctab  FPKM $$@.tmp && mv  $$@.tmp $$@

$(call lib2quant_folder,$(1))$(3).transcripts.raw.$(quant_method).tsv:  $(call lib2quant_folder,$(1))$(3).$(quant_method).gtf
	$$(call p_info,Warning! Stringtie produces FPKMS and does not provide counts. Generating pseudo counts $$@.)
	$$(call stringtie2counts_t,$(1),$$(@D),$(call lib2quant_folder,$(1))$(1)/t_data.ctab,$(quant_method),$$@)

$(call lib2quant_folder,$(1))$(3).genes.raw.$(quant_method).tsv: $(call lib2quant_folder,$(1))$(3).$(quant_method).gtf
	$$(call p_info,Warning! Stringtie produces FPKMS and does not provide counts. Generating pseudo counts $$@.)
	$$(call stringtie2counts_g,$(1),$$(@D),$(call lib2quant_folder,$(1))$(1)/t_data.ctab,$(quant_method),$$@)

endef

# dynamic rules
$(foreach l,$(se),$(eval $(call make-stringtie-quant-rule,$(l),$(quant_method),$(l).se,gene)))	
$(foreach l,$(pe),$(eval $(call make-stringtie-quant-rule,$(l),$(quant_method),$(l).pe,gene)))

# stringtie* specific rules
$(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).tsv)
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/%xons.raw.$(quant_method).tsv: 
	$(call p_info,Warning! Stringtie produces FPKMS and does not provide counts. Generating empty file $@.)
	@$(call empty_file,$@)

# Run Cuffmerge on all your assemblies to create a single merged transcriptome annotation
$(name)/$(mapper)/$(quant_method)/$(name).merged.gtf:  $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.cuff.gtf) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.cuff.gtf) $(name)/$(mapper)/$(quant_method)/$(name).assemblies.txt
	$(call run_cuffmerge,$(quant_method),$@)

# 
$(name)/$(mapper)/$(quant_method)/$(name).assemblies.txt: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.cuff.gtf) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.cuff.gtf)
	rm -f $@
	touch $@.tmp && for f in $^; do echo $$f >> $@.tmp ; done && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/transcripts.fpkm.$(quant_method).$(quant_norm_tool).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).transcripts.fpkm.$(quant_method).$(quant_norm_tool).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).transcripts.fpkm.$(quant_method).$(quant_norm_tool).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/genes.fpkm.$(quant_method).$(quant_norm_tool).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).genes.fpkm.$(quant_method).$(quant_norm_tool).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).genes.fpkm.$(quant_method).$(quant_norm_tool).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/%nlib.stringtie.stringtie.tsv:
	$(call p_info, Warning! Unable to generate  nlib file $@ with $(quant_method).)
	@$(call empty_file,$@)

$(name)/$(mapper)/$(quant_method)/exons.fpkm.$(quant_method).$(quant_norm_tool).tsv: 
	$(call p_info, Warning! Exon quantification file $@ with $(quant_method) will be empty.)
	@$(call empty_file,$@)

$(name)/$(mapper)/$(quant_method)/exons.raw.$(quant_method).tsv: 
	$(call p_info, Warning (WIP)! Exon quantification file $@ with $(quant_method) will be empty.)
	@$(call empty_file,$@)
endif

##############################################################
# Exon quantification - DexSeq
ifndef dexseq_params
dexseq_params=
endif

# 'no', the exons that can not be assigned to a single gene are ignored otherwise the genes are merged
ifndef dexseq_index_params
dexseq_prepare_annotation_params=
endif
#dexseq_prepare_annotation_params=-r yes
################################
# if PE then add option -p yes
# sam/bam (-f bam) file needs to be sorted by read name or chr (-r name)
# strand specific (== htseq)
# 1 - bam
# 2 - output
# 3 - lib
# 4 - gtf
define run_dexseq=
	samtools view -F 4 $(1) | python $(IRAP_DIR)/Rlibs/DEXSeq/python_scripts/dexseq_count.py  $(dexseq_params) $(if $(call is_pe_lib,$(3)),-p yes) $(call htseq_strand_params,$(3)) $(4) - $(2).tmp && \
	grep -E "(_ambiguous|_empty|_lowaqual|_notaligned)" $(2).tmp > $(2).stats && \
	grep -v -E "(_ambiguous|_empty|_lowaqual|_notaligned)" $(2).tmp > $(2).tmp2 && \
	mv $(2).tmp2 $(2) && rm -f $(2).tmp
endef


%.gtf.DEXSeq.gff: %.gtf
	python $(IRAP_DIR)/Rlibs/DEXSeq/python_scripts/dexseq_prepare_annotation.py $(dexseq_prepare_annotation_params)  $< $@.tmp && mv $@.tmp $@

ifeq ($(strip $(exon_quant)),y)
ifeq ($(strip $(exon_quant_method)),dexseq) 

# add the generation of the flatten annotation to stage0 iff dexseq is selected
exon_length=$(name)/data/$(gtf_file_basename).dexseq.lengths.Rdata
SETUP_DATA_FILES+=$(exon_length)

## htseq bam file needs to be sorted by name
# $1 - lib
# $2 - bam file prefix (includes .se|.pe)
# $3 - gtf file
define make-dexseq-quant-rule=
$(call lib2quant_folder,$(1))$(2).exons.raw.dexseq.tsv: $(call lib2bam_folder,$(1))$(2).hits.byname.bam $(3)
	mkdir -p $$(@D) && $$(call run_dexseq,$$<,$$@,$(1),$(3))
endef

$(name)/$(mapper)/$(quant_method)/exons.raw.dexseq.tsv: $(foreach p, $(pe),$(call lib2quant_folder,$(p))$(p).pe.exons.raw.$(exon_quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.exons.raw.$(exon_quant_method).tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@


$(foreach l,$(se),$(eval $(call make-dexseq-quant-rule,$(l),$(l).se, $(gtf_file_abspath).DEXSeq.gff)))	
$(foreach l,$(pe),$(eval $(call make-dexseq-quant-rule,$(l),$(l).pe, $(gtf_file_abspath).DEXSeq.gff)))

endif
endif
#######################################################
# RSEM

ifeq ($(quant_method),rsem)

# For now RSEM only works with the STAR alignments
ifneq ($(mapper),star)
$(error RSEM needs to be ran with STAR)
endif

star_map_params+=
# --quantMode TranscriptomeSAM will be enabled in the star code because transcript quantification is enabled

# extra parameters passed when preparing the reference
rsem_prepare_ref_params?=
# extra parameters passed to rsem-calculate-expression 
rsem_params?=
rsem_params+=--alignments  --estimate-rspd  --calc-ci --no-bam-output --seed 12321 -p $(max_threads) --ci-memory $(max_mem)

##--bam
rsem_reference_name=$(subst .gtf,rsem,$(gtf_file_abspath))/rsemref
rsem_reference_target=$(rsem_reference_name)/rsem.irap
# --strandedness non|forward|reverse
# --paired-end
# Add the reference preparation to STAGE0
SETUP_DATA_FILES+=$(rsem_reference_target)

$(rsem_reference_target): $(reference_abspath)  $(gtf_file_abspath)
	mkdir -p $(@D) &&	irap_wrapper.sh rsem rsem-prepare-reference  --gtf $(gtf_file_abspath) $(reference_abspath)  $(rsem_reference_name) $(rsem_prepare_ref_params) && touch $@

# 1 bam
# 2 output prefix
# 3 PE/SE info
#  --sam/--bam [--paired-end] input reference_name sample_name
# TODO: stranded data: --forward-prob 0
define run_rsem=
	irap_wrapper.sh rsem rsem-calculate-expression $(rsem_params)   $(3) $(1) $(rsem_reference_name) $(2)
endef
## 
#sample_name.genes.results
#sample.name.isoforms.results
# 1 lib
# 2 bam prefix
# 3 paired-end/se options
define make-rsem-quant-rule=
$(call lib2quant_folder,$(1))$(2).genes.results $(call lib2quant_folder,$(1))$(2).isoforms.results: $(call lib2bam_folder,$(1))$(2).hits.bam.trans.bam 
	(mkdir -p $$(@D) && $$(call run_rsem,$$<,$$(@D)/$(2),$(3),$(1))) || (rm -rf $$@ && exit 1)

$(call lib2quant_folder,$(1))$(2).genes.tpm.rsem.rsem.tsv: $(call lib2quant_folder,$(1))$(2).genes.results
	cut -f 1,6 $$< |tail -n +2 > $$@.tmp && mv $$@.tmp $$@

$(call lib2quant_folder,$(1))$(2).transcripts.raw.rsem.tsv: $(call lib2quant_folder,$(1))$(2).isoforms.results $(call lib2quant_folder,$(1))$(2).genes.results
	cut -f 1,5 $$< |tail -n +2 > $$@.tmp && mv $$@.tmp $$@

endef

# quantification
$(foreach l,$(se),$(eval $(call make-rsem-quant-rule,$(l),$(l).se,)))
$(foreach l,$(pe),$(eval $(call make-rsem-quant-rule,$(l),$(l).pe,--paired-end)))


%.genes.raw.rsem.tsv: %.genes.results
	cut -f 1,5 $< |tail -n +2 > $@.tmp && mv $@.tmp $@

%.genes.fpkm.rsem.rsem.tsv: %.genes.results
	cut -f 1,7 $< |tail -n +2 > $@.tmp && mv $@.tmp $@


%.transcripts.fpkm.rsem.rsem.tsv: %.isoforms.results %.genes.results
	cut -f 1,7 $< |tail -n +2 > $@.tmp && mv $@.tmp $@

%.transcripts.tpm.rsem.rsem.tsv: %.isoforms.results %.genes.results
	cut -f 1,6 $< |tail -n +2 > $@.tmp && mv $@.tmp $@


$(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/transcripts.fpkm.$(quant_method).$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.fpkm.$(quant_method).$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.transcripts.fpkm.$(quant_method).$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@


$(name)/$(mapper)/$(quant_method)/genes.fpkm.$(quant_method).$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).genes.fpkm.$(quant_method).$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).genes.fpkm.$(quant_method).$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@

endif



########################################################
# kallisto
ifeq ($(mapper),kallisto)
SETUP_DATA_FILES+=$(mapTrans2gene)
endif

ifeq ($(quant_method),kallisto)

# 
ifneq ($(mapper),none)
$(warning kallisto does not need the reads to be aligned)
endif

# transcripts
$(call file_exists,$(user_trans_abspath))


# Note: currently if kallisto is ran in quant mode no bam file is kept
kallisto_index_params?=
# --fusion --bias
kallisto_quant_params?=


#
kallisto_index_name=$(trans_abspath)_kallisto/kallisto_index
kallisto_index=$(trans_abspath)_kallisto/kallisto_index.irap

# Add the reference preparation to STAGE0
SETUP_DATA_FILES+=$(kallisto_index)

# code now in irap_map.mk
$(kallisto_index): $(trans_abspath)
	$(call run_kallisto_index)




# out_directory - 1
# fastq files - 2
# single-end? - 4
# read lenth - 3
# lib - 5
define run_kallisto_quant=
 irap_wrapper.sh kallisto kallisto quant $(kallisto_quant_params) -i $(kallisto_index_name) $(if $(4),,--single) $(call kallisto_strand_params,$(5)) -l $(3) -s 1  -t $(max_threads) -o $(1)  $(2)
endef
# 

# abundance.tsv
# 1 lib
# 2 file prefix (lib.[pe,se])
define make-kallisto-quant-rule=


$(call lib2quant_folder,$(1))$(1)/$(1).abundance.tsv: $(call libname2fastq,$(1)) $(kallisto_index)
	(mkdir -p $$(@D) && $$(call run_kallisto_quant,$$(@D),$(call libname2fastq,$(1)),$($(1)_rs),$(call is_pe_lib,$(1)),$(1)) && mv $$(@D)/abundance.tsv $$@)


$(call lib2quant_folder,$(1))$(2).transcripts.raw.kallisto.tsv: $(call lib2quant_folder,$(1))$(1)/$(1).abundance.tsv
	cut -f 1,4 $$< |tail -n +2 > $$@.tmp && mv $$@.tmp $$@

$(call lib2quant_folder,$(1))$(2).transcripts.tpm.kallisto.kallisto.tsv: $(call lib2quant_folder,$(1))$(1)/$(1).abundance.tsv
	cut -f 1,5 $$< |tail -n +2 > $$@.tmp && mv $$@.tmp $$@


$(call lib2quant_folder,$(1))$(2).genes.raw.kallisto.tsv: $(call lib2quant_folder,$(1))$(2).transcripts.raw.kallisto.tsv $(mapTrans2gene)
	libTSVAggrTransByGene -i $$<  -m $(mapTrans2gene) -o $$@.tmp && mv  $$@.tmp $$@


endef


# quantification
$(foreach l,$(se),$(eval $(call make-kallisto-quant-rule,$(l),$(l).se)))
$(foreach l,$(pe),$(eval $(call make-kallisto-quant-rule,$(l),$(l).pe)))



$(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@


endif


########################################################
# Salmon - Quasi mapping mode

ifeq ($(quant_method),salmon)

# 
ifneq ($(mapper),none)
$(warning salmon does not need the reads to be aligned)
endif

# transcripts
$(call file_exists,$(user_trans_abspath))

# K - minimum matching
# the user should adjust the K to 1/3 of the read length
# use 21 as default (assuming minum read length of 30)
ifndef salmon_index_params
salmon_index_params=--type quasi -k 21
endif

ifndef salmon_quant_params
salmon_quant_params=
endif


salmon_index_name=$(trans_abspath)_salmon/salmon_index
salmon_index=$(trans_abspath)_salmon/salmon_index.irap

# Add the reference preparation to STAGE0
SETUP_DATA_FILES+=$(salmon_index)


$(salmon_index): $(trans_abspath)
	$(call run_salmon_index)
##	mkdir -p $(@D) && irap_w
##	mkdir -p $(@D) && irap_wrapper.sh salmon salmon index $(salmon_index_params) -i $(salmon_index_name) -t $< && touch $@



# LibType
#	Paired-end 	Single-end
#-fr-unstranded 	-l IU 	-l U
#-fr-firststrand 	-l ISR 	-l SR
#-fr-secondstrand 	-l ISF 	-l SF
irap_salmon_lib_type1=$(if $(findstring $(1),first),SR,$(if $(findstring $(1),second),SF,U))
irap_salmon_lib_type=$(if $(call is_pe_lib,$(1)),I,)$(call irap_salmon_lib_type1,$(1))

# out_directory - 1
# read lenth - 2 - ignored
# lib - 3
# fastq files - 4
define run_salmon_quant_se=
 irap_wrapper.sh salmon salmon quant -i $(salmon_index_name) $(salmon_quant_params)  --libType $(call irap_salmon_lib_type,$(3)) -r $(4) -o $(1) -p $(max_threads)
endef
define run_salmon_quant_pe=
 irap_wrapper.sh salmon salmon quant -i $(salmon_index_name) $(salmon_quant_params)  --libType  $(call irap_salmon_lib_type,$(3)) -1 $(word 1,$(4)) -2 $(word 2,$(4))   -o $(1) -p $(max_threads)
endef
# 

# abundance.tsv
# 1 lib
# 2 file prefix (lib.[pe,se])
# 3 pe|se
define make-salmon-quant-rule=

$(call lib2quant_folder,$(1))$(1)/$(1).abundance.tsv: $(call libname2fastq,$(1)) $(salmon_index)
	(mkdir -p $$(@D) && $$(call run_salmon_quant_$(3),$$(@D),$($(1)_rs),$(1),$(call libname2fastq,$(1))) && mv $$(@D)/quant.sf $$@)


$(call lib2quant_folder,$(1))$(2).transcripts.raw.salmon.tsv: $(call lib2quant_folder,$(1))$(1)/$(1).abundance.tsv
	cut -f 1,5 $$< |tail -n +2 > $$@.tmp && mv $$@.tmp $$@

$(call lib2quant_folder,$(1))$(2).transcripts.tpm.salmon.salmon.tsv: $(call lib2quant_folder,$(1))$(1)/$(1).abundance.tsv
	cut -f 1,4 $$< |tail -n +2 > $$@.tmp && mv $$@.tmp $$@


$(call lib2quant_folder,$(1))$(2).genes.raw.salmon.tsv: $(call lib2quant_folder,$(1))$(2).transcripts.raw.salmon.tsv $(mapTrans2gene)
	libTSVAggrTransByGene -i $$< -m $(mapTrans2gene) -o $$@.tmp && mv  $$@.tmp $$@

endef

# quantification
$(foreach l,$(se),$(eval $(call make-salmon-quant-rule,$(l),$(l).se,se)))
$(foreach l,$(pe),$(eval $(call make-salmon-quant-rule,$(l),$(l).pe,pe)))


$(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/genes.raw.$(quant_method).tsv: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).tsv) $(foreach s,$(se),$(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).tsv) 
	( $(call pass_args_stdin,$(call merge_tsv,$(quant_method)),$@,$^) )  > $@.tmp && mv $@.tmp $@


endif


#####################################################################
# relative isoform usage (RIU) and dominant transcript
ifeq ($(transcript_quant),y)

SETUP_DATA_FILES+=$(mapTrans2gene)


$(name)/$(mapper)/$(quant_method)/%.transcripts.riu.$(quant_method).irap.$(expr_ext): $(name)/$(mapper)/$(quant_method)/%.transcripts.raw.$(quant_method).$(expr_ext) $(mapTrans2gene)
	irap_transcript_gene_rel_expr --ifile $< --$(expr_format) --map $(mapTrans2gene) --cores $(max_threads)  --gene_col 1 --trans_col 2  --out $@ || (rm -f $@ && exit 1)

$(name)/$(mapper)/$(quant_method)/transcripts.riu.$(quant_method).irap.tsv: $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.riu.$(quant_method).irap.tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.riu.$(quant_method).irap.tsv)
	( $(call pass_args_stdin,irap_merge_tsv.sh,$@,$^) ) > $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/transcripts.riu.$(quant_method).irap.mtx.gz: $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.riu.$(quant_method).irap.mtx.gz) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.riu.$(quant_method).irap.mtx.gz)
	( $(call pass_args_stdin,irap_merge_mtx,$@,-o $@ --in "$^") ) || ( rm -f $@ || exit 1 )


$(name)/$(mapper)/$(quant_method)/%.transcripts.dt.$(quant_method).irap.$(expr_ext): $(name)/$(mapper)/$(quant_method)/%.transcripts.riu.$(quant_method).irap.$(expr_ext)  $(mapTrans2gene)
	irap_riu2dominant --fc $(dt_fc) -i $< --$(expr_format) --map $(mapTrans2gene) --cores $(max_threads) --gene_col 1 --trans_col 2  --out $@ || (rm -f $@ && exit 1)

# dominant transcript
$(name)/$(mapper)/$(quant_method)/transcripts.dt.$(quant_method).irap.$(expr_ext):$(name)/$(mapper)/$(quant_method)/transcripts.riu.$(quant_method).irap.$(expr_ext) $(mapTrans2gene)
	irap_riu2dominant --fc $(dt_fc) -i $<  --$(expr_format) --map $(mapTrans2gene) --cores $(max_threads) --gene_col 1 --trans_col 2  --out $@ || (rm -f $@ && exit 1)

# include the raw counts
ifeq ($(transcript_expr),n)

define transcripts_quant_file=
endef
define transcripts_quant_files=
endef
else

STAGE3_S_OFILES+= $(name)/$(mapper)/$(quant_method)/transcripts.raw.$(quant_method).$(expr_ext)

# do not get the dominant transcript
ifeq ($(dt_fc),n)

trans_file_target=riu
else
trans_file_target=dt
# include riu files are part of the set of files produced in stage3 
STAGE3_S_OFILES+= $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.riu.$(quant_method).irap.tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.riu.$(quant_method).irap.tsv)
endif


STAGE3_S_TARGETS+= $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.$(trans_file_target).$(quant_method).irap.$(expr_ext)) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.$(trans_file_target).$(quant_method).irap.$(expr_ext)) 


# riu files are always produced
STAGE3_S_OFILES+= $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.riu.$(quant_method).irap.$(expr_ext)) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.riu.$(quant_method).irap.$(expr_ext)) 


# include the raw counts
STAGE3_S_OFILES+= $(name)/$(mapper)/$(quant_method)/transcripts.riu.$(quant_method).irap.$(expr_ext)


# useful functions
define transcripts_quant_file=
$(if $(filter y,$(transcript_quant)),$(name)/$(mapper)/$(quant_method)/transcripts.$(trans_file_target).$(quant_method).irap.$(expr_ext) $(name)/$(mapper)/$(quant_method)/transcripts.$(trans_file_target).$(quant_method).irap.$(expr_ext),)
endef

define transcripts_quant_files=
$(if $(filter y,$(transcript_quant)),$(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.$(trans_file_target).$(quant_method).irap.$(expr_ext)) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.$(trans_file_target).$(quant_method).irap.$(expr_ext)),)
endef


#STAGE3_S_TARGETS+= $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).tsv) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).tsv)

# Atlas specific - get fpkms and tpms
ifdef atlas_run
# 
STAGE3_S_TARGETS+= $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.fpkm.$(quant_method).irap.$(expr_ext) $(call lib2quant_folder,$(p))$(p).pe.transcripts.tpm.$(quant_method).irap.$(expr_ext)) $(foreach s,$(se),  $(call lib2quant_folder,$(s))$(s).se.transcripts.fpkm.$(quant_method).irap.$(expr_ext) $(call lib2quant_folder,$(s))$(s).se.transcripts.tpm.$(quant_method).irap.$(expr_ext))
# include fpkm/TPMs 
STAGE3_S_OFILES+= $(foreach p,$(pe), $(call lib2quant_folder,$(p))$(p).pe.transcripts.fpkm.$(quant_method).irap.$(expr_ext) $(call lib2quant_folder,$(p))$(p).pe.transcripts.tpm.$(quant_method).irap.$(expr_ext)) $(foreach s,$(se),  $(call lib2quant_folder,$(s))$(s).se.transcripts.fpkm.$(quant_method).irap.$(expr_ext) $(call lib2quant_folder,$(s))$(s).se.transcripts.tpm.$(quant_method).irap.$(expr_ext))
endif

endif
else
# no transcript quantification
define transcripts_quant_file=
endef
define transcripts_quant_files=
endef
endif


##################################
## collect QC stats
$(name)/$(mapper)/$(quant_method)/%.genes.raw.$(quant_method).quant_qc.tsv: $(name)/$(mapper)/$(quant_method)/%.genes.raw.$(quant_method).$(expr_ext)
	irap_quant_qc --ifile $<  --feature gene --$(expr_format) --metric raw --gtf $(gtf_file_abspath) --out $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/%.exons.raw.$(exon_quant_method).quant_qc.tsv: $(name)/$(mapper)/$(quant_method)/%.exons.raw.$(exon_quant_method).$(expr_ext)
	irap_quant_qc --ifile $< --feature exon --metric raw --$(expr_format) --gtf $(gtf_file_abspath) --out $@.tmp && mv $@.tmp $@

$(name)/$(mapper)/$(quant_method)/%.transcripts.raw.$(quant_method).quant_qc.tsv: $(name)/$(mapper)/$(quant_method)/%.transcripts.raw.$(quant_method).$(expr_ext)
	irap_quant_qc --ifle $< --feature transcript --$(expr_format) --metric raw --gtf $(gtf_file_abspath) --out $@.tmp && mv $@.tmp $@

#####################################
##
ifeq ($(isl_mode),y)
STAGE3_S_TARGETS+= $(a_quant_qc_stats)
STAGE3_S_OFILES+= $(a_quant_qc_stats)
else
STAGE4_TARGETS+=$(a_quant_qc_stats)
endif


