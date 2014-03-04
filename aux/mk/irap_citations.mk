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
#    $Id: scripts/irap Nuno Fonseca Fri Jan 4 18:20:42 2013$
# =========================================================
# Citations to the software used and version used.
# Provide users an easy way to acknowledge the software they use (besides IRAP ;-))


# load the version´s file generated during installation
$(call file_exists,$(IRAP_DIR)/aux/mk/irap_versions.mk)
include $(IRAP_DIR)/aux/mk/irap_versions.mk

############################################################
# Map IRAP options to programs
qc_on=FASTQC FASTX bowtie1
qc_off=empty

quant_norm_none=empty
de_none=empty
quant_none=empty
gse_none=empty
# tophat1 tophat2 smalt gsnap soapsplice bwa1 bwa2 bowtie1 bowtie2 gem star
mapper_tophat1=tophat1
mapper_tophat2=tophat2
mapper_smalt=SMALT
mapper_gsnap=GSNAP
mapper_soapsplice=SOAPsplice
mapper_bwa1=bwa1
mapper_bwa2=bwa2
mapper_bowtie1=bowtie1
mapper_bowtie2=bowtie2
mapper_gem=GEM
mapper_star=STAR
mapper_osa=OSA

# Quant. methods
quant_basic=empty
quant_htseq1=htseq
quant_htseq2=htseq
quant_cufflinks1=cufflinks1
quant_cufflinks1_nd=cufflinks1
quant_cufflinks2=cufflinks2
quant_cufflinks2_nd=cufflinks2
quant_flux_cap=FLUX_CAPACITOR
quant_scripture=scripture
quant_bitseq=BitSeq
quant_ireckon=Ireckon

# quant_norm_method cufflinks1 cufflinks2 cufflinks1_nd cufflinks2_nd deseq flux_cap edger 
quant_norm_cufflinks1=cufflinks1
quant_norm_cufflinks1_nd=cufflinks1
quant_norm_cufflinks2=cufflinks2
quant_norm_cufflinks2_nd=cufflinks2
quant_norm_flux_cap=FLUX_CAPACITOR
quant_norm_deseq=deseq

# DE
de_cuffdiff1=cuffdiff1
de_cuffdiff2=cuffdiff2
de_cuffdiff1_nd=cuffdiff1
de_cuffdiff2_nd=cuffdiff2
de_deseq=DESeq
de_dexseq=DEXSeq
de_bayseq=baySeq

# GSE
gse_piano=piano
# TODO: analysis/reports

################################################################
# key_citation
# key_bibtex
# key_version (load from another file created during installation)
empty_citation=
empty_bibtex=

# first and foremost, IRAP :)
IRAP_citation="IRAP 2012-2013"
IRAP_bibtex=irap.bib
IRAP_version=$(version)

#TODO: update
STAR_citation="Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, and Thomas R. Gingeras. STAR: ultrafast universal RNA-seq aligner Bioinformatics first published online October 25, 2012"
STAR_bibtex=STAR.bib

#### 
tophat1_citation=Cole Trapnell, Lior Pachter, and Steven L. Salzberg. TopHat: discovering splice junctions with RNA-Seq Bioinformatics (2009) 25(9): 1105-1111
tophat1_bibtex=tophat1.bib


#### 
tophat2_citation=Kim, Daehwan, et al. "TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions." Genome Biol 14.4 (2013): R36.
tophat2_bibtex=tophat2.bib

#### 
SMALT_citation=http://www.sanger.ac.uk/resources/software/smalt/
SMALT_bibtex= 

#### 
GSNAP_citation=Thomas D. Wu and Serban Nacu.Fast and SNP-tolerant detection of complex variants and splicing in short reads Bioinformatics (2010) 26(7): 873-881 
GSNAP_bibtex=gsnap.bib

#### 
SOAPsplice_citation=Huang, Songbo, Jinbo Zhang, Ruiqiang Li, Wenqian Zhang, Zengquan He, Tak-Wah Lam, Zhiyu Peng, and Siu-Ming Yiu. "SOAPsplice: genome-wide ab initio detection of splice junctions from RNA-Seq data." Frontiers in genetics 2 (2011).
SOAPsplice_bibtex=soapsplice.bib

#### 
bwa1_citation=Li, Heng, and Richard Durbin. "Fast and accurate short read alignment with Burrows–Wheeler transform." Bioinformatics 25, no. 14 (2009): 1754-1760.
bwa1_bibtex=bwa1.bib

#### 
bwa2_citation=Li, Heng, and Richard Durbin. "Fast and accurate long-read alignment with Burrows–Wheeler transform." Bioinformatics 26, no. 5 (2010): 589-595.
bwa2_bibtex=bwa2.bib

#### 
bowtie1_citation=Langmead, Ben, Cole Trapnell, Mihai Pop, and Steven L. Salzberg. "Ultrafast and memory-efficient alignment of short DNA sequences to the human genome." Genome Biol 10, no. 3 (2009): R25.
bowtie1_bibtex=bowtie1.bib

#### 
bowtie2_citation=Langmead, Ben, and Steven L. Salzberg. "Fast gapped-read alignment with Bowtie 2." Nature methods 9, no. 4 (2012): 357-359.
bowtie2_bibtex=bowtie2.bib

#### 
osa_citation=Hu, Jun, Huanying Ge, Matt Newman, and Kejun Liu. "OSA: a fast and accurate alignment tool for RNA-Seq." Bioinformatics 28, no. 14 (2012): 1933-1934.
osa_bibtex=osa.bib


#### 
GEM_citation=Marco-Sola, Santiago, Michael Sammeth, Roderic Guigó, and Paolo Ribeca. "The GEM mapper: fast, accurate and versatile alignment by filtration." Nature Methods 9, no. 12 (2012): 1185-1188.
GEM_bibtex=gem.bib

#### 
FASTQC_citation=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
FASTQC_bibtex=

#### 
FASTX_citation=http://hannonlab.cshl.edu/fastx_toolkit/
FASTX_bibtex=

#### 
cufflinks1_citation=Trapnell, Cole, Brian A. Williams, Geo Pertea, Ali Mortazavi, Gordon Kwan, Marijke J. Van Baren, Steven L. Salzberg, Barbara J. Wold, and Lior Pachter. "Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation." Nature biotechnology 28, no. 5 (2010): 511-515.
cufflinks1_bibtex=cufflinks1.bib

####
cufflinks2_citation=Trapnell, Cole, et al. "Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks." Nature protocols 7.3 (2012): 562-578.
cufflinks2_bibtex=cufflinks2.bib

####
htseq_citation=Simon Anders, Paul Theodor Pyl, Wolfgang Huber: HTSeq – A Python framework to work with high-throughput sequencing data. bioRxiv, 20 Feb 2014, http://doi.org/10.1101/002824.
htseq_bibtex=htseq.bib

####
FLUX_CAPACITOR_citation=Montgomery, Stephen B., Micha Sammeth, Maria Gutierrez-Arcelus, Radoslaw P. Lach, Catherine Ingle, James Nisbett, Roderic Guigo, and Emmanouil T. Dermitzakis. "Transcriptome genetics using second generation sequencing in a Caucasian population." Nature 464, no. 7289 (2010): 773-777.
FLUX_CAPACITOR_bibtex=flux_capacitor.bib


####
BitSeq_citation=Glaus, Peter, Antti Honkela, and Magnus Rattray. "Identifying differentially expressed transcripts from RNA-seq data with biological variation." Bioinformatics 28, no. 13 (2012): 1721-1728.
BitSeq_bibtex=bitseq.bib

####
Ireckon_citation=Mezlini, Aziz M., Eric JM Smith, Marc Fiume, Orion Buske, Gleb Savich, Sohrab Shah, Sam Aparicion, Derek Chiang, Anna Goldenberg, and Michael Brudno. "iReckon: Simultaneous isoform discovery and abundance estimation from RNA-seq data." Genome Research (2012).
Ireckon_bibtex=ireckon.bib


####
cuffdiff1_citation=Trapnell, Cole, Brian A. Williams, Geo Pertea, Ali Mortazavi, Gordon Kwan, Marijke J. Van Baren, Steven L. Salzberg, Barbara J. Wold, and Lior Pachter. "Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation." Nature biotechnology 28, no. 5 (2010): 511-515.
cuffdiff1_bibtex=cuffdiff1.bib

####
cuffdiff2_citation=Trapnell, Cole, David G. Hendrickson, Martin Sauvageau, Loyal Goff, John L. Rinn, and Lior Pachter. "Differential analysis of gene regulation at transcript resolution with RNA-seq." Nature Biotechnology (2012).
cuffdiff2_bibtex=cuffdiff2.bib

####
DESeq_citation=Anders, Simon, and Wolfgang Huber. "Differential expression analysis for sequence count data." Genome Biol 11, no. 10 (2010): R106.
DESeq_bibtex=deseq.bib

####
edgeR_citation=Robinson, Mark D., Davis J. McCarthy, and Gordon K. Smyth. "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." Bioinformatics 26, no. 1 (2010): 139-140.
edgeR_bibtex=edger.bib

####
mmseq_citation=Turro, Ernest, Shu-Yi Su, Angela Gonçalves, L. J. Coin, Sylvia Richardson, and Alex Lewin. "Haplotype and isoform specific expression estimation using multi-mapping RNA-seq reads." Genome Biol 12, no. 2 (2011): R13.
mmseq_bibtex=

####
DEXSeq_citation=Anders, Simon, Alejandro Reyes, and Wolfgang Huber. "Detecting differential usage of exons from RNA-seq data." Genome Research (2012).
DEXSeq_bibtex=dexseq.bib

####
baySeq_citation=Hardcastle, Thomas J., and Krystyna A. Kelly. "baySeq: empirical Bayesian methods for identifying differential expression in sequence count data." BMC bioinformatics 11, no. 1 (2010): 422.
baySeq_bibtex=bayseq.bib


piano_citation=Väremo, Leif, Jens Nielsen, and Intawat Nookaew. "Enriching the gene set analysis of genome-wide data by incorporating directionality of gene expression and combining statistical hypotheses and methods". Nucleic acids research 41.8 (2013): 4378-4391.
piano_bibtex=piano.bib

#################################################################
phony_targets+= show_citations

# prog 1
# citation 2
# bibtex_file 3
# prog version 4
define pprint_prog_info=

echo "* $(1) ; Version:  $(4) ; Citation: $(2)"

endef

# prog_name - 1
define prog_info=
$(if $($(1)_citation),$(call pprint_prog_info,$(1),$($(1)_citation),$($(1)_bibtex),$($(1)_version)),)
endef

# Scientific software
progs_used=$(qc_$(qual_filtering))  $(mapper_$(mapper)) $(quant_$(quant_method))  $(quant_norm_$(quant_norm_method))  $(de_$(de_method)) $(gse_$(gse_tool)) IRAP

silent_targets+= show_citations citations_file

show_citations: 
	$(info Summary of software used:)
	$(foreach p,$(sort $(progs_used)),$(call prog_info,$(p)))

citations_file: $(name)/report/software.txt

#.PHONY: $(name)/report/software.txt
# it is only dependent on the configuration file
$(name)/report/software.txt: $(conf)
	rm -f $(name)/report/software.txt
	touch $@.tmp && \
	( $(foreach p,$(sort $(progs_used)),$(call prog_info,$(p))) ) >> $@.tmp &&\
	mv $@.tmp $@



