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
mapper_HISAT2=HISAT2

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
quant_nurd=NURD
quant_salmon=salmon
quant_kallisto=kallisto
quant_stringtie=stringtie
quant_stringtie_nd=stringtie
quant_rsem=rsem

# quant_norm_method cufflinks1 cufflinks2 cufflinks1_nd cufflinks2_nd deseq flux_cap edger 
quant_norm_cufflinks1=cufflinks1
quant_norm_cufflinks1_nd=cufflinks1
quant_norm_cufflinks2=cufflinks2
quant_norm_cufflinks2_nd=cufflinks2
quant_norm_flux_cap=FLUX_CAPACITOR
quant_norm_deseq=DESeq
quant_norm_irap=IRAP

# DE
de_cuffdiff1=cuffdiff1
de_cuffdiff2=cuffdiff2
de_cuffdiff1_nd=cuffdiff1
de_cuffdiff2_nd=cuffdiff2
de_deseq=DESeq
de_deseq2=DESeq2
de_dexseq=DEXSeq
de_bayseq=baySeq
de_voom=baySeq
de_ebseq=EBSeq

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
IRAP_citation="Nuno A. Fonseca, Robert Petryszak, John Marioni, Alvis Brazma. iRAP - an integrated RNA-seq Analysis Pipeline. bioRxiv 2014. doi: http://dx.doi.org/10.1101/005991"
IRAP_bibtex=irap.bib
IRAP_version=$(version)
IRAP_analysis="Pipeline"

#TODO: update
STAR_citation="Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, and Thomas R. Gingeras. STAR: ultrafast universal RNA-seq aligner Bioinformatics first published online October 25, 2012"
STAR_bibtex=STAR.bib
STAR_analysis="Reads alignment"

HISAT2_citation="Kim, Daehwan, Ben Langmead, and Steven L. Salzberg. "HISAT: a fast spliced aligner with low memory requirements." Nature methods 12.4 (2015): 357-360."
HISAT2_bibtex=HISAT2.bib
HISAT_analysis="Reads alignment"

#### 
tophat1_citation=Cole Trapnell, Lior Pachter, and Steven L. Salzberg. TopHat: discovering splice junctions with RNA-Seq Bioinformatics (2009) 25(9): 1105-1111
tophat1_bibtex=tophat1.bib
tophat1_analysis="Reads alignment"

#### 
tophat2_citation=Kim, Daehwan, et al. "TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions." Genome Biol 14.4 (2013): R36.
tophat2_bibtex=tophat2.bib
tophat2_analysis="Reads alignment"

#### 
SMALT_citation=http://www.sanger.ac.uk/resources/software/smalt/
SMALT_bibtex= 
SMALT_analysis="Reads alignment"



#### 
GSNAP_citation=Thomas D. Wu and Serban Nacu.Fast and SNP-tolerant detection of complex variants and splicing in short reads Bioinformatics (2010) 26(7): 873-881 
GSNAP_bibtex=gsnap.bib
GSNAP_analysis="Reads alignment"

#### 
SOAPsplice_citation=Huang, Songbo, Jinbo Zhang, Ruiqiang Li, Wenqian Zhang, Zengquan He, Tak-Wah Lam, Zhiyu Peng, and Siu-Ming Yiu. "SOAPsplice: genome-wide ab initio detection of splice junctions from RNA-Seq data." Frontiers in genetics 2 (2011).
SOAPsplice_bibtex=soapsplice.bib
SOAPsplice_analysis="Reads alignment"

#### 
bwa1_citation=Li, Heng, and Richard Durbin. "Fast and accurate short read alignment with Burrows–Wheeler transform." Bioinformatics 25, no. 14 (2009): 1754-1760.
bwa1_bibtex=bwa1.bib
bwa1_analysis="Reads alignment"

#### 
bwa2_citation=Li, Heng, and Richard Durbin. "Fast and accurate long-read alignment with Burrows–Wheeler transform." Bioinformatics 26, no. 5 (2010): 589-595.
bwa2_bibtex=bwa2.bib
bwa2_analysis="Reads alignment"

#### 
bowtie1_citation=Langmead, Ben, Cole Trapnell, Mihai Pop, and Steven L. Salzberg. "Ultrafast and memory-efficient alignment of short DNA sequences to the human genome." Genome Biol 10, no. 3 (2009): R25.
bowtie1_bibtex=bowtie1.bib
bowtie1_analysis="Reads alignment"

#### 
bowtie2_citation=Langmead, Ben, and Steven L. Salzberg. "Fast gapped-read alignment with Bowtie 2." Nature methods 9, no. 4 (2012): 357-359.
bowtie2_bibtex=bowtie2.bib
bowtie2_analysis="Reads alignment"

#### 
osa_citation=Hu, Jun, Huanying Ge, Matt Newman, and Kejun Liu. "OSA: a fast and accurate alignment tool for RNA-Seq." Bioinformatics 28, no. 14 (2012): 1933-1934.
osa_bibtex=osa.bib
osa_analysis="Reads alignment"

#### 
GEM_citation=Marco-Sola, Santiago, Michael Sammeth, Roderic Guigó, and Paolo Ribeca. "The GEM mapper: fast, accurate and versatile alignment by filtration." Nature Methods 9, no. 12 (2012): 1185-1188.
GEM_bibtex=gem.bib
GEM_analysis="Reads alignment"

#### 
FASTQC_citation=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
FASTQC_bibtex=
FASTQC_analysis="QC"

#### 
FASTX_citation=http://hannonlab.cshl.edu/fastx_toolkit/
FASTX_bibtex=
FASTX_analysis="Filtering and trimming"

#### 
cufflinks1_citation=Trapnell, Cole, Brian A. Williams, Geo Pertea, Ali Mortazavi, Gordon Kwan, Marijke J. Van Baren, Steven L. Salzberg, Barbara J. Wold, and Lior Pachter. "Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation." Nature biotechnology 28, no. 5 (2010): 511-515.
cufflinks1_bibtex=cufflinks1.bib
cufflinks1_analysis="Gene and transcript quantification"

####
cufflinks2_citation=Trapnell, Cole, et al. "Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks." Nature protocols 7.3 (2012): 562-578.
cufflinks2_bibtex=cufflinks2.bib
cufflinks2_analysis="Gene and transcript quantification"

#
rsem_citation=Li, Bo, and Colin N. Dewey. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC bioinformatics 12.1 (2011): 1.
rsem_bibtex=rsem.bib
rsem_analysis="Gene and transcript quantification"

####
htseq_citation=Simon Anders, Paul Theodor Pyl, Wolfgang Huber: HTSeq – A Python framework to work with high-throughput sequencing data. bioRxiv, 20 Feb 2014, http://doi.org/10.1101/002824.
htseq_bibtex=htseq.bib
htseq_analysis="Gene, transcript and exon quantification"

####
FLUX_CAPACITOR_citation=Montgomery, Stephen B., Micha Sammeth, Maria Gutierrez-Arcelus, Radoslaw P. Lach, Catherine Ingle, James Nisbett, Roderic Guigo, and Emmanouil T. Dermitzakis. "Transcriptome genetics using second generation sequencing in a Caucasian population." Nature 464, no. 7289 (2010): 773-777.
FLUX_CAPACITOR_bibtex=flux_capacitor.bib
FLUX_CAPACITOR_analysis="Gene and transcript quantification"

####
BitSeq_citation=Glaus, Peter, Antti Honkela, and Magnus Rattray. "Identifying differentially expressed transcripts from RNA-seq data with biological variation." Bioinformatics 28, no. 13 (2012): 1721-1728.
BitSeq_bibtex=bitseq.bib
BitSeq_analysis="Differential transcript expression"

####
Ireckon_citation=Mezlini, Aziz M., Eric JM Smith, Marc Fiume, Orion Buske, Gleb Savich, Sohrab Shah, Sam Aparicion, Derek Chiang, Anna Goldenberg, and Michael Brudno. "iReckon: Simultaneous isoform discovery and abundance estimation from RNA-seq data." Genome Research (2012).
Ireckon_bibtex=ireckon.bib
Ireckon_analysis="Transcript quantification"

####
cuffdiff1_citation=Trapnell, Cole, Brian A. Williams, Geo Pertea, Ali Mortazavi, Gordon Kwan, Marijke J. Van Baren, Steven L. Salzberg, Barbara J. Wold, and Lior Pachter. "Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation." Nature biotechnology 28, no. 5 (2010): 511-515.
cuffdiff1_bibtex=cuffdiff1.bib
cuffdiff1_analysis="Differential gene expression"

####
cuffdiff2_citation=Trapnell, Cole, David G. Hendrickson, Martin Sauvageau, Loyal Goff, John L. Rinn, and Lior Pachter. "Differential analysis of gene regulation at transcript resolution with RNA-seq." Nature Biotechnology (2012).
cuffdiff2_bibtex=cuffdiff2.bib
cuffdiff2_analysis="Differential gene expression"

#### 
stringtie_citation=Mihaela Pertea, Geo M Pertea, Corina M Antonescu, Tsung-Cheng Chang, Joshua T Mendell and Steven L Salzberg. "StringTie enables improved reconstruction of a transcriptome from RNA-seq reads", Nature Biotechnology (2015).
stringtie_bibtex=stringtie.bib
stringtie_analysis="Gene/transcript quantification"

####
DESeq_citation=Anders, Simon, and Wolfgang Huber. "Differential expression analysis for sequence count data." Genome Biol 11, no. 10 (2010): R106.
DESeq_bibtex=deseq.bib
DESeq_analysis="Differential gene expression"


####
DESeq2_citation=Anders, Love MI, Huber W and Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology, 15, pp. 550.
DESeq2_bibtex=deseq2.bib
DESeq2_analysis="Differential gene expression"

EBSeq_citation=Leng, N., J.A. Dawson, J.A. Thomson, V. Ruotti, A.I. Rissman, B.M.G. Smits, J.D. Haag, M.N. Gould, R.M. Stewart, and C. Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments, Bioinformatics, 2013.
EBSeq_bibtex=ebseq.bib
EBSeq_analysis="Differential gene and transcript expression"

####
edgeR_citation=Robinson, Mark D., Davis J. McCarthy, and Gordon K. Smyth. "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." Bioinformatics 26, no. 1 (2010): 139-140.
edgeR_bibtex=edger.bib
edgeR_analysis="Differential gene expression"

####
#mmseq_citation=Turro, Ernest, Shu-Yi Su, Angela Gonçalves, L. J. Coin, Sylvia Richardson, and Alex Lewin. "Haplotype and isoform specific expression estimation using multi-mapping RNA-seq reads." Genome Biol 12, no. 2 (2011): R13.
#mmseq_bibtex=


####
DEXSeq_citation=Anders, Simon, Alejandro Reyes, and Wolfgang Huber. "Detecting differential usage of exons from RNA-seq data." Genome Research (2012).
DEXSeq_bibtex=dexseq.bib
DEXSeq_analysis="Differential exon usage"

####
baySeq_citation=Hardcastle, Thomas J., and Krystyna A. Kelly. "baySeq: empirical Bayesian methods for identifying differential expression in sequence count data." BMC bioinformatics 11, no. 1 (2010): 422.
baySeq_bibtex=bayseq.bib
baySeq_analysis="Differential gene expression"

kallisto_citation=Nicolas L Bray, Harold Pimentel, Páll Melsted, Lior Pachter. "Near-optimal probabilistic RNA-seq quantification." Nature Biotechnology (2016).
kallisto_bibtex=kallisto.bib
kallisto_analysis="Gene/transcript quantification"

salmon_citation=Patro, Rob, Geet Duggal, and Carl Kingsford. "Salmon: Accurate, versatile and ultrafast quantification from rna-seq data using lightweight-alignment." bioRxiv (2015): 021592
salmon_bibtex=salmon.bib
salmon_analysis="Gene/transcript quantification"


piano_citation=Väremo, Leif, Jens Nielsen, and Intawat Nookaew. "Enriching the gene set analysis of genome-wide data by incorporating directionality of gene expression and combining statistical hypotheses and methods". Nucleic acids research 41.8 (2013): 4378-4391.
piano_bibtex=piano.bib
piano_analysis="Enrichment analysis" 
#################################################################
phony_targets+= show_citations

# prog 1
# citation 2
# bibtex_file 3
# prog version 4
define pprint_prog_info=

echo "$(5)	$(1)	$(4)	$(2)"

endef

# prog_name - 1
define prog_info=
$(if $($(1)_citation),$(call pprint_prog_info,$(1),$($(1)_citation),$($(1)_bibtex),$($(1)_version),$($(1)_analysis)))
endef

# Scientific software
progs_used=$(qc_$(qual_filtering))  $(mapper_$(mapper)) $(quant_$(quant_method))  $(quant_norm_$(quant_norm_tool))  $(de_$(de_method)) $(gse_$(gse_tool)) IRAP

silent_targets+= show_citations citations_file

show_citations: 
	$(info Summary of software used:)
	echo "Analysis	Software	Version	Citation" > $@.tmp && \
	$(foreach p,$(sort $(progs_used)),$(call prog_info,$(p)))

citations_file: $(name)/report/software.tsv

%versions.tsv: $(conf)
	echo "Analysis	Software	Version	Citation" > $@.tmp && \
	( $(foreach p,$(sort $(progs_used)),$(call prog_info,$(p))) ) >> $@.tmp &&\
	mv $@.tmp $@

# it is only dependent on the configuration file
# software.tsv is a tab separated value file
$(name)/report/software.tsv: $(conf)
	rm -f $@
	touch $@.tmp && \
	echo "Analysis	Software	Version	Citation" > $@.tmp && \
	( $(foreach p,$(sort $(progs_used)),$(call prog_info,$(p))) ) >> $@.tmp &&\
	mv $@.tmp $@



