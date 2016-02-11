#!/usr/bin/env Rscript
# =========================================================
# Copyright 2012-2016,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
# R script to install all R packages required by iRAP
repo<-"http://www.stats.bris.ac.uk/R/"

print(.libPaths())

######
# upgrade all installed packages before starting the installation
try(update.packages(repos=repo,instlib=.libPaths()[1],ask=FALSE))

#biocLite("BiocUpgrade",ask=FALSE, suppressUpdates=TRUE)
#try(remove.packages("BiocInstaller"))
source("http://bioconductor.org/biocLite.R")
#biocLite("BiocInstaller",ask=FALSE, suppressUpdates=TRUE)
# not available in older versions of bioconductor
try(biocUpdatePackages(pkgs = c("BiocGenerics","Biobase", "IRanges","AnnotationDbi"),ask=FALSE,instlib=.libPaths()[1]))

source("http://bioconductor.org/biocLite.R")

#biocValid()
#################
cat("-----------------------------------------------------\n")
cat("Installing packages:\n")

packages2install<-c("intervals","gclus",'R2HTML',"agricolae",
                    "optparse","brew","reshape","gtools","gdata","caTools",
                    "sfsmisc","gplots","lattice","data.table",
                    'edgeR',
                    'DESeq','DESeq2','DEXSeq','baySeq',
                    'limma','marray','igraph')
for (p in packages2install ) {
  cat("PACKAGE:",p,"\n")
  biocLite(p,ask=FALSE, suppressUpdates=TRUE)
}


biocLite("piano",ask=FALSE, suppressUpdates=TRUE)
#biocLite("org.Hs.eg.db",ask=FALSE, suppressUpdates=TRUE)
biocLite('RCurl',ask=FALSE, suppressUpdates=TRUE)
# http://bioconductor.org/packages/2.13/data/annotation/src/contrib/GO.db_2.10.1.tar.gz fails to install
#biocLite('GO.db',ask=FALSE, suppressUpdates=TRUE)
#biocLite("topGO",ask=FALSE, suppressUpdates=TRUE)
#biocLite("biomaRt",ask=FALSE, suppressUpdates=TRUE)

#biocLite('goseq',ask=FALSE, suppressUpdates=TRUE)

species2db<-matrix(c('org.Ag.eg.db','Anopheles',
                     'org.At.tair.db','Arabidopsis',
                     'org.Bt.eg.db','Bovine',
                     'org.Ce.eg.db','Worm',
                     'org.Cf.eg.db','Canine',
                     'org.Dm.eg.db','Fly',
                     'org.Dr.eg.db','Zebrafish',
                     'org.EcK12.eg.db','E coli strain K12',
                     'org.Gg.eg.db','Chicken',
                     'org.Hs.eg.db','Human',
                     'org.Mm.eg.db','Mouse',
                     'org.Mmu.eg.db','Rhesus',
                     'org.Pf.plasmo.db','Malaria',
                     'org.Pt.eg.db','Chimp',
                     'org.Rn.eg.db','Rat',
                     'org.Sc.sgd.db','Yeast',
                     'org.Sco.eg.db','Streptomyces coelicolor',
                     'org.Ss.eg.db','Pig',
                     'org.Tgondii.eg.db','Toxoplasma gondii',
                     'org.Xl.eg.db','Xenopus'),byrow=T,ncol=2)
colnames(species2db)<-c("db","species")
for (p in species2db[,'db']) {
  biocLite(p,ask=FALSE,suppressUpdates=TRUE)
}

q(status=0)
