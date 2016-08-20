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

#.Library.site<-c()
#.Library <- c()
#.libPaths(.libPaths()[1])
message("Using library: ", .libPaths()[1])
print(.libPaths())

# Check if version is ok
version <- getRversion()
currentVersion <- sprintf("%d.%d", version$major, version$minor)
message("R version:",version)
if ( version$major < 3 || (version$major>=3 && version$minor<2) ) {
  cat("ERROR: R version should be 3.2 or above\n")
  q(status=1)
}

##########################################################
# install all packages to the iRAP folder
lib.vars <- ls(environment(.libPaths),all=TRUE)
if ( length(lib.vars)!=1 && lib.vars[1]!=".lib.loc") {
  cat("ERROR: .lib.loc not found...unable to proceed.\n")
  q(status=1)
}

# Keep only the iRAP path (that should be the first from .libPaths())
assign(".lib.loc",.libPaths()[1],envir=environment(.libPaths))

message("Installing packages to: ", .libPaths())
message("_____________________________________________________")
source("http://bioconductor.org/biocLite.R")


######
# Update all packages
# update.packages(repos=repo,instlib=.libPaths()[1],ask=FALSE)


#biocValid()
#################
cat("-----------------------------------------------------\n")
cat("Installing packages:\n")

# Install an older version of the package (required by DEseq2 and DEXseq)
biocLite("Rcpp",ask=FALSE, suppressUpdates=FALSE)
install.packages("http://cran.rstudio.com/src/contrib/00Archive/RcppArmadillo/RcppArmadillo_0.6.600.4.0.tar.gz")

packages2install<-c("intervals","gclus",'R2HTML',"agricolae","bit64",
                    "optparse","brew","reshape","gtools","gdata","caTools",
                    "sfsmisc","gplots","lattice","data.table",
                    'edgeR','piano','RCurl','GO.db',
                    'DESeq','xtable','DESeq2','DEXSeq','baySeq',
                    'limma','marray','igraph')
for (p in packages2install ) {
  cat("PACKAGE:",p,"\n")
  biocLite(p,ask=FALSE, suppressUpdates=c("^RcppArmadillo"))
}

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
  biocLite(p,ask=FALSE,, suppressUpdates=c("^RcppArmadillo"))
}

q(status=0)
