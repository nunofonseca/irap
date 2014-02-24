#!/usr/bin/env Rscript
## /*******************************************************************************************
##  *
##  * Copyright (c) 2012 Nuno Fonseca. All rights reserved.
##  * This code is freely available for academic purposes.
##  * If you intend to use it for commercial purposes then please contact the author first.
##
##  * Author: Nuno A. Fonseca 
##  * Date: 2012-03-24
##  * $Id: irap.txt Nuno Fonseca Wed Jan 16 23:26:59 2013$
##  *
##  ******************************************************************************************/
suppressPackageStartupMessages(library("optparse"))

IRAP.DIR <- Sys.getenv(c("IRAP_DIR"))
if ( IRAP.DIR == "" ) {
  cat("ERROR: environment variable IRAP_DIR is not set\n")
  q(status=1)
}

# specify our desired options in a list
#
source(paste(IRAP.DIR,"aux/R","irap_utils.R",sep="/"))
pdebug.enabled <- TRUE

html.template <- get.path2template("mappers_comp_sum")




# compare summary info of multiple mappers
# mappers_comp_sum.R "label1,lebel2,..." tsv_file1  tsv_file.... 
usage <- "mappers_comp_sum.R --labels 'label1,label2,...' --tsv 'tsvfile1,tsvfile2,...' [options]"
filenames <- c("tsv_file","annotation") ;#filenames that must exist (if defined)
option_list <- list(
                    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
                    make_option(c("-q", "--quietly"), action="store_false", dest="verbose", help="Print little output"),
                    make_option(c("-i", "--tsv"), type="character", dest="tsv_files", default=NULL,help="Comma or space separated list of TSV file names"),
                    make_option(c("-l", "--labels"), type="character", dest="labels", default=NULL,help="Comma or space separated names to associate to each TSV file"),
                    make_option(c("-o", "--out"), type="character",default=NULL,help="Output file name prefix. An HTML and TSV.gz file will be created."),
                    
                    make_option(c("-t", "--title"), type="character", dest="title", default="",help="Report title"),
                    make_option(c("--css"), type="character",default="",help="Path to the irap.css file")
)

# check multiple options values
multiple.options = list()
mandatory <- c("tsv_files","out","labels")
opt <- myParseArgs(usage = usage, option_list=option_list,filenames.exist=filenames,multiple.options=multiple.options,mandatory=mandatory)


opt$tsv_files <- strsplit(mytrim(opt$tsv_files),"[ ,]+")[[1]];
opt$labels <- strsplit(mytrim(opt$labels),"[ ,]+")[[1]];
o.html.file <-paste(opt$out,".html",sep="")
pdebug("TSV FILES=",opt$tsv_files)
pdebug("LABELS=",opt$labels)
pdebug("out=",opt$out)

if ( length(opt$labels)<=1 ) {
  pdebug("Nothing to compare...")
  q(status=0)
}
library(brew)
#

# load each TSV FILE
pdebug("Loading TSV files...")
all <- load.matrices(opt$tsv_files)
pdebug("Loading TSV files...done.")

# Average and std

# One matrix with the average and the other with the sd
aggr.med <- matrix(nrow=nrow(all[[1]]),ncol=0)
aggr.sd <- matrix(nrow=nrow(all[[1]]),ncol=0)
for ( m in names(all)) {
  rownames(all[[m]]) <- all[[m]][,1]
  #remove first col
  all[[m]]<-all[[m]][,-1]
  med <- apply(all[[m]],1,FUN=median)
  sd <- apply(all[[m]],1,FUN=sd)
  aggr.sd <- cbind(aggr.sd,sd)
  aggr.med <- cbind(aggr.med,med)
}
#  rownames(aggr.all[[m]]) <- rownames((all[[m]]))
colnames(aggr.med) <- opt$labels
colnames(aggr.sd) <- opt$labels

get.ptable <- function(aggr.med,aggr.sd) {
  
  total.reads=aggr.med["Reads mapped",]+aggr.med["Unmapped",]
  denom <- matrix(rep(total.reads,nrow(aggr.med)),byrow=T,nrow=nrow(aggr.med))
  med <- aggr.med/denom*100
  sd  <- aggr.sd/denom*100
  ci.u <- med+sd
  ci.l <- med-sd
  ci.l[ci.l<0] <- 0
  return(list(med=med,sd=sd,u=ci.u,l=ci.l))
}

ci.l <- aggr.med-aggr.sd
ci.l[ci.l<0] <- 0
ci.u <- aggr.med+aggr.sd

aggr <- list(med=aggr.med,sd=aggr.sd,l=ci.l,u=ci.u)
paggr <- get.ptable(aggr.med,aggr.sd)

title <- opt$title
file.prefix <- opt$out
#save.image()
pdebug("Generating html ",o.html.file)
brew.wrapper(html.template,o.html.file)


q(status=0)

scripts/mappers_comp_sum.R --tsv "test_files/th2_align_overall_comparison.png.tsv  test_files/th1_align_overall_comparison.png.tsv "  --labels "TH2 TH1" -o 2del




