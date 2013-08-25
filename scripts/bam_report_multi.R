#!/bin/env Rscript
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



library(R2HTML)

#library(myUtils)
library(emBAM)
options("cores"=multicore:::detectCores())

# Parallelize
#
args <- commandArgs(trailingOnly=TRUE)
#}
# optional
html.dir <- args[1]
bam.files <- args[-1]

# duplicated code (see irap_report_mapping)
bam.compare.plots <- function(html.dir,bam.files) {

  data <- lapply(bam.files,em.bam.counts.df,TRUE)
  df <- as.data.frame(data[[1]][,1])
  rownames(df)<-rownames(data[[1]])

  if(length(data)>1) {
    for (i in c(2:length(data))) {
      df <- cbind(df,data[[i]]$Count)
    }
  }
  colnames(df) <- gsub(".*/","",gsub("(\\.pe|\\.se).*","",bam.files))

  ########
  #
  legend <- c("Alignments","Reads mapped","Unmapped","Uniquely mapped","Spliced reads")
  a.data <- df[c("Alignments","Primary","Unmapped","Uniquely mapped reads","Spliced reads"),]
  rownames(a.data) <- legend
  file <- "align_overall_comparison.png"
  filepath <- paste(html.dir,file,sep="/")
  cols <- rainbow(length(rownames(a.data)))
  cex <- 0.95

  n <- length(colnames(a.data))
  width=(480*round(0.5+n/6/2,0))
  png(filepath,width=width,height=480)
  bp <- barplot(as.matrix(a.data),ylab="Number of alignments",las=2 , cex.axis=cex, cex.names=cex*0.90,col=cols,beside=T)

  par( xpd=NA )  
  legend( "top", legend=rownames(a.data),fill=cols,horiz=TRUE, cex=0.90)

  # save TSV file
  data2file <- as.matrix(a.data)
  data2file <- cbind(rownames(data2file),data2file)
  write.tsv(data2file,paste(filepath,".tsv",sep=""))

  # Add %
  #
  dev.off()
}

##################################################
system(paste("mkdir -p ",html.dir,sep=""))
bam.compare.plots(html.dir,bam.files)


###################################################
