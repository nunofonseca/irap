#!/usr/bin/env Rscript
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

###############################################################
suppressPackageStartupMessages(library("optparse"))

IRAP.DIR <- Sys.getenv(c("IRAP_DIR"))
if ( IRAP.DIR == "" ) {
  cat("ERROR: environment variable IRAP_DIR is not set\n")
  q(status=1)
}
#
# specify our desired options in a list
#
source(paste(IRAP.DIR,"aux/R","irap_utils.R",sep="/"))
pdebug.enabled <- TRUE
#######################
usage <- "irap_cont_report.R --reads num.reads --bam file --out out_file_prefix"
filenames <- c("bam_file")

option_list <- list(
    make_option(c("-i", "--bam"), type="character", dest="bam_file", help="BAM file"),
    make_option(c("-o", "--out"), type="character",default=NULL,help="Output file name prefix."),
    make_option(c("-r", "--reads"), type="numeric",dest="reads",default=0,help="Number of reads checked"),
  make_option(c("--debug"),action="store_true",dest="debug",default=FALSE,help="Debug mode")
)

# check multiple options values
##multiple.options = list(
##  metric=c('rpkm','tpm','fpkm','fpkm-uq','uq-fpkm')
##  )
multiple.options <- NULL
mandatory <- c("bam_file","out")

#pinfo("saved")
opt <- myParseArgs(usage = usage, option_list=option_list,filenames.exist=filenames,multiple.options=multiple.options,mandatory=mandatory)

pdebug.enabled <- opt$debug


## Stats from the BAM
p <- pinfo("Generating ",opt$out,".tsv")
cmd <- paste0("samtools view ",opt$bam_file," |cut -f 3|sort |uniq -c|awk 'BEGIN {OFS=\"\t\";} {print $2,$1;}' > ",opt$out,".tmp")
status <- system(cmd)
if ( status != 0 ) {
    q(status=status)
}

## generic stats
cmd <- paste0("irapBAM2stats bam=",opt$bam_file)
status <- system(cmd)
if ( status != 0 ) {
    q(status=status)
}


# read the tables
x <- read.table(paste0(opt$out,".tmp"),sep="\t")
x$class <- gsub(":.*","",x$V1)
y <- aggregate(x$V2,by=list(class=x$class),FUN=sum)

# read the tables
z <- read.table(paste0(opt$bam_file,".stats.csv"),sep=",")
z <- rbind(z,c("ReadsUnmapped",opt$reads))
y$class <- paste("class:",y$class,sep="")
x$class <-  paste("species:",x$V1,sep="")
x <- x[,c("class","V2")]

colnames(z) <- c("label","N")
colnames(y) <- c("label","N")
colnames(x) <- c("label","N")
z <- rbind(z,y)
z <- rbind(z,x)

status=system(paste0("rm ",opt$out,".tmp "))
pinfo("Creating ",opt$out,"...")
write.tsv(z,opt$out,header=T)
pinfo("All done!")
# 
q(status=status)

# TODO: validate args
#if (!is.na(libs.info)) {
#  libs.info.v <- strsplit(libs.info,",")[[1]];
#} else {
#
#}

if (is.na(min.count)) {
  min.count <- 0
}

pinfo("starting...")

source(paste(IRAP.DIR,"aux/R","irap_misc.R",sep="/"))
library(R2HTML)

###############################################
pinfo("reading table...")

table <- read.table(data.file,sep="\t",header=T)
rownames(table) <- as.character(table[,1])
table <- table[,-1]

pinfo("reading table...done.")

system(paste("mkdir -p ",out.dir));

HTMLInitFile( outdir=out.dir, filename="index", Title=paste("",sep=""), useGrid=FALSE)
pinfo("HTML file opened")

# group
if (!is.na(groups)) {
  v <- strsplit(groups,",")[[1]];
  HTML("<H1></H1>")
  #
  groups <- unique(v)
  pinfo("groups:",groups)
  df<-data.frame(matrix(ncol=length(groups),nrow=nrow(table)))
  colnames(df)<-groups
  rownames(df) <- rownames(table)
  for ( g in groups ) {
    df[,g] <- rowSums(table[,v==g])
  }
  print(head(df))
  deseq.htseq.summary.report(df,groups,min.count,"grouped_libs",html.dir=out.dir)
}

HTML("<H1></H1>")
deseq.htseq.summary.report(table,colnames(table),min.count,"all_libs",html.dir=out.dir)

HTMLEndFile()

q()
