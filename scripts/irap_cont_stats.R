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

# Setup files

reads_per_chr_file <- paste0(opt$out,".tmp")
generic_stats_file <- paste0(opt$bam_file,".stats.csv")

## Reads per chr from the BAM
pinfo("Generating ",opt$out,".tsv")
cmd <- paste0("samtools view ",opt$bam_file," |cut -f 3|sort |uniq -c|awk 'BEGIN {OFS=\"\t\";} {print $2,$1;}' > ", reads_per_chr_file)
status <- system(cmd)
if ( status != 0 ) {
    q(status=status)
}

## generic stats
pinfo("irapBAM2stats bam=",opt$bam_file)
cmd <- paste0("irapBAM2stats bam=",opt$bam_file)
status <- system(cmd)
if ( status != 0 ) {
    q(status=status)
}

## Read the tables to generate a contamination summary

# The generic status should be populated whatever

out_col_names <- c("label","N")
contamination_summary <- read.table(generic_stats_file, sep=",", col.names = out_col_names)
contamination_summary <- rbind(contamination_summary,c("ReadsUnmapped",opt$reads))

# The reads per chr may not be populated if 0 reads mapped to chromosomes in
# the contamination indices

if ( file.size(reads_per_chr_file) > 0){
    reads_per_chr <- read.tsv(reads_per_chr_file, header=FALSE)
    reads_per_chr$class <- gsub(":.*","",reads_per_chr$V1)
    reads_per_class <- aggregate(reads_per_chr$V2,by=list(class=reads_per_chr$class),FUN=sum)

    reads_per_class$class <- paste("class:",reads_per_class$class,sep="")
    reads_per_chr$class <-  paste("species:",reads_per_chr$V1,sep="")
    reads_per_chr <- reads_per_chr[,c("class","V2")]

    colnames(reads_per_chr) <- colnames(reads_per_class) <- out_col_names

    # Add the results into the summary
    contamination_summary <- do.call(rbind, list(contamination_summary, reads_per_class, reads_per_chr))
}else{
    pinfo("There were no mappings to contamination indices")
}

status=system(paste0("rm ", reads_per_chr_file ))
pinfo("Creating ",opt$out,"...")
write.tsv(contamination_summary, opt$out,header=T)
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
