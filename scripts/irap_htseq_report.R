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
#    $Id$
# =========================================================

###############################################################



IRAP.DIR <- Sys.getenv(c("IRAP_DIR"))
if ( IRAP.DIR == "" ) {
  cat("ERROR: environment variable IRAP_DIR is not set\n")
  q(status=1)
}
source(paste(IRAP.DIR,"aux/R","irap_utils.R",sep="/"))
#source(paste(IRAP.DIR,"scripts","irap_misc.R",sep="/"))
             
args <- commandArgs(trailingOnly=TRUE)
out.dir <- args[1]
data.file <- args[2]
min.count <- args[3]
groups <- args[4]
lib.info <- args[5]
if ( is.na(data.file) ) {
  perror("Usage: irap_htseq_report.R out_dir htseq_count_file [min_count] [groups] [lib.info]")
  q(status=1)
}

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
