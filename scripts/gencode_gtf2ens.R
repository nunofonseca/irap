#!/usr/bin/env Rscript
# =========================================================
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
#    $Id: irap.txt Nuno Fonseca Sun Jan 13 14:02:59 2013$
# =========================================================
args <- commandArgs(trailingOnly=TRUE)

usage<-"gencode_gtf2ens gtf_file\n";
if (length(args)!=1) {
  cat("ERROR!");
  cat(usage);
  q(status=1);
}


gtf_file<-args[1]

######################################################
load.gencode.gtf2ens <- function(gtf.file) {
  
  gtf<-read.table(gtf.file,sep="\t",header=F,quote="\"")
  cnames <- c("seqid","source","feature","start","end","score","strand","frame","attributes")
  colnames(gtf)<-cnames[0:ncol(gtf)]
  gtf.attributes.names<-c("gene_type")
  gene.type <- gsub(";.*","",gsub(".* gene_type ","",gtf$attributes))
  head(gene.type)
  gtf$source <- gene.type
  return(gtf)
}

gtf <- load.gencode.gtf2ens(gtf_file)
write.table(gtf,sep='\t',row.names=F,col.names=F,quote=F)
q(status=0)
