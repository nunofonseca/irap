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
# tsv2bedGraph.R test_se_2l/tophat1/htseq1/SE1.se.genes.raw.htseq1.tsv gene test_se_2l/data/chr19.gff3 > a.bg
args <- commandArgs(trailingOnly=TRUE)

usage<-"gff32csv gff3_file\n";
if (length(args)!=1) {
  cat("ERROR!");
  cat(usage);
  q(status=1);
}

# TODO: add option to filter by source (e.g protein coding)
# validate arguments
gff3_file<-args[1]

# load the gff3 file
gff3<-read.table(gff3_file,sep="\t",header=F)
all_names<-c('seqname','source','feature','start','end','score','strand','frame','attributes','comments')
colnames(gff3)<-all_names[0:ncol(gff3)]
gff3$attributes <- as.character(gff3$attributes)
# convert the start/end to numbers
gff3$start <- as.integer(gff3$start)
gff3$end <- as.integer(gff3$end)

# add the attributes
tags<-c("ID","Name","Alias","Parent","Target","Gap","Derives","Note","DBxref","Ontology_term","Is_circular")

#tags <- c("ID")
for (tag in tags ) {
  pattern <- paste("",tag,"=([^;]+)",sep="")
  m<-regexec(pattern,as.vector(gff3$attributes))
#m<-regexec("^ID=([^;]+);",as.vector(gff3$attributes))
  for ( i in c(1:length(m)) ) { if ( m[[i]][1]==-1 ) { m[[i]]=NA; } }
  x<-regmatches(as.vector(gff3$attributes),m)
  for ( i in c(1:length(x)) ) { if ( length(x[[i]])==0 ) { x[[i]]=c(NA,NA); } }
  vals<-matrix(unlist(x),ncol=2,byrow=T,)[,2]
  gff3[,tag]<-vals
}


## # very slow...
## get_attrs <- function(v,attr) {
##   att<-strsplit(as.character(v),split=";")
##   for ( a in att[[1]] ) {
##     s<-strsplit(a,split="=")
##     var<-s[[1]][1]
##     value<-s[[1]][2]
##     if (var==attr) {
##       return(value)
##     }
##     #print(paste(">",var,"=",value,"",sep=""))
##     #v[var] <- value
##   }
##   NA
## }
## for (n in c(1:nrow(gff3)) ) {
##   gff3[n,"ID"] <- get_attrs(gff3[n,"attributes"],"ID")
## }

write.table(gff3,sep=',',row.names=F,col.names=T,quote=F)
q(status=0)
