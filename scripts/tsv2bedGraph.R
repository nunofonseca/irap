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

usage<-"tsv2bedGraph tsv_file feature@(gene,exon,CDS)  gff3_file.csv [chr_file]\n";
if (length(args)<3||length(args)>4) {
  cat("ERROR!");
  cat(usage);
  q(status=1);
}
# chr file: should contain the name of all chromossomes (can be used to filter)
#    chr_name\tlength


# TODO: add option to filter by source (e.g protein coding)
# validate arguments
tsv_file <-args[1]
feature<-args[2]
gff3_file<-args[3]
chr_file <- args[4]
chrs <- c()
#feature<-"gene"
#gff3_file<-"test_se_2l/data/chr19.gff3"
#tsv_file <-"test_se_2l/tophat1/htseq1/SE2.genes.htseq1.tsv"

# load the tsv file (2 cols)
tsv<-read.table(tsv_file,sep="\t",header=F)
# first column should have the ID
colnames(tsv) <- c('ID','Value')

# load the gff3 file in CSV format (it can take a while to do it, hence do it only once)
gff3<-read.table(gff3_file,sep=",",header=T)
## gff3<-read.table(gff3_file,sep="\t",header=F)
## all_names<-c('seqname','source','feature','start','end','score','strand','frame','attributes','comments')
## colnames(gff3)<-all_names[0:ncol(gff3)]
## # add the attributes
## tags<-c("ID","Name","Alias","Parent","Target","Gap","Derives","Note","DBxref","Ontology_term","Is_circular")
## for (tag in tags ) {
##   gff3[,tag]<-rep(x=NA,nrow(gff3))
## }
## for (n in c(1:nrow(gff3)) ) {  
##   att<-strsplit(as.character(gff3$attributes[n]),split=";")
##   # look for the Id
##   for ( a in att ) {
##     s<-strsplit(a,split="=")
##     #print(s)
##     gff3[n,s[[1]][1]] <- s[[1]][2]
##   }
## }

sel<-gff3[gff3$feature==feature,]

#
sel<-sel[,c('ID','seqname','start','end')]
sel$ID <- as.character(sel$ID)
tsv$ID <- as.character(tsv$ID)
t1<-merge(sel,tsv,by="ID")
t1<-t1[order(t1$start),]
track_label<- "default"
# filter
if ( ! is.null(chr_file) ) {
  if ( ! is.null(chr <- file) ) {
    chr.table<-read.table(chr_file,sep="\t",header=F)
    chrs <- as.character(chr.table[,1])
    #print(chrs)
  }
  selection <- (t1$seqname %in% chrs)
  t1 <- t1[selection,]
}
# bedgraph does not support overlapping features
# 
# bedgraph
#It is possible to have different values for the same coord
cat(paste("track type=bedGraph name=",track_label," description=",track_label," visibility=full\n",sep=""))
#chr star end val
write.table(t1[,-1],sep='\t',row.names=F,col.names=F,quote=F)
q(status=0)

