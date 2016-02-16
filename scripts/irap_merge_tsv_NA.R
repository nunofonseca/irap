#!/usr/bin/env Rscript

IRAP.DIR <- Sys.getenv(c("IRAP_DIR"))
if ( IRAP.DIR == "" ) {
  cat("ERROR: environment variable IRAP_DIR is not set\n")
  q(status=1)
}
source(paste(IRAP.DIR,"aux/R","irap_utils.R",sep="/"))

#args <- commandArgs(trailingOnly=TRUE)
args <- importArgsfromStdin()
#pwarning(args)
if (length(args)<1) {
  cat("ERROR! usage: irap_merge_tsv_NA.R  file1.tsv [file2.tsv file3.tsv ...]\n")
  q(status=1);
}

files <- args
f1<-files[1]
if (!file.exists(f1)) {
  perror("Unable to find file ",f1)
  q(status=1)
}
files<-files[-1]
t1<-read.table(f1,sep="\t",as.is=c(T,F),comment.char="",quote="\"")
t1 <- fix.cufflinks.fpkms(t1)
colnames(t1)<-c("Gene",basename(f1))
for ( f in  files ) {
        if (!file.exists(f)) {
          perror("Unable to find file ",f)
          q(status=1)
        }
	t2<-read.table(f,sep="\t",as.is=c(T,F),comment.char="",quote="\"")
        t2 <- fix.cufflinks.fpkms(t2)
        colnames(t2) <- c("Gene",basename(f))
	m<-merge(t1,t2,by="Gene",all=TRUE)
        t1 <- NULL
        t2 <- NULL
        t1 <- m
        m  <- NULL        
}
#change NAs to 0
t1[is.na(t1)]<-0
write.table(t1,sep="\t",quote=F,row.names=F,col.names=T)
q(status=0)
