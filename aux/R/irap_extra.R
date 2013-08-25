# functions under devel & testing
IRAP.DIR <- Sys.getenv(c("IRAP_DIR"))
if ( IRAP.DIR == "" ) {
  cat("ERROR: environment variable IRAP_DIR is not set\n")
  q(status=1)
}

transcripts.per.gene <- function(mat) {
  if ( c("gene_id","transcript_id") %in% colnames(mat)) {
    pwarning("transcripts.per.gene: row names incorrect",colnames(mat))
    return(NULL)
  }
  xu<-unique(mat[,c("gene_id","transcript_id")])
  qntrans<-table(xu$gene_id)
  qntrans<-as.data.frame(qntrans)
  names(qntrans)<-c("id","ntrans")
  qntrans
}
# info about #transcripts per gene
load.transbygene <- function() {
  file <- paste(IRAP.DIR,"/homosapiens.trans.per.gene.tsv",sep="")
  x<-read.tsv(file)
  rownames(x) <- x[,"id"]
  x
}
#plot(databox$quant.gene.nerror[,1],gtrans[rownames(databox$quant.gene.nerror),2],method="spearman")
