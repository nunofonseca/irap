#!/usr/bin/env Rscript3
## /*******************************************************************************************
##  *
##  * Copyright (c) 2012 Nuno Fonseca. All rights reserved.
##  * This code is freely available for academic purposes.
##  * If you intend to use it for commercial purposes then please contact the author first.
##
##  * Author: Nuno A. Fonseca 
##  * Date: 2012-03-24
##  * $Id: irap.txt Nuno Fonseca Sun Feb 3 15:37:04 2013$
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
source(paste(IRAP.DIR,"","emBAM.R",sep="/"))

pdebug.enabled <- FALSE


# 
usage <- "bam_report.R -b|--bam  bam_file [--htmldir -d  directory] [--fastq_files | -f fastq_file(s)] [-c|--cores num_cores]"
option_list <- list(
  make_option(c("-b", "--bam"), type="character", dest="bam_file", default=NULL,help="BAM file name"),
  make_option(c("-d", "--htmldir"), type="character", dest="htmldir", default="",help="htmldir ([default %default]"),
  make_option(c("-f", "--fastq"), type="character", dest="fastq_files",default="",help="fastq files ([default %default]"),
  make_option(c("-c", "--cores"), type="character",default="3",dest="num_cores",help="Number of cores to use ([default %default])"),
  make_option(c("--debug"),action="store_true",dest="debug",default=FALSE,help="Debug mode")
)

mandatory <- c("bam_file")
filenames <- c("bam_file") ;#filenames that must exist (if defined)
opt <- myParseArgs(usage = usage, option_list=option_list,filenames.exist=filenames,mandatory=mandatory)
# ensure that the path include / in the end
opt$htmldir <- paste(gsub("/$","",opt$htmldir),"/",sep="")
pdebug.enabled <- opt$debug

suppressPackageStartupMessages(library(R2HTML))
#suppressPackageStartupMessages(library(emBAM))

bam.file <- opt$bam_file
# optional
html.dir <- opt$htmldir
# optional (if PE then provide the two)
fastq.files <- opt$fastq_files

tryCatch(num.cores <- as.integer(as.numeric(opt$num_cores)),warning=
         function(w) {
           perror("Invalid number of cores ",opt$num_cores)
           q(status=3)    
       }
)
if (num.cores<1) {
  perror("Invalid number of cores ",opt$num_cores)
  q(status=3)    
}

irap.assert(num.cores>0)

if ( num.cores>parallel:::detectCores()) {
  num.cores <- parallel:::detectCores()
  pwarning("The number of cores to use exceeds the cores available. Reducing the limit to ",parallel:::detectCores())
}

options("cores"=num.cores)

# TODO: check arguments
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol 

plot.b2b.barplot <- function(minus,plus,groups,minus.labels=NULL,bottom.labels=c("-","+"),plus.labels=NULL,cex=0.8,show.var=FALSE) {

  plus <- as.numeric(plus)
  minus <- as.numeric(minus)
  
  par(mar=c(0.5, 5, 0.5, 1))
  plot.new()
  max.val <- max(max(abs(plus)),max(abs(minus))) 

  num.ticks <- 5
  if (is.na(max.val) ) {
    cat("bam_report: Unable to generate b2b plot\n")
    return(NULL)
  }
  ticks <- seq(-max.val, max.val,max.val/num.ticks)
  # this should be done differently
  # if the numbers in the vector are integers then round to 0 decimal numbers
  if (prod(is.wholenumber(minus))==0 && prod(is.wholenumber(plus))==0 ) {
    ticks<-round(ticks,0)
  } else {
    ticks<-round(ticks,2)
  }

  height <- max(3,length(groups))
  y.neg <- height*-0.2
  bar.h <- height/((length(groups)+1))*0.7
  y <- seq(1,length(groups),1)
  
  plot.window(xlim=c(-max.val, max.val), ylim=c(y.neg, height+(y.neg*-1)))

  # vertical line
  lines(rep(0, 2), c(0, height), col="gray")
  #segments(-max.val, y, max.val, y, lty="dotted",col="light grey")
  rect(-minus, y-bar.h/2, 0, y+bar.h/2, col="grey")
  rect(0, y-bar.h/2, plus, y+bar.h/2, col="dark grey")
  if (show.var) {
      variation <- plus-minus
      plusvar <- variation
      plusvar[plusvar<=0] <- 0
      minusvar <- variation
      minusvar[minusvar>=0] <- 0
      
      rect(minusvar, y-bar.h/2, 0, y+bar.h/2, col="light grey")  
      rect(0, y-bar.h/2, plusvar, y+bar.h/2, col="light grey")
  }
  mtext(groups, at=y, adj=1, side=2, las=2,cex=cex)
  par(cex.axis=cex, mex=0.5)
  axis(1, at=ticks, labels=abs(ticks), pos=0)
  
#  tw <- y.neg*max(strwidth(bottom.labels[1]),strwidth(bottom.labels[2]))
#  rect(-tw, -1-h, 0, -1+h, col="dark grey")
#  rect(0, -1-h, tw, -1+h, col="light grey")
  text(0, y.neg, bottom.labels[1], pos=2)
  text(0, y.neg, bottom.labels[2], pos=4)

  # print the labels
  for ( i in c(1:length(plus.labels)) ) {
    text(x=plus[i],y=i,label=plus.labels[i],pos=4,cex=cex)
  }
  for ( i in c(1:length(minus.labels)) ) {
    text(x=-minus[i],y=i,label=minus.labels[i],pos=2,cex=cex)
  }
  # print the variation
  box("inner", col="grey")
}

bam.validate.counts <- function(bam.file) {


}
###############################################

#
bam.readsperseq.plot <- function(bam.readsperseq) {
  
  #par(mfrow=c(1,2))
  colors <- rainbow(nrow(bam.readsperseq))
  plot(bam.readsperseq[,"RefLen"],bam.readsperseq[,"Mapped"],
       xlab="Reference Length",ylab="Reads Mapped",col=colors)
  text(bam.readsperseq[,"RefLen"],bam.readsperseq[,"Mapped"],
       labels=as.character(bam.readsperseq[,"Reference"]),
       col=colors,
       cex=0.7,
       pos=1,
       adj=c(0,-1)
       )
}

############################################################

plot.pe.alignments <- function(df) {

  first <- c(df["Paired 1st mapped",1],df["Paired 1st properly",1])
  second <- c(df["Paired 2nd mapped",1],df["Paired 2nd properly",1])

  first.labels <- c("",paste("~",round(first[2]*100/first[1],0),"%",sep=""))
  second.labels <- c("",paste("~",round(second[2]*100/second[1],0),"%",sep=""))
  plot.b2b.barplot(minus=first,plus=second,groups=c("Mapped","OK"),
                   minus.labels=first.labels,plus.labels=second.labels,
                   bottom.labels=c("1st Mate","2nd Mate"))
}
#
plot.mapped.reads <- function(df) {
  count.label <- "Count"
  count.label.p <- "Count plus"
  count.label.m <- "Count minus"

  data2 <- matrix(as.integer(c(df["Primary",count.label],df["Uniquely mapped reads",count.label],df["Multimap reads",count.label],df["Spliced reads",count.label],df["Unmapped",count.label],
                               df["Primary",count.label.p],df["Uniquely mapped reads",count.label.p],df["Multimap reads",count.label.p],df["Spliced reads",count.label.p],df["Unmapped",count.label.p],
                               df["Primary",count.label.m],df["Uniquely mapped reads",count.label.m],df["Multimap reads",count.label.m],df["Spliced reads",count.label.m],df["Unmapped",count.label.m])),
                  ncol=3,byrow=F)
  colnames(data2) <- c("Both strands","+","-")
  rownames(data2) <- c("Reads Mapped","Uniquely Mapped","Multimaps","Reads Spliced","Unmapped")

  data.pmv <- data2[,2:3]
  Variation<- data.pmv[,1]-data.pmv[,2]
  colnames(data.pmv) <- c("Plus","Minus")

  groups <- names(data.pmv[,1])
  plus  <- data.pmv[,"Plus"]
  minus  <- data.pmv[,"Minus"]

    # print the %
  plus.lab <- c("",paste("~",round(plus["Uniquely Mapped"]*100/plus["Reads Mapped"],0),"%",sep=""),
           paste("~",round(plus["Multimaps"]*100/plus["Reads Mapped"],0),"%",sep=""),
           paste("~",round(plus["Reads Spliced"]*100/plus["Reads Mapped"],0),"%",sep=""),
           paste("~",round(plus["Unmapped"]*100/(plus["Reads Mapped"]+plus["Unmapped"]),0),"%",sep=""))

  minus.lab <- c("",paste("~",round(minus["Uniquely Mapped"]*100/minus["Reads Mapped"],0),"%",sep=""),
           paste("~",round(minus["Multimaps"]*100/minus["Reads Mapped"],0),"%",sep=""),
           paste("~",round(minus["Reads Spliced"]*100/minus["Reads Mapped"],0),"%",sep=""),
           paste("~",round(minus["Unmapped"]*100/(minus["Reads Mapped"]+minus["Unmapped"]),0),"%",sep=""))

  
  plot.b2b.barplot(plus=plus,minus=minus,groups=groups,minus.labels=minus.lab,plus.labels=plus.lab)
}

#
plot.mapped.hits <- function(df) {
    count.label <- "Count"
    df[is.na(df)] <- 0
    data3 <- matrix(as.integer(c(df["Alignments (perfect)",count.label],df["Alignments (1-difference)",count.label],df["Alignments (2-difference)",count.label],
                                 as.integer(df["Alignments",count.label])-sum(as.integer(c(df["Alignments (perfect)",count.label],df["Alignments (1-difference)",count.label],df["Alignments (2-difference)",count.label]))))),
                    ncol=1,byrow=F)
    
    colnames(data3) <- c("Mapping")
    rownames(data3) <- c("Alignments (perfect)","Alignments (1-difference)","Alignments (2-difference)","Alignments (>2-difference)")
    bp<-barplot(data3,beside=FALSE, ylab="Number of alignments" ) 
    cex <- 0.8

    y <- 0
    legend( bp[1], y+data3["Alignments (perfect)",1]/2, paste("0-diff ~",round(data3["Alignments (perfect)",1]*100/sum(data3[,1]),0),"%",sep=""), bty="n", xjust=0.5, yjust=0.5,cex=cex)
    y <- y+data3["Alignments (perfect)",1]
    legend( bp[1], y+data3["Alignments (1-difference)",1]/2 ,paste("1-diff ~",round(data3["Alignments (1-difference)",1]*100/sum(data3[,1]),0),"%",sep=""), bty="n", xjust=0.5, yjust=0.5,cex=cex)
    y <- y+data3["Alignments (1-difference)",1]
    legend( bp[1], y+data3["Alignments (2-difference)",1]/2 ,paste("2-diff ~",round(data3["Alignments (2-difference)",1]*100/sum(data3[,1]),0),"%",sep=""), bty="n", xjust=0.5, yjust=0.5,cex=cex)
    
    y <- y+data3["Alignments (2-difference)",1]
    if ( sum(data3[,1])>0 && round(data3["Alignments (>2-difference)",1]*100/sum(data3[,1]),0)<5 ) {
      y.just=0; y <- y +data3["Alignments (>2-difference)",1]/2
    } else {
      y.just <- 0.5;
    }
    legend( bp[1], y+data3["Alignments (>2-difference)",1]/2 ,paste(">2-diff ~",round(data3["Alignments (>2-difference)",1]*100/sum(data3[,1]),0),"%",sep=""), bty="n", xjust=0.5, yjust=y.just,cex=cex)
    
 }

#
# some (most :() mappers don't include the unaligned reads in the BAM file
# therefore it is necessary to pass the initial number of reads
bam.report <- function(bam.file,html.dir,total.num.reads=NA) {
  bam.flavour <- em.bam.flavour(bam.file)
  
  HTMLInitFile( outdir=html.dir, filename="index", Title="BAM report" )
  HTML( "<h1>BAM Information</h1>" )
  HTML( paste("Filename:", bam.file) )
  HTML( paste("BAM flavour:", bam.flavour) )


  # pprint header information
  # Mapper
  # CL

  HTML( "<h1>Summary</h1>")
  # TODO: uncomment
  bam.sum <- em.bam.sum.stats.table(bam.file)
  HTML(bam.sum,align="center", big.mark='.')
  #plot.bam.sum.stats(bam.sum)
  #

  #
  ################
  # TODO: reenable cache
  df <- em.bam.counts.df(bam.file,use.cache=TRUE,num.reads=total.num.reads)
  #save.image("debug.Rdata")
  count.label <- "Count"
  count.label.p <- "Count plus"
  count.label.m <- "Count minus"
  data <- matrix(as.integer(c(df["All entries",count.label],df["Valid entries",count.label],df["Duplicate",count.label],df["Alignments",count.label],df["Alignments (spliced)",count.label],df["Paired",count.label],df["Paired ok",count.label],df["Primary",count.label])),
                   ncol=1,byrow=F)
  rownames(data) <- c("Entries","Valid","Duplicate","Alignments","Spliced","Paired","Pair mapped","Reads")
  colnames(data) <- c("Summary")


  Unmapped <- df["Unmapped","Count"]
  data2file <- rbind(data,Unmapped)
  #print(data)
  data2file <- cbind(rownames(data2file),data2file[,1])
  colnames(data2file) <- c("Labels",bam.file)
  write.tsv(data2file,paste(bam.file,".tsv",sep=""))
  
  ##################################################################################################
  my.html.plot(filename="alignments_overall.png",html.dir=html.dir, caption="Alignments", to.plot=function() {
    # remove paired?
    data <- data[-nrow(data),]
    barplot(data,beside=T, ylab="Alignments",cex.names=0.8,las = 3)
  });
  
  my.html.plot(filename="mapped_reads_bp.png",html.dir=html.dir, caption="Mapped Reads (Details)", to.plot=function() {
    plot.mapped.reads(df)                 
  });
  my.html.plot(filename="mapped_reads_hits.png",html.dir=html.dir, caption="Mapped Reads (Alignments)", to.plot=function() {
    plot.mapped.hits(df)                 
  });

  # if (PE)
  if (data["Paired",1]>0) {
      HTML( "<h1>PE details</h1>")
      #
      my.html.plot(filename="pe_alignments.png",html.dir=html.dir, caption="Pair-end Alignments", to.plot=function() { plot.pe.alignments(df) })
      # insert distribution
      #l<-em.bam.pe.isize(bam.file,query=em.query("Paired"))
      #if (!is.null(l)) {
      #  ins <- abs(l)
      #  my.html.plot(filename="pe_ins_size.png",html.dir=html.dir, caption="Insert size", to.plot=function() { boxplot(ins,horizontal=FALSE,cex.axis=0.9,ylab="bp",xlab="Estimated Insert Size") } )
      #} else {
      #  pdebug("No isize")
      #}
      # insert length per reference
      #
  }

  ########################################
  HTML( "<h1>Reads per sequence</h1>")
  #bai.file <- em.bam.index.file(bam.file)
  bam.readsperseq <- em.bam.readsperseq.table(bam.file)
  img.file <- paste(html.dir,"/readsperseq_sum.png",sep="")
  png(img.file)#,width=300,height=150)
  bam.readsperseq.plot(bam.readsperseq)
  dev.off()
  HTMLInsertGraph("readsperseq_sum.png", Caption = "Mapped Reads/Reference Length")#,WidthHTML=width)
  ## img.file <- paste(html.dir,"/readsperseq_density.png",sep="")
  ## png(img.file)#,width=300,height=150)
  ## plot(density(log(bam.readsperseq[,"Mapped"])),main="")
  ## dev.off()
  ## HTMLInsertGraph("readsperseq_density.png", Caption = "Mapped Reads Distribution (per reference)",ylab="")#,WidthHTML=width)
  #if (length(bam.readsperseq[,"Mapped"])>2) {
  #  my.html.plot(filename="readsperseq_density.png",html.dir=html.dir, caption="Mapped Reads Distribution (per reference)", to.plot=function() {
  #               plot(density(log(bam.readsperseq[,"Mapped"])),main="")                 
  #             });
  #}
  ## plot the table
  ##
  HTML( "<h3>Alignments details</h3>")  
  HTML(df,align="center", big.mark='.')

  HTMLEndFile()
}


#  my.html.plot(file="readsperseq_density.png",html.dir=html.dir, caption="Mapped Reads Distribution (per reference)", to.plot=function() { plot(density(log,seq(1,300))) } )
my.html.plot <- function(filename=NULL,html.dir=NULL,rel.dir=NULL,caption=NULL,
                           bg="white",
                           width=400,height=400,to.plot=NULL) {

  if (is.null(html.dir)) {
    html.dir <- ""
  }
  if ( is.null(rel.dir)) {
    rel.dir <- ""
  }
                                        # automatic filename
  if ( is.null(filename)) {
    time.label <- format(Sys.time(), "%d%m%Y%H%M%S")
    file <- paste("graph_",time.label,".png",sep="")
  }
  plot.filename <- paste(html.dir,filename,sep="")
  png(filename=plot.filename,width=width, height=height,
      bg=bg)
  to.plot()
  dev.off()
  HTMLInsertGraph(paste(rel.dir,filename, sep=""), Caption = caption)  
}
#################################################################
#
#
bam_report.init <- function(bam.file,html.dir=NULL) {
  
  if (!file.exists(bam.file)) {
    pwarning("Bam file ",bam.file," not found.")
    exit(1)
  }
  
  if (is.null(html.dir) || is.na(html.dir)) {
    html.dir <- paste(bam.file,"_html",sep="")
  }

  html.dir <- paste(sub("/$","",html.dir),"/",sep="")
  if (!file.exists(html.dir)) {
    system(paste("mkdir -p ",html.dir))
    pinfo("Created directory ",html.dir)
  }
  html.dir
}

##################################################
html.dir <- bam_report.init(bam.file,html.dir=html.dir)

if(!is.na(fastq.files) && !is.null(fastq.files) ) {
  files <- strsplit(fastq.files," +")[[1]];
  num.reads <- 0
  for (f in files) {
    # get the number of reads
    #cmd <- paste("num_reads.sh ",f,sep="")
    cmd <- paste("grep -c '^[^ ]'  ",f,sep="")
    r<-system(cmd,intern=T)
    num.reads <- num.reads + as.integer(r)/4
  }
  pinfo("Num reads=",num.reads)
  bam.report(bam.file,html.dir,num.reads)
} else {
  bam.report(bam.file,html.dir)
}
q(status=0)
###################################################
