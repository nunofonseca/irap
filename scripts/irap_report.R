#!/usr/bin/env Rscript
# =========================================================
# Copyright 2012-2013,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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

#
#library(emBAM)
source("/home/nf/Research/myR/src/myUtils/R/myUtils.R")
source("/home/nf/Research/myR/src/emBAM/R/emBAM.R")

args <- commandArgs(trailingOnly=TRUE)
conf_file <- args[1]

################################################################################################
# Useful functions

pinfo <- function(...) {
  cat(paste("[INFO] ",...,"\n",sep=""))
}

pwarning <- function(...) {
  cat(paste("[WARNING] ",...,"\n",sep=""))
}
pmissing <- function(...) {
  cat(paste("[MISSING FILE] ",...,"\n",sep=""))
}

###########################
mytrim <- function(s) {
  sub("^[[:blank:]]+", "",sub("[[:blank:]]+$", "",s))
}

myread.list <- function (s) {
  if(!length(s)) { return(c()) }
  strsplit(mytrim(as.character(s))," +")[[1]]
}

myread.string <- function (s) {
  mytrim(as.character(s))
}

myderef <- function(s) {
  v <- myread.list(s)
  r <- c()
  for ( i in v) {
    #pinfo(i,"=",get(i))
    r <- append(as.character(get(i)),r)
  }
  r
}
################################################################################################
# File functions
file.required <- function(file,fatal=TRUE) {
  if (!file.exists(file)) {
    if (fatal) {
      cat(paste("ERROR: ",file," not found\n",sep=""))
      q()
    }
    pmissing(file)
    return(NULL)
  }
  file
}

#################################################################################################
conf.get.value <- function(conf,val) {
  idx <- grep(paste("^",val,"$",sep=""),conf[,1],ignore.case=TRUE,value=FALSE,perl=TRUE)
  if (is.null(idx)) { return(NULL) }
  return(myread.list(conf[idx,2]))
}

get.libname <- function(s) {
  sub("/[12]$","",s)
}
libname2filename <- function(s) {
  sub("/","_",s)
}
is.pe.lib <- function(s) {
  return(length(grep("/",s))>0)
}
pe.lib2rawfilename <- function(s) {
  id <- sub(".*/","",s)
  lib <- get.libname(s)
  files <- myderef(lib)
  myread.list(files)[as.numeric(id)]
}

# 
path2qual.report <- function(dir.name,file.prefix) {
  root.dir <- dir.name
  #print(root.dir)
  #print(file.prefix)
  #pinfo(">>>>>",dir.name,"----",file.prefix)
  filt.qual.dir <- paste(root.dir,file.prefix,"_fastqc/",sep="")
  filt.qual.zip <- paste(root.dir,file.prefix,".fastqc.zip",sep="")
  filt.qual.index <- paste(filt.qual.dir,"fastqc_report.html",sep="")
  filt.qual.plot.base.quality <- paste(filt.qual.dir,"Images/per_base_quality.png",sep="")
  filt.qual.plot.seq.quality<- paste(filt.qual.dir,"Images/per_sequence_quality.png",sep="")
  #print(filt.qual.dir)
  #print("------------------")
  # check if dir exists
  if (!file.exists(filt.qual.dir)) {
    pinfo("Report directory ",filt.qual.dir," not found. Looking for zip file...")
    if (!file.exists(filt.qual.zip)) {
      pmissing(filt.qual.zip)
      return(NULL)
    }
    # unzip file
    unzip.cmd <- paste("unzip -d ",root.dir," ",filt.qual.zip,sep="")
    system(unzip.cmd)
    if (!file.exists(filt.qual.dir)) {
      pwarning("Unzip file ",filt.qual.zip," apparently failed.")
      return(NULL)
    }
  }
  files <- c(filt.qual.index,filt.qual.plot.base.quality,filt.qual.plot.seq.quality)
  names(files) <- c("index","base.quality","seq.quality")
  #print(files)
  return(files)
}

# bp.data
# colnames=filenames/libname
# rownames=filtering step
# quality.details.matrix(bp.data,"/home/nf/Research/Projects/WIP/EBI1/aehts2/examples/marioni2008/report/riq/raw_data/")
quality.details.matrix <- function(bp.data,dir.name) {
  data <- bp.data
  # extra row
  data<-rbind(data,rep("",ncol(data)))
  rownames(data)[length(rownames(data))]<-"Raw data"
  # another row with the links for the qual plots
  data<-rbind(data,rep("",ncol(data)))
  rownames(data)[length(rownames(data))]<-"Quality plots"
  id <- 1
  #print(bp.data)
  for ( lib in colnames(bp.data) ) {
    lib <- mytrim(lib)
    # Handle the special case (filename) of Paired-end libraries
    if (is.pe.lib(lib)) {
      #pinfo("PE")
      filtered.file <- paste(libname2filename(lib),".f",sep="")
      #pinfo("filtered.file ",filtered.file)
      unfiltered.prefix <- sub("^.*/","",sub(".(fastq|fq).*$","",pe.lib2rawfilename(lib)))
      #unfilt.qual.plot <- file.required(paste(dir.name,unfiltered.prefix,".qual.png",sep=""),fatal=FALSE)
      files.unfilt <- path2qual.report(paste(dir.name,"raw_data/",sep=""),unfiltered.prefix)
      pinfo(files.unfilt)
    } else {
      pinfo("SE")
      filtered.file <- paste(lib,".f",sep="")
      #filt.qual.plot <- file.required(paste(dir.name,lib,".f.qual.png",sep=""),fatal=FALSE)
      #unfilt.qual.plot <- file.required(paste(dir.name,sub(".fastq.*","",myderef(lib)),".qual.png",sep=""),fatal=FALSE)
      #pinfo(myderef(lib),"--->",sub("^.*/","",sub(".(fastq|fq).*$","",myderef(lib))))
      files.unfilt <- path2qual.report(paste(dir.name,"raw_data/",sep=""),sub("^.*/","",sub(".(fastq|fq).*$","",myderef(lib))))
    }
    #
    #trace(path2qual.report)
    #debug(path2qual.report)
    files.filt   <- path2qual.report(paste(dir.name,"raw_data/",sep=""),filtered.file)
    # add a link to the the detailed report
    data["Passed",lib] <- paste(bp.data["Passed",lib],"&nbsp;<a class='qrl' href='",files.filt["index"],"' target='_top'>QR</a>",sep="")
    data["Raw data",lib] <- paste("<a class='qrl' href='",files.unfilt["index"],"' target='_top'>QR</a>",sep="")
    # Add a div to show the plots before and after filtering
    idon  <- paste("qual",id,"on",sep="")
    idoff <- paste("qual",id,"off",sep="")
    width <- 450
    html  <- paste("<a id='",idoff,"' href='javascript:toggle(\"",idon,"\",\"",idoff,"\");'>Quality plots</a>",
                  "<div id='",idon,"' style='display: none'><table class='plots' style='text-align: center'><tr><td>Raw data</td><td>Filtered</td></tr>",
                   "<tr><td><img  width='",width,"' src='",files.unfilt["seq.quality"],"' alt='Quality plot'></td><td><img  width='",width,"' src='",files.filt["seq.quality"],"' alt='Quality plot'></td></tr>",
                   "<tr><td><img width='",width,"' src='",files.unfilt["base.quality"],"' alt='Per base quality plot'></td><td><img width='",width,"' src='",files.filt["base.quality"],"' alt='Per base quality plot'></td></tr>",
                   "</table></div>",sep="")
    data["Quality plots",lib] <- html
    
    id <- id+1
  }
  data
}

########################################################3
# Report
# 1)Table with the number of reads filtered per file
# 2)Plot (showing all files)
# 3)Quality plots (base and sequence level) for each file
# 4)Link to the detailed report generated by fastqc
# return a matrix with the data
filtering.stats <- function(pe.v,se.v) {

  pinfo("Generating quality report")
  ######################################
  # reads info
  df <- raw.summary.dataframe(pe.v,se.v)
  ########################
  # Plot 
  legend <- c("Passed", "Pairing", "Ns", "Contamination","Quality")
  cols <- rainbow(length(legend))
  bp.data <- t(as.matrix(df[,c("Reads4","Reads3","Reads2","Reads1","Reads0")]))
  colnames(bp.data) <- as.character(df[,c("Lib")])
  rownames(bp.data) <- legend
# TODO: reduce font with the number of cols
  n <- length(colnames(bp.data))
  cex <- NULL
  if ( n>10 ) { cex <- 0.6 } else { cex=0.5}
  
  file <- "read_filtering_plot.png"
  filepath <- paste(report_dir,file,sep="/")
  
  width=(480*round(0.5+n/6/2,0))
  png(filepath,width=width,height=480)
                                        #bp.data
  bp <- barplot( bp.data, beside = FALSE, col=cols, ylab="Number of reads", las=2 , cex.axis=cex, cex.names=0.8)
                                        #par( xpd=NA )
  par( xpd=NA )
  
  legend( "bottomright", legend=legend, fill=cols,horiz=TRUE, cex=0.7)
                                        # print the %
  for( i in c(1:length(colnames(bp.data))) ) {
    total <- sum(bp.data[,i])
    f <- c()
    y <- c()
    f[1] <- round(bp.data[1,i]*100/total)
    y[1] <- bp.data[1,i]
    f[2] <- round(bp.data[2,i]*100/total)
    y[2] <- sum(bp.data[c(1:2),i])  
    f[3] <- round(bp.data[3,i]*100/total)
    y[3] <- sum(bp.data[c(1:3),i])
    f[4] <- round(bp.data[4,i]*100/total)
    y[4] <- sum(bp.data[c(1:4),i])
    f[5] <- round(bp.data[5,i]*100/total)
    y[5] <- sum(bp.data[c(1:5),i])
    for( r in c(1:5) ){
      if (f[r]>0.01) {
        #pinfo(f[r])
        legend( bp[i], y[r]-bp.data[r,i]/2, paste("~",f[r], "%", sep=""), bty="n", xjust=0.5, yjust=0.5)      
      }
    }
  }
  dev.off()

  HTML("<H2>Filtering Report</H2>")
  HTMLInsertGraph(file, Caption = "Quality Filtering",WidthHTML=width)

  #
  #
  #df
  
  dir.name <- paste(name,"/report/riq/",sep="")
  # change table to have links to the detailed quality reports
  tab <- quality.details.matrix(bp.data,dir.name)
  #print(tab)
  # Change full pathnames to relative 
  tab <- gsub(paste(name,"/report/",sep=""),"",tab)
  HTML(tab,align="center", big.mark='.')
}

totnumreads2filtered <- function(v) {
  v2 <- v
  v2[1,2] <- v[1,2]-v[1,3]
  v2[1,3] <- v[1,3]-v[1,4]
  v2[1,4] <- v[1,4]-v[1,5]
  v2[1,5] <- v[1,5]-v[1,6]
  v2
}
# stats
# lib, file, csv data
raw.summary.dataframe <- function(pe,se) {  
  filt_stats <- matrix(ncol=11,nrow=0, dimnames = list(NULL,c("Lib","File", "length","qual","ins","sd", "Reads0","Reads1","Reads2","Reads3","Reads4")))
  df <- data.frame(filt_stats)
  i <- 1
  for (lib in pe ) {
    linfo <- as.numeric(c(myderef(paste(lib,"_rs",sep="")),myderef(paste(lib,"_qual",sep="")),myderef(paste(lib,"_ins",sep="")),myderef(paste(lib,"_sd",sep=""))))
    
    files <- myread.list(myderef(lib))
                                        #libfiles[lib] <- files
    files <- sub("(.gz|.bzip2)","",sub("(.fastq|.fq)",".f",files,ignore.case =TRUE),ignore.case=TRUE,perl=TRUE)
                                        #filteredlibfiles[lib] <- files
    ## sing.file <- sub("_1.f.fastq","_1.sing.fastq",files[1])
    ## d <- read.csv(paste(name,"/report/riq/",sing.file,".csv",sep=""),header=FALSE)
    ## df[i,] <- append(append(c(lib,sing.file),linfo),d[1,2:6])
    ## i <- i+1
    
    # use / because some files use a _ or other characters in the filename
    d <- totnumreads2filtered(read.csv(file.required(paste(name,"/report/riq/",files[1],".csv",sep="")),header=FALSE))
    df[i,] <- append(append(c(paste(lib,"/1",sep=""),files[1]),linfo),d[1,2:6])
    i <- i+1
    d <- totnumreads2filtered(read.csv(file.required(paste(name,"/report/riq/",files[2],".csv",sep="")),header=FALSE))
    df[i,] <- append(append(c(paste(lib,"/2",sep=""),files[2]),linfo),d[1,2:6])
    i <- i+1
  }
  for (lib in se ) {
    linfo <- as.numeric(c(myderef(paste(lib,"_rs",sep="")),myderef(paste(lib,"_qual",sep="")),0,0))
    files <- strsplit(as.character(myderef(lib))," ")[[1]]
                                        #libfiles[lib] <- files
    files <- sub("(.gz|.bzip2)","",sub("(.fastq|.fq)",".f",files,ignore.case =TRUE),ignore.case=TRUE,perl=TRUE)
                                        #filteredlibfiles[lib] <- files
    d <- totnumreads2filtered(read.csv(file.required(paste(name,"/report/riq/",files[1],".csv",sep="")),header=FALSE))

    df[i,] <- append(append(c(lib,files[1]),linfo),d[1,2:6])
    i <- i+1
  }
  df
}

alignments.report <- function(conf,report.dir) {

  mapper  <- conf.get.value(conf,"mapper")
  dir.prefix <- paste("bam",mapper,sep="/")
  libs.pe <- conf.get.value(conf,"pe")
  libs.se <- conf.get.value(conf,"se")
  libs <- unique(append(libs.pe,libs.se))
  nlibs <- length(libs)
  # Plot for all libraries
  if (length(libs.pe)>0) { libs.pe <- paste("<a href='",dir.prefix,"/",libs.pe,"/index.html'>Align. Details</a>",sep="") }
  if (length(libs.se)>0) { libs.se <- paste("<a href='",dir.prefix,"/",libs.se,"/index.html'>Align. Details</a>",sep="") }
  urls <- append(libs.pe,libs.se)
  m <- matrix(c(urls),nrow=1,byrow=T)
  colnames(m) <- libs
  HTML("<H2>Alignments</H2>")

  overall.plot <- paste(dir.prefix,"/align_overall_comparison.png",sep="")
  html <- paste("<img src='",overall.plot,"' alt=''>",sep="")
  HTML(html)

  ###############
  # Overall table
  #                          lib1   lib2     ... libk
  # Reads Mapped           value (%)                      | average (%)
  # Reads not mapped       value (%) ...                  |
  # Reads uniquely mapped
  # Reads spiced
  #  
  bam.dir <- paste(name,mapper,sep="/")
  bam.files <- c()
  pe <- conf.get.value(conf,"pe")
  se <- conf.get.value(conf,"se")
  if (! is.null(pe)) { bam.files <- paste(pe,".pe.hits.bam",sep="") }
  if (! is.null(se)) { bam.files <- append(bam.files,paste(se,".se.hits.bam",sep=""))}  
  bam.files <- paste(bam.dir,"/",bam.files,sep="")
  
  data <- lapply(bam.files,em.bam.counts.df,TRUE)
  df <- as.data.frame(data[[1]][,1])
  rownames(df)<-rownames(data[[1]])
  if(length(data)>1) {
    for (i in c(2:length(data))) {
      df <- cbind(df,data[[i]]$Count)
    }
  }
  colnames(df) <- gsub(".*/","",gsub("(\\.pe|\\.se).*","",bam.files))

  n.reads <- df[c("Unmapped","Primary","Uniquely mapped reads","Multimap reads","Spliced reads"),]
  table.rnames <- c("Unmapped","Mapped","Uniquely mapped","Multimaps","Spliced")
  rownames(n.reads) <- table.rnames

  tot.reads <- apply(n.reads[c(1,2),],MARGIN=2,FUN=sum)

  #%
  p.reads <- n.reads
  p.reads[1,] <- round(p.reads[1,]/tot.reads*100,2)
  p.reads[2,] <- round(p.reads[2,]/tot.reads*100,2)
  p.reads["Uniquely mapped",] <- round(n.reads["Uniquely mapped",]/as.numeric(n.reads["Mapped",])*100,2)
  p.reads["Multimaps",] <- round(n.reads["Multimaps",]/n.reads["Mapped",]*100,2)
  p.reads["Spliced",] <- round(n.reads["Spliced",]/n.reads["Mapped",]*100,2)
  # average
  p.reads<-cbind(p.reads,apply(p.reads,MARGIN=1,FUN=median))
  colnames(p.reads) <- append(colnames(n.reads),"Median")

  p.reads <- rbind(p.reads,append(as.vector(m),NA))
  n.reads <- rbind(n.reads,append(as.vector(m),NA))
  rownames(p.reads) <- append(table.rnames," ")
  rownames(n.reads) <- append(table.rnames," ")

  # add a link for details
  HTML(data.frame(p.reads),align="center", big.mark='.',caption="Percentage of reads")
  HTML(data.frame(n.reads),align="center", big.mark='.')
  
  #############################3
  #
  #HTML("<H3>Alignments: detailed reports</H2>")
  # Table with the library names and links to the detailed reports
  #HTML(data.frame(m),align="center", big.mark='.',row.names=FALSE)

}

###################################################
# WIP
ge.report <- function(conf,report.dir) {
  assembler  <- conf.get.value(conf,"assembler")
  mapper  <- conf.get.value(conf,"mapper")
  name <- conf.get.value(conf,"name")
  gene_class_file<-tools:::file_path_as_absolute(paste(name,"data/gene_class.txt",sep="/"))

  HTML("<H2>GE</H2>")
  if ( assembler=="cufflinks1" || assembler=="cufflinks2" ) {
    dir <- paste(name,mapper,assembler,sep="/")
    reportdir <- paste(name,"report",mapper,assembler,sep="/")
    system(paste("mkdir -p ",reportdir,sep=""))
    system(paste("cd ",dir,"; cufflinks_plots.R ",gene_class_file,sep=""))
    system(paste("cp ",dir,"/ge.html ",reportdir,"/index.html",sep=""))
    system(paste("cp ",dir,"/*.png ",reportdir,sep=""))
    HTML(paste("<a href='",mapper,"/",assembler,"/index.html'>",assembler," report</a>",sep=""))
  }
  # DE SEQ
}

#################################################################################################
if ( is.na(conf_file) ) {
  cat("ERROR: Missing conf file (arg1)\n")
  q()
}
qual_filtering <- "yes"
conditions <- NULL
contrasts <- NULL
pe <-
se <- 

pinfo("conf_file=",conf_file)

pinfo("Loading conf_file...")
conf.table <- read.delim(conf_file,sep="=",comment.char="#",header=FALSE,stringsAsFactors=FALSE)
pinfo("Loading conf_file...done")
pinfo("Configuration:")
for ( i in 1:nrow(conf.table) ) {
  assign(as.character(conf.table[i,1]),mytrim(conf.table[i,2]))
  pinfo("          ",conf.table[i,1],"=",conf.table[i,2])
}
#################
# Default values
# TODO: fix this... the default values should come from irap
if ( sum(as.character(conf.table[,1])=="mapper")==0 ) {
  entry <- c("mapper","tophat1")
  conf.table[nrow(conf.table)+1,] <- entry
}
if ( sum(as.character(conf.table[,1])=="assembler")==0 ) {
  entry <- c("assembler","cufflinks1")
  conf.table[nrow(conf.table)+1,] <- entry
}

rownames(conf.table) <- conf.table[,1]
############################################333
#print(conf.table)
report_dir <- paste(name,"/report/",sep="")
css_file <- "irap.css"

###########################
# Check if directory exists
file.required(report_dir)

#####################
#          1
# General information
#
library(R2HTML)

HTMLInitFile( outdir=report_dir, filename="index", Title=paste("Report ",name,sep=""), CSSFile=css_file, useGrid=FALSE)
# Javascript code
HTML( '<script language="javascript"> 
function toggle(showHideDiv, switchTextDiv) {
	var ele = document.getElementById(showHideDiv);
	var text = document.getElementById(switchTextDiv);
	if(ele.style.display == "block") {
    		ele.style.display = "none";
		text.innerHTML = "show";
  	}
	else {
		ele.style.display = "block";
		text.innerHTML = "hide";
	}
} 
</script>
')

HTML( "<h1>Project Information</h1>" )
HTML( paste("Project name:", name) )
HTML( paste("Species:", species) )

if ( is.null(contrasts) ) {
  # default
  if (!is.null(conditions)) {
      contrasts <- "conditions"
  } else {
    contrasts <- NULL
  }
}

if (is.null(contrasts)) {
  HTML( paste("Contrasts:", "-") )
  has.contrasts <- 0
} else {
  HTML( paste("Contrasts:", contrasts) )
  contrasts.v <- myread.list(contrasts)
  has.contrasts <- 1
}

if (is.null(qual_filtering) || qual_filtering!="yes") {
  HTML( paste("Reads quality filtering:", "Disabled") )
  qual_filt <- 0
} else {
  HTML( paste("Reads quality filtering:", "Enabled") )
  qual_filt <- 1
}
######
HTML( paste("<h2>Data Summary</h2>") )

if (! is.null(pe) ) {
  HTML( paste("Pair-end libraries:",pe) )
  pe.v <- myread.list(pe)
} else {
  pe.v <- c()
}
if (! is.null(se) ) {
  HTML( paste("Single-end libraries:",se) )
  se.v <- myread.list(se)
} else {
  se.v <- c()
}

# Summary
# if has.cond then show the data organized by condition
# otherwise show the data by libraries (se, pe)
# generate the table
# HTML(iris, file=tmpfic)
# load the files from [name]/report/riq/[lib].f.csv
lib2cond <- list()
lib2type <- list()
libfiles <- list()
filedata <- list()
libdata  <- list()
filteredlibfiles <- list()
#cond.v <- myread.list(conditions)
if (has.contrasts) {
  for (c in contrasts.v) {
    pinfo(c)
    contrast.v <- myread.list(c)
    for ( l in contrast.v) {
      for (lib in myderef(c) ) {
        lib2cond[lib] <- c
      }
    }
  }
}

##############################################################################
# Step 1 - Filtering

for (lib in pe.v ) {
    lib2type[lib] <- "PE"
}
for (lib in se.v ) {
    lib2type[lib] <- "SE"
}

pinfo(qual_filt)
if (!qual_filt) {

} else {
  pinfo("QUAL filtering")
  t <- filtering.stats(pe.v,se.v)
  print(t)
}

#############################################################################################
# Alignment stats
#
alignments.report(conf.table,report_dir)

###########################################
# Expression stats
ge.report(conf.table,report_dir)


#############################################################################################
#print(lib2cond)
#print(libfiles)
#print(filedata)
#print(libdata)
    

HTMLEndFile()
q()

