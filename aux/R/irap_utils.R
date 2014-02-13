# =========================================================
# Copyright 2012-2013,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
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
#    $Id: aux/R/irap_utils.R 0.1.1 Nuno Fonseca Fri Dec 21 01:07:37 2012$
# =========================================================
# TODO: create an object for the conf. information

irap_version<-"0.3.3"

IRAP.DIR <- Sys.getenv(c("IRAP_DIR"))
if ( IRAP.DIR == "" ) {
  cat("ERROR: environment variable IRAP_DIR is not set\n")
  q(status=1)
}

# PATH to the CSS file used in all HTML reports
irap_css_file <- paste(IRAP.DIR,"/aux/css/irap.css",sep="")
irap_menu_css_file <- paste(IRAP.DIR,"/aux/css/menu.css",sep="")
irap.css <- "irap.css"
############################################################
# Useful functions
pinfo <- function(...) {
  cat(paste("[INFO] ",...,"\n",sep=""))
}

pwarning <- function(...) {
  cat(paste("[WARNING] ",...,"\n",sep=""),file=stderr())
}
pmissing <- function(...) {
  cat(paste("[MISSING FILE] ",...,"\n",sep=""))
}

perror <- function(...) {
  cat(paste("[ERROR] ",...,"\n",sep=""),file=stderr())
}

##################################################
# debugging
pdebug.enabled <- FALSE
pdebug.stage <- NULL

pdebug <- function(...) {
  if ( pdebug.enabled) {
    cat(paste("[DEBUG] ",...,"\n",sep=""),file=stderr())
  }
}

pdebug.rdata <- function() {
# reset some variables in order to properly debug an Rdata object data
  IRAP.DIR <- Sys.getenv(c("IRAP_DIR"))
  if ( IRAP.DIR == "" ) {
    IRAP.DIR="/home/nf/Research/Projects/WIP/IRAP/irap_install"
    Sys.setenv(IRAP_DIR=IRAP.DIR)
  }
  assign("IRAP.DIR",Sys.getenv(c("IRAP_DIR")),envir = .GlobalEnv)
  pinfo("debug stage: ",pdebug.stage)
}

pdebug.enable <- function() {
  assign("pdebug.enabled",TRUE,envir = .GlobalEnv)
}
pdebug.disable <- function() {
  assign("pdebug.enabled",FALSE,envir = .GlobalEnv)
}

pdebug.status <- function(status=NULL) {
  if (!is.null(status)) {
    if ( status==TRUE ) {
      pdebug.enable()
    } else {
      pdebug.disable()
    }
  } 
  pdebug.enabled
}

pdebug.save.state <- function(file,stage="undef") {  
  if (pdebug.enabled) {
    assign("pdebug.stage",stage,envir = .GlobalEnv)
    save.image(file=paste(file,".Rdata",sep=""))
  }
}
irap.assert <- function(assertion,...) {
  if ( assertion==FALSE ) {
    perror("Internal error: Assertion failed!")
    perror(...)
    save.image("debug.Rdata")
    q(status=3)
  }
}
######################################################
# Formats
formats.cols <- list(
  gff3=c("seqid","source","type","start","end","score","strand","phase","attributes","comments"),
  gtf=c("seqid","source","feature","start","end","score","strand","frame","attributes")
  )

# TODO: improve error handling
load.gtf <- function(gtf.file,feature=NULL,selected.attr=NULL) {

  gtf<-read.table(gtf.file,sep="\t",header=F,quote="\"")
  cnames <- formats.cols$gtf
  colnames(gtf)<-cnames[0:ncol(gtf)]
  #feature <- "CDS"
  #selected.attr <- c("gene_id","transcript_id")
  if (!is.null(feature)) {
    gtf<- gtf[gtf$feature==feature,]
  }
  gtf$attributes <- as.character(gtf$attributes)
  gtf.attributes.names<-c("gene_id","transcript_id","exon_number","gene_name","gene_biotype","transcript_name","protein_id")
  if ( !is.null(selected.attr) ) {
    gtf.attributes.names<- gtf.attributes.names[gtf.attributes.names %in% selected.attr]
  }
  for (tag in gtf.attributes.names ) {
    pattern <- paste("",tag," ([^;]+)",sep="")
    m<-regexec(pattern,as.vector(gtf$attributes))
    for ( i in c(1:length(m)) ) { if ( m[[i]][1]==-1 ) { m[[i]]=NA; } }
    x<-regmatches(as.vector(gtf$attributes),m)
    for ( i in c(1:length(x)) ) { if ( length(x[[i]])==0 ) { x[[i]]=c(NA,NA); } }
    vals<-matrix(unlist(x),ncol=2,byrow=T,)[,2]
    gtf[,tag]<-vals
  }
  gtf
}

# Given a gtf file returns a vector with the length of the genes
get.gene.length.from.gtf.file <- function(gtf.file,filter.biotype=NULL,length.mode="union.exons") {
  gtf <- load.gtf(gtf.file)
  get.gene.length.from.gtf(gtf,filter.biotype,length.mode)
}
#
# Given a matrix obtained from a gtf file returns a vector with the length of the genes
get.gene.length.from.gtf <- function(gtf,filter.biotype=NULL,length.mode="union.exons") {
  # TODO: validate gtf
  # protein coding
  if ( !is.null(filter.biotype) ) {
    gtf <- gtf[gtf$gene_biotype==filter.biotype,]
  }
  # compute the length for each exon
  gtf <- gtf[gtf$feature=="exon",]
  gene.ids <- unique(gtf$gene_id)
  glen <- unlist(lapply(gene.ids,get.gene.length,gtf,mode=length.mode))
  names(glen) <- gene.ids
  glen  
}


get.gene.length <- function(gene.id,gtf.data,mode="sum.exons",lim=+Inf,do.plot=FALSE) {
  library("intervals")
  #gtf.data <- gtf2
  i <- Intervals(gtf.data[gtf.data$gene_id==gene.id  & gtf.data$feature=="exon",c("start","end")])
  ne <- length(size(i))
  #
  res <- c()
  if ( mode != "sum.exons" ) {
    # a) reduce -> for gene length
    i <- reduce(i)
    si <- size(i)
    res <- sum(si[si<lim])+length(si)
  } else {
    # b) interval_union
    # size(i)    
    si <- size(i)    
    res <- sum(si[si<lim])+ne 
  }
  if ( do.plot ) 
    plot(i)
  res
}

# Given a matrix obtained from a gtf file returns a vector with the length of the transcripts
get.transcript.length.from.gtf <- function(gtf,filter.biotype=NULL) {
  # TODO: validate gtf
  # protein coding
  if ( !is.null(filter.biotype) ) {
    gtf <- gtf[gtf$gene_biotype==filter.biotype,]
  }
  # compute the length for each transcript
  gtf <- gtf[gtf$feature=="exon",]
  transcript.ids <- unique(gtf$transcript_id)
  tlen <- unlist(lapply(transcript.ids,get.transcript.length,gtf))
  names(tlen) <- transcript.ids
  tlen  
}

#
get.transcript.length <-  function(transcript.id,gtf.data) {
  library("intervals")
  i <- Intervals(gtf.data[gtf.data$transcript_id==transcript.id & gtf.data$feature=="exon",c("start","end")])
  sum(size(i))+length(i)/2
}

# Given a matrix obtained from a gtf file returns a matrix with the length of the  exons of a given gene/transcript
get.exon.length.from.gtf <- function(gtf,filter.biotype=NULL) {
  # TODO: validate gtf
  # protein coding
  if ( !is.null(filter.biotype) ) {
    gtf <- gtf[gtf$gene_biotype==filter.biotype,]
  }
  # compute the length for each exon
  gtf <- gtf[gtf$feature=="exon",]  
  elen <- abs(gtf$start-gtf$end)+1
  gtf$elength <- elen
  gtf[,c("gene_id","transcript_id","exon_number","elength")]
}


load.gff3 <- function(file,type="gene") {
  # load the gff3 file
  gff3<-read.table(file,sep="\t",header=F,quote="\"")
  cnames <- formats.cols$gff3
  colnames(gff3)<-cnames[0:ncol(gff3)]
  gff3$attributes <- as.character(gff3$attributes)

  if ( !is.na(type) ) {
    gff3 <- gff3[gff3$type==type,]
  }
  # convert the start/end to numbers
  gff3$start <- as.integer(gff3$start)
  gff3$end <- as.integer(gff3$end)

  # add the attributes
  tags<-c("ID","Name","Alias","Parent","Target","Gap","Derives","Note","DBxref","Ontology_term","Is_circular")
  pdebug("tags")
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

  gff3
}

######################################################
# RPKMs
# deprecated
counts2RPKMs <- function(count.matrix,annot.table=NULL) {

  if (is.null(annot.table)) {
    rpkms <- count.matrix*NA
  } else {
    len <- annot.table$len
    names(len) <- annot.table$ID
    len <- len[rownames(count.matrix)]
    rpkm <- function(counts,len) {
      counts / (len/1e3) / (sum(counts)/1e6)
    }
    rpkms <- round(apply(count.matrix,2,rpkm,len),2)
  }
  rpkms
}

# RPKMs
countstable2rpkms <- function(table,lens,exitonerror=TRUE) {
  # check if there missing features
  missing.feat <- (!rownames(table) %in% names(lens))

  if ( sum(missing.feat) ) {
    perror("Length of ",paste(rownames(table)[missing.feat],sep=",")," not found.")
    if (exitonerror) { q(status=1) }
    return(NULL)
  }
  v.compute.rpkm <-  function(l,lens) {
    #(l*1e6)/(sum(l)*lens[names(l)]/1000)
    10^9*l/(sum(l)*lens[names(l)])
  }
  round(apply(table,2,v.compute.rpkm,lens),2)
}


######################################################
# Annot
# extract a named vector of all terms
goID2term <- function(goid) {
  suppressPackageStartupMessages(library("GO.db"))
  Term(goid)
}

species2db.m<-matrix(c('org.Ag.eg.db','Anopheles',
'org.At.tair.db','Arabidopsis',
'org.Bt.eg.db','Bovine',
'org.Ce.eg.db','Worm',
'org.Cf.eg.db','Canine',
'org.Dm.eg.db','Fly',
'org.Dr.eg.db','Zebrafish',
'org.EcK12.eg.db','E coli strain K12',
'org.Gg.eg.db','Chicken',
'org.Hs.eg.db','Human',
'org.Mm.eg.db','Mouse',
'org.Mmu.eg.db','Rhesus',
'org.Pf.plasmo.db','Malaria',
'org.Pt.eg.db','Chimp',
'org.Rn.eg.db','Rat',
'org.Sc.sgd.db','Yeast',
'org.Sco.eg.db','Streptomyces coelicolor',
'org.Ss.eg.db','Pig',
'org.Tgondii.eg.db','Toxoplasma gondii',
'org.Xl.eg.db','Xenopus'),byrow=T,ncol=2)
colnames(species2db.m)<-c("db","species")

# Collect gene annotation from AnnotationDbi
#  o Consistent across analyses
#  o Reliably accessible; obtain with biocLite
#  o Careful secondary curation at Bioconductor
#

species2dbs <- function(species.name) {
  suppressPackageStartupMessages(library("GO.db"))

  if (is.na(species.name) || is.null(species.name)) {
    return(NA);
  }
  # lower case
  # normalize species_name (remove trailing .*,_-> )
  n.species <- tolower(gsub("_"," ",gsub("\\.[0-9]+","",species.name)))
  pdebug("nspecies=",n.species)
  go.db <- NA;
  pfam.db <- NA;
  lgn.db <- NA; # long gene name summary(org.Hs.egGENENAME)
  symbol.db <- NA; # to map name to iid
  ens.db <- NA
  kegg.db <- NA
  # refseq.dn <- NA; #
  # PA  pathways
  # Kegg
  # human (support alias)
  if ( regexpr("^human",n.species,ignore.case=T,perl=T)!=-1 ||
       regexpr("^homo.*",n.species,ignore.case=T,perl=T)!=-1 ||
      n.species=="hs" ) {
    species.name <- "Homo Sapiens"
    library("org.Hs.eg.db")
    go.db <- org.Hs.egGO
    pfam.db <- org.Hs.egPFAM
    lgn.db <- org.Hs.egGENENAME
    symbol.db <- org.Hs.egSYMBOL
    ensembl.db <- org.Hs.egENSEMBL
    kegg.db <- org.Hs.egPATH
  } else if ( regexpr("^mus.*",n.species,ignore.case=T,perl=T)!=-1 ||
       regexpr("^mouse",n.species,ignore.case=T,perl=T)!=-1 ||
       n.species == "mus musculus" )
    {
      species <- "mus musculus"
      library("org.Mm.eg.db")
      go.db <- org.Mm.egGO
      pfam.db <- org.Mm.egPFAM
      lgn.db <- org.Mm.egGENENAME
      symbol.db <- org.Mm.egSYMBOL
      ensembl.db <- org.Mm.egENSEMBL
      kegg.db <- org.Mm.egPATH
    } else  if (
      regexpr("^ecoli.*",n.species,ignore.case=T,perl=T)!=-1 ||
      regexpr("^e.coli.*",n.species,ignore.case=T,perl=T)!=-1 ||      
      n.species == "ecK" ) {
      species <- "EcK12"
      library("org.EcK12.eg.db")
      go.db <- org.EcK12.egGO
      # how to get the entrez ids from the ensembl gene id?
      pfam.db <- NA
      lgn.db <- org.EcK12.egGENENAME
      symbol.db <- org.EcK12.egSYMBOL
      ensembl.db <- NA
      kegg.db <- org.EcK12.egPATH
    } else  if (
      regexpr("^chicken.*",n.species,ignore.case=T,perl=T)!=-1 ||
      regexpr("^Gg.*",n.species,ignore.case=T,perl=T)!=-1 
      ) {
      species <- "Chicken"
      library('org.Gg.eg.db')
      go.db <- org.Gg.egGO
      pfam.db <- org.Gg.egPFAM
      lgn.db <- org.Gg.egGENENAME
      symbol.db <- org.Gg.egSYMBOL
      ensembl.db <- org.Gg.egENSEMBL
      kegg.db <- org.Gg.egPATH
    } else  if (
      regexpr("^pig.*",n.species,ignore.case=T,perl=T)!=-1 ||
      regexpr("^Ss.*",n.species,ignore.case=T,perl=T)!=-1 
      ) {
      species <- "Pig"
      library('org.Ss.eg.db')
      go.db <- org.Ss.egGO
      pfam.db <- NA
      lgn.db <- org.Ss.egGENENAME
      symbol.db <- org.Ss.egSYMBOL
      ensembl.db <- NA
      kegg.db <- org.Ss.egPATH
    } else {
      pwarning("species2dbs: Unable to get annotation for species ",n.species,".")
      return(NA);
    }
  # other species
  # Add as needed
  return(list("n.species"=n.species,
              "species"=species,
              "go.db"=go.db,
              "pfam.db"=pfam.db,
              "lgn.db"= lgn.db,
              "symbol.db"=symbol.db,
              "kegg.db"=kegg.db,
              "ens.db"=ensembl.db
              ))
}

annot.expand.fields <- function(annot.table) {
  chr <- unlist(strsplit(x=as.character(annot.table$locus),split=":.*"))
  get.gene.start <- function(s) {
    l <- unlist(strsplit(x=s,split=":",fixed=T))
    loc <- unlist(strsplit(x=l[2],split="..",fixed=T))
    as.numeric(loc[1])
  }
  get.gene.len <- function(s) {
    l <- unlist(strsplit(x=s,split=":",fixed=T))
    loc <- unlist(strsplit(x=l[2],split="..",fixed=T))
    abs(as.numeric(loc[1])-as.numeric(loc[2]))
  }
  # expand annot info
  # TODO: add length of longest transcript
  len <- sapply(as.character(annot.table$locus),get.gene.len)
  annot.table$len<-len
  annot.table$chr<-chr
  annot.table$start<-sapply(as.character(annot.table$locus),get.gene.start)
  return(annot.table)
}
# entrez gene identifier
annot.get.egi <- function(gid,dbs) {
  library("GO.db")
  db <- dbs$symbol.db
  egi <- revmap(db)[[gid]]
  #pdebug(gid,"---",egi)
  if (is.null(egi)) {
    db <- dbs$ens.db
    egi<-revmap(db)[[gid]]
  }
  #pick the first  (ensembl mapping may return mul. entries)
  if (is.null(egi) || is.na(egi)) {
    return(NA)
  }
  return(egi[1]) 
}
#egi - entrez gene identifier
annot.get.go <- function(egi,dbs) {
  if (is.null(egi)||is.na(egi)) { return(NA) }
  db <- dbs$go.db
  go <- db[[egi]]
  if (is.null(go)||is.na(go)) {
    return(NA);
  }
  l <- unlist(go)  
  go <- paste(gsub("\"","",grep("^GO:",l,value=T)),collapse=";")
  #print(go)
  return(go)
}
annot.get.go.term <- function(egi,dbs) {
  library("GO.db")

  if (is.null(egi)||is.na(egi)) { return(NA) }
  db <- dbs$go.db
  go <- db[[egi]]
  if (is.null(go)||is.na(go)) {
    return(NA);
  }
  l <- unlist(go)  
  go <- gsub("\"","",grep("^GO:",l,value=T))
  terms <- paste(goID2term(go),collapse=";")
  return(terms)
}

annot.get.pfam <- function(egi,dbs) {

  if (is.na(egi)||is.na(egi)) { return(NA) }
  db <- dbs$pfam.db
  val <- db[[egi]]
  if (is.null(val)||is.na(val)) {
    return(NA);
  }
  val <- unlist(val)
  return(paste(val,collapse=";"))
}
annot.get.lname <- function(egi,dbs) {

  if (is.na(egi)||is.na(egi)) { return(NA) }
  db <- dbs$lgn.db
  val <- db[[egi]]
  if (is.null(val)||is.na(val)) {
    return(NA);
  }
  return(val)
}
annot.get.kegg <- function(egi,dbs) {

  if (is.na(egi)||is.na(egi)) { return(NA) }
  db <- dbs$kegg.db
  val <- db[[egi]]
  if (is.null(val)||is.na(val)) {
    return(NA);
  }
  return(paste(unlist(val),collapse=";"))
}
#org.Hs.egPATH
#Each KEGG pathway has a name and identifier

######################################################
# HTML Templates

irap.ctrs.Figure <- 0
irap.ctrs.Table  <- 0

reset.irap.ctr <- function() {
  assign("irap.ctrs.Table", 0, env=globalenv())
  assign("irap.ctrs.Figure", 0, env=globalenv())
  return(TRUE)
}

irap.ctr <- function(c.name,add.label=TRUE) {
  if ( !c.name %in% c("Figure","Table") ) {
    return(NULL);
  }
  if ( c.name == "Figure" ) {
    x <- irap.ctrs.Figure
    var <- "irap.ctrs.Figure"
    lab <- "Figure"
  } else {
    x <- irap.ctrs.Table
    var <- "irap.ctrs.Table"
    lab <- "Table"
  }
  x <- x+1
  assign(var, x, env=globalenv())
  if (add.label) return(paste(lab," ",x,sep=""))
  return(x)
}


irap.output.html <- function(html,info.msg=NULL) {
  if ( !is.null(info.msg) && info.msg!="") {
    cat(irap.info2html(info.msg))
  }
  cat(html)
}
    

# Return a message in HTML format
irap.info2html.msgs.counts <- 0
irap.info2html <- function(msg) {
  assign("irap.info2html.msgs.counts",irap.info2html.msgs.counts+1,envir = .GlobalEnv)
  return(paste("<div class=\"info_wrapper\"><a href=\"javascript:void(0)\" onclick=\"document.getElementById('light",irap.info2html.msgs.counts,"').style.display='block';document.getElementById('fade",irap.info2html.msgs.counts,"').style.display='block'\" class=\"button_info\">?</a></div>
               <div id=\"light",irap.info2html.msgs.counts,"\" class=\"info_content\">",
               "<H2>Information</H2>",
               "<p>",
               msg,
               "</p><a href=\"javascript:void(0)\" onclick=\"document.getElementById('light",irap.info2html.msgs.counts,"').style.display='none';document.getElementById('fade",irap.info2html.msgs.counts,"').style.display='none'\" class=\"button\">X</a></div><div id=\"fade",irap.info2html.msgs.counts,"\" class=\"black_overlay\"></div>
               ",
               sep=""))
}

get.path2template <- function(template) {
  p<-paste(IRAP.DIR,"aux/html",template,sep="/")
  if ( length(grep(".html",p,fixed=T))==0 ) {
    p <- paste(p,".html",sep="")
  }
  p
}

# generates a png file and a PDF file
# return the HTML to display the image
# with (optionally) links to the pdf
# and data (TSV) file
# returns a list with the filenames for png, pdf, data, and
# the html code
# (png.file=,pdf.file=,data.file=,html=)
# to.plot - function to make the plot
gen.plot2report <- function(filename=NULL,
                     dir=NULL,
                     bg="white",
                     width=400,
                     height=400,
                     size.unit="px",
                     html=TRUE,
                     ps=TRUE,
                     data.table=NA,
                     caption="",
                     to.plot=NULL) {

  if ( is.null(dir) ) {
    if ( grepl("^/",filename) ) {
      # abspath
      dir <- ""
    } else {
      dir <- "."
    }
  }
  if ( dir != "" ) {
    dir <- paste(dir,"/",sep="")
  }  
  # automatic filename  
  if ( is.null(filename)) {
    time.label <- format(Sys.time(), "%d%m%Y%H%M%S")
    filename <- paste("graph_",time.label,".png",sep="")
  }
  ps.file <- NA
  data.file <- NA
  png.file <- paste(dir,filename,sep="")
  # width and height is in pixels
  # convert to in
  if ( size.unit == "px" ) {
    # width: 600px=20cm=7.8in
    width.in <- width*7.8/600
    height.in <- height*7.8/600 
  } else {
    width.in <- width
    height.in <- height
  }
  # PNG
  png(filename=png.file,width=width.in, height=height.in,bg=bg,res=300,units="in")
  # Catch the errors
  err <- try(to.plot())
  #print(err)
  dev.off()
  if ( class(err)=="try-error") {
    pwarning("Failed to generate plot ",png.file)
    return(NULL);
  }
  # PS
  if (ps==TRUE) {
    ps.file <- paste(png.file,".eps",sep="")
    #ps(file=ps.file)
    postscript(file=ps.file)
    err <- try(to.plot())
    dev.off()
    if ( class(err)=="try-error") {
      pwarning("Failed to generate plot ",png.file)
      return(NULL);
    }
  }
  # DATA
  if ( sum(!is.na(data.table)) > 0 ) {
    data.file <- paste(png.file,".tsv",sep="")
    write.tsv(data.table,data.file)
  }
  png.file
  # HTML
  if ( html == TRUE ) {
    html <- paste("<DIV class=\"dplot\"><IMG src='",basename(png.file),"' border =0 width=",width," height=",height,"><BR/>",sep="")
    files <- c()
    if (ps==TRUE) {
      files <- append(files,basename(ps.file))
      names(files)[length(files)] <- "PS"
    }
    if (!is.na(data.file)) {
      files <- append(files,basename(ps.file))
      names(files)[length(files)] <- "TSV"
    }
    #
    html <- paste(html,html.download.bar(files),"</DIV>",sep="")
    html <- paste(html,"<p class='caption'>",caption,"</p>",sep="")
  } else {
    HTML=FALSE
  }
  return(list(png.file=png.file,html=html,data.file=data.file,ps.file=ps.file))
}

html.download.bar <- function(files) {

  if( length(files) ==0 ) { return("") }
  html <- "<DIV class=\"downbar\"> Download: "
  
  for ( f in names(files) ) {
    html <- paste(html,"&nbsp;<a class='download' href='",basename(files[f]),"'>",f,"</a>",sep="")
  }
  html <- paste(html,"</DIV>",sep="")
  return(html)
}

# deprecated
gen.plot <- function(filename=NULL,
                     dir=".",
                     bg="white",
                     caption="",
                     width=400,height=400,
                     to.plot=NULL) {

  # automatic filename
  if ( is.null(filename)) {
    time.label <- format(Sys.time(), "%d%m%Y%H%M%S")
    file <- paste("graph_",time.label,".png",sep="")
  }
  plot.filename <- paste(dir,filename,sep="/")
  png(filename=plot.filename,width=width, height=height,bg=bg)
  to.plot()
  dev.off()
  pdf(file=paste(plot.filename,".pdf",sep=""),width=width, height=height)
  to.plot()
  dev.off()
  plot.filename
}

######################################################
HTML.pprint.conf <- function(conf,names=NULL) {
  import.conf.variables(conf)
  HTML( "<h1>Project Information</h1>" )
  if ( ! is.null(names) ) {
    name <- paste(names,sep=" ",collapse=",")
  }
  HTML( paste("Project name:", name) )
  HTML( paste("Species:", species) )

  contrasts <- conf.get.value(conf,"contrasts")
  if ( is.null(contrasts) ) {
     # default (deprecated)
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
    HTML (paste("Contrasts:", paste(contrasts,collapse=' ', sep=' '),sep=' '))
    
    #contrasts.v <- myread.list(contrasts)
    contrasts.v <- contrasts
    has.contrasts <- 1
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
      HTML("<UL>" )
      for (c in contrasts.v) {
        pinfo(c)
        contrast.v <- myread.list(c)
        pinfo(paste(c,"-->",contrast.v,sep="//"))
        html <- ""
        for ( l in contrast.v) {
          visible <- paste("<LI>",l,"=",myderef(l),sep=" ",collapse=" ")
          html <- "<UL>"
          for (lib in myread.list(myderef(l)) ) {
            pinfo(lib)
            html <- paste(html,"<LI>",lib,"=",myderef(lib),sep=" ",collapse=" ")
            lib2cond[lib] <- c
          }
          html <- paste(html,"</UL>",sep=" ",collapse=" ")
        }
        pinfo(visible)
        pinfo(html)
        HTML(HTML.toogle.div(visible,html))
      }
      HTML("</UL>")
    }
    qual_filtering<- conf.get.value(conf,"qual_filtering")
    if ( is.null(qual_filtering) || qual_filtering!="on") {
      HTML( paste("Reads quality filtering:", "Disabled") )
    } else {
      HTML( paste("Reads quality filtering:", "Enabled") )
    }
######
    HTML( paste("Data") )
    HTML( paste("<UL>") )
    pe <- conf.get.value(conf,"pe")
    if (! is.null(pe) ) {
      HTML( paste("<LI>Pair-end libraries:",paste(pe,collapse=" "), sep=" ", collapse=" ") )
      pe.v <- myread.list(pe)
    } else {
      pe.v <- c()
    }
    se <- conf.get.value(conf,"se")
    if (! is.null(se) ) {
      HTML( paste("<LI>Single-end libraries:",paste(se,collapse=" "), sep=" ") )
      se.v <- myread.list(se)
    } else {
      se.v <- c()
    }
    
    HTML( paste("</UL>") )
  }
}

irap.HTMLInit <- function(outdir,filename,title) {
  HTMLInitFile(outdir=outdir,
               filename=filename,
               Title=title,
               useLaTeX=FALSE,
               CSSFile=basename(irap_css_file),               
               useGrid=FALSE)
  system(paste("cp ",irap_css_file," ",outdir,sep=""));
  system(paste("cp ",irap_menu_css_file," ",outdir,sep=""));
  # Javascript code

  HTML( '<script language="javascript"> 
function toggle(showHideDiv, switchTextDiv) {
	var ele = document.getElementById(showHideDiv);
	var text = document.getElementById(switchTextDiv);
	if(ele.style.display == "block") {
    		ele.style.display = "none";
		text.innerHTML = "+";
  	}
	else {
		ele.style.display = "block";
		text.innerHTML = "-";
	}
} 
</script>
')

}

irap.HTMLEnd <- function() {
  # IRAP footer
  HTML("<HR>")
  HTML(paste("IRAP ",irap_version,sep=""))
  HTMLEndFile()
}
############### 
div.ids <- 0
HTML.toogle.div <- function(visible.text, hidden.text) {
  assign("div.ids",div.ids+1,envir = .GlobalEnv)
  idon  <- paste("t",div.ids,"on",sep="")
  idoff <- paste("t",div.ids,"off",sep="")
  paste(
        "<a class=\"button\" id='",idoff,"' href='javascript:toggle(\"",idon,"\",\"",idoff,"\");'>",visible.text,"</a>",
        "<div class=\"button\" id='",idon,"' style='display: none'>",hidden.text,"</div>",sep="")
}

pprint.fieldname <- function(x) {
  return(gsub("\\."," ",x))
}

######################################################
#Mem
my.ls <- function (pos = 1, sorted = FALSE, envir = as.environment(pos))
{
    .result <- sapply(ls(envir = envir, all.names = TRUE),
                      function(..x) object.size(eval(as.symbol(..x),
                                                     envir = envir)))
    if (sorted) {
      .result <- rev(sort(.result))
    }
    .ls <- as.data.frame(rbind(as.matrix(.result), `**Total` = sum(.result)))
    names(.ls) <- "Size"
    .ls$Size <- formatC(.ls$Size, big.mark = ",", digits = 0,
                        format = "f")
    .ls$Mode <- c(unlist(lapply(rownames(.ls)[-nrow(.ls)], function(x)
                                mode(eval(as.symbol(x),
                                          envir = envir)))), "-------")
    .ls
}
#######################################################
# tries to fix a known bug in cufflinks (some gene entries are not summed up)
# http://seqanswers.com/forums/showthread.php?t=5224
fix.cufflinks.fpkms <- function(t,id.col=1,value.col=2) {
  f<-aggregate(t[,value.col],by=list(Gene=as.character(t[,id.col])),sum)
    if ( nrow(t) !=  nrow(f) ) {
    pwarning("Removed duplicates duplicates:",nrow(t),"--->",nrow(f),"")
  }
  f
}

#######################################################
#
write.tsv <- function(x,file,header=TRUE,rownames.label=NULL) {
  #
  if (!is.null(x)) {
    if (!is.matrix(x)) {
      x <- as.matrix(x)
    }
    x <- apply(x,c(1,2),gsub,pattern="'",replacement="")
    if ( !is.null(rownames.label) ) {
      y <- cbind(rownames(x),x)
      colnames(y) <- append(rownames.label,colnames(x))
      x <- y
    }
  }
  write.table(x,file,sep="\t",row.names=F,col.names=header,quote=F)
  return(1)
}

read.tsv <- function(file,header=T) {
  read.csv(file,sep="\t",header,quote="\"")
}

#
# Return true if the matrix x
# corresponds to a "Empty" file produced by IRAP
# 
is.irap.empty <- function(x) {

  if (typeof(x)!="list") return (TRUE);
  if (ncol(x)==1 && colnames(x)[1]=="Empty.Empty") return (TRUE)
  return(FALSE)
}

my.factor2numeric <- function(v) {
  as.numeric(levels(v)[v])
}
# Given two tables
# merge them using the two given columns
# handle name clashes
mergeAnnot <- function(table,annot,table.field,annot.field="ID") {
  # Annotation file (generated from the GTF file)
  # gene_id locus@(chr:start...end) Description GO
  m<-merge(table,annot,by.x=table.field,by.y=annot.field,all.x=TRUE,sort=F)
  # handle col name clashes
  # rename the initial columns
  newNames <- gsub(".x$","",colnames(m))
  pdebug("Cols=",newNames)
  colnames(m) <- newNames
  m
}

matrices.ok<- function(matrices) {
  # assumption: all the tables have the same col. names, the same number of columns and rows
  n.rows <- nrow(matrices[[1]])
  n.cols <- ncol(matrices[[1]])
  labels <- names(matrices)
  for (l in labels) {
    if (nrow(matrices[[l]])!=n.rows) {
      stop(paste("Error: number of rows differs (",l,"/",labels[1],")\n",sep=""))
      return(FALSE);
    }
    if (ncol(matrices[[l]])!=n.cols) {        
      stop(paste("Error: number of cols differs (",l,"/",labels[1],")\n",sep=""))
      return(FALSE);
    }
  }
  return(TRUE)
}


load.matrices <- function(files,check.col.number=TRUE) {

  if (! is.vector(files)) {
    stop("files must be a vector', 'load.tsv.files'.")
  }

  matrices <- list()
  i <- 1
  for (f in files) {
    if (!file.exists(f)) {
      stop(paste("file does not exist ',",f,", 'load.tsv.files'.",sep=""))
    }
    pinfo("Found ",f,", loading...")
    matrices[[f]] <- read.tsv(f)
    i <- i+1
  }

  if (check.col.number) {
    if (!matrices.ok(matrices)) {
      exit(1)
    }
  }
  return(matrices)
}
#
listOflists2matrix <- function(l) {
  # names
  m <- as.matrix(as.matrix(l[[1]]),ncol=1)
  # merge
  for ( n in names(l)[-1] ) {
    m <- merge(m,as.matrix(l[[n]],ncol=1),by=0)
    rownames(m) <- m[,1]
    m <- m[,-1]
  }
  colnames(m) <- names(l)
  m
}
matrices.list2matrix <- function(matrices) {
  aggr.table <- matrix(ncol=length(names(matrices)),nrow=nrow(matrices[[1]]))
  colnames(aggr.table) <- names(matrices)
  rownames(aggr.table) <- matrices[[1]][,1]
  aggr.table[,1] <- rowSums(matrices[[1]][,-1])
  for (i in seq(2,length(matrices))) {
    m <- matrices[[i]]
    row.names <- m[,1]
    sum <- rowSums(m[,-1])
    aggr.table[,i] <- sum[match(row.names,rownames(aggr.table))]
  }
  aggr.table
}

matrix.sum.cols <- function(matrix,labels) {

  sum.mat <- matrix(ncol=length(unique(labels)),nrow=nrow(matrix))
  colnames(sum.mat) <- unique(labels)
  rownames(sum.mat) <- rownames(matrix)
  for (l in unique(labels)) {
    sum <- rowSums(matrix[,labels==l])
    sum.mat[,l] <- sum
  }
  sum.mat
}

merge.matrices <- function(matrices) {

  col.names <-   unlist(lapply(names(matrices),paste,colnames(matrices[[1]])[-1],sep=""))
  aggr.table <- matrix(ncol=length(col.names),nrow=nrow(matrices[[1]]))
  colnames(aggr.table) <- col.names
  rownames(aggr.table) <- matrices[[1]][,1]
  i1 <- 1
  for (m.name in names(matrices) ) {
      m <- matrices[[m.name]]
      head(m)
      for (i2 in seq(2,ncol(m))) {
        pinfo(i2)
        aggr.table[,i1] <- m[,i2]
        i1 <- i1+1
      }
    }
  aggr.table
}

# Given the map of a file to a group and a vector with colnames
# produce the vector with the respective group names
map.conds2cols <- function(label2group,cols) {
  deseq.conds <- cols
  i <- 1
  for (c in cols) {    
    deseq.conds[i] <- label2group[[deseq.conds[i]]]
    i <- i+1
  }
  deseq.conds
}


# compute the average values and sd
data2groups <- function(data.matrix,conds) {
  groups <- unique(conds)
  avg.matrix <- data.matrix[,c(1:length(groups))]
  colnames(avg.matrix) <- groups
  sd.matrix <- avg.matrix
  for ( g in groups ) {
    sel.cols <- conds==g
    pdebug("group=",g," sel.cols=",sum(sel.cols))
    irap.assert(sum(sel.cols)>0,"Columns selected")
    if (sum(sel.cols)==1) {
      # single element
      pwarning("group ",g," with a single column, unable to compute sd")
      sd.matrix[,g] <- 0
      avg.matrix[,g] <- data.matrix[,sel.cols]
    } else {
      sd.matrix[,g] <- apply(data.matrix[,sel.cols],1,FUN=sd)
      avg.matrix[,g] <- apply(data.matrix[,sel.cols],1,FUN=mean)
    }
  }
  list(mean=avg.matrix,sd=sd.matrix)
}

rowVariance <- function(x) {

  sqrt=function(x) { x * x }
  n <- rowSums(!is.na(x))
  n[n<=1] <- NA
  var <- rowSums(sqrt(x - rowMeans(x))/(n-1))
  return(var)
}
#test
rep.v<-function(v,n) {
  return(unlist(lapply(v,replicate,n=n)))
#  n.v <- c()
#  for ( e in v ) {
#    n.v <- append(n.v,replicate(n,e))
#  }
#  n.v
}

# Load the annotation (tsv file or cached Rdata file) Try loading
# first a R representation of the data to speedup things...loading the
# TSV file can take a while
load.annot <- function(file) {
  if ( is.null(file) || !is.character(file) ) {
    return(NULL)
  }
  pdebug("Loading annotation=",file)
  cached.annot <- paste(file,".Rdata",sep="")
  if ( file.exists(cached.annot) ) {
    load(cached.annot)
    if ( exists("gene.annot")) {
      annot.table <- gene.annot
    }
  } else {
    annot.table <- read.tsv(file)
    save(list=c("annot.table"),file=cached.annot)
  }
  annot.table <- annot.expand.fields(annot.table)
  pdebug("Loading annotation (done)")
  return(annot.table)
}
# load a file with a quant. matrix
# returns NULL in case of failure
quant.load <- function(f,clean.cuff=FALSE) {
  tsv.data <- NULL
  tryCatch(tsv.data <- read.tsv(f),error=function(x) NULL)
  if ( !is.null(tsv.data) && ncol(tsv.data)>1 ) {
    rownames(tsv.data) <- as.character(tsv.data[,1])
    tsv.data <- tsv.data[,-1]
    if (clean.cuff) {
      sel<-grep("^CUFF.*",rownames(tsv.data),perl=T,invert=TRUE)
      tsv.data <- tsv.data[sel,]
    }
    return(tsv.data)
  }
  NULL
}
# 1=ok
# 0 error
quant.check.quant.matrix.OK <- function(filename) {
  
  # ncol==0 for empty irap tsv files
  tsv.data <- quant.load(filename,clean.cuff=T)
  if ( !is.null(tsv.data) && ncol(tsv.data)>1 ) {
    #pinfo(f,"--->",paste(apply(tsv.data,2,sum),sep=","))
    if ( sum(as.numeric(is.na(tsv.data)))> 0 || sum(rowSums(tsv.data)>0) == nrow(tsv.data) ) {
      pwarning("File ",filename," with only 0s or NAs!")
      return(FALSE)
    } else {
      # look for columns with 0s
      col.sum <- apply(tsv.data,2,sum)
      if ( 0 %in% col.sum || sum(col.sum<nrow(tsv.data)*0.02)>0 ) {
        pwarning("File ",filename," with a column with only 0s!")
        return(FALSE)
      }
      return(TRUE)
    }
  }
  return(FALSE)
}
################################################
#

quant.heatmap <- function(data,nlim=50,do.cor=FALSE,density.info="histogram",...) {
  suppressPackageStartupMessages(library("gplots"))
  var <- rowVariance(data)
  # select a subset of rows
  select <- order(var,decreasing=TRUE)[c(1:nlim)]
  colors <- topo.colors(nlim)
  data <- data[select,]            
  if ( do.cor ) {
    data <- cor(data)
  }
  # sort the genes by variability
  heatmap.2(data,col = colors, scale = "none",cexCol=0.9,cexRow=0.6,keysize=1,density.info=density.info,trace="none",...)
}

# general function to plot an heatmap
gen.heatmap <- function(data,colors=NULL,ncolors=NULL,do.cor=FALSE,trace="none",cexCol=0.9,cexRow=0.9,keysize=1,density.info="histogram",...) {
  suppressPackageStartupMessages(library("gplots"))
  if ( is.null(ncolors) ) {
    ncolors <- nrow(data)*ncol(data)
  }
  if ( is.null(colors) ) {
    colors <- topo.colors(ncolors)
    if ( do.cor ) {
      data <- cor(data)
    }
  } 
  # sort the genes by variability
  heatmap.2(data,col = colors, scale = "none",cexCol=cexCol,cexRow=cexRow,keysize=keysize,density.info=density.info,trace=trace,...)
}

plot.panel.label <- function(label,add=")",cex=1.5,col="black") {
  if ( par("ylog") ) {
    at=10^par("usr")[4]
  } else {
    at=par("usr")[4]
  }
  mtext(paste(label,add,sep=""), side=2, at=at,line=2.5,cex=cex,las=2,col=col)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol 
plot.b2b.barplot <- function(minus,plus,groups,minus.labels=NULL,bottom.labels=c("-","+"),plus.labels=NULL,cex=0.8,show.var=FALSE) {
  
  plus <- as.numeric(plus)
  minus <- as.numeric(minus)
  
  par(mar=c(0.5, 5, 0.5, 1))
  plot.new()
  max.val <- max(max(abs(plus)),max(abs(minus))) 

  num.ticks <- 5
  if (is.na(max.val) ) {
    cat("Unable to generate b2b plot\n")
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

################################################
# File functions
################################################

file.required <- function(file,fatal=TRUE) {
  if (!file.exists(file)) {
    if (fatal) {
      cat(paste("ERROR: ",file," not found\n",sep=""))
      q(status=1)
    }
    pmissing(file)
    return(NULL)
  }
  file
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
#################################################################################################

myParseArgs <- function(usage,option_list,filenames.exist=NULL,multiple.options=NULL,mandatory=NULL) {

  # get command line options, if help option encountered print help and exit,
  # otherwise if options not found on command line then set defaults,
  parser <- OptionParser(usage = usage, option_list=option_list)
  opt <- parse_args(parser)
  
  for ( m in mandatory ) {
    if ( is.null(opt[[m]]) ) {
        perror("Parameter ",m," needs to be defined")
        q(status=1)
    }
  }
  
  for ( p in filenames.exist ) {
    if (! is.null(opt[[p]]) ) {
      if (! file.exists(opt[[p]]) ) {
        perror("File ",opt[[p]]," not found")
        q(status=1)
      }
    }
  }

  for ( op in names(multiple.options) ) {
    if ( ! opt[[op]] %in% multiple.options[[op]] ) {
      perror("Invalid value ",opt[[op]]," for option ",op)
      q(status=1)
    }
  }
  return(opt)
}
#
#################################################################################################
conf.get.value <- function(conf,val) {
  #pinfo(">>>>>>>>>>>",val)
  idx <- grep(paste("^",val,"$",sep=""),conf[,1],ignore.case=TRUE,value=FALSE,perl=TRUE)
  
  #pinfo("AAAAAAAAAAAAAAA:",val,idx,is.null(idx),typeof(idx))
  if (is.null(idx) || length(idx)==0) {
    #
    if ( length(grep(".*_strand$",val,ignore.case=TRUE,value=FALSE,perl=TRUE))>0 ) {
      return(NA)
    }
    return(NULL)
  }
  return(myread.list(conf[idx,2]))
}

conf.set.value <- function(var,value,conf) {
  if ( ! var %in% rownames(conf) ) {
    conf <- rbind(c(var,value),conf)
    rownames(conf)[nrow(conf)] <- var
  } else {
    pinfo(conf[var,])
    conf[var,] <- c(var,value)
  }
  return(conf)
}

conf.is.qc.enabled <- function(conf) {
  qc.enabled <- as.logical(conf.get.value(conf,"qc.enabled"))
  return(qc.enabled)
}

conf.get.default.value <- function(var) {
  cmd <- paste("grep -i '^def_",var,"' ",IRAP.DIR,"/scripts/irap | cut -f 2 -d= | tail -n 1",sep="")
  #pinfo(cmd)
  val <- system(cmd,intern=TRUE)
  #pinfo("default ",var,"=",val)
  if (is.null(val)) { return(NULL) }
  return(val)
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


pe.lib2rawfilename2 <- function(s,conf) {
  id <- sub(".*/","",s)
  lib <- get.libname(s)
  files <- conf.get.value(conf,lib)
  return(files[as.numeric(id)])
}
#
load.configuration.file <- function(conf_file) {

  # 
  pinfo("Loading configuration file...")
  conf.table <- read.delim(conf_file,sep="=",comment.char="#",header=FALSE,stringsAsFactors=FALSE)
  pinfo("Loading conf_file...done")
  pinfo("Configuration:")
  for ( i in 1:nrow(conf.table) ) {
    assign(as.character(conf.table[i,1]),mytrim(conf.table[i,2]))
    pinfo("",conf.table[i,1],"=",conf.table[i,2])
  }
  pinfo("Loading default parameters...")
  #################
  # Default values
  vars <- c("mapper","quant_method","de_method","qual_filtering","qual_filtering")
  for ( v in vars ) {
    if ( sum(as.character(conf.table[,1])==v)==0 ) {
      # get default value from IRAP main file
      entry <- c(v,conf.get.default.value(v))
      conf.table[nrow(conf.table)+1,] <- entry
    }
  }
  rownames(conf.table) <- conf.table[,1]
  pinfo("Loading default parameters...done.")
  conf.table[nrow(conf.table)+1,] <- c("toplevel.dir",paste(dirname(conf_file),conf.table["name",2],sep="/"))
  #alias
  if ( ! "qc" %in%  rownames(conf.table) ) {
    conf.table[nrow(conf.table)+1,] <- c("qc",conf.table["qual_filtering",2])
  }
  conf.table[nrow(conf.table)+1,] <- c("qc.enabled",(conf.get.value(conf.table,"qc")=="on"||conf.get.value(conf.table,"qc")=="yes"))
  #
  rownames(conf.table) <- conf.table[,1]
  
  conf.table
}

# imports to global env. all the variables loaded from the configuration
# file and stored in conf.table
import.conf.variables <- function(conf.table) {
  for ( i in 1:nrow(conf) ) {
    assign(as.character(conf[i,1]),mytrim(conf[i,2]),envir = .GlobalEnv)
  }
}

irap.conf2list<- function(conf.table) {
  irap.conf <- list()
  for ( i in 1:nrow(conf) ) {
    irap.conf[[as.character(conf[i,1])]] <- mytrim(conf[i,2])
  }
  return(irap.conf)
}

#####################################################3
#
# filter by the source column in the gtf file
source.filt.groups.def <- list(
  "all"=NA,
  "protein coding"=NA,
  "pseudogenes"=NA,
  "xRNA"=NA
)
# Ids of the genes on each group
source.filt.groups <- list("all"=TRUE,
                    "pseudogenes"=NA,
                    "protein coding"=NA,
                    "xRNA"=NA
                    )



# Initialization
init.source.filter <- function(table) {
  if ( "source" %in% colnames(table) ) {
    source.filt.groups.def <- list()
    source.filt.groups <- list()
    sources <- levels(table$source)
    #pdebug(sources)
    # all
    source.filt.groups.def$"all" <- sources
    source.filt.groups$"all" <- table[,"ID"]

    protein_coding <- c("protein_coding")
    source.filt.groups.def$"protein coding" <- protein_coding
    protein_coding.filt <- sapply(table$source,`%in%`,protein_coding)
    source.filt.groups$"protein coding" <- table[protein_coding.filt,"ID"]

    pseudogenes <- grep("pseudogene",sources,value=T)
    source.filt.groups.def$pseudogenes <- pseudogenes
    pseudogenes.filt <- sapply(table$source,`%in%`,pseudogenes)
    source.filt.groups$pseudogenes <- table[pseudogenes.filt,"ID"]
    if (length(pseudogenes.filt)>0) {
      pdebug(sum(pseudogenes.filt))
    }
    
    xRNA <- grep("RNA$",sources,value=T)
    source.filt.groups.def$"xRNA" <- table[xRNA,"ID"]
    xRNA.filt <- sapply(table$source,`%in%`,xRNA)
    source.filt.groups$"xRNA" <- table[xRNA.filt,"ID"]
    # update the global variable
    assign("source.filt.groups.def",source.filt.groups.def, envir = .GlobalEnv)
    assign("source.filt.groups",source.filt.groups, envir = .GlobalEnv)
  }
  T
}
#
apply.source.filter <- function(table,filter.name) {

  #ids<-annot.table[sapply(annot.table$source,`%in%`, source.filt.groups[[filter.name]]),"ID"]  
  ids <- source.filt.groups[[filter.name]]
  if (is.matrix(table) || is.data.frame(table) ) {
    table <- table[rownames(table) %in% ids,]
  } else {
    if (is.vector(table)) {
      table <- table[names(table) %in% ids]
    }
  }
  #aggr.data.annot <- mergeAnnot(aggr.data,annot.table[,c("ID","source")],table.field="row.names",annot.field="ID")
  return(table)
}
#
get.source.filename <- function(file.prefix,type) {
  if ( type =="all" || type=="All" ) {
    return(gsub(" ","_",paste(file.prefix,".html",sep="")))
  }
  gsub(" ","_",paste(file.prefix,"_",type,".html",sep=""))
}
# outprefix -> menu html
get.source.filter.menu <- function(file.prefix) {
   get.source.menu.entry <- function(name,file.prefix) {
     paste("<a href='",basename(get.source.filename(file.prefix,name)),"'>",name,"</a>",sep="")
   }
   sources.menu <- paste(sapply(names(source.filt.groups.def),get.source.menu.entry,file.prefix),collapse="|")
   sources.menu
 }

get.source.filter.names <- function() {
  names(source.filt.groups.def)
}

get.source.filter.size <- function(filter.name) {
  length(source.filt.groups[[filter.name]])
}
