# =========================================================
# Copyright 2012-2017,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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

irap_version<-"0.8.5d3"


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
  cat(paste0("[INFO] ",...,"\n"))
}

pwarning <- function(...) {
  cat(paste0("[WARNING] ",...,"\n"),file=stderr())
}
pmissing <- function(...) {
  cat(paste0("[MISSING FILE] ",...,"\n"))
}

perror <- function(...) {
  cat(paste0("[ERROR] ",...,"\n"),file=stderr())
}

capitalize <- function(s) {
  if ( is.na(s) || is.null(s) ) { return(s); }
  return(paste0(toupper(substring(s, 1,1)),substring(s,2)))
}
##################################################
# debugging
pdebug.enabled <- FALSE
pdebug.stage <- NULL

pdebug <- function(...) {
  if ( pdebug.enabled) {
    cat(paste0("[DEBUG] ",...,"\n"),file=stderr())
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
    save.image(file=paste0(file,".Rdata"))
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
attributes.cols <- list(
  gencode=c("gene_id","transcript_id","exon_number","gene_name","gene_type","gene_status","transcript type","level","transcript_name","havana_gene","exon_id","gene_biotype"),
  ensembl=c("gene_id","transcript_id","exon_number","gene_name","gene_biotype","transcript_name","protein_id","exon_id","exonic_part_number")
  )
# TODO: improve error handling
load.gtf <- function(gtf.file,feature=NULL,selected.attr=NULL,gtf.format="auto") {
  suppressPackageStartupMessages(library(parallel))
  gtf<-qload.tsv(f=gtf.file,header=F,comment.char = "#")
  cnames <- formats.cols$gtf
  colnames(gtf)<-cnames[0:ncol(gtf)]

  if (!is.null(feature)) {
    gtf<- gtf[gtf$feature %in% feature,,drop=FALSE]
  }
  gtf$attributes <- as.character(gtf$attributes)
  # detect format
  if ( gtf.format=="auto" ) {
    if (sum(grepl("level",head(gtf$attributes,50))) >0) {
      gtf.format <- "gencode"
    } else {
      gtf.format <- "ensembl"
    }
    cat("GTF attributes ",gtf.format,"\n")
  }
  gtf.attributes.names<-attributes.cols[[gtf.format]]

  if ( !is.null(selected.attr) ) {
      selected.attr.i <- append(selected.attr,"gene_biotype")
      gtf.attributes.names<- gtf.attributes.names[gtf.attributes.names %in% selected.attr.i]
  }

  num.attr <- length(gtf.attributes.names)
  attr2vec <- function(s,gtf.attributes.names) {
    a<-strsplit(mytrim(s),split=";([ ]?)+")
    var.value <- function(s) {
      rm <- regexec(pattern="^([^ ]+) (.+)$",s)
      a2 <- regmatches(s,rm)
      return(a2[[1]][-1])
    }
    a2 <- unlist(lapply(a[[1]],var.value))
    #a2 <- unlist(strsplit(a[[1]]," ",fixed=T))
    m <- matrix(a2,nrow=2,byrow=F)
    x <- m[2,,drop=FALSE]
    names(x) <- m[1,]
    r <- rep(NA,num.attr)
    names(r) <- gtf.attributes.names
    attr.int <- names(x)[names(x)%in% gtf.attributes.names]
    r[attr.int] <- x[attr.int]
    r <- gsub("^\"","",gsub("\"$","",r))
    return(r)
  }
  attr <- mclapply(gtf$attributes,attr2vec,gtf.attributes.names)
  print(length(attr))
  if ( length(attr)!=0 ) {
  #attr2vec(gtf$attributes[1])
    vals<-matrix(unlist(attr),ncol=length(gtf.attributes.names),byrow=T)
    colnames(vals) <- gtf.attributes.names
    gtf <- cbind(gtf,vals)
  }
  # try to determine the biotype column (use source by default...)
  biotypes <- gtf[,biotype.column(gtf)]
  pinfo("biotype col:",biotype.column(gtf))
  gtf$biotype <- biotypes
  #print(head(gtf))
  return(gtf[,! colnames(gtf) %in% c("attributes")])
}

# deprecated
load.gencode.gtf <- function(gtf.file,feature=NULL,selected.attr=NULL) {
  return(load.gtf(gtf.file,feature=feature,selected.attr=selected.attr,gtf.format="encode"))
}

# Given a gtf file returns a vector with the length of the genes
get.gene.length.from.gtf.file <- function(gtf.file,filter.biotype=NULL,length.mode="union.exons") {
  gtf <- load.gtf(gtf.file)
  get.gene.length.from.gtf(gtf,filter.biotype,length.mode)
}


get.gene.length.from.gtf <- function(gtf,filter.biotype=NULL,length.mode="union.exons") {
  library(data.table)
  dt <- as.data.table(gtf)
  if ( !is.null(filter.biotype) ) {
    dt<-subset(dt,gene_biotype==filter.biotype)
  }

  dt<-subset(dt,feature=="exon",c("feature","gene_id","start","end","seqid"),key=c('seqid'))

  pfun <- function(chr,dt) {
      dt<-subset(dt,seqid==chr,key=c('gene_id'))
      #setkey(dt,gene_id)
      gene.ids <- unique(as.character(dt$gene_id))
      v <- rep(0,times=length(gene.ids))
      names(v) <- gene.ids
      #pinfo(chr,":",length(gene.ids))
      for ( gene.id in gene.ids ) {
          #print(gene.id)
          v[gene.id] <- get.gene.length(gene.id,gtf.data=subset(dt,gene_id==gene.id,drop=FALSE),mode=length.mode,exons.only=TRUE)
      }
      return(v)
  }
  chrs <- as.character(unique(dt$seqid))
  #pinfo("Chrs",chrs)
  x <- mclapply(chrs,FUN=pfun,dt=dt,mc.allow.recursive = FALSE)
  #print(length(x))
  return(unlist(x))
}


get.gene.length <- function(gene.id,gtf.data,mode="sum.exons",lim=+Inf,do.plot=FALSE,exons.only=FALSE) {
  library(data.table)

  if (!exons.only) {
    dt <- subset(gtf.data,feature=="exon")
  } else {
    dt <- gtf.data
  }

  start.end <- as.matrix(subset(dt,gene_id==gene.id,sel=c("start","end"),drop=FALSE))


  suppressPackageStartupMessages(library("intervals"))
  library(methods)

  i <- Intervals(start.end)
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
  return(res)
}
#
addGeneFeature2gtf <- function(gtf) {
  #
  #save.image()
  #gtf <- gtf2
  if (!is.data.frame(gtf)) {
    perror("Internal error: expected a data.frame and got ",typeof(gtf))
    q(status=1)
  }
  gtf2 <- gtf
  levels(gtf2$feature) <- append(levels(gtf2$feature),"gene")
  gtf2 <- gtf2[FALSE,]
  genes <- unique(as.character(gtf$gene_id))
  #g <- genes[1]
  for (g in genes ) {
    cat(".")
    sel <- as.character(gtf$gene_id)==g
    #print(gtf[sel,])
    gene.sel <- gtf[sel,]
    gtf <- gtf[!sel,]
    #print(gene.sel)
    ###print(gene.sel[,"start"])
    if (is.factor(gene.sel[,"start"]) ) {
      start <- min(min(my.factor2numeric(gene.sel[,"start"])),min(my.factor2numeric(gene.sel[,"end"])))
      pdebug("factor")
    } else {
      start <- min(min(gene.sel$start),min(gene.sel$end))
      pdebug("not factor")
    }
    if (is.factor(gene.sel[,"end"]) ) {
      end <- max(max(my.factor2numeric(gene.sel[,"end"])),max(my.factor2numeric(gene.sel[,"start"])))
    } else {
      end <- max(max(gene.sel[,"start"]),max(gene.sel[,"end"]))
    }
    new.entry <- c(gene.sel[1,])
    new.entry$feature <- "gene"
    new.entry$start <- start
    new.entry$end <- end
    # reset values
    vals2reset <- c("transcript_id","exon_number","transcript_name")
    for (val  in vals2reset) {
      if ( val %in% names(new.entry) ) {
        new.entry[val] <- ""
      }
    }
    #gtf2 <- rbind(gtf2,new.entry)
    #pinfo(new.entry)
    gtf2[nrow(gtf2)+1,] <- new.entry
  }
  #pinfo("GTF2:",nrow(gtf2))
  return(gtf2)
}


## # Given a matrix obtained from a gtf file returns a vector with the length of the transcripts
## get.transcript.length.from.gtf <- function(gtf,filter.biotype=NULL) {
##   # TODO: validate gtf
##   # protein coding
##   suppressPackageStartupMessages(library(parallel))
##   library(RSQLite)
##   gtf.con <- dbConnect(RSQLite::SQLite(), ":memory:")
## # create a table for the GTF data
##   dbWriteTable(gtf.con,"GTF",gtf)

##   bio.sql <- ""
##   if ( !is.null(filter.biotype) ) {
##     bio.sql <- paste(" AND gene_biotype='",filter.biotype,sep="")
##   }
##   # compute the length for each transcript
##   gene.ids <- unique(dbGetQuery(gtf.con,paste("select gene_id from GTF where feature='exon' ",bio.sql,sep=""))[,1])
  
##   get.gene.transcripts.length <- function(gene.id,gtf.con,bio.sql) {
##     gtf <-  dbGetQuery(gtf.con,paste("select * from GTF where feature='exon' and gene_id='",gene.id,"' ",bio.sql,sep=""))
##     transcript.ids <- unique(gtf$transcript_id)
##     tlen <- unlist(mclapply(transcript.ids,get.transcript.length,gtf))
##     names(tlen) <- transcript.ids
##     return(tlen)
##   }
##   tlen <- unlist(mclapply(gene.ids,get.gene.transcripts.length,gtf.con,bio.sql=bio.sql))
##   dbDisconnect(gtf.con)
##   return(tlen)
## }

get.transcript.length.from.gtf <- function(gtf,filter.biotype=NULL) {
  # TODO: validate gtf
  # protein coding
    suppressPackageStartupMessages(library(parallel))
    library(data.table)
    dt <- as.data.table(gtf)
    
    if ( !is.null(filter.biotype)) {
        dt <- dt[dt$gene_biotype==filter.biotype,]
    }
    dt<-subset(dt,dt$feature=="exon",c("feature","transcript_id","start","end"))
    dt$elen <- abs(dt$end-dt$start)+1
    setkey(dt,"transcript_id")

    tids <- unique(dt$transcript_id)
    tlen <- rep(0,times=length(tids))
    names(tlen) <- tids
    len <- aggregate(dt$elen,by=list(tid=dt$transcript_id),FUN=sum)
    tlen <- len[,2]
    names(tlen) <- len[,1]
    dt <- NULL
    return(tlen)
}

## # Given a matrix obtained from a gtf file returns a vector with the length of the transcripts
## get.transcript.length.from.gtf <- function(gtf,filter.biotype=NULL) {
##   # TODO: validate gtf
##   # protein coding
##   suppressPackageStartupMessages(library(parallel))
##   if ( !is.null(filter.biotype) ) {
##     gtf <- gtf[gtf$gene_biotype==filter.biotype,,drop=FALSE]
##   }
##   # compute the length for each transcript
##   gtf <- gtf[gtf$feature=="exon",,drop=FALSE]
##   transcript.ids <- unique(gtf$transcript_id)
##   tlen <- unlist(mclapply(transcript.ids,get.transcript.length,gtf))
##   names(tlen) <- transcript.ids
##   tlen  
## }


# Sum the length of the exons
get.transcript.length <-  function(transcript.id,gtf.data) {
  s <- subset(gtf.data,transcript_id==transcript.id,c("start","end"),drop=FALSE)
  len <- sum(abs(s$start-s$end))+nrow(s)
  return(len)
}

# Given a matrix obtained from a gtf file returns a matrix with the length of the  exons of a given gene/transcript
get.exon.length.from.gtf <- function(gtf,filter.biotype=NULL) {
  # TODO: validate gtf
  # protein coding
  if ( !is.null(filter.biotype) ) {
    gtf <- gtf[gtf$gene_biotype==filter.biotype,,drop=FALSE]
  }
  # compute the length for each exon
  gtf <- gtf[gtf$feature=="exon",,drop=FALSE]  
  elen <- abs(gtf$start-gtf$end)+1
  gtf$elength <- elen
  return(gtf[,c("exon_id","gene_id","transcript_id","exon_number","elength")])
}


load.gff3 <- function(file,type="gene",attrs=TRUE,selected.attrs=NULL) {
  # load the gff3 file
  gff3<-try(read.table(file,sep="\t",header=F,quote="\"",comment.char ="",
                       nrows=5000000))
#                       colClasses=c('character','character','character','integer','integer','character','character','character','character','character')))
  if(class(gff3)=='try-error') {
    pwarning("Loading failed")
    gff3=as.data.frame(matrix(nrow=0,ncol=length(formats.cols$gff3)))
  }
  #print(head(gff3))
  cnames <- formats.cols$gff3[seq(1,ncol(gff3))]
  colnames(gff3)<-cnames[0:ncol(gff3)]
  gff3$attributes <- as.character(gff3$attributes)

  if ( !is.na(type) ) {
    gff3 <- gff3[gff3$type==type,,drop=FALSE]
  }
  # convert the start/end to numbers
  gff3$start <- as.integer(gff3$start)
  gff3$end <- as.integer(gff3$end)

  if(attrs ) {
    # add the attributes
    if (is.null(selected.attrs)) {
      tags<-c("ID","Name","Alias","Parent","Target","Gap","Derives","Note","DBxref","Ontology_term","Is_circular")
    } else {
      tags <- selected.attrs
    }
    pdebug("tags")
    if (nrow(gff3)>0 ) {
      for (tag in tags ) {
        pattern <- paste0("",tag,"=([^;]+)")
        m<-regexec(pattern,as.vector(gff3$attributes))
                                        #m<-regexec("^ID=([^;]+);",as.vector(gff3$attributes))
        for ( i in c(1:length(m)) ) { if ( m[[i]][1]==-1 ) { m[[i]]=NA; } }
        x<-regmatches(as.vector(gff3$attributes),m)
        for ( i in c(1:length(x)) ) { if ( length(x[[i]])==0 ) { x[[i]]=c(NA,NA); } }
        vals<-matrix(unlist(x),ncol=2,byrow=T,)[,2]
        gff3[,tag]<-vals
      }
    } else {
      m2 <- as.data.frame(matrix(nrow=0,ncol=length(tags)))
      colnames(m2) <- tags
      gff3 <- cbind(gff3,m2)
    }
  }
  return(gff3)
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

############################################################
# RPKM/FPKM 
# TUQ=True then it returns UQ-FPKM = UQ(FPKM) 
# UQ=TUR then it returns FPKM-UQ (as described in PCAWG)
countstable2rpkms <- function(table,lens,mass.labels=NULL,exitonerror=TRUE,UQ=FALSE,TUQ=FALSE) {
  # check if there missing features
    missing.feat <- (!rownames(table) %in% names(lens))
    
    if ( sum(missing.feat) ) {
        message("ERROR: Length of ",paste(rownames(table)[missing.feat],sep=",")," not found.")
        if (exitonerror) { q(status=1) }
        return(NULL)
    }
    v.compute.rpkm <-  function(l,lens,mass.labels=NULL,UQ=FALSE) {
        #(l*1e6)/(sum(l)*lens[names(l)]/1000)        
        if ( is.null(mass.labels) ) {
            mass.vals <- l
        } else {
            mass.vals <- l[mass.labels]
        }
        # UQ?
        if ( UQ ) {# 
            mass.vals.no.zero <- mass.vals[mass.vals>0]
            x <- summary(mass.vals.no.zero)[5]
            #mass.vals <- mass.vals.no.zero[mass.vals.no.zero>=x]
            mass.vals <- x*length(mass.labels)
            #message("UQ:",x,"-",sum(mass.vals.no.zero),"--->",sum(mass.vals))
        }
        tot.mass <- sum(mass.vals)     
        if ( tot.mass == 0 ) {
            print("Tot.mass==0!!!!")
            tot.mass <- 1
        }
        l <- as.integer(l)
        #message("Tot.mass:",tot.mass)
        tot.mass <- as.integer(tot.mass)
        expr <- (10^9/tot.mass*l/lens)
        if (TUQ ) {            
            mass.vals.no.zero <- expr[expr>0]
            x <- summary(mass.vals.no.zero)[5]
            expr <- round(expr/x,2)
            #message("TUQ:",x)
        }
        #return((10^9*l)/(tot.mass*lens))
        return(expr)
    }
    # 
    if ( is.vector(table) ) {
        stopifnot(sum(!names(table)%in%names(lens))==0)
        # fix ordering & and convert to numeric
        lens <- as.numeric(lens[names(table)])
        res <- round(v.compute.rpkm(table,lens,mass.labels,UQ),2)
        names(res) <- names(table)
        return(res)
    } else {
        stopifnot(sum(!rownames(table)%in%names(lens))==0)
        # fix ordering & and convert to numeric
        lens <- as.numeric(lens[rownames(table)])
        res <- round(apply(table,2,v.compute.rpkm,lens,mass.labels,UQ),2)
        rownames(res) <- rownames(table)
        return(res)
    }
}

####################################################################
# TPMs/GPM
# RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome
# BMC Bioinformatics 2011, 12:323  doi:10.1186/1471-2105-12-323
# http://arxiv.org/abs/1104.3889 - TPM=(FPKM/SUM(FPKM))*10^6
# TODO: call the FPKM function
countstable2tpm <- function(table,lens,exitonerror=TRUE) {
  # check if there missing features
  missing.feat <- (!rownames(table) %in% names(lens))

  if ( sum(missing.feat) ) {
    perror("Length of ",paste(rownames(table)[missing.feat],sep=",")," not found.")
    if (exitonerror) { q(status=1) }
    return(NULL)
  }
  v.compute.tpm <-  function(l,lens) {
    ta <- l/lens[names(l)]
    tap <- ta/sum(ta)
    return(tap*10^6)
  }
  if ( is.vector(table) ) {
    return(round(v.compute.tpm(table,lens),2))
  } else {
    return(round(apply(table,2,v.compute.tpm,lens),2))
  }
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
# Assumption: gene ids are ensembl ids
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
    species <- "Homo Sapiens"
    library("org.Hs.eg.db")
    go.db <- org.Hs.egGO
    #pfam.db <- org.Hs.egPFAM
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
      lgn.db <- org.EcK12.egGENENAME
      symbol.db <- org.EcK12.egSYMBOL
      # how to get the entrez ids from the ensembl gene id?
      # this will not until we can map the ensembl ids
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
      lgn.db <- org.Ss.egGENENAME
      symbol.db <- org.Ss.egSYMBOL
      ensembl.db <-org.Ss.egENSEMBL
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
    loc <- unlist(strsplit(x=l[2],split="(\\.\\.|\\-)"))
    as.numeric(loc[1])
  }
  get.gene.len <- function(s) {
    l <- unlist(strsplit(x=s,split=":",fixed=T))
    loc <- unlist(strsplit(x=l[2],split="(\\.\\.|\\-)"))
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
# gets the internal ids associated to ensgids
# 
annot.get.egi <- function(gids,dbs) {
  library("GO.db")
  id2egi<-ls(dbs$symbol.db)
  names(id2egi) <- ls(revmap(dbs$symbol.db))
  egi <- id2egi[gids]
  #pdebug(gid,"---",egi)
  if (is.null(egi) || sum(!is.na(egi))==0) {
      id2egi<-ls(dbs$ens.db)
      names(id2egi) <- ls(revmap(dbs$ens.db))
      egi <- id2egi[gids]      
  }
  #pick the first  (ensembl mapping may return mul. entries)
  if (is.null(egi) || sum(!is.na(egi))==0 ) {
    return(NA)
  }
  return(egi) 
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
  if (add.label) return(paste0(lab," ",x))
  return(x)
}

brew.wrapper <- function(...) {
  x <- try(brew(...,envir=parent.frame()))
  if ( inherits(x,"try-error") ) {
    perror("Failed to generate HTML.",attr(x,"condition"))
    if (pdebug.enabled) {
      code <- brew(...,envir=parent.frame(),run=FALSE,parseCode=FALSE)
      cat(paste(unlist(code),collapse="\n",sep="\n"))
    }
    q(status=2)
  }
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
    p <- paste0(p,".html")
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
                     res=300,
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
    dir <- paste0(dir,"/")
  }  
  # automatic filename  
  if ( is.null(filename)) {
    time.label <- format(Sys.time(), "%d%m%Y%H%M%S")
    filename <- paste0("graph_",time.label,".png")
  }
  ps.file <- NA
  data.file <- NA
  png.file <- paste0(dir,filename,".png")
  png.file <- sub(".png.png$",".png",png.file)
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
  png(filename=png.file,width=width.in, height=height.in,bg=bg,res=res,units="in")
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
    ps.file <- paste0(png.file,".eps")
    #ps(file=ps.file)
    postscript(file=ps.file,fonts=c("serif", "sans"))
    err <- try(to.plot())
    dev.off()
    if ( class(err)=="try-error") {
      pwarning("Failed to generate plot ",png.file)
      return(NULL);
    }
  }
  # DATA
  if ( sum(!is.na(data.table)) > 0 ) {
    data.file <- paste0(png.file,".tsv")
    write.tsv(data.table,data.file)
  }
  png.file
  # HTML
  if ( html == TRUE ) {
    html <- paste0("<DIV class=\"dplot\"><IMG src='",basename(png.file),"' border =0 width=",width," height=",height,"><BR/>")
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
    html <- paste0(html,html.download.bar(files),"</DIV>")
    html <- paste0(html,"<p class='caption'>",caption,"</p>")
  } else {
    HTML=FALSE
  }
  return(list(png.file=png.file,html=html,data.file=data.file,ps.file=ps.file))
}

html.download.bar <- function(files) {

  if( length(files) ==0 ) { return("") }
  html <- "<DIV class=\"downbar\"> Download: "
  
  for ( f in names(files) ) {
    html <- paste0(html,"&nbsp;<a class='download' href='",basename(files[f]),"'>",f,"</a>")
  }
  html <- paste0(html,"</DIV>")
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
    file <- paste0("graph_",time.label,".png")
  }
  plot.filename <- paste(dir,filename,sep="/")
  png(filename=plot.filename,width=width, height=height,bg=bg)
  to.plot()
  dev.off()
  pdf(file=paste0(plot.filename,".pdf"),width=width, height=height)
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
  HTML( paste0("Project name: ", name) )
  HTML( paste0("Species: ", species) )

  contrasts <- conf.get.value(conf,"contrasts")
  if ( is.null(contrasts) ) {
    # default (deprecated)
    #conditions <- conf.get.value(conf,"conditions")
    #if (!is.null(conditions)) {
    #  contrasts <- "conditions"
    #} else {
      contrasts <- NULL
    #}
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
  HTML(paste0("IRAP ",irap_version))
  HTMLEndFile()
}
############### 
div.ids <- 0
HTML.toogle.div <- function(visible.text, hidden.text) {
  assign("div.ids",div.ids+1,envir = .GlobalEnv)
  idon  <- paste0("t",div.ids,"on")
  idoff <- paste0("t",div.ids,"off")
  paste0(
        "<a class=\"button\" id='",idoff,"' href='javascript:toggle(\"",idon,"\",\"",idoff,"\");'>",visible.text,"</a>",
        "<div class=\"button\" id='",idon,"' style='display: none'>",hidden.text,"</div>")
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
write.tsv <- function(x,file,header=TRUE,rownames.label=NULL,fix=TRUE) {
  #
  if (!is.null(x)) {
    if (!is.matrix(x)) {
      x <- as.matrix(x)
    }
    if (fix) {
      x <- apply(x,c(1,2),gsub,pattern="'",replacement="")
    }
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
#read.table(file,sep = "\t", header=header, quote = "\"",comment.char="",check.names=FALSE)
  return(qload.tsv(file,header))
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
  #pinfo("merge ",table.field," ",annot.field)
  #pinfo(colnames(table))
  #pinfo(colnames(annot))
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
  avg.matrix <- data.frame(data.matrix[,c(1:length(groups))])
  colnames(avg.matrix) <- groups
  rownames(avg.matrix) <- rownames(data.matrix)
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
  #if ( file.exists(cached.annot) ) {
  #  load(cached.annot)
  #  if ( exists("gene.annot")) {
  #    annot.table <- gene.annot
  #  }
  #} else {
  #
  if (!file.exists(file) ) {
    perror("File ",file," not found")
    quit(status=1)
  }
  annot.table <- tryCatch(read.tsv(file,header=TRUE),error=function(x) return(NULL))
  if ( is.null(annot.table) ) {
    pdebug("Loading annotation (failed)")
    return(NULL)
  }       
  #  save(list=c("annot.table"),file=cached.annot)
  #}
  pinfo("Load annot...")
  if ( nrow(annot.table) == 0  ) {
    pdebug("Loading annotation (done) - empty file")
    annot.table <- NULL
  } else {
    annot.table <- annot.expand.fields(annot.table)
    pdebug("Loading annotation (done)")
  }
  return(annot.table)
}
# keep backwards compatibility by using read.table when data.table is not 
# available
qload.tsv <- function(f,header=NULL,comment.char="",nrows=-1L) {
  tsv.data <- NULL
  if (require("data.table",quietly=TRUE,character.only=TRUE) &&
      compareVersion(as.character(packageVersion("data.table")),"1.9.6")>=0) {
    library("data.table")
    if ( sum(grep(".gz$",f)) ) {
      f <- paste("zcat ",f,sep="")
    } else {
      f <- paste("cat ",f,sep="")
    }
    # not optimal, but faster than read.table
    if ( comment.char!="") {
      f <- paste(f," | grep -v \"^",comment.char,"\"",sep="")
    }
    if (is.null(header)) header="auto"
    tryCatch(tsv.data <- fread(input=f,sep = "\t", nrows=nrows, header=header,check.names=FALSE,data.table=FALSE),error=function(x) NULL)
  } else 
    tryCatch(tsv.data <- read.table(f,sep = "\t", header=header, comment.char=comment.char, quote = "\"", nrows=nrows, check.names=FALSE),error=function(x) NULL)
  return(tsv.data)
}

# load a file with a quant. matrix
# returns NULL in case of failure
quant.load <- function(f,clean.cuff=FALSE) {
  tsv.data <- NULL

  tsv.data <- qload.tsv(f,header=FALSE)
  if ( !is.null(tsv.data) && ncol(tsv.data)>1 ) {
    if ( sum(grepl("(Gene|Exon|Transcript|ID|feature)",tsv.data[1,1],ignore.case=T))!=0 || tsv.data[1,1]=="") {
      # reload to include the header
      tsv.data <- qload.tsv(f,header=TRUE)
    }
  }
  if(is.null(tsv.data)) return(NULL);
  rownames(tsv.data) <- as.character(tsv.data[,1])
  tsv.data <- tsv.data[,-1,drop=FALSE]
  if (clean.cuff) {
    sel<-grep("^CUFF.*",rownames(tsv.data),perl=T,invert=TRUE)
    tsv.data <- tsv.data[sel,,drop=FALSE]
  }
  return(tsv.data)
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

quant.heatmap <- function(data,nlim=100,key.xlabel="Expression",do.cor=FALSE,density.info="histogram",cexRow=0.6,cexCol=0.9,...) {
  suppressPackageStartupMessages(library("gplots"))
  #data <- log2(quant.data.group$mean+1)
  var <- rowVariance(data)
  ## if (is.na(nlim)) {
  ##   ncol <- max(100,nrow(data))
  ## } else {
  ##   ncol <- min(nlim,ncol(data))
  ## }
  if (is.na(nlim) && do.cor==FALSE) {
    #ncol <- nlim
    if (is.na(nlim) || (nlim==nrow(data) && nlim>500)) {    
      nlim <- 500
    }
  }
  # select a subset of rows
  if (!is.na(nlim)) {
    select <- order(var,decreasing=TRUE)[c(1:nlim)]
    data <- data[select,,drop=FALSE]
  }
  if ( do.cor ) {
    data <- cor(data)
    key.xlabel <- paste0(key.xlabel," correlation")
  }
  colors <- topo.colors(ncol(data))

  # sort the genes by variability
  suppressWarnings(irap.heatmap.2(x=as.matrix(data),col = colors, scale = "none",cexCol=cexCol,cexRow=cexRow,keysize=1,density.info=density.info,trace="none",key.xlabel=key.xlabel,...))
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

plot.panel.label <- function(label,add=")",x=2.5,cex=1.5,col="black",adj=2) {
  if ( par("ylog") ) {
    at=10^par("usr")[4]
  } else {
    at=par("usr")[4]
  }
  mtext(paste0(label,add), side=2, at=at,line=x,cex=cex,las=2,col=col,adj=adj,padj=0)
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
      cat(paste0("ERROR: ",file," not found\n"))
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

myParseArgs <- function(usage,option_list,filenames.exist=NULL,multiple.options=NULL,mandatory=NULL,...) {

  # get command line options, if help option encountered print help and exit,
  # otherwise if options not found on command line then set defaults,
  parser <- OptionParser(usage = usage, option_list=option_list)
  opt <- parse_args(parser,...)
  
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
#############################################################
conf.get.value <- function(conf,val,get.default=FALSE) {
  #pinfo(">>>>>>>>>>>",val)
  idx <- grep(paste0("^",val,"$"),conf[,1],ignore.case=TRUE,value=FALSE,perl=TRUE)
  
  #pinfo("AAAAAAAAAAAAAAA:",val,idx,is.null(idx),typeof(idx))
  if (is.null(idx) || length(idx)==0) {
    #
    if ( length(grep(".*_strand$",val,ignore.case=TRUE,value=FALSE,perl=TRUE))>0 ) {
      return(NA)
    }
    if ( get.default ) {
      return(conf.get.default.value(val))
    }
    return(NULL)
  }
  v <- myread.list(conf[idx,2])
  if ( sum(is.null(v))>0 && get.default ) {
    v <- conf.get.default.value(val)
  }
  return(v)
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
  cmd <- paste0("grep -i '^def_",var,"' ",IRAP.DIR,"/{scripts/irap,aux/mk/*.mk} | cut -f 2 -d= |head -n 1")
  #pinfo(cmd)
  val <- system(cmd,intern=TRUE)
  if ( is.null(val)) { return(NULL); }
  if ( length(val)==0 ) { return(NULL); }
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

tech.replicates2list <- function(val) {
  if (is.null(val)) {
    return(val)
  }
  l <- list()  
  for ( s in strsplit(val,";")[[1]]) {
    l[[paste(",",s,",",sep="")]] <- strsplit(s,",")[[1]]
  }
  return(l)
}
groups2matrix <- function(exp.conf) {

  groups.val <- conf.get.value(exp.conf,"groups")
  if ( is.null(groups.val) ) {
    # get the groups used in the contrasts
    contrasts.val <- conf.get.value(exp.conf,"contrasts")
    if ( is.null(contrasts.val) ) {
      return(NULL)
    }
    # gather all groups from the contrasts
    groups.val <- unique(unlist(lapply(contrasts.val,FUN=conf.get.value,conf=exp.conf)))    
  }

  libs.ids <- libraries2ids(exp.conf)
  libs <- names(libs.ids)
  m <- matrix(ncol=length(libs),nrow=length(groups.val),dimnames=list(groups.val,libs))
  #
  for ( g in groups.val ) {
    group.def <- conf.get.value(exp.conf,g)
    sel <- colnames(m) %in% group.def
    m[g,] <- libs.ids[colnames(m)]
    m[g,!sel] <- NA
  }  
  m
}
libraries2ids <- function(exp.conf) {
                                        # all libraries
  libs <- sort(unique(append(conf.get.value(exp.conf,"se"),conf.get.value(exp.conf,"pe"))))
     # technical replicates
  tech.replicates.l <-   tech.replicates2list(conf.get.value(exp.conf,"technical.replicates"))
  if ( is.null(tech.replicates.l)) {
    libs.ids <- seq(1,length(libs))
    names(libs.ids) <- libs    
  } else {
    x <- seq(1,length(names(tech.replicates.l)))  
    xtimes <- unlist(lapply(tech.replicates.l,length))
                                        #print(xtimes)
    libs.ids <- rep(x,times=xtimes)
    names(libs.ids) <- unlist(tech.replicates.l)    
  }
  return(libs.ids)
}

irap.get.palette <- function(n) {
  suppressPackageStartupMessages(library(gplots))
  suppressPackageStartupMessages(library(RColorBrewer))
  #
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  return(getPalette(n))
}

ids2colours.mapping <- function(libs.ids) {
  #mapping <- gsub("..$","",irap.get.palette(max(libs.ids)))
  mapping <- irap.get.palette(max(libs.ids))
  names(mapping) <- seq(1,length(mapping))  
  return(mapping)
}
libids2colours <- function(libs.ids) {

  all.lib.colours <- ids2colours.mapping(libs.ids)
  libs.colours <- sapply(libs.ids,FUN=table.get.colour,palette=all.lib.colours)
  names(libs.colours) <- names(libs.ids)
  return(libs.colours)
}

table.get.colour <- function(val,palette) {
  if (is.na(val) || is.null(NA) ) { return(NA); }
  return(palette[as.numeric(as.character(val))])
}

#
#table <- groups.m.html
#table.col <- groups.m.colours
#table <- contrasts.m.html
#table.col <- contrasts.m.colours
table.set.bgcolor <- function(table,table.col) {
  if ( is.matrix(table) || is.data.frame(table) ) {
    table.v <- as.vector(apply(table,c(1,2),as.character))
    table.col.v <- as.vector(table.col)
  } else {
    table.v <- table
    table.col.v <- table.col
  }
  irap.assert(length(table.v)==length(table.col.v))
  # html
  html.v <- paste0("<div style=\"border: 0px; padding: 0px; background-color: ",table.col.v,";\">",table.v,"</div>")
  # exclude NA
  html.v[is.na(table.v)] <- NA
  # rebuild the matrix
  if ( is.matrix(table) || is.data.frame(table) ) {
    newm <- matrix(html.v,ncol=ncol(table),nrow=nrow(table),byrow=F)
    rownames(newm) <- rownames(table)
    colnames(newm) <- colnames(table)
  } else {
    newm <- html.v
    names(newm) <- names(table)
  }
  return(newm)
}

vars2table.html <- function(conf,vars,table.label="",colnames2add=NULL) {
  values <- sapply(vars,conf.get.value,conf=conf,get.default=TRUE)
  #vars.with.vals <- sapply(vars[!unlist(lapply(values,is.null))],capitalize)
  vars.with.vals <- vars[!unlist(lapply(values,is.null))]
  #pinfo(length(append(vars.with.vals,unlist(values))))
  vals <- matrix(append(vars.with.vals,unlist(values)),byrow=F,ncol=2,dimnames=NULL)
  #vals[,1] <- gsub("_"," ",vals[,1])
  if ( !is.null(colnames2add)) {
    vals <- rbind(colnames2add,vals)
    rownames(vals) <- NULL
  }
  vals <- as.matrix(vals)
  xt <- xtable(vals,
               caption="",
               label=table.label)
  
  html <- print.xtable(print.results=FALSE,
                       xt,
                       type="html",
                       sanitize.text.function=function(x) {x},
                       include.colnames=FALSE,
                       include.rownames=FALSE,
                       html.table.attributes = " class='expinfo'")
  return(html)
}

lib.info <- function(exp.conf,lib) {
  # check if we have any data for the library
  if ( is.null(conf.get.value(exp.conf,lib)) ) {
    perror("Unable to find library ",lib," in configuration file.")
    q(status=1)
  }

  info <- conf.get.value(exp.conf,paste0(lib,"_info"))
  if ( is.null(info) || is.na(info)) {
    info <- ""
  }
  v <- c(lib,
         paste(conf.get.value(exp.conf,lib),collapse=",",sep=","),
         ifelse(lib %in% conf.get.value(exp.conf,"pe"),"PE","SE"),
         conf.get.value(exp.conf,paste(lib,"_rs",sep="")),
         conf.get.value(exp.conf,paste(lib,"_qual",sep="")),
         info)
  names(v) <- c("Name","File(s)","Type","Read Size","Quality encoding","Notes/Info")
  return(v)
}

contrasts2matrix <- function(exp.conf) {

   # get the groups used in the contrasts
  contrasts.val <- conf.get.value(exp.conf,"contrasts")
  if ( is.null(contrasts.val) ) {
    return(NULL)
  }
  # gather all groups from the contrasts  
  groups.val <- unique(unlist(lapply(contrasts.val,FUN=conf.get.value,conf=exp.conf)))    
  m <- as.data.frame(matrix(ncol=length(groups.val)+1,byrow=T,nrow=length(contrasts.val),dimnames=list(contrasts.val,append("Contrast",groups.val))))
  #
  groups.ids <- seq(1,length(groups.val))
  names(groups.ids) <- groups.val
  #v <- contrasts.val[2]
  for ( v in contrasts.val ) {
    cont.def <- conf.get.value(exp.conf,v)
    sel <- names(groups.ids) %in% cont.def
    m[v,] <- append(paste(cont.def,collapse=" - "),groups.ids[colnames(m)[-1]])
    m[v,!append(TRUE,sel)] <- NA
  }  
  return(m)
}

#
load.configuration.file <- function(conf_file) {

  # 
  pinfo("Loading configuration file...")
  conf.table <- read.delim(conf_file,sep="=",comment.char="#",header=FALSE,stringsAsFactors=FALSE)
  pinfo("Loading conf_file...done")
  pinfo("Configuration:")
  line2include <- c()
  for ( i in 1:nrow(conf.table) ) {
    if ( mytrim(as.character(conf.table[i,1])) != "" ) {
      assign(as.character(conf.table[i,1]),mytrim(conf.table[i,2]))
      pinfo("",conf.table[i,1],"=",conf.table[i,2])
      line2include <- append(line2include,i)
    } 
  }
  # exclude empty lines
  conf.table <- conf.table[line2include,,drop=FALSE]
  pinfo("Loading default parameters...")
  ###################
  # handle duplicates
  conf.table <- conf.table[!duplicated(conf.table[,1]),]
  #################
  # Default values
  vars <- c("mapper","quant_method","de_method","qual_filtering","gse_tool")
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
  "tRNA"=NA,
  "rRNA"=NA,
  "lincRNA"=NA,
  "retained intron"=NA
  )
# Ids of the genes on each group
source.filt.groups <- list("all"=TRUE,
                           "pseudogenes"=NA,
                           "protein coding"=NA,
                           "tRNA"=NA,
                           "rRNA"=NA,
                           "retained intron"=NA,
                           "lincRNA"=NA
                    )

source.filt.query <- list("all"=TRUE,
                          "pseudogenes"=c("pseudogene"),
                          "protein coding"=c("protein_coding"),
                          "tRNA"=c("tRNA"),
                          "retained intron"=c("retained_intron"),
                          "lincRNA"=c("lincRNA")
                          )
source.column <- NULL
#                          "rRNA"=c("rRNA"),
biotype.column <- function(table) {
  if ( sum("source" %in% colnames(table))>0 ) {
    sources <- unique(as.character(table$source))
    if ( sum(grepl("HAVANA",sources,ignore.case=FALSE))>0 ) {
      pinfo("cols:",colnames(table))
      if ( sum("gene_type" %in% colnames(table))>0 ) {
        return("gene_type")
      }
      # 
      return("gene_biotype")        
    } else {
      if ( sum(grepl("(havana|ena|WormBase|irgsp|VectorBase|bgi|ensembl|jgi|tair|JCVI|brad|cmer|FlyBase|ibsc|oge|ensembl_genomes|GFBGP|cirad|pgsb|PomBase|iwgsc|gramene
excepds|pythiumgenomedatabase|Ensembl_Fungi)",sources,ignore.case=FALSE))>0 ) {
        # more recent ensembl files do not have the biotype in the source column
        return("gene_biotype")
      }
    }
  }
  return("source")
}

# Initialization
init.source.filter <- function(table,id.col="ID") {
  # define the source column
  source.column <- "biotype"
  pinfo("biotype column=",source.column)
  pinfo(colnames(table))
  if ( source.column %in% colnames(table) ) {
    source.filt.groups.def <- list()
    source.filt.groups <- list()
    # slower
    sources <- unique(as.character(table[,source.column]))
    
    pinfo("#Sources/biotypes:",length(sources))
    pinfo(sources)
    # all
    source.filt.groups.def$"all" <- sources
    source.filt.groups$"all" <- as.character(table[,id.col])

    for (f in names(source.filt.query) ) {
      if (f!="all") {
        q <- source.filt.query[f]
        qr <- grep(q,sources,value=T)
        #pinfo("Filter: ",f," ",q,"->",paste(qr,collapse="|"))
        source.filt.groups.def[[f]] <- qr
        re <- paste0("(",paste(qr,collapse="|"),")")
        filt <- NULL
        if (length(qr)!=0) 
          filt <- sapply(table[,source.column],grepl,pattern=re)
        source.filt.groups[[f]] <- as.character(table[filt,id.col])
        pinfo("Filter ",f,": found ",length(source.filt.groups[[f]]))
      }
    }
    ## pseudogenes <- grep("pseudogene",sources,value=T)
    ## source.filt.groups.def$pseudogenes <- pseudogenes
    ## pseudogenes.filt <- sapply(table$source,grepl,grepl,pattern=paste("(",paste(pseudogenes,collapse="|"),")",sep=""))
    ## source.filt.groups$pseudogenes <- as.character(table[pseudogenes.filt,id.col])
    ## if (length(pseudogenes.filt)>0) {
    ##   pdebug(sum(pseudogenes.filt))
    ## }
    
    #xRNA <- grep("RNA$",sources,value=T)
    #source.filt.groups.def$"xRNA" <- as.character(table[xRNA,id.col])
    #xRNA.filt <- sapply(table$source,grepl,grepl,perl=TRUE,pattern=paste("(",paste(xRNA,collapse="|"),")",sep=""))
    #source.filt.groups$"xRNA" <- as.character(table[xRNA.filt,id.col])

    ## tRNA <- grep("tRNA$",sources,value=T)
    ## source.filt.groups.def$"tRNA" <- as.character(table[tRNA,id.col])
    ## tRNA.filt <- sapply(table$source,grepl,grepl,perl=TRUE,pattern=paste("(",paste(tRNA,collapse="|"),")",sep=""))
    ## source.filt.groups$"tRNA" <- as.character(table[tRNA.filt,id.col])

    # update the global variable
    assign("source.filt.groups.def",source.filt.groups.def, envir = .GlobalEnv)
    assign("source.filt.groups",source.filt.groups, envir = .GlobalEnv)
  }
  return(T)
}
#
apply.source.filter <- function(table,filter.name) {

  ids <- source.filt.groups[[filter.name]]
  if (is.matrix(table) || is.data.frame(table) ) {
    table <- table[rownames(table) %in% ids,]
  } else {
    if (is.vector(table)) {
      table <- table[names(table) %in% ids]
    }
  }
  return(table)
}
#
get.source.filename <- function(file.prefix,type) {
  if ( type =="all" || type=="All" ) {
    return(gsub(" ","_",paste0(file.prefix,".html")))
  }
  gsub(" ","_",paste0(file.prefix,"_",type,".html"))
}
# outprefix -> menu html
get.source.filter.menu <- function(file.prefix) {
   get.source.menu.entry <- function(name,file.prefix) {
     paste0("<a href='",basename(get.source.filename(file.prefix,name)),"'>",name,"</a>")
   }
   sources.menu <- paste(sapply(sort(names(source.filt.groups.def)),get.source.menu.entry,file.prefix),collapse="|")
   sources.menu
 }

get.source.filter.names <- function() {
  names(source.filt.groups.def)
}

get.source.filter.size <- function(filter.name) {
  length(source.filt.groups[[filter.name]])
}


boxplot.n <- function(significance=NULL,comparisons=NULL,show.n=TRUE,adj.n=0,...) {

  if (!is.null(significance) ) {
    significance[is.na(significance)] <- ""
    significance[significance!=""] <- paste("p=",significance[significance!=""],sep="")
  }
  bp <- boxplot(...)
  if (show.n) {
    text(seq(1,ncol(bp$stats)),bp$stats[5,]+strheight("T")*2+adj.n,
         paste("n=",bp$n,""),
         col="darkgrey")
  }
  if (!is.null(comparisons)) {
    bars <- unique(unlist(comparisons))
    min.height <- bp$stats[5,bars]+strheight("T")*2+adj.n
  # adjust height to avoid overlaps
    heights.l <- rep(NA,length(bars))
    min.h <- max(min.height)
    for (c in seq(1,length(comparisons)) ) {
      heights.l[c] <- min.h+c*strheight("T")*2+adj.n
    }
    for (c in seq(1,length(comparisons)) ) {
      x1 <- comparisons[[c]]
      x2 <- c(heights.l[c],heights.l[c])
      txt <- significance[c]
      segments(x1[1],x2[1],x1[2],x2[2],col="gray")
      segments(x1[1],x2[1],x1[1],x2[1]-strheight("T")/2,col="gray")
      segments(x1[2],x2[2],x1[2],x2[2]-strheight("T")/2,col="gray")
      text((x1[1]+x1[2])/2,x2[2]+strheight("T")/2,txt,col="gray")
    }
  }

  bp
}

###############################################
# 
fisherNetworkPlot <- function (gsaRes,
                               adjusted = FALSE, significance = 0.001,
                               overlap = 1,
                               long.names=NULL,
                               label = "names", cexLabel = 0.9, 
                               ncharLabel = 15, cexLegend = 1,
                               nodeSize = c(10, 40), edgeWidth = c(1,15),
                               edgeColor = NULL, scoreColors = NULL, main) 
{
  library(igraph)
  library(marray)
  
  if (overlap <= 0) 
    stop("argument overlap has to be at least 1")
  if (length(nodeSize) != 2) 
    stop("argument nodeSize has to have length 2")
  if (length(edgeWidth) != 2) 
    stop("argument edgeWidth has to have length 2")
  if (class(adjusted) != "logical") 
    stop("argument adjusted has to be TRUE or FALSE")
  if (missing(main))  {
    main <- paste("Exact Fisher-test (significance <=",significance,")",sep="")
  } else {    
    if (class(main) != "character") 
      stop("argument main has to be a character string")
  }                              #
  if (!all(c("pvalues", "p.adj", "gsc","effect.size","sSizes") %in% names(gsaRes)))
    stop("gsaRes needs to contain pvalues, p.adj, gsc, effect.size and sSizes")

  #
  geneSetNames <- rownames(gsaRes$resTab)
  gsc <- gsaRes$gsc[geneSetNames]

  # number of genes selected (DE) in each gene set
  sSizes <- gsaRes$sSizes[geneSetNames]
  effect.size <- gsaRes$effect.size[geneSetNames]
  # use p-values or adjusted p-values
  if (adjusted) {
    pValues <- gsaRes$p.adj[geneSetNames]
  } else {
    pValues <- gsaRes$pvalues[geneSetNames]
  }
                                        #significance=1

                                        # which genes (index) have a significant p-value
  indSignificant <- which(abs(pValues) <= significance)
  
  if (length(indSignificant) < 1) {
    stop("less than one gene sets were selected, cannot plot (tip: adjust the significance cutoff)")
    return(NULL)
  }
  pSignificant <- pValues[indSignificant]
  # overlap  genes between gene sets
  overlapMat <- matrix(nrow = length(indSignificant), ncol = length(indSignificant))
  for (i in 1:nrow(overlapMat)) {
    for (j in i:ncol(overlapMat)) {
      tmp <- sum(gsc[[indSignificant[i]]] %in% gsc[[indSignificant[j]]])
      overlapMat[i, j] <- tmp
            overlapMat[j, i] <- tmp
    }
  }
  # size = number of genes in a gene set
  gsSize <- diag(overlapMat)
  overlapMat[overlapMat < overlap] <- 0    
  adjMat <- overlapMat > 0
  tmp <- adjMat
  diag(tmp) <- 0
  if (all(!tmp))  {
    warning("no overlap between gene sets found, try to decrease argument overlap or increase argument significance...fixing it for now")
    diag(tmp) <- 1
    return(NULL)
  }
  g <- graph.adjacency(tmp, mode = "undirected", diag = FALSE)
  
  edgeOverlap <- rep(NA, ecount(g))
  if (ecount(g)>0) {
    for (iEdge in 1:ecount(g)) {
      tmp <- get.edge(g, iEdge)
      edgeOverlap[iEdge] <- overlapMat[tmp[1], tmp[2]]
    }
      # edge width
    eWidth <- (edgeOverlap - min(edgeOverlap))/(max(edgeOverlap) - 
                                                min(edgeOverlap)) * (edgeWidth[2] - edgeWidth[1]) + edgeWidth[1]
    
    if (is.null(edgeColor)) 
      edgeColor = c("gray80", "gray70", "gray60")
    
    tmp <- seq(min(edgeOverlap), max(edgeOverlap), length.out = length(edgeColor) +  1)
    eColor <- rep(edgeColor[1], ecount(g))
    for (i in 2:length(edgeColor)) {
      eColor[edgeOverlap > tmp[i]] <- edgeColor[i]
    }
  } else {
    eColor=NULL
    eWidth=10
    tmp <- NULL
  }
                                        # Node size - proportional to the size of the subsets of gene sets
                                        #nodeSize <- c(3,20)
  vSize <- (gsSize - min(gsSize))/(max(gsSize) - min(gsSize)) * 
    (nodeSize[2] - nodeSize[1]) + nodeSize[1]

  if ( length(vSize)==1 && is.na(vSize) ) {
    # max(gsSize=min)
    vSize=50
  }
  
  effect.size2 <- effect.size[names(gsc)[indSignificant]]
  eSize <- ( - min(effect.size2))/(max(effect.size2) - min(effect.size2)) * 
    (nodeSize[2] - nodeSize[1]) + nodeSize[1]
  
  if (is.null(scoreColors)) {
    gradColors <- maPalette(low="white",mid="yellow",high="red",k=100)
    vColor <- rep(NA,length(effect.size2))
    tmp <- effect.size2>1
    if ( sum(tmp) > 0 )
      vColor[tmp] <- gradColors[floor(effect.size2[tmp] * 99/max(effect.size2[tmp],na.rm=T)) + 1]
    tmp <- effect.size2<1
    if (sum(tmp) > 0 ) 
      vColor[tmp] <- gradColors[1]
  } else {
      if (length(scoreColors) < 2) 
        stop("argument scoreColors should contain at least two colors")
      tmp <- scoreColors
      gradColors <- colorRampPalette(tmp, interpolate = "linear")(100)
      vColor <- gradColors[floor(effect.size2 * 99/max(effect.size2, 
                                                      na.rm = TRUE)) + 1]
    }
  vColor[is.na(vColor)] <- "#CCCCCC" 
  
                                        #label <- "numbers"
  tmp <- try(label <- match.arg(label, c("names", "numbers", 
                                         "namesAndLong"), several.ok = FALSE), 
             silent = TRUE)
  if (class(tmp) == "try-error") {
    stop("argument label has to be set to either 'names' or 'numbers'")    
  }
  tmp <- names(gsc)[indSignificant]
                                        #ncharLabel <- 10
  for (i in 1:length(tmp)) {
    if (nchar(tmp[i]) > ncharLabel) 
        tmp[i] <- paste(substr(tmp[i], 1, ncharLabel), "...", 
                        sep = "")
  }
                                        # Labels
  if (label == "names")  {
    vLabels <- tmp
  } else {
    if (label == "numbers") 
        vLabels <- 1:length(indSignificant)
    else {
      if (label == "namesAndLongNames") 
        vLabels <- tmp
      else {
        if (label == "namesAndSizes") 
          vLabels <- paste(tmp, " (", gsSize, ")", sep = "")
        }
    }
  }
  
                                        # keep code from piano
  if ( ecount(g)>0 ) {
    eWeight <- rep(NA, ecount(g))
      for (iEdge in 1:ecount(g)) {
        tmp <- get.edge(g, iEdge)
        tmp1 <- gsSize[tmp[1]]
        tmp2 <- gsSize[tmp[2]]
        if (tmp1 > tmp2) {
          eWeight[iEdge] <- 1/(tmp1 * sum(adjMat[, tmp[1]]))
        }
        else {
          eWeight[iEdge] <- 1/(tmp2 * sum(adjMat[, tmp[2]]))
        }
      }
    if (min(eWeight) != max(eWeight)) {
      eWeight <- (eWeight - min(eWeight))/(max(eWeight) - 
                                           min(eWeight)) * (100 - 0) + 0
    } else {
      eWeight <- rep(50, length(eWeight))
    }
  } else {
    eWeight <- NULL
  }
  lay <- layout.fruchterman.reingold(g, area = vcount(g)^4, 
                                     repulserad = vcount(g)^5,
                                     weights = eWeight)
  opar <- par()
  par("mar"=c(3,3,3,3))
  if (!is.null(long.names)) {
    layout(matrix(c(1, 0, 3, 2, 2, 3), ncol = 2),
           widths = c(1,7),
           heights = c(2, 3, 2) )
  } else {
    layout(matrix(c(1, 3, 2, 2 ), ncol = 2), widths = c(1,5),
           heights = c(1, 1 ))
  }
  #
  #
  #layout.show(3)
  #dev.off()
  
  maColorBar(seq(1, max(effect.size2, na.rm = TRUE),
                 by = max(effect.size2, na.rm = TRUE)/100),
             horizontal=FALSE, gradColors,
             k = 7, main = "Observed/\nExpected",
             cex.main = cexLegend)
               
  plot(g, vertex.size = vSize, vertex.color = vColor,
       vertex.label = as.character(vLabels), 
       vertex.label.family = "sans", vertex.label.color = "black", 
       vertex.label.cex = 1, edge.width = eWidth, edge.color = eColor,
       layout = lay, main = main)
  #long.names <- vLabels
  #long.names[1] <- "asdadajhljhasd fusdhflasdjfh sdfjsdh flasjkfh asdjlkf has"
  if (!is.null(long.names) ) {
    m <- min(ncharLabel,length(vLabels))
    labels <- paste(vLabels[seq(1,m)],long.names[seq(1,m)],sep=" - ")

    par("mar"=c(2,2,0,2))
    #par("mar"=c(2,2,2,2))
    frame()
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="grey90")
    legend("topleft",labels,bty='n')
  }
  # interactive
  #tkplot
  par(mar=opar$mar)
  par(bg=opar$bg)
  return(g)
  }


# reimplement Piano's code
runGSAhyper2 <- function (genes, pvalues, pcutoff, universe,
                         gsc, gsSizeLim = c(1, Inf),
                         adjMethod = "fdr") 
{
    
  if (length(gsSizeLim) != 2) 
    stop("argument gsSizeLim should be a vector of length 2")
  if (missing(genes)) {
    stop("argument genes is required")
  }
  else {
    genes <- as.vector(as.matrix(genes))
    if (class(genes) != "character") 
      stop("argument genes should be a character vector")
    if (length(unique(genes)) != length(genes)) 
      stop("argument genes should contain no duplicated entries")
  }
  if (missing(pvalues)) {
    pvalues <- rep(0, length(genes))
  }
  else {
    pvalues <- as.vector(as.matrix(pvalues))
    if (class(pvalues) != "numeric") 
      stop("argument pvalues should be a numeric vector")
    if (length(pvalues) != length(genes)) 
      stop("argument pvalues should be the same length as argument genes")
    if (max(pvalues) > 1 | min(pvalues) < 0) 
            stop("pvalues need to lie between 0 and 1")
  }
  if (missing(pcutoff)) {
    if (all(pvalues %in% c(0, 1))) {
      pcutoff <- 0
    }
    else {
      pcutoff <- 0.05
    }
  }
  else {
    if (length(pcutoff) != 1 & class(pcutoff) != "numeric") 
      stop("argument pcutoff should be a numeric of length 1")
    if (max(pcutoff) > 1 | min(pcutoff) < 0) 
      stop("argument pcutoff needs to lie between 0 and 1")
  }
  if (missing(gsc)) {
    stop("argument gsc needs to be given")
  }
  else {
    if (class(gsc) != "GSC") 
      stop("argument gsc should be of class GSC, as returned by the loadGSC function")
  }
  if (missing(universe)) {
    if (!all(pvalues == 0)) {
      universe <- genes
      message("Using all genes in argument genes as universe.")
    }
    else {
            universe <- unique(unlist(gsc$gsc))
            message("Using all genes present in argument gsc as universe.")
          }
  }
  else {
    if (class(universe) != "character") 
      stop("argument universe should be a character vector")
    if (!all(pvalues == 0)) 
      stop("if universe is given, genes should be only the genes of interest, i.e. pvalues should all be set to 0.")
  }
  if (!all(unique(unlist(gsc$gsc)) %in% universe)) 
    warning("there are genes in gsc that are not in the universe, these will be removed before analysis")
  if (!all(genes %in% universe)) {
    warning("not all genes given by argument genes are present in universe, these will be added to universe")
    universe <- c(universe, genes[!genes %in% universe])
  }
  if (length(unique(universe)) != length(universe)) 
    stop("argument universe should contain no duplicated entries")
  tmp <- try(adjMethod <- match.arg(adjMethod, c("holm", "hochberg", 
                                                 "hommel", "bonferroni", "BH", "BY", "fdr", "none"), several.ok = FALSE), 
             silent = TRUE)
  if (class(tmp) == "try-error") {
    stop("argument adjMethod set to unknown method")
  }
  pvalues[pvalues == 0] <- -1e-10
  goi <- genes[pvalues < pcutoff]
  if (length(goi) < 1) 
    stop("no genes selected due to too strict pcutoff")
  bg <- universe[!universe %in% goi]
  gsc <- gsc$gsc
  #
  library(parallel)
  gs_sizes <- mclapply(gsc,length)  
  gsc <- gsc[gs_sizes >=gsSizeLim[1] &  gs_sizes <=gsSizeLim[2]]
  message(paste("Analyzing the overrepresentation of ", length(goi), 
                " genes of interest in ", length(gsc), " gene sets, using a background of ", 
                length(bg), " non-interesting genes.", sep = ""))

  do.test <- function(index,gsc,universe,goi,bg) {
    gs <- gsc[[index]]
    nogs <- universe[!universe %in% gs]
    ctab <- rbind(c(sum(goi %in% gs), sum(goi %in% nogs)), 
                  c(sum(bg %in% gs), sum(bg %in% nogs)))
    p <- fisher.test(ctab, alternative = "greater")$p.value
    rownames(ctab) <- c("Significant", "Non-significant")
    colnames(ctab) <- c("Genes in gene set", "Genes not in gene set")
    return(c(p, NA, sum(goi %in% gs), sum(bg %in% 
                                             gs),
             sum(goi %in% nogs), sum(bg %in% nogs)))           
  }
  a<-  mclapply(names(gsc),do.test,gsc,universe,goi,bg)
  names(a) <- names(gsc)
  resTab <- matrix(unlist(a),byrow=T,ncol=6)
  rownames(resTab) <- names(gsc)
  colnames(resTab) <- c("p-value", "Adjusted p-value", "Significant (in gene set)", 
                        "Non-significant (in gene set)", "Significant (not in gene set)", 
                        "Non-significant (not in gene set)")
  
  #hist(resTab[,"p-value"])

  res <- list()
  res$pvalues <- resTab[,"p-value"]
  res$p.adj <- p.adjust(res$pvalues, method = adjMethod)
  resTab[, 2] <- res$p.adj
  res$resTab <- resTab
  #res$contingencyTable <- contTabList
  res$gsc <- gsc
  return(res)
}

###################################################
#

importArgsfromStdin <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if ( length(args)== 0 ) {
    return(args)
  }
  if ( args[1]=="-stdin" ) {
    f <- file("stdin")
    open(f)
    line <- readLines(f,n=1)
    ##write(line, stderr())
    # process line
    args <- commandline.args2vector(line)
    close(f)
  }
  return(args)
}

commandline.args2vector <- function(line) {
  args <- mytrim(strsplit(mytrim(line),split="[ ]+")[[1]])
  args <- sub("^\"","",sub("\"$","",args))
  # where do the args names are (-)
  i <- rev(grep("^-",args))
  new.args <- c()
  #print(i)
  if ( length(i)==0 ) {
    return(args)
  }
  for ( idx in seq(1,length(i))) {
    if ( idx==1 && i[1]==length(args) ) {
      # last option does not have a value
      new.args <- args[i[1]]
    } else {
      if (idx==1) {
        last <- length(args)
      } else {
        last <- i[idx-1]-1
      }
      cur.l <- i[idx]+1
      #print(cur.l)
      #print(last)
      sel.i <- seq(cur.l,last)
      if (cur.l>last ) {
        # no value
        new.args <- append(args[cur.l-1],new.args)
      } else {
        sel <- paste(args[sel.i],sep=" ",collapse=" ")
        new.args <- append(append(args[cur.l-1],sel),new.args)
      }
      #print(new.args)
    }
    #print(new.args)        
    #print("------------------------------")
  }
  #new.args
  return(new.args)
}

commandline.args2vector.tests <- function() {
  sum(length(commandline.args2vector(""))==0)==1
  sum((commandline.args2vector("-1")==c("-1")))==1
  sum((commandline.args2vector("-1 a")==c("-1","a")))==2
  sum((commandline.args2vector("-1 a b")==c("-1","a b")))==2
  sum((commandline.args2vector("-1 a b -2")==c("-1","a b","-2")))==3
  sum((commandline.args2vector("-1 a b -2 -3")==c("-1","a b","-2","-3")))==4
  sum(commandline.args2vector("-0 -1 a b -2")==c("-0","-1","a b","-2"))==4
  sum(commandline.args2vector("a b")==c("a","b"))==2
}
## quantile_norm(quantile_norm(df),quantile_norm(df)$means)
quantile_norm_vect <- function(v,qn_values) {
  lv <- length(v)
  lqn <- length(qn_values)
  if ( lv!=lqn ) {
    perror("length of vector v (",lv,") is different from qn_means' length (",lqn,")")
  }
  
  p <- rank(v,ties.method="max")
  return(qn_values[p])
}

quantile_norm <- function(df,means=NULL){
  if ( ! is.data.frame(df) ) {
    perror("Expected a data frame")
  }  
  #
  if ( ! is.null(means)) {
    l1 <- nrow(df)
    l2 <- length(means)
    if ( l1 != l2 ) {
      pwarning("Number of rows in data frame (",l1,") does not match with the length of quantile normalized means vector (",l2,")")
      if ( l1 > l2 ) {
        perror("Unable to proceed")
      }
      # l1 <l2
      offset <- l2-l1
      means <- means[append(2,seq(offset+2,l2))]
      #pinfo(length(means),"==",l1)
    }  
  }

  # increasing
  ranks <- apply(df,2,rank,ties.method="max")
  # sort: increasing
  if (is.null(means) ) {
    means <- apply(data.frame(apply(df, 2, sort)), 1, mean, na.rm=T)
  }
  df_qn<- apply(ranks, 2, quantile_norm_vect, means)
  rownames(df_qn) <- rownames(df)
  return(list(qn=df_qn,means=means))
}
