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
#    $Id: irap.txt Nuno Fonseca Tue Jan 29 23:13:16 2013$
# =========================================================

##
#
process.cmdline.args <- function(cmd) {
  args <- commandArgs(trailingOnly=TRUE)
    
  usage <- paste(cmd," --tsv file --min min.reads --contrasts contrast.def --labels label1;label2;... --annotation tsv.file [--annot-genes-only --feat gene|transcript] --out outprefix ",sep="")
  option_list <- list(
    make_option(c("--independent-filtering"),action="store_true",dest="indfilter",default=FALSE,help="Use independent filtering (DESeq2 only) [default %default]"),
    make_option(c("-m", "--min"), type="character", dest="min_count", default=NULL,help="exclude genes with counts < min"),
    make_option(c("--contrasts"), type="character", dest="contrasts", default=NULL,help="contrasts definition."),
    make_option(c("--tech-rep"), type="character", dest="tech_replicates", default=NULL,help="technical replicates. E.g., 'L1R1,L1R2;L2R1,L2R2,L2R3'"),
    make_option(c("--labels"), type="character", dest="labels", default=NULL,help="labels/contrasts names/aka factor values"),
    make_option(c("-a", "--annotation"), type="character",default=NULL,help="Annotation file (TSV format)"), 
    make_option(c("--annot-genes-only"), action="store_true",default=FALSE,dest="only.annot.genes",help="Only use genes in the DE analysis that appear in the annot. file."),
    make_option(c("--feature"), type="character",default="gene",dest="feature",help="Type of feature: gene, transcript (default %default%)."),
    make_option(c("--out"), type="character", dest="out", default=NULL,help="Output file prefix"),
    make_option(c("-i", "--tsv"), type="character", dest="tsv_file", default=NULL,help="TSV file name"),
    make_option(c("--debug"),action="store_true",dest="debug",default=FALSE,help="Debug mode")
  )

  # check multiple options values
  filenames <- c("tsv_file")
  multiple.options = list()
  mandatory <- c("tsv_file","out","labels","contrasts")
  opt <- myParseArgs(usage = usage, option_list=option_list,filenames.exist=filenames,multiple.options=multiple.options,mandatory=mandatory)

  if (opt$debug) {
    pdebug.enable()
  }
  #
  # back compat.
  # ;, separate list of groups : ex. file1,file2;file3,file4
  opt$labels.v <- strsplit(opt$labels,",")[[1]];
  # list with the files/cols that belong to a specific group
  opt$groups.def <- strsplit(opt$contrasts,";")[[1]]
  if (length(opt$labels.v)!=length(opt$groups.def)) {
    perror("labels=",opt$labels.v)
    perror("contrast=",opt$groups.def)
    perror("Contrast and labels mismatch!\n");
    q(status=1);
  }
  # label2group
  opt$label2group <- list();
  opt$contrast.l <- sapply(opt$groups.def,strsplit,",")
  names(opt$contrast.l) <- opt$labels.v
  for ( l in names(opt$contrast.l)) {
    for ( x in opt$contrast.l[[l]] ) {
      opt$label2group[[x]] <- l
    }
  }
  opt$tech.replicates.l <- NULL
  # tech_replicates  - technical replicates
  if ( !is.null(opt$tech_replicates) ) {
    tech.replicates.def <- strsplit(opt$tech_replicates,";")[[1]]
    opt$tech.replicates.l <- list()
    pinfo("Data with ",length(tech.replicates.def)," technical replicates...")
    for ( s in tech.replicates.def) {
      opt$tech.replicates.l[[paste(",",s,",",sep="")]] <- strsplit(s,",")[[1]]
    }
  }
  #
  opt$annot <- NULL
  if ( is.null(opt$annotation)) {
    pinfo("Annot.file not provided\n")
    if (opt$only.annot.genes) {
      perror("One annotation file must be provided when enabling --annot-genes-only")
      q(status=2)
    }
  } else {
    pinfo("Annot. file=",opt$annotation)
    is.empty <- FALSE
    annot <- NULL
    annot <- load.annot(as.character(opt$annotation))
    if ( is.null(annot) ) {
      is.empty <- TRUE
    } else if ( nrow(annot) ==0 ) {
      is.empty <- TRUE
    }
    if (is.empty==TRUE) {
      pinfo("Empty annotation file.")
      opt$annot <- NULL
      opt$annotation <- NULL
    } else {
      opt$annot <- annot
    }
    
  }
  
  ######################
  pinfo(" Matrix/counts=",opt$tsv_file)
  pinfo(" min_count=",opt$min_count)
  pinfo(" Labels=",opt$labels.v)
  pinfo(" Contrasts=",opt$contrast.l)
  pinfo(" Output prefix=",opt$out)
  ############################
  opt
}

####################################################
## Given the map of a file to a group and a vector with colnames
## produce the vector with the respective group names
## change to one liner ...
map.conds2cols <- function(label2group,cols) {
  conds <- cols
  i <- 1
  for (c in cols) {    
    conds[i] <- label2group[[conds[i]]]
    i <- i+1
  }
  names(conds) <- cols
  conds
}
# Handle technical replicates by merging
# the multiple libraries into a single one
handle.tech.replicates <- function(data,opt) {

  if ( is.null(opt$tech.replicates.l) ) {
    return
  }
  all <- unlist(opt$tech.replicates.l)
  processed <- c()
  nc<-ncol(data)
  for (n in colnames(data)) {
    if ( n %in% all & ! n %in% processed ) {
      x<-all[all==n]
      cols2aggr <- strsplit(gsub(",[0-9]+$","",gsub("^,","",names(x))),split=",")[[1]]
      #
      #pinfo("cols2aggr",names(x),"--->",cols2aggr)
      irap.assert(length(cols2aggr)>0,"empty technical replicates ",n,"...",x)      
      cols <- names(data)
      if ( sum(cols2aggr %in% cols)!=length(cols2aggr) ) {
        perror("technical replicates: some columns were not found in the file ")
        missing.cols <- cols2aggr[!cols2aggr %in% cols]
        perror("Columns not found: ",paste(missing.cols,sep=","))
        q(status=1)   
      }
      # sum
      aggr <- rowSums(data[,cols2aggr,drop=FALSE])
      cols <- cols[!cols %in% cols2aggr]
      data <- data[,cols,drop=FALSE]
      data[,n] <- aggr           
      #print(head(data))
      pinfo("Technical replicates: ",paste(cols2aggr,sep=",",collapse=","))
      processed <- append(processed,cols2aggr)
    } else {
      processed <- append(processed,n)
    }
  }
  nc2<-ncol(data)
  if ( nc != nc2) {
    pinfo("Technical replicates: matrix was reduced from ",nc," to ",nc2," columns") 
  }
  data
}
# filter the gene count matrix
# optionally filter by minimum number of reads
# optionally filter by list of genes provided
filter.read.counts.table <- function(data,opt) {
  data.f <- data

  # filter out column names not in the contrasts
  cols.used <- unique(unlist(opt$contrast.l))
  cols.names <- colnames(data)
  if ( sum(! cols.used %in% cols.names)!=0 ) {
    perror("Some columns were not found in the file ",opt$tsv_file)
    #
    missing.cols <- cols.used[!cols.used %in% cols.names]
    perror("Columns not found: ",paste(missing.cols,sep=","))
    q(status=1)
  }
  cols.names <- cols.names[!is.na(match(cols.names,cols.used))]
  data.f <- data.f[,cols.names,drop=FALSE]
  if (opt$min_count>0) {
    pinfo("Filtering out genes with low counts (<=",opt$min_count,")...")
    rows.sel <- apply(data,1,max)>opt$min_count ;# filter out the rows with the maximum number of reads under the given threshold
    data.f <- data.f[rows.sel,,drop=FALSE]
    pinfo("Filtering out genes with low counts (<=",opt$min_count,")...done.")
  }
  if(opt$only.annot.genes) {
    pinfo("Filtering out genes not in annotation file...")
    data.f <- data.f[rownames(data.f) %in% opt$annot$ID,,drop=FALSE]
    pinfo("Filtering out genes not in annotation file...done.")
  }
  data.f
}

save.de2tsv <- function(de.table,ofile.prefix) {
  ofile <- paste(opt$out,"de.tsv",sep="/")
  pinfo("Writing to TSV ",ofile)
  write.table(de.table,ofile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
