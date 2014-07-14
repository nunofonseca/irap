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
#    $Id$
# =========================================================
# 
#source(paste(IRAP.DIR,"aux/R","irap_compare.R",sep="/"))
source(paste(IRAP.DIR,"aux/R","irap_utils.R",sep="/"))
source(paste(IRAP.DIR,"aux/R","irap_misc.R",sep="/"))

#install.packages("sfsmisc")
#library(sfsmisc)

# pprint irap names
pprint_irap_names <- function(label) {

  label <- sub("cufflinks1_nd","cufflinks1",label)
  label <- sub("cufflinks2_nd","cufflinks2",label)
  label <- sub("htseq1","htseq-u",label)
  label <- sub("htseq2","htseq-ine",label)
  return(sub("#"," x ",sub("_-_"," x ",label)))
}

ens2geneNames <- function(v,annot.table) {
  ens2gene <- function(ensid,annot.table) {
    name <- as.character(annot.table$Name[annot.table$ID==ensid])
    if ( length(name)==0 || is.null(name) || name=="") {
      name <- ensid
    }
    name
  }
  sapply(v,ens2gene,annot.table)
}

##################################################################
plot.selection.auto <- function(file.prefix,cor.data,counts.data,de.base.mean,de.fold.change,de.pval,fdr,annot.data=NULL) {
  NGENES <- 6
  # pipe1<-p$row
  # pipe2=p$col
  # expr.data<-counts.data
  plot.sel <- function(file=NULL,pipe1,pipe2,cor.data,expr.data,de.base.mean,de.fold.change,de.pval,fdr,bp.log="") {    
    if ( !is.null(file) ) {
      # High resolution
      png(filename=paste(file,".png",collapse="",sep=""),width=16,height=4,units="in",res=300)
    }
    pipe1.label <- pprint_irap_names(sub("X"," x ",pipe1))
    pipe2.label <- pprint_irap_names(sub("X"," x ",pipe2))
    #layout(matrix(c(1,2,3,4),ncol=2,byrow=TRUE))
    #dev.off()
    par(mfrow=c(2,2),xpd=F)
    par(mar=c(4,5,4,3))
    par(mgp=c(2,1,0))

    # boxplot
    par(bty="l")
    boxplot(list("Multiple pipelines"=as.vector(cor.data)) ,outline=F,ylab="Spearman correlation",xlab="Multiple pipelines",xatx="n")
    plot.panel.label("A",add="",x=2)
    # xlab=paste("P(",pipe1,",",pipe2,")=",round(cor.data[pipe1,pipe2],2),sep=""))
    # show the correlation aggregated by mapper?
    if ( is.vector(cor.data) ) {
      points(cor.data[paste(pipe1," X ",pipe2,sep="")],col="purple",pch=0,cex=1.2)
    } else {
      points(cor.data[pipe1,pipe2],col="purple",pch=0,cex=1.2)
    }
    #text(1,cor.data[pipe1,pipe2],labels=paste(pipe1.label," VS ",pipe2.label,sep=""),cex=0.8,pos=1)
    # expression scatter
    df <- data.frame(listOflists2matrix(gene.expr[c(pipe1,pipe2)]))
    colnames(df) <- c(pipe1.label,pipe2.label)
    # FIX!
    par(bty="l")
    par(mar=c(4,4,4,3))
    cor.scatterplot(df,pipe1.label,pipe2.label,value.name="",panel=NULL,print.cor=FALSE)
    plot.panel.label("B",add="",x=2)    
    # heatmaps
    # gen.heatmap(as.matrix(counts[,c(comp1,comp2)]),ncolors=2,dendrogram="none",main=paste("Counts",sep=""),density.info="none",Rowv=FALSE,Colv=FALSE,cexRow=0.65,cexCol=0.65)
    # fold change (MA plot)
    #
    de.key <-  names(de.pval)[names(de.pval) %in% c(paste(pipe1,"X",pipe2,sep=""),paste(pipe2,"X",pipe1,sep=""))]
    de.data <- data.frame(genes=names(de.base.mean[[de.key]]),baseMean=as.vector(de.base.mean[[de.key]]),log2FoldChange=as.vector(de.fold.change[[de.key]]),padj=as.vector(de.pval[[de.key]]))
    rownames(de.data) <- de.data$genes
    # exclude NA
    de.data <- de.data[!is.na(de.data$log2FoldChange) & !is.infinite(de.data$log2FoldChange),]
    #de.data <- de.data[!is.na(de.data$log2FoldChange),]
    #fold.change.plot(de.data,fdr=fdr,y.lab=paste("log2(",pipe1.label,"/",pipe2.label,")",sep=""))
    fold.change.plot(de.data,fdr=fdr,y.lab=paste("log2(fold change)",sep=""),x.lab="log(mean expression)",ylim=c(min(de.data$log2FoldChange),max(de.data$log2FoldChange)))
    plot.panel.label("C",add="",x=1)    
    # highlight (replot) the selected genes
    p.vals <- de.pval[[de.key]]
    genes<- head(sort(p.vals[names(p.vals) %in% rownames(de.data)]) ,NGENES)
    genes.d <- de.data[names(genes),]
    genes.d
    # check if fold change is always in the same direction
    if (sum(genes.d$log2FoldChange>0)==nrow(genes.d)) {
      pinfo("Adding a negative fold change example...if possible")
      de.sorted <- rownames(de.data)[order(de.data$padj,decreasing=F)]
      fc.sorted <- de.data[de.sorted,"log2FoldChange"]
      names(fc.sorted) <- de.sorted
      sel.gene <- head(fc.sorted[fc.sorted<0],1)
      if (! is.null(sel.gene) ) {
        genes.d[nrow(genes.d),] <- de.data[names(sel.gene),]
        genes <- append(genes[1:length(genes)-1],p.vals[names(sel.gene)])
        pinfo("Added ",names(sel.gene))
        print(de.data[names(sel.gene),])
      }
    } else {
      if (sum(genes.d$log2FoldChange<0)==nrow(genes.d)) {
        pinfo("Adding a positive fold change example...if possible")
        de.sorted <- rownames(de.data)[order(de.data$padj,decreasing=F)]
        fc.sorted <- de.data[de.sorted,"log2FoldChange"]
        names(fc.sorted) <- de.sorted
        sel.gene <- head(fc.sorted[fc.sorted>0],1)
        if (! is.null(sel.gene) ) {
          genes.d[nrow(genes.d),] <- de.data[names(sel.gene),]
          genes <- append(genes[1:length(genes)-1],p.vals[names(sel.gene)])
          pinfo("Added ",names(sel.gene))
          print(de.data[names(sel.gene),])
        }
      }
    }
    # get the gene names from the annot.table (if possible)
    colors <- rainbow(nrow(genes.d))
    points(log(genes.d$baseMean),genes.d$log2FoldChange,pch=0,cex=1.2,col=colors)
    # plot a barplot showing the expression values of the selected genes with the respective p.value (de and ttest)
    bcolors <- rainbow(2)
    par(mar=c(7,5,4,4),xpd=TRUE,bty="n")
    log <- bp.log
    sel.genes <- t(as.matrix(df[names(genes),]))
    colnames(sel.genes) <- ens2geneNames(colnames(sel.genes),annot.table)
    if ( sum(sel.genes==0)>0) { log <- "" }
    #col.axis="darkviolet",
    par(mgp=c(4,1,0))
    
    bp <- barplot(height=sel.genes,log=log,ylab="Counts",las=2,beside=T,col=bcolors,cex.names=0.8,
            col.axis="darkviolet",
            yaxt="n")
    aY <- axTicks(2)
    if ( max(aY) <= 9999 ) {
      axis(2, at=aY, label= round(aY,0),las=1,cex.axis=1)
    } else {
      axis(2, at=aY, label= axTexpr(2, aY),las=1,cex.axis=0.9)
    }
    for (n in seq(1,ncol(sel.genes)) ) {
      axis(1,at=(colSums(bp)/2)[n],colnames(sel.genes)[n],col.axis=colors[n],las=2,cex.axis=0.8,lwd.ticks=0)
    }
    legend("topright", inset=c(0,-0.35), legend=c(pipe1.label,pipe2.label), pch=20, col=bcolors, title="", bty='n')
    plot.panel.label("D",add="",x=2)
    par(mgp=c(3,1,0))
    if ( !is.null(file) ) {
      dev.off()
    }
  }
  # bp$stats
  # pick 2 values close to the 1st quantile, median and 3rd quantile
  get.comparison <- function(data.m,val,to.avoid=NULL) {
    # remove the cols/rows to.avoid
    if ( is.matrix(data.m)) {
      if(!is.null(to.avoid)) {
        data.m <- data.m[!rownames(data.m)%in%to.avoid,]
        data.m <- data.m[,!colnames(data.m)%in%to.avoid]
      }
      if ( nrow(data.m)==0 ) {
        return(NULL)
      }
      data.v <- data.m[upper.tri(data.m)]
      data.v <- data.v[!is.na(data.v)]
      m <- as.vector(data.v)
      closest.val <- m[which(abs(m-val)==min(abs(m-val)))]
      print(closest.val)
      if (length(closest.val) > 1 ) {
        closest.val <- closest.val[1]
      }                                        #
      pos <- which(data.m==closest.val,arr.ind=T)
      print(pos)
      row <- rownames(data.m)[pos[1,"row"]]
      col <- colnames(data.m)[pos[1,"col"]]
      list(row=row,col=col)
    } else {
      if (is.vector(data.m)) {
        m <- data.m
        if(!is.null(to.avoid)) {
          m <- m[! names(m) %in% to.avoid]                    
        }
        closest.val <- m[which(abs(m-val)==min(abs(m-val)))]
        print(closest.val)
        if (length(closest.val) > 1 ) {
          closest.val <- closest.val[1]
        }                                        #
        return(names(closest.val))
      } else {
        return(NULL)
      }
    }
  }
  npipes2sel <- c(1,2,3,4,5,6,7,8,9,10,11,12,13)
  sel.pipelines <- NULL
  if ( is.vector(cor.data) ) {
    bp <- boxplot(cor.data,plot=F)
  } else {
    bp <- boxplot(cor.data[upper.tri(cor.data)],plot=F)
  }
  # get the two pipelines
  for (sel.pipes in npipes2sel ) {
    p<-get.comparison(cor.data,bp$stats[3],sel.pipelines)
    if (!is.null(p)) {
      if (length(p) == 1 ) {
        sel.pipelines <- unique(append(sel.pipelines,p))
        p <- strsplit(p,split=" X ")[[1]]
        pipe1 <- p[1]
        pipe2 <- p[2]
      } else {
        sel.pipelines <- unique(append(sel.pipelines,append(p$row,p$col)))
        pipe1 <- p$row
        pipe2 <- p$col
      }
      pinfo("plotting to file with prefix ",file.prefix)
      file <- NULL
      if ( !is.null(file.prefix) ) {
        file <- paste(file.prefix,"_median",sel.pipes,sep="")
      }
      gen.plot2report(filename=paste(paste(file,".png",sep="")),
                      html=F,
                      width=600,
                      height=600,
                      to.plot=function() { 
                        plot.sel(NULL,pipe1,pipe2,cor.data,counts.data,de.base.mean,de.fold.change,de.pval,fdr)
                      })
    }
  }
  # worst
  sel.pipelines <- NULL
  for (sel.pipes in npipes2sel ) {
    p<-get.comparison(cor.data,bp$stats[2],sel.pipelines)
    if (!is.null(p)) {
      if (length(p) == 1 ) {
        sel.pipelines <- unique(append(sel.pipelines,p))
        p <- strsplit(p,split=" X ")[[1]]
        pipe1 <- p[1]
        pipe2 <- p[2]
      } else {
        sel.pipelines <- unique(append(sel.pipelines,append(p$row,p$col)))
        pipe1 <- p$row
        pipe2 <- p$col
      }
      file <- NULL
      if ( !is.null(file.prefix) ) {
        file <- paste(file.prefix,"_low",sel.pipes,sep="")
      }
      gen.plot2report(filename=paste(paste(file,".png",sep="")),
                      html=F,
                      width=600,
                      height=600,
                      to.plot=function() { plot.sel(NULL,pipe1,pipe2,cor.data,counts.data,de.base.mean,de.fold.change,de.pval,fdr) }
                      )
    }
  }
  # higher correlation
  sel.pipelines <- NULL
  for (sel.pipes in npipes2sel ) {
    p<-get.comparison(cor.data,bp$stats[4],sel.pipelines)
    if (!is.null(p)) {
      if (length(p) == 1 ) {
        sel.pipelines <- unique(append(sel.pipelines,p))
        p <- strsplit(p,split=" X ")[[1]]
        pipe1 <- p[1]
        pipe2 <- p[2]
      } else {
        sel.pipelines <- unique(append(sel.pipelines,append(p$row,p$col)))
        pipe1 <- p$row
        pipe2 <- p$col
      }
      file <- NULL
      if ( !is.null(file.prefix) ) {
        file <- paste(file.prefix,"_high",sel.pipes,sep="")
      }
      gen.plot2report(filename=paste(paste(file,".png",sep="")),
                      html=F,
                      width=600,
                      height=600,
                      to.plot=function() { 
                        plot.sel(NULL,pipe1,pipe2,cor.data,counts.data,de.base.mean,de.fold.change,de.pval,fdr)
                      })
    }
  }
}
#########################################################
selected.examples <- function(comparisons=4,genes=5) {

  abs.max <- function(x) {
    if (is.null(x)) { return(NULL) }
    x <- x[!is.na(x)]
    sel <- x[abs(x)==max(abs(x),na.rm=T)]
    names(sel) <- NULL
    sel[1]
  }
  k <- comparisons
  G <- genes
  # Sort the comparisons with highest fold changes
  s <- sort(unlist(lapply(log2.fold.change.deseq,abs.max)),decreasing=T)
  # remove duplicated comparisons
  sn <- names(s)
  keep <- rep(T,length(s))
  for ( i in c(1:length(s))) {
    n.splitted <- strsplit(sn[i],"X")
    inv.n <- paste(n.splitted[[1]][2],"X",n.splitted[[1]][1],sep="")
    if( inv.n %in% sn[1:i] ) { keep[i] <- F }
  }  
  # Select the K comparisons with highest fold changes
  sn <- sn[keep]
  sn <- sn[c(1:k)]       
  # For each comparison pick G genes with the highest p-values and counts
  sel.genes <- list()
  for ( c in sn ) {
    print(c)
    sel.genes[[c]] <- head(sort(p.val.deseq[[c]]),G)
  }
  return(sel.genes)
}

# Now we have a comparison, and G genes we sim.flux.150.nde.30.pe.conf sim.flux.200.nde.30.pe.conf sim.flux.150.nde.60.se.conf
#   a) plot a scatter plot with correlation coefficient
#   b) heatmap of all genes
#   b) fold change with p-values and selected genes highlighted
#   c) plot a barplot showing the expression values of the selected genes with the respective p.value (de and ttest
plot.selection <- function(comp1,comp2,genes,counts,baseMean,fold.change,p.values,p.values2,fdr) {

  par(mfrow=c(1,3),xpd=F)
  # expression scatter
  cor.scatterplot(counts,comp1,comp2,paste(comp1,"/",comp2,sep=""))
  # heatmaps
  # gen.heatmap(as.matrix(counts[,c(comp1,comp2)]),ncolors=2,dendrogram="none",main=paste("Counts",sep=""),density.info="none",Rowv=FALSE,Colv=FALSE,cexRow=0.65,cexCol=0.65)
  # fold change (MA plot)
  #
  de.data <- data.frame(baseMean=baseMean,log2FoldChange=fold.change,padj=p.values)
  # exclude NA
  de.data <- de.data[!is.na(de.data$log2FoldChange),]
  fold.change.plot(de.data,fdr=fdr,y.lab=paste("log2(",comp1,"/",comp2,")",sep=""))
  # highlight (replot) the selected genes
  genes.d <- de.data[names(genes),]
  points(genes.d$baseMean,genes.d$log2FoldChange,pch=0,cex=0.8,col="black")
  #points(genes.d$baseMean,genes.d$log2FoldChange,pch=22,cex=0.8,col="black")
  # plot a barplot showing the expression values of the selected genes with the respective p.value (de and ttest)
  bcolors <- rainbow(2)
  barplot2(t(as.matrix(counts[names(genes),c(comp1,comp2)])),log="y",ylab="Read Counts",las=2,beside=T,col=bcolors,cex.names=0.8)
  legend("topright", inset=c(-0.2,0), legend=c(comp1,comp2), pch=22, col=bcolors, title="", bty='n')
}


irap.compare.raw.quant <- function(tsv1.data,tsv2.data,annot.table) {
# fix colnames of simulated quantification (remove the suffix .true)
  colnames(tsv1.data) <- gsub(".true","",colnames(tsv1.data))
  colnames(tsv2.data) <- gsub(".true","",colnames(tsv2.data))

  if ( nrow(tsv1.data)!=nrow(tsv2.data) ) {
    pinfo("Processing simulated data? #rows:",nrow(tsv1.data),"!=",nrow(tsv2.data))
  }
  # change the column names
  colnames(tsv2.data) <- paste(colnames(tsv2.data),".y",sep="")
  # merge (consider only the rows that appear on both files)
  i <- intersect(rownames(tsv1.data),rownames(tsv2.data))
  aggr.data <- cbind(tsv1.data[i,], tsv2.data[i,])
  #head(aggr.data)

  labels.v<-c("A","B")   
#
  tsv1.data <-tsv1.data[i,]
  tsv2.data <-tsv2.data[i,]

  ###################################
  # Compute the paired t-test for each gene
  # to check if the means are the sam
  ttest.wrapper <- function(row,m1,m2) {
    x <- tryCatch(t.test(as.numeric(m1[row,]),as.numeric(m2[row,]),paired=T)$p.value,error=function(x) 1.0)
    if ( is.na(x) ) { return(1.0) }
    return(x)
  }
################################################
# 
#N1 <- colSums(tsv1.data+1)
#N2 <- colSums(tsv2.data+1)
# transform the data by adding one
# normalize by library size and get the value in millions
# get the value in millions to use the matrix with deseq
# Note: quant.nerr may contain NAs (add 1)
#       prop.table(matrix,2) ==< by column
  tsv1.data.n <- round(prop.table(as.matrix(tsv1.data+1),2)*10^6,0)
  tsv2.data.n <- round(prop.table(as.matrix(tsv2.data+1),2)*10^6,0)

###############################
  pinfo("Checking for sign. different values using DESeq...")
  suppressPackageStartupMessages(library(DESeq))
  conds <- append(rep(labels.v[1],ncol(tsv1.data)),rep(labels.v[2],ncol(tsv2.data)))
# TODO: aggregate the cols from the same sample?
#rows.sel <- apply(aggr.data[,-1],1,max)>min_count ;# filter out the rows with the maximum number of reads under the given threshold
#cds <- newCountDataSet(aggr.data[rows.sel,-1],conds)
  pinfo(conds)
  cds <- newCountDataSet(aggr.data,conds)
  cds <- estimateSizeFactors(cds)

  if ( sum(is.na(sizeFactors(cds))) > 0 ) {
    save.image("err.debug.Rdata")
    perror("Unable to estimate size factors")
    q(status=1)
  }

# from the manual first computes for each gene an empirical dispersion
#value (a.k.a. a raw SCV value), then fits by regression a
#dispersion-mean relationship and finally chooses for each gene a
#dispersion parameter that will be used in subsequent tests from the
#empirical and the fitted value according to the 'sharingMode'
#argument.
# pooled Use the samples from all conditions with replicates
#to estimate a single pooled empirical dispersion value, called
#"pooled", and assign it to all samples.
# blind - ‘blind’ - Ignore the
#sample labels and compute a gene's empirical dispersion value as if
#all samples were replicates of a single condition. This can be done
#even if there are no biological replicates. This method can lead to
#loss of power; see the vignette for details. The single estimated
#dispersion condition is called "blind" and used for all samples.
#Use the samples from all conditions with replicates
#to estimate a single pooled empirical dispersion value
# sharing: After the empirical dispersion values have been computed for
# each gene, a dispersion-mean relationship is fitted for sharing
# information across genes in order to reduce variability of the
# dispersion estimates. After that, for each gene, we have two values:
# the empirical value (derived only from this gene's data), and the
# fitted value (i.e., the dispersion value typical for genes with an
# average expression similar to those of this gene).

# SharingMode="maximum"
  result <- try(cds <- estimateDispersions(cds,method="pooled",fitType="local",sharingMode="maximum"));
  if(class(result) == "try-error") {
   # pooled-CR: crossed factors
    print("Ooops, it seems that you need to manually tune DEseq  to estimate the dispersion.")
    q("no",status=1)
  }
  t<-counts(cds,normalized=TRUE)
  de <- nbinomTest(cds,labels.v[1],labels.v[2])
  rownames(de) <- de$id
  pinfo("Checking for sign. different values using DESeq...done.")
##################################################################
#
  pdebug.save.state("compare_raw_quant","stage3")
  if (!init.source.filter(annot.table)) {
    perror("Internal error while initializing gene filter.")
    q(status=1)
  }
#head(annot.table)
# for debugging
####################################
  filter.name <- "protein coding"
  pinfo("Filter: ",filter.name)

  tsv1.data.f <- apply.source.filter(tsv1.data,filter.name)
  tsv2.data.f <- apply.source.filter(tsv2.data,filter.name)
  
  tsv1.data.n.f <- apply.source.filter(tsv1.data.n,filter.name)
  tsv2.data.n.f <- apply.source.filter(tsv2.data.n,filter.name)
  
##############################################
  comp.data <- list()
  comp.data[["labels"]] <- labels.v
  comp.data[["filter"]] <- filter.name
  comp.data[["pearson"]] <- cor(apply(tsv1.data.f,1,sum),apply(tsv2.data.f,1,sum),method="pearson")
  comp.data[["pearson.n"]] <- cor(apply(tsv1.data.n.f,1,sum),apply(tsv2.data.n.f,1,sum),method="pearson")
  comp.data[["spearman"]] <- cor(apply(tsv1.data.f,1,sum),apply(tsv2.data.f,1,sum),method="spearman")
  comp.data[["spearman.n"]] <- cor(apply(tsv1.data.n.f,1,sum),apply(tsv2.data.n.f,1,sum),method="spearman")
  comp.data[["de"]] <- apply.source.filter(de,filter.name)
  comp.data[["normalized.counts"]] <- apply.source.filter(t,filter.name)
  return(comp.data)
}

compute.cor <- function(expr.list,cor.method="spearman") {

  comps <- c()
  pipelines <- names(expr.list)
  for ( n1 in seq(1,length(pipelines)-1) ) {
    for ( n2 in seq(n1+1,length(pipelines)) ) {
      i <- intersect(names(expr.list[[n1]]),names(expr.list[[n2]]))
      
      cor <- cor.test(x=expr.list[[n1]][i],y=expr.list[[n2]][i],method=cor.method,exact=T)
      s <- cor$estimate
      names(s) <- paste(pipelines[n1]," X ",pipelines[n2],sep="")
      comps <- append(comps,s)
    }
  }
  return(comps)
}
#debug
## setwd("/home/nf/Research/WorkingDocs/RNA_comp/ibm/")
## load("/home/nf/Research/WorkingDocs/RNA_comp/ibm/ibm.html_data.Rdata")
## oprefix="test"
## plot.selection.auto(oprefix,p.cor,gene.nexpr,base.mean.deseq,log2.fold.change.deseq,p.val.deseq,fdr)
## p.cor
## p.cor
## spearman.corr
## pearson.corr
## gene.nexpr
## base.mean.deseeq
## names(log2.fold.change.deseq)
## # comparisons
## names(p.val.deseq)
## fdr
## ls()
