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
source(paste(IRAP.DIR,"aux/R","deseq_shared.R",sep="/"))

#install.packages("sfsmisc")
library(sfsmisc)

# pprint irap names
pprint_irap_names <- function(label) {

  label <- sub("cufflinks1_nd","cufflinks1",label)
  label <- sub("cufflinks2_nd","cufflinks2",label)
  label <- sub("htseq1","htseq-u",label)
  label <- sub("htseq2","htseq-ine",label)
  label
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
    par(mfrow=c(2,2),xpd=F)
    par(mar=c(4,4,4,4))
    # boxplot
    par(bty="l")
    boxplot(list("Multiple pipelines"=as.vector(cor.data)) ,outline=F,ylab="Pearson correlation",xlab="Multiple pipelines",xatx="n")
    plot.panel.label("A")
    # xlab=paste("P(",pipe1,",",pipe2,")=",round(cor.data[pipe1,pipe2],2),sep=""))
    # show the correlation aggregated by mapper?
    points(cor.data[pipe1,pipe2],col="purple",pch=20,cex=1.2)    
    # expression scatter
    df <- data.frame(listOflists2matrix(gene.expr[c(pipe1,pipe2)]))
    colnames(df) <- c(pipe1.label,pipe2.label)
    # FIX!
    par(bty="l")
    par(mar=c(4,4,4,3))
    cor.scatterplot(df,pipe1.label,pipe2.label,value.name="",panel="B")
    # heatmaps
    # gen.heatmap(as.matrix(counts[,c(comp1,comp2)]),ncolors=2,dendrogram="none",main=paste("Counts",sep=""),density.info="none",Rowv=FALSE,Colv=FALSE,cexRow=0.65,cexCol=0.65)
    # fold change (MA plot)
    #
    de.key <-  names(de.pval)[names(de.pval) %in% c(paste(pipe1,"X",pipe2,sep=""),paste(pipe2,"X",pipe1,sep=""))]
    de.data <- data.frame(genes=names(de.base.mean[[de.key]]),baseMean=as.vector(de.base.mean[[de.key]]),log2FoldChange=as.vector(de.fold.change[[de.key]]),padj=as.vector(de.pval[[de.key]]))
    rownames(de.data) <- de.data$genes
    # exclude NA
    de.data <- de.data[!is.na(de.data$log2FoldChange) & !is.infinite(de.data$log2FoldChange),]
    #fold.change.plot(de.data,fdr=fdr,y.lab=paste("log2(",pipe1.label,"/",pipe2.label,")",sep=""))
    fold.change.plot(de.data,fdr=fdr,y.lab=paste("log2(fold change)",sep=""),x.lab="log(mean expression)")
    plot.panel.label("C")    
    # highlight (replot) the selected genes
    p.vals <- de.pval[[de.key]]
    genes<- head(sort(p.vals[names(p.vals) %in% rownames(de.data)]) ,NGENES)
    genes.d <- de.data[names(genes),]
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
      if (sum(genes.d$log2FoldChange>0)==nrow(genes.d)) {
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
    points(genes.d$baseMean,genes.d$log2FoldChange,pch=20,cex=0.8,col="purple")
    # plot a barplot showing the expression values of the selected genes with the respective p.value (de and ttest)
    bcolors <- rainbow(2)
    par(mar=c(7,4,4,4),xpd=TRUE,bty="n")
    log <- bp.log
    sel.genes <- t(as.matrix(df[names(genes),]))
    colnames(sel.genes) <- ens2geneNames(colnames(sel.genes),annot.table)
    if ( sum(sel.genes==0)>0) { log <- "" }
    barplot(height=sel.genes,log=log,ylab="Counts",las=2,beside=T,col=bcolors,cex.names=0.8,
            col.axis="darkviolet",
            yaxt="n")
    aY <- axTicks(2)
    if ( max(aY) <= 9999 ) {
      axis(2, at=aY, label= round(aY,0),las=1,cex.axis=1)
    } else {
      axis(2, at=aY, label= axTexpr(2, aY),las=1,cex.axis=0.9)
    }
    legend("topright", inset=c(0,-0.4), legend=c(pipe1.label,pipe2.label), pch=20, col=bcolors, title="", bty='n')
    plot.panel.label("D")    
    if ( !is.null(file) ) {
      dev.off()
    }
  }
  # bp$stats
  # pick 2 values close to the 1st quantile, median and 3rd quantile
  get.comparison <- function(data.m,val,to.avoid=NULL) {
    # remove the cols/rows to.avoid
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
    }
    #
    pos <- which(data.m==closest.val,arr.ind=T)
    print(pos)
    row <- rownames(data.m)[pos[1,"row"]]
    col <- colnames(data.m)[pos[1,"col"]]
    list(row=row,col=col)
  }
  npipes2sel <- c(1,2,3,4,5)
  sel.pipelines <- NULL
  bp <- boxplot(cor.data[upper.tri(cor.data)],plot=F)
  # get the two pipelines
  for (sel.pipes in npipes2sel ) {
    p<-get.comparison(cor.data,bp$stats[3],sel.pipelines)
    if (!is.null(p)) {      
      sel.pipelines <- unique(append(sel.pipelines,append(p$row,p$col)))
      pinfo("plotting to file with prefix",file.prefix)
      file <- NULL
      if ( !is.null(file.prefix) ) {
        file <- paste(file.prefix,"_median",sel.pipes,sep="")
      }
      gen.plot2report(filename=paste(paste(file,".png",sep="")),
                      html=F,
                      width=600,
                      height=600,
                      to.plot=function() { 
                        plot.sel(NULL,p$row,p$col,cor.data,counts.data,de.base.mean,de.fold.change,de.pval,fdr)
                      })
    }
  }
  # worst
  sel.pipelines <- NULL
  for (sel.pipes in npipes2sel ) {
    p<-get.comparison(cor.data,bp$stats[2],sel.pipelines)
    if (!is.null(p)) {
      sel.pipelines <- unique(append(sel.pipelines,append(p$row,p$col)))
      file <- NULL
      if ( !is.null(file.prefix) ) {
        file <- paste(file.prefix,"_low",sel.pipes,sep="")
      }
      gen.plot2report(filename=paste(paste(file,".png",sep="")),
                      html=F,
                      width=600,
                      height=600,
                      to.plot=function() { plot.sel(NULL,p$row,p$col,cor.data,counts.data,de.base.mean,de.fold.change,de.pval,fdr) }
                      )
    }
  }
  # higher correlation
  sel.pipelines <- NULL
  for (sel.pipes in npipes2sel ) {
    p<-get.comparison(cor.data,bp$stats[4],sel.pipelines)
    if (!is.null(p)) {
      sel.pipelines <- unique(append(sel.pipelines,append(p$row,p$col)))
      file <- NULL
      if ( !is.null(file.prefix) ) {
        file <- paste(file.prefix,"_high",sel.pipes,sep="")
      }
      gen.plot2report(filename=paste(paste(file,".png",sep="")),
                      html=F,
                      width=600,
                      height=600,
                      to.plot=function() { 
                        plot.sel(NULL,p$row,p$col,cor.data,counts.data,de.base.mean,de.fold.change,de.pval,fdr)
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
  points(genes.d$baseMean,genes.d$log2FoldChange,pch=20,cex=0.8,col="black")
  # plot a barplot showing the expression values of the selected genes with the respective p.value (de and ttest)
  bcolors <- rainbow(2)
  barplot2(t(as.matrix(counts[names(genes),c(comp1,comp2)])),log="y",ylab="Read Counts",las=2,beside=T,col=bcolors,cex.names=0.8)
  legend("topright", inset=c(-0.2,0), legend=c(comp1,comp2), pch=20, col=bcolors, title="", bty='n')
}

