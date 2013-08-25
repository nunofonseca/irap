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
#    $Id: irap.txt Nuno Fonseca Sun Jan 27 01:52:00 2013$
# =========================================================


# plotMA
fold.change.plot <- function(res,y.lab="log2(Fold change)",x.lab="log(Mean)",fdr=0.1,log="",ylim=NA) {

  if (!(is.data.frame(res) && all(c("baseMean", "log2FoldChange") %in% colnames(res))))
    stop("'res' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")

  linecol = "#ff000080"

  x = subset(res, baseMean!=0)
  py = x$log2FoldChange

  if(missing(ylim))
    ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1

  if (missing(x.lab))
    x.lab <- paste(x.lab," [FDR=",fdr,";DE genes=",sum(res$padj < fdr)," (",sum(res$padj[res$log2FoldChange>0]<fdr),"/",sum(res$padj[res$log2FoldChange<0]<fdr),")]",sep="")
  vals<-round(res$padj*100)
  #colors<-append(colorRampPalette(c("darkred","orange"))(10),rep("darkgrey",90))
  colors<-append(colorRampPalette(c("red"))(10),rep("darkgrey",90))
  
  plot(
       log(res$baseMean),
       res$log2FoldChange,
       log=log,pch=20,cex=0.3,
       xlab=x.lab,ylab=y.lab,
       col = ifelse(res$padj < fdr,"red","darkgrey") )
  abline(h=0, lwd=3, col=linecol)
}
#       col = ifelse(res$padj < fdr,"red","black") )
#  plot(x$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
#       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
#       cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
cor.scatterplot <- function(df=NULL,colA=NULL,colB=NULL,x=NULL,y=NULL,
                            value.name="",smooth="",log="xy",panel=NULL) {
  library(sfsmisc)
  del <- 0.5
  if (!is.null(df)) {
    x <- df[,colA]
    y <- df[,colB]
  }
  s<-round(cor(x,y,method="spearman"),2)
  p<-round(cor(x,y,method="pearson"),2)
  if ( !is.null(log) ) {
    x<-x+del
    y<-y+del
  } else {
    log=""
  }
  if ( value.name =="") {
    main <- ""
  } else {
    main <- paste(value.name," (pearson r=",p,",spearman Rs=",s,")",sep="")
  }
  pch <- "."
  type <- "p"
  col <- "black"
  cex <- 1.1
  if (log!="") {
    if ( log=="xy" ) {
      plot(x,y,log=log,type=type,pch=pch,xlab=colA,ylab=colB,main=main,cex=cex,col=col,
           yaxt="n",xaxt="n")
    } else {
      if ( log=="x" ) {
        plot(x,y,log=log,type=type,pch=pch,xlab=colA,ylab=colB,main=main,cex=cex,col=col,
             xaxt="n")
      } else {
      plot(x,y,log=log,type=type,pch=pch,xlab=colA,ylab=colB,main=main,cex=cex,col=col,
           yaxt="n")
      }
    }
    if ( length(grep("x",log))>0 ) {
      aX <- axTicks(1); axis(1, at=aX, label= axTexpr(1, aX))
    }
    if ( length(grep("y",log))>0 ) {
      aY <- axTicks(2); axis(2, at=aY, label= axTexpr(2, aY))
    }
  } else {
    plot(x,y,log=log,type=type,pch=pch,xlab=colA,ylab=colB,main=main,cex=cex,col=col)
  }

  if ( !is.null(panel) ) {
    plot.panel.label(panel)
  }
  if ( value.name =="" ) {
    if ( par("ylog") ) {
      Sx <- min(x)+ceiling(del)
      Sy <- 10^par("usr")[4]
    } else {
      Sx <- min(x)+ceiling(del)
      Sy <- par("usr")[4]
    }
    text(paste("Rs=",s,sep=""), x=Sx,y=Sy,cex=1.3,pos=1,col="grey30")    
  }  
  if (smooth=="lowess") {    
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) {
        lines(stats::lowess(x[ok], y[ok]), 
            col = "red")
      }
  } else {
    abline(0,1,col="red")
  }
}

log.pairs <- function(x,...) {
   pairs(log(x), lower.panel=panel.smooth, upper.panel=panel.cor,pch=".",logit="xy")
}
# useful for pairs()
panel.cor.spearman <- function(...) {
  panel.cor(...,cor.method="spearman")
}

panel.cor <- function(x, y, cor.method="pearson", digits=2, prefix="", cex.cor=1.0,logit=NULL,...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    if(!is.null(logit)) {
      if ( logit=="xy" || logit=="x") { x <- exp(x) }
      if ( logit=="xy" || logit=="y") { y <- exp(y) }
    }
    # pearson correlation
    r <- cor(x,y,method="pearson")
    txt <- round(r,digits)
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)

    text(0.5, 0.5, txt,cex=cex.cor)
}

my.html.plot <- function(filename=NULL,html.dir=NULL,rel.dir=NULL,caption=NULL,
                           bg="white",
                           width=400,height=400,to.plot=NULL) {
  suppressPackageStartupMessages(library(R2HTML))
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
  plot.filename <- paste(html.dir,filename,sep="/")
  png(filename=plot.filename,width=width, height=height,
      bg=bg)
  to.plot()
  dev.off()
  pdf(file=paste(plot.filename,".pdf",sep=""),width=width, height=height)
  to.plot()
  dev.off()
  HTMLInsertGraph(paste(rel.dir,filename, sep=""), Caption = caption,WidthHTML=width)  
}

