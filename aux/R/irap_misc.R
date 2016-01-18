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
fold.change.plot <- function(res,y.lab="log2(Fold change)",x.lab="log(Mean)",fdr=0.1,log="",ylim=NA,main=NULL) {

  if (!(is.data.frame(res) && all(c("baseMean", "log2FoldChange") %in% colnames(res))))
    stop("'res' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")

  linecol = "#ff000080"
  linecol = makeTransparent("grey10")

  x = subset(res, baseMean!=0)
  py = x$log2FoldChange

  if(missing(ylim))
    ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1

  if (missing(x.lab))
    x.lab <- paste(x.lab," [FDR=",fdr,";DE genes=",sum(res$padj < fdr)," (",sum(res$padj[res$log2FoldChange>0]<fdr),"/",sum(res$padj[res$log2FoldChange<0]<fdr),")]",sep="")
  vals<-round(res$padj*100)
  #colors<-append(colorRampPalette(c("darkred","orange"))(10),rep("darkgrey",90))
  colors<-append(colorRampPalette(c("red"))(10),rep("darkgrey",90))
  xvals <- log(res$baseMean)
  plot(
    xvals,
    res$log2FoldChange,
    log=log,pch=20,cex=0.3,
    xlab=x.lab,ylab=y.lab,
    main=main,
    col = ifelse(res$padj <= fdr,"red","darkgrey") )  
  abline(h=0, lwd=3, col=linecol)

  # show the fdr
  if ( par("ylog") ) {
    Sx <- 10^(par("usr")[1]+(par("usr")[2]-10^par("usr"))*0.2)
    Sy <- 10^par("usr")[4]
  } else {
    Sx <- par("usr")[1]+(par("usr")[2]-par("usr")[1])*0.2
    Sy <- par("usr")[4]
  }
  text(paste("FDR=",fdr,sep=""), x=Sx,y=Sy,cex=1.1,pos=1,col="red")
}
#       col = ifelse(res$padj < fdr,"red","black") )
#  plot(x$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
#       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
#       cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
cor.scatterplot <- function(df=NULL,colA=NULL,colB=NULL,x=NULL,y=NULL,
                            value.name="",smooth="",log="xy",panel=NULL,
                            print.cor=TRUE) {
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
  if ( value.name =="" && print.cor ) {
    if ( par("ylog") ) {
      Sx <- 10^(par("usr")[1]+(par("usr")[2]-10^par("usr")[1])*0.2)
        #min(x)+ceiling(del)
      Sy <- 10^par("usr")[4]
    } else {
      Sx <- min(x)+ceiling(del)
      Sy <- par("usr")[4]
    }
    #pinfo("Sx",Sx)
    text(paste("Rs=",s,"\nRp=",p,sep=""), x=Sx,y=Sy,cex=1.3,pos=1,col="grey30")
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
    x1 <- par("xaxt")
    x2 <- par("yaxt")
    par(usr = c(0, 1, 0, 1),    xaxt='n',yaxt='n')
    if(!is.null(logit)) {
      if ( logit=="xy" || logit=="x") { x <- exp(x) }
      if ( logit=="xy" || logit=="y") { y <- exp(y) }
    }
    # pearson correlation
    r <- cor(x,y,method=cor.method)
    txt <- round(r,digits)
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)

    text(0.5, 0.5, txt,cex=cex.cor)
    par(usr = usr,    xaxt="s",yaxt="s")
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

makeTransparent<-function(someColor, alpha=100) {
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


# wrapper to the consensus function found in the agricolae package
# to plot coloured leafs 
irap.consensus <- function (leafsColours=NA,horiz=FALSE,text.cex=0.8,...) 
{

  colourTreeLabel <- function(n,cols) {
    if (is.leaf(n)) {
      a <- attributes(n)    
      labCol <- cols[a$label]
      attr(n, "nodePar") <- list(lab.col=labCol)
    }
    n
  }

  con<-consensus(...)
  n <- dendrapply(as.dendrogram(con$dendrogram), colourTreeLabel, leafsColours)
  # replot
  plot(n, cex=1, xlab="",yaxt='n', ann=FALSE,horiz=horiz)
  text(label=con$table.dend$percentage, con$table.dend$xaxis,con$table.dend$height, cex = text.cex, col = "red")
  con$n <- n
  return(con)
}

# rewrite to add the option of changing the key labels 
irap.heatmap.2 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                            key.xlabel = "Expression",
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    srtRow = NULL, srtCol = NULL, adjRow = c(0, NA), adjCol = c(NA, 
        0), offsetRow = 0.5, offsetCol = 0.5, key = TRUE, keysize = 1.5, 
    density.info = c("histogram", "density", "none"), denscol = tracecol, 
    symkey = min(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, 
    main = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, 
    lwid = NULL,...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    if (is.null(srtCol)) 
        axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
            offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
            padj = adjCol[2])
    else {
        if (is.numeric(srtCol)) {
            if (missing(adjCol) || is.null(adjCol)) 
                adjCol = c(1, NA)
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
                strheight("M"), labels = labCol, adj = adjCol, 
                cex = cexCol, srt = srtCol)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtCol ignored.")
    }
    if (is.null(srtRow)) {
        axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, 
            tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
    }
    else {
        if (is.numeric(srtRow)) {
            xpd.orig <- par("xpd")
            par(xpd = NA)
            ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                line = -0.5, tick = 0)
            text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
                y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
                srt = srtRow)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtRow ignored.")
    }
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, key.xlabel, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram",cex.main=0.8)
            par(cex = 0.5)
            mtext(side = 2, "Frequency", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}

