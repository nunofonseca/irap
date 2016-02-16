#!/usr/bin/env Rscript
# =========================================================
# Copyright 2012-2016,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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

###############################################################

IRAP.DIR <- Sys.getenv(c("IRAP_DIR"))
if ( IRAP.DIR == "" ) {
  cat("ERROR: environment variable IRAP_DIR is not set\n")
  q(status=1)
}
source(paste(IRAP.DIR,"aux/R","irap_utils.R",sep="/"))
             
args <- commandArgs(trailingOnly=TRUE)
conf.file <- args[1]
status.file <- args[2]
out.file <- args[3]
out.dir <- ""
# other variables
stages <- c("stage1","stage2","stage3","stage4")
print.html <- FALSE
pinfo("Usage: irap_report_status conf.file status.file.csv [file.html]")

# validate input
if ( is.na(status.file) || is.na(conf.file) ) {
  perror("Missing arguments")
  q(status=1)
}
file.required(status.file)
file.required(conf.file)

if ( is.na(out.file) ) {
  pinfo("No output file given, printing status to stdout")
  print.html <- FALSE
} else {
  pinfo("Storing status report in ",out.file)
  print.html <- TRUE
  out.file <- sub(".htm[l]*$","",out.file)
  out.dir <- dirname(out.file)
  out.file <- basename(out.file)
  # validate input
  if ( ! file.exists(out.dir) ) {
    perror("Directory ",out.dir," not found")
    q(status=1)
  }
}

conf <- load.configuration.file(conf.file)
import.conf.variables(conf)

data<-read.table(status.file,sep=",",comment.char="",quote="\"")
names(data) <- c("exp","stage","qc","mapper","quant","de","exp.files","comp.files","comp.perc")


# 
HTML.pprint.table <- function(table,caption="") {
  css.class <- table
  css.class[table==100] <- paste("<span class='done100'>",table[table==100],"</SPAN>",sep="")
  css.class[table<100 & table>=50] <- paste("<span class='donegt50'>",table[table<100 & table>=50],"</SPAN>",sep="")
  css.class[table>0 & table<50] <- paste("<span class='donelt50'>",table[table>0 & table<50],"</SPAN>",sep="")
  css.class[table==0] <- paste("<span class='done0'>",table[table==0],"</SPAN>",sep="")
  HTML(css.class,row.names=TRUE,caption=caption,align="left")
}

exps <- as.character(unique(data$exp))

if ( print.html ) {
  library(R2HTML)
  irap.HTMLInit(outdir=out.dir,filename=out.file, title=paste("",sep=""))
#  HTMLInitFile( outdir=out.dir, filename=out.file, title=paste("",sep=""), CSSFile=basename(irap_css_file), useGrid=FALSE)
  system(paste("cp ",irap_css_file," ",out.dir,sep=""));
  pinfo("HTML file opened")
  HTML.pprint.conf(conf,exps)
  HTML( "<h1>Completion Status</h1>" )
}

########################################1
pinfo("overall status")

# summary.table per stage
sum0 <-aggregate(x=data[,c("exp.files","comp.files","comp.perc")], by=list(exp=data$exp,stage=data$stage),FUN=sum)
sum0[,"comp.perc"] <- round(sum0[,"comp.files"]*100/sum0[,"exp.files"])

sum0.perc <- sum0[,c("exp","stage","comp.perc")]
sum0.perc.table <- xtabs(comp.perc~exp+stage, data=sum0.perc)


print(sum0.perc.table)
if ( print.html ) {
  HTML.pprint.table(sum0.perc.table,"Overall Status")
}
##############################
# stage1

pe.libs <- conf.get.value(conf,"pe")
se.libs <- conf.get.value(conf,"se")
libs <- append(pe.libs,se.libs)
#pinfo(libs)
df <- data.frame(matrix(data=0,ncol=length(libs),nrow=length(exps)))
colnames(df) <- libs
rownames(df) <- exps
for ( n in exps ) {
  for (l in se.libs) {    
    path2file <- paste(n,"/data/",l,".f.fastq",sep="")
    #pinfo(path2file)
    if ( file.exists(path2file) ) {
      df[n,l] <- 100
    } else {
      pwarning("File ",path2file," not found")
    }
  }
}
print(df)
if ( print.html ) {
  HTML.pprint.table(df,"Stage1 Status")
}

##############################
# mappers (stage=2)
pinfo("Stage2 status")

stage2.data <- data[data[,"stage"]=="stage2",]
sum2 <-aggregate(x=stage2.data[,c("exp.files","comp.files","comp.perc")], by=list(exp=stage2.data$exp,stage=stage2.data$stage,mapper=stage2.data$mapper),FUN=sum)
sum2.perc <- sum2[,c("exp","mapper","comp.perc")]
sum2.perc$mapper<-factor(sum2.perc$mapper)
sum2.perc.table <- xtabs(comp.perc~exp+mapper, data=sum2.perc)

print(sum2.perc.table)
if ( print.html ) {
  HTML.pprint.table(sum2.perc.table,"Stage2 Status")
}

# quant (stage=3)
pinfo("Stage3 status")
stage3.data <- data[data[,"stage"]=="stage3",]
sum3 <-aggregate(x=stage3.data[,c("exp.files","comp.files","comp.perc")], by=list(exp=stage3.data$exp,stage=stage3.data$stage,quant=stage3.data$quant),FUN=sum)
sum3.perc <- sum3[,c("exp","quant","comp.perc")]
sum3.perc$quant<-factor(sum3.perc$quant)
sum3.perc.table <- xtabs(comp.perc~exp+quant, data=sum3.perc)

print(sum3.perc.table)
if ( print.html ) {
  HTML.pprint.table(sum3.perc.table,"Stage3 Status")
}

# quant (stage=4)
pinfo("Stage4 status")
stage4.data <- data[data[,"stage"]=="stage4",]
stage4.data <- stage4.data[!stage4.data[,"de"]=="none",]
if ( nrow(stage4.data) == 0 ) {
  pinfo("DE disabled?")
  if ( print.html ) {
    HTML("DE disabled")
  }  
} else {
  sum4 <-aggregate(x=stage4.data[,c("exp.files","comp.files","comp.perc")], by=list(exp=stage4.data$exp,stage=stage4.data$stage,de=stage4.data$de),FUN=sum)
  sum4.perc <- sum4[,c("exp","de","comp.perc")]
  sum4.perc$de<-factor(sum4.perc$de)
  sum4.perc.table <- xtabs(comp.perc~exp+de, data=sum4.perc)
  print(sum4.perc.table)
  if ( print.html ) {
    HTML.pprint.table(sum0.perc.table,"Stage4 Status")
  }
}

if ( print.html ) {
  irap.HTMLEnd()
}

q(status=0)
