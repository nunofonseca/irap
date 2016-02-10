#!/usr/bin/env Rscript
# args = package names
args <- commandArgs(trailingOnly=TRUE)

for ( p in args ) {
  v <- try(packageVersion(p))
  cat(paste(p,"_VERSION=",v,"\n",sep=""))
}
q(status=0)
