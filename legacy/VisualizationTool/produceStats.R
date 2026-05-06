#!/usr/bin/env Rscript

#------------------------------------------------------------
# author: Willy Bruhn
# Script to visualize some statistics of a VARUS-run.
#
#
#------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Missing parameters for produceStats.R", call.=FALSE)
}

file=args[1]
outPath=args[2]
# yLIM=200
# if(length(args)>=3){
#   yLIM=args[3]
# }


dist = read.csv2(file)

pdf(paste(outPath,"stats.pdf",sep=""))
barplot(dist$totalObservations, ylim = c(0,10000))
#barplot(dist$totalObservations)
dev.off()
