#!/usr/bin/Rscript

# Version 2023-02-28
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

### This script evaluates standart sex differences in depth per site across a given scaffold

## Export variables and load libraries
rm(list=ls())

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

## Get the average normalised depth for homogametes
HOMGAM <- list.files(OUTDIR2, pattern="HOMGAM", full.names=T)
HOMGAM_inds <- length(HOMGAM)
HOMGAM_depth <- read.delim(HOMGAM[1], sep="\t", header=F)
for(i in 2:HOMGAM_inds) {
 HOMGAM_depth <- cbind(HOMGAM_depth[,1:2], (HOMGAM_depth[,3] + read.delim(HOMGAM[i], sep="\t", header=F)[,3]))
}
HOMGAM_depth[,3] <- HOMGAM_depth[,3]/HOMGAM_inds

## Get the average normalised depth for heterogametes
HETGAM <- list.files(OUTDIR2, pattern="HETGAM", full.names=T)
HETGAM_inds <- length(HETGAM)
HETGAM_depth <- read.delim(HETGAM[1], sep="\t", header=F)
for(i in 2:HETGAM_inds) {
  HETGAM_depth <- cbind(HETGAM_depth[,1:2], HETGAM_depth[,3] + read.delim(HETGAM[i], sep="\t", header=F)[,3])
}
HETGAM_depth[,3] <- HETGAM_depth[,3]/HETGAM_inds

## Calculate the normalised depth ratio per site between sexes and correct divisions by zero
SEX_DIFF <- cbind(HOMGAM_depth[,1:2], HETGAM_depth[,3]/HOMGAM_depth[,3])
SEX_DIFF[which(SEX_DIFF[,3] == "Inf"),3] <- 1

## Write per site bed file
write.table(SEX_DIFF, file=paste(OUTDIR2, "/", PROJECT, "_", SCAFFOLD, "_sex_diff.txt", sep=""), quote=FALSE, sep='\t', row.names = F, col.names = F)
q()
