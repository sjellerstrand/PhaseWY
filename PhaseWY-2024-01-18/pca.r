#!/usr/bin/Rscript

### Plot Pricipal component 1 and 2
## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(viridis)
library(gridExtra)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
assign(args[[i]][1], args[[i]][2])
}

## read in data
pca <- read.delim(paste(VCF_OUT, ".eigenvec", sep=""), sep=" ", head=F)
eigenval <- scan(paste(VCF_OUT, ".eigenval", sep=""))
pca <- pca[,-1]

## set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

## sort out the individual sample locations and geographic regions
info <- read.delim(INDS, header=T)
#order <- read.delim(paste(METADATA, "/order.txt", sep=""), header=F)
geo <- matrix(NA, nrow(pca), 5)
for(i in 1:nrow(pca)) {
geo[i,1:5] <-as.matrix(as.character(c(pca[i,1], info[which(pca[i,1] == info[,1]),c(1, 3, 4, 5)])))
}
geo <- as.data.frame(geo)
#geo[,3] <- factor(geo[,3], levels=rev(order[,1]))
colnames(geo) <- c("ID", "ID", "Sex","Library", "PCRfree")
pca <- cbind(pca, geo)

## calculate the cumulative sum of the percentage variance explained
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)
cumsum(pve$pve)

## Create plots
pca <- as_tibble(data.frame(pca, geo[,4]))
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") + ylab("Percentage variance explained") + theme_light()
b <- ggplot(pca, aes(PC1, PC2), col) + geom_point(size = 3, aes(color=Sex)) + coord_equal() + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + scale_color_viridis(discrete = T)

## Write plots as .jpeg
jpeg(paste(VCF_OUT, "_PCA.jpeg", sep=""), width=1000, height=1000, quality=100)
grid.arrange(a, b, nrow=2, ncol=1)
dev.off()
q()
