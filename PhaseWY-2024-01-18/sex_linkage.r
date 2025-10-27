#!/usr/bin/Rscript

# Version 2023-10-11
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

## Export variables and load libraries
rm(list=ls())
library(vcfR)
library(data.table)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

options(scipen=999)

THRESHOLD <- as.numeric(THRESHOLD)
WINDOW <- as.numeric(WINDOW)
STEP <- as.numeric(STEP)
SCAFFOLD_LENGTH <- as.numeric(SCAFFOLD_LENGTH)

### Import data
vcf1 <- read.vcfR(paste(OUTDIR6, "/", PROJECT, "_", SCAFFOLD, "_phased_all_variants_mac2.vcf.gz", sep=""), verbose = F)
if(nrow(extract.gt(vcf1)) == 1) {
  vcf2 <- cbind(t(as.matrix(getFIX(vcf1))), as.matrix(extract.gt(vcf1, element="GT")))
} else {
  vcf2 <- cbind(getFIX(vcf1), extract.gt(vcf1, element="GT"))
}
rm(vcf1)
sample_info <- read.table(INDS, sep='\t', head=T)

# Find heterogametic individuals
heterogametes <- sample_info[which(sample_info$Sex == "HETGAM"),1]
homogametes <- sample_info[which(sample_info$Sex == "HOMGAM"),1]

# Seperate haplotypes
haplotypes <- as.data.frame(t(cbind(matrix(vcf2[,1]), matrix(vcf2[,2]), matrix(NA, nrow(vcf2), (ncol(vcf2)-7)*2))))
rownames(haplotypes)[c(1,2)] <- c("Scaffold", "Pos")
for(i in 1:(ncol(vcf2)-7)) {
  rownames(haplotypes)[c(2+(i*2-1),2+(i*2))] <- c(paste(colnames(vcf2)[7+i], "_left", sep="") ,paste(colnames(vcf2)[7+i], "_right", sep=""))
  haplotypes[2+(i*2-1),] <- unlist(strsplit(vcf2[,7+i], "|", fixed=T))[c(T,F)]
  haplotypes[2+(i*2),] <- unlist(strsplit(vcf2[,7+i], "|", fixed=T))[c(F,T)]
}

# Determine sliding window parameters
if(SCAFFOLD_LENGTH == 0) { ## Kommer det ens gå att ladda in isåfall?  ####################
  ###### Kill script and write info about linkage etc. Fundera på om det behöver sättas en högre gräns än 0, och om dem ska behandlas annorlunda ###################
  quit() ### FUnkar detta?
} else if(SCAFFOLD_LENGTH < WINDOW*4) {
  WINDOW0 <- WINDOW
  STEP0 <- STEP
  WINDOW <- ceiling(SCAFFOLD_LENGTH/4)
  STEP <- ceiling((STEP0*SCAFFOLD_LENGTH)/(4*WINDOW0))
  writeLines(paste("Scaffold is ", SCAFFOLD_LENGTH, " bp long, which is few relative to the defined window size of ", WINDOW0, ". The initial window size has been reduced to ", WINDOW, " with a sliding step of ", STEP, ".", sep=""))
}
if(SCAFFOLD_LENGTH == 1) {
  startpos <- 1
  endpos <- 1
} else {
  startpos <- seq(1, SCAFFOLD_LENGTH - WINDOW, STEP)
  endpos <- startpos + WINDOW - 1
}
if(endpos[length(endpos)] < SCAFFOLD_LENGTH) {
  startpos <- c(startpos, SCAFFOLD_LENGTH - WINDOW + 1)
  endpos <- c(endpos, SCAFFOLD_LENGTH)
}
WIND_OVERLAP <- (WINDOW/STEP) - 1
positions <- as.numeric(haplotypes[2,])
haplotype_names <- rownames(haplotypes)[3:nrow(haplotypes)]

### Evaluate windows for sex linkage with K-means
sexlink <- matrix(NA, 1, (17+length(heterogametes)))
colnames(sexlink) <- c("Start pos [bp]", "End pos [bp]", "Sex-linked status", "Number of variants in window", "Ratio of total heterogametic individuals in smallest cluster [#inds/#inds]", "Ratio of haplotypes in smallest cluster [#haps/#haps]", "Smallest cluster [#inds]", "Largest cluster [#inds]", "TotSS", "Smallest cluster WithinSS", "Largest cluster WithinSS", "BetweenSS", "Window length [bp]", "Variants/bp", "Window mid position [bp]", "Phase information available", "Border change", heterogametes)
i <- startpos[1]
stop_search <- NA
phase_switch <- NA
border2 <- NA
while(is.na(stop_search)) {
  
  # Reset some parameters
  i0 <- i
  status <- NA
  phase <- matrix("Unknown", 1, length(heterogametes))
  colnames(phase) <- colnames(sexlink)[18:ncol(sexlink)]

  # Evaluate K-means
  snps <- as.data.frame(haplotypes[3:nrow(haplotypes), which(between(positions, startpos[i], endpos[i]))], row.names=haplotype_names)
  distmat <- dist(snps)
  if(length(which(distmat > 0)) == 0) {
    status <- "Autosomal: no variation in window"
  } else {
    km <- kmeans(distmat, 2, nstart=10)
    
    # Find smallest cluster
    if(km$size[1] == km$size[2]) {
      min <- 1
      max <- 2
    } else {
      min <- which(km$size == min(km$size))
      max <- which(km$size == max(km$size))
    }
    min_cluster_temp <- strsplit(names(km$cluster)[which(km$cluster == min)], "_", fixed=T)
    min_cluster <- rbind(sapply(lapply(min_cluster_temp, head, -1), paste, collapse="_"), sapply(min_cluster_temp, tail, 1))
    
    # Are any homogametes clustered in the smallest cluster?
    for(j in 1:length(homogametes)) {
      if(length(which(min_cluster[1,] == homogametes[j])) > 0) {
        status <- "Autosomal: homogamete in cluster"
        break
      }
    }
    
    # Count and evaluate heterogamete haplotypes in the smallest cluster
    heterogametes_count <- 0
    for(j in 1:length(heterogametes)) {
      if(length(which(min_cluster[1,] == heterogametes[j])) > 0) {
        
        # Do both haploypes of any heterogametes occur in cluster?
        if(length(which(min_cluster[1,] == heterogametes[j])) > 1 && is.na(status)) {
          status <- "Autosomal: both haplotypes of heterogamete in cluster"
        }
        heterogametes_count <- heterogametes_count + 1
      }
    }
    
    # Does the fraction of heterogametes in the smallest cluster meet the set threshold for it to be classed as sex-linked?
    if(heterogametes_count/length(heterogametes) >= THRESHOLD && is.na(status)) {
      status <- "Sex-linked"
    } else if(is.na(status)) {
      status <- "Autosomal: too few heterogametes in cluster"
    }
  }
  # If the current window is sex-linked, do the following:
  if(status == "Sex-linked") {
    # Find phase
    for(j in 1:ncol(min_cluster)) {
      phase[1,which(colnames(phase) == min_cluster[1,j])] <- min_cluster[2,j]
    }
  }

  # Register data for window
  sexlink[i,1] <- startpos[i]
  sexlink[i,2] <- endpos[i]
  sexlink[i,3] <- status
  sexlink[i,4] <- ncol(snps)
  if(status == "Autosomal: no variation in window") {
    sexlink[i,5:12] <- NA
  } else {
    sexlink[i,5] <- heterogametes_count/length(heterogametes)
    sexlink[i,6] <- km$size[min]/(km$size[min] + km$size[max])
    sexlink[i,7] <- km$size[min]
    sexlink[i,8] <- km$size[max]
    sexlink[i,9] <- km$totss
    sexlink[i,10] <- km$withinss[min]
    sexlink[i,11] <- km$withinss[max]
    sexlink[i,12] <- km$betweenss
  }
  sexlink[i,13] <-  endpos[i] - startpos[i] + 1
  sexlink[i,14] <- as.numeric(sexlink[i,4]) / as.numeric(sexlink[i,13])
  sexlink[i,15] <- as.numeric(startpos[i] + (as.numeric(sexlink[i,13])-1)/2)
  if(length(which(is.na(phase))) == 0 && length(which(phase == "Unknown")) == 0) {
    sexlink[i,16] <- "Yes"
  } else {
    sexlink[i,16] <- "No"
  }
  if(i > 1 && (all(sexlink[i-1, 18:ncol(sexlink)] != phase[1,]) || !is.na(border2))) {
    if(all(sexlink[i-1, 18:ncol(sexlink)] != phase[1,])) {
      if(i-WIND_OVERLAP < 0) {
        border1 <- 1
      } else {
        border1 <- floor(i-WIND_OVERLAP)
      }
      sexlink[border1:i, 17] <- "Yes"
      border2 <- ceiling(i+WIND_OVERLAP)
    }
    else if(!is.na(border2) && i <= border2) {
      sexlink[i, 17] <- "Yes"
      if(i == border2) {
        border2 <- NA
      }
    }
  }
  sexlink[i,18:ncol(sexlink)] <- phase[1,]
  
  # Continue loop?
  if(endpos[i] == endpos[length(endpos)]) {
    stop_search <- "Stop"
  } else if(i != i0) {
    i <- i0
  } else {
    sexlink <- rbind(sexlink, matrix(NA, 1, (17+length(heterogametes))))
    i <- i0+1
  }
}

# Summarize coherent regions
sexlink[is.na(sexlink[,17]), 17] <- "No"
regions <- matrix(NA, 1, (4+length(heterogametes)))
colnames(regions) <- c("chrom", "chromStart", "chromEnd", "Status", heterogametes)

# If there is any sex-linkage
if(length(which(sexlink[,3] == "Sex-linked")) > 0) {
  stop_search <- NA
  i <- 1
  
  while(is.na(stop_search)) {
    i0 <- i
    beginreg2 <- NA
    endreg2 <- NA
    status3 <- NA
    
    while(sexlink[i0, 17] == sexlink[i, 17] && i != nrow(sexlink)) {
      i <- i + 1
    }
    
    # Define sex-linkage status of region
    if(sexlink[i0, 17] == "Yes") {
      status3 <- "Unknown"
    } else if(sexlink[i-1,3] == "Sex-linked") {
      status3 <- "Sex-linked"
    } else {
      status3 <- "Pseudoautosomal"
    }
    
    # Define start of region
    if(i0 == 1) {
      beginreg2 <- 0
    } else if(status3 == "Unknown") {
      beginreg2 <- as.numeric(regions[nrow(regions)-1,3]) + 1 - 1
    } else {
      beginreg2 <- as.numeric(sexlink[i0,1]) - 1
    }
    
    # Define end of region
    if(i == nrow(sexlink) && sexlink[i0, 17] == sexlink[i, 17]) {
      endreg2 <- SCAFFOLD_LENGTH
      if(status3 == "Unknown") {
        regions[nrow(regions),] <- c(SCAFFOLD, beginreg2, endreg2, status3, rep("Unknown", length(heterogametes)))
      } else {
        regions[nrow(regions),] <- c(SCAFFOLD, beginreg2, endreg2, status3, sexlink[i-1,18:ncol(sexlink)])
      }
    } else if(i == nrow(sexlink) && sexlink[i0, 17] != sexlink[i, 17]) {
      if(status3 == "Unknown") {
        endreg2 <- as.numeric(sexlink[i,1]) - 1
        regions[nrow(regions),] <- c(SCAFFOLD, beginreg2, endreg2, status3, rep("Unknown", length(heterogametes)))
        regions <- rbind(regions, c(SCAFFOLD, as.numeric(endreg2)+1, SCAFFOLD_LENGTH, sexlink[i,3], sexlink[i,18:ncol(sexlink)]))
      } else {
        endreg2 <- as.numeric(sexlink[i-1,2])
        regions[nrow(regions),] <- c(SCAFFOLD, beginreg2, endreg2, status3, sexlink[i-1,18:ncol(sexlink)])
        if(status3 == "Unknown") {
          regions[nrow(regions),] <- c(SCAFFOLD, beginreg2, endreg2, status3, rep("Unknown", length(heterogametes)))
        }
      }
    } else if(status3 == "Unknown") {
      endreg2 <- as.numeric(sexlink[i,1])-1
      regions[nrow(regions),] <- c(SCAFFOLD, beginreg2, endreg2, status3, rep("Unknown", length(heterogametes)))
    } else {
      endreg2 <- as.numeric(sexlink[i-1,2])
      regions[nrow(regions),] <- c(SCAFFOLD, beginreg2, endreg2, status3, sexlink[i-1,18:ncol(sexlink)])
    }

    # Continue loop?
    if(i == nrow(sexlink)) {
      stop_search <- "Stop"
    } else {
      regions <- rbind(regions, matrix(NA, 1, (4+length(heterogametes))))
    }
  }
# If there is no evidence of sex-linkage, write whole scaffold as autosomal
} else {
  regions[1,] <- c(SCAFFOLD, 1, SCAFFOLD_LENGTH, "Autosomal", rep("Unknown", length(heterogametes)))
}
if(is.matrix(regions)) {
  regions <- as.data.frame(regions)
} else {
  regions <- t(as.data.frame(regions))
}

colnames(regions)[1] <- paste("#", colnames(regions)[1], sep="")

# Write phase and sex-linkage information from search as bed files
sexlink <- cbind("#chrom"=SCAFFOLD, sexlink)
write.table(sexlink, file=paste(OUTDIR6, "/", PROJECT, "_", SCAFFOLD, "_phase_windows.txt", sep=""), quote=FALSE, sep='\t', row.names = F, col.names = T)
write.table(regions, file=paste(OUTDIR6, "/", PROJECT, "_", SCAFFOLD, "_phase_info.bed", sep=""), quote=FALSE, sep='\t', row.names = F, col.names = T)
if(length(which(regions[,4] == "Sex-linked")) > 0) {
  for(IND in colnames(regions)[5:length(colnames(regions))]) {
    het_left <- regions[which(regions[,which(colnames(regions)==IND)] == "left"), 1:3]
    het_right <- regions[which(regions[,which(colnames(regions)==IND)] == "right"), 1:3]
    write.table(het_left, file=paste(OUTDIR6, "/", PROJECT, "_", SCAFFOLD, "_", IND, "_het_left.bed", sep=""), quote=FALSE, sep='\t', row.names = F, col.names = T)
    write.table(het_right, file=paste(OUTDIR6, "/", PROJECT, "_", SCAFFOLD, "_", IND, "_het_right.bed", sep=""), quote=FALSE, sep='\t', row.names = F, col.names = T)

  }
}
quit()
