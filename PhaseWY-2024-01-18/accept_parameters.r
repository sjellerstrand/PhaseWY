#!/usr/bin/Rscript

# Version 2023-10-16
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand


### This script goes through the input parameters and aborts the pipeline if any non-accepted values are provided.

## Export variables and load libraries
rm(list=ls())

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

## Make numeric
MIN_SCAFFOLD_SIZE <- as.numeric(MIN_SCAFFOLD_SIZE)
MIN_DP <- as.numeric(MIN_DP)
MIN_MEAN <- as.numeric(MIN_MEAN)
MAX_MEAN <- as.numeric(MAX_MEAN)
MISSING <- as.numeric(MISSING)
THRESHOLD <- as.numeric(THRESHOLD)
WINDOW <- as.numeric(WINDOW)
STEP <- as.numeric(STEP)

## If something is wrong, abort program and print error-message
ABORT <- "FALSE"
OUTPUT <- NA

### Check if phasing should be performed
if(PHASE != "No" && PHASE != "Yes") {
  OUTPUT <- c(OUTPUT, paste('Should phasing be performed? Alternatively, only output callable sites (step 1). "Yes" or "No" are accepted. Currently you have specified: PHASE=', PHASE, sep=""))
  ABORT="TRUE"
}

### Check if analysis is only performed on scaffolds annotated with coding region
if(CDS_SCAFFOLDS_ONLY != "No" && CDS_SCAFFOLDS_ONLY != "Yes") {
  OUTPUT <- c(OUTPUT, paste('Should analysis only be performed on scaffolds with annotated coding regions?, "No" or "Yes" are accepted. Currently you have specified: CDS_SCAFFOLDS_ONLY=', CDS_SCAFFOLDS_ONLY, sep=""))
  ABORT="TRUE"
}

### Check the minimum scaffold length
if(is.na(MIN_SCAFFOLD_SIZE) || MIN_SCAFFOLD_SIZE < 0 || round(MIN_SCAFFOLD_SIZE) != MIN_SCAFFOLD_SIZE) {
  OUTPUT <- c(OUTPUT, paste('The minimum scaffolds length for analysis must be set as an integer with a positive value. Currently you have specified: MIN_SCAFFOLD_SIZE=', MIN_SCAFFOLD_SIZE, sep=""))
  ABORT="TRUE"
}

### Check sex-chromosome system. "ZW", "XY", or "NA" are accepted.
if(SEX_SYSTEM != "ZW" && SEX_SYSTEM != "XY" && SEX_SYSTEM != "NA") {
  OUTPUT <- c(OUTPUT, paste('No valid sex-determination system. "ZW", "XY", or "NA" are accepted. Currently you have specified: SEX_SYSTEM=', SEX_SYSTEM, sep=""))
  ABORT="TRUE"
}

### Check if there are non-diploids
if(NONDIPLOIDS != "No" && NONDIPLOIDS != "Yes") {
  OUTPUT <- c(OUTPUT, paste('Are there any non-diploids (i.e. haploids)? "Yes" or "No" are accepted. Currently you have specified: NONDIPLOIDS=', NONDIPLOIDS, sep=""))
  ABORT="TRUE"
}

### Check the minimum accepted depth
if(is.na(MIN_DP) || is.na(MIN_DP) || MIN_DP < 0) {
  OUTPUT <- c(OUTPUT, paste('The minimum accepted depth must be a number of at least 0. Currently you have specified: MIN_DP=', MIN_DP, sep=""))
  ABORT="TRUE"
}

### Check the minimum accepted mean depth
if(!is.na(MIN_DP) && (is.na(MIN_MEAN) ||  is.na(MIN_MEAN) || MIN_MEAN < MIN_DP)) {
  OUTPUT <- c(OUTPUT, paste('The maximum accepted mean depth should be a number of at least 0 and at least as large as MIN_DP. Currently you have specified: MIN_MEAN=', MIN_MEAN, sep=""))
  ABORT="TRUE"
}

### Check the maximum accepted mean depth
if(!is.na(MIN_MEAN) && (is.na(MAX_MEAN) || is.na(MAX_MEAN) || MAX_MEAN <= 0 || MAX_MEAN < MIN_MEAN)) {
  OUTPUT <- c(OUTPUT, paste('The maximum accepted mean depth should be a number larger than 0 and larger than MIN_MEAN. Currently you have specified: MAX_MEAN=', MAX_MEAN, sep=""))
  ABORT="TRUE"
}

### Check the minimum accepted missingness
if(is.na(MISSING) || is.na(MISSING) || MISSING < 0 || MISSING > 1) {
  OUTPUT <- c(OUTPUT, paste('The minimum accepted missingness should be a ratio of at least 0 and no larger than 1. Currently you have specified: MISSING=', MISSING, sep=""))
  ABORT="TRUE"
}

### Check the ratio of heterogametes needed in smallest cluster for it to be called sex linked.
if((SEX_SYSTEM == "ZW" || SEX_SYSTEM == "XY") && (is.na(THRESHOLD) || THRESHOLD > 1 || THRESHOLD <= 0)) {
  OUTPUT <- c(OUTPUT, paste('When a sex chromsome system is specified, the threshold for identying sex-linkage must be set as a ratio larger than 0 and no larger than 1. Currently you have specified: THRESHOLD=', THRESHOLD, sep=""))
  ABORT="TRUE"
}

### Check the number of variants used per window
if((SEX_SYSTEM == "ZW" || SEX_SYSTEM == "XY") && (is.na(WINDOW) || WINDOW == 0 || round(WINDOW) != WINDOW)) {
  OUTPUT <- c(OUTPUT, paste('When a sex chromsome system is specified, the window size for identying sex-linkage must be set as an integer larger than 0. Currently you have specified: WINDOW=', WINDOW, sep=""))
  ABORT="TRUE"
}

### Check the number of variants for each sliding window step
if((SEX_SYSTEM == "ZW" || SEX_SYSTEM == "XY") && (!is.na(WINDOW) && WINDOW > 0 && round(WINDOW) == WINDOW) && (is.na(STEP) || STEP == 0 || round(STEP) != STEP || STEP > WINDOW)) {
  OUTPUT <- c(OUTPUT, paste('When a sex chromsome system is specified, the sliding window step size for identying sex-linkage must be set as an integer larger than larger than 0 and not larger than the window. Currently you have specified: STEP=', STEP, sep=""))
  ABORT="TRUE"
}

### Check if a list of scaffolds to phase is provided
if(PHASE_SCAFFOLD_LIST != "NA" && file.exists(PHASE_SCAFFOLD_LIST) == FALSE) {
    OUTPUT <- c(OUTPUT, paste('List of scaffolds to be phased. "NA" or path to list accepted. Currently you have specified neither NA, nor a path to an existing file: PHASE_SCAFFOLD_LIST=', PHASE_SCAFFOLD_LIST, sep=""))
    ABORT="TRUE"
} else if(PHASE_SCAFFOLD_LIST != "NA" && PHASE == "No") {
    OUTPUT <- c(OUTPUT, paste('List of scaffolds to be phased is only accepted when PHASE=Yes. Currently you have specified: PHASE=No', sep=""))
    ABORT="TRUE"
}

### Check if analysis is only performed on scaffolds in the provided list "PHASE_SCAFFOLD_LIST"
if(LIST_ONLY != "No" && LIST_ONLY != "Yes") {
  OUTPUT <- c(OUTPUT, paste('Should analysis only be performed on listed scaffolds given by "PHASE_SCAFFOLD_LIST"? "Yes" or "No" are accepted. Currently you have specified: LIST_ONLY=', LIST_ONLY, sep=""))
  ABORT="TRUE"
} else if(PHASE_SCAFFOLD_LIST == "NA" && LIST_ONLY == "Yes") {
  OUTPUT <- c(OUTPUT, paste('Should analysis only be performed on listed scaffolds given by "PHASE_SCAFFOLD_LIST"? "Yes" or "No" are accepted. Currently you have specified: LIST_ONLY=Yes, but not provided a list of scaffolds.', sep=""))
  ABORT="TRUE"
}

### Check if a file with CDS regions is provided
if(CDS != "NA" && file.exists(CDS) == FALSE) {
    OUTPUT <- c(OUTPUT, paste('Bed file with CDS regions. "NA" or path to bed file accepted. Currently you have specified neither NA, nor a path to an existing file: CDS=', CDS, sep=""))
    ABORT="TRUE"
}

### Check if a file with exonic regions is provided
if(EXON != "NA" && file.exists(EXON) == FALSE) {
    OUTPUT <- c(OUTPUT, paste('Bed file with exonic regions. "NA" or path to bed file accepted. Currently you have specified neither NA, nor a path to an existing file: EXON=', EXON, sep=""))
    ABORT="TRUE"
}

### Check if a file with annotated features is provided
if(GENES != "NA" && file.exists(GENES) == FALSE) {
    OUTPUT <- c(OUTPUT, paste('Bed file with annotated features. "NA" or path to bed file accepted. Currently you have specified neither NA, nor a path to an existing file: GENES=', GENES, sep=""))
    ABORT="TRUE"
}

### Check if a file with annotated features is provided
if(REPEATS != "NA" && file.exists(REPEATS) == FALSE) {
    OUTPUT <- c(OUTPUT, paste('Bed file with repeat regions. "NA" or path to bed file accepted. Currently you have specified neither NA, nor a path to an existing file: REPEATS=', REPEATS, sep=""))
    ABORT="TRUE"
}

### Check if program should abort
OUTPUT <- c(OUTPUT, ABORT)

### Write output
OUTPUT <- OUTPUT[-1]
writeLines(OUTPUT)

## Exit program
q()