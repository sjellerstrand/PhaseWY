#!/bin/bash -l

# Version 2023-10-16
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

# Step 0. Setup anaysis

## Make sure all parameters are accepted. If something is wrong, abort program and print error-message
PARAMETERS=$(Rscript $FUNCTIONS/accept_parameters.r --no-save --args PHASE=$PHASE CDS_SCAFFOLDS_ONLY=$CDS_SCAFFOLDS_ONLY \
MIN_SCAFFOLD_SIZE=$MIN_SCAFFOLD_SIZE SEX_SYSTEM=$SEX_SYSTEM NONDIPLOIDS=$NONDIPLOIDS MIN_DP=$MIN_DP \
MIN_MEAN=$MIN_MEAN MAX_MEAN=$MAX_MEAN MISSING=$MISSING THRESHOLD=$THRESHOLD WINDOW=$WINDOW STEP=$STEP \
PHASE_SCAFFOLD_LIST=$PHASE_SCAFFOLD_LIST LIST_ONLY=$LIST_ONLY CDS=$CDS EXON=$EXON GENES=$GENES REPEATS=$REPEATS);

echo "$PARAMETERS" | head -n-1;
ABORT_PROGRAM=$(echo "$PARAMETERS" | tail -n1);
if [ $ABORT_PROGRAM == TRUE ]; then
  echo -e "\nExiting program\n"
  exit;
fi;

## Create folders
mkdir $OUTDIR/A8_PhaseWY;
OUTDIR=$OUTDIR/A8_PhaseWY;
mkdir $OUTDIR/0_file_info \
$OUTDIR/1_target_region \
$OUTDIR/2_sex_depth_difference \
$OUTDIR/3_scaffold_input \
$OUTDIR/4_whatshap \
$OUTDIR/5_shapeit4 \
$OUTDIR/6_sex_linkage \
$OUTDIR/7_concat_vcfs \
$OUTDIR/final_output \
$OUTDIR/Scaffold_progress;

## Setup file info
OUTDIR0=$OUTDIR/0_file_info;
if [ $SEX_SYSTEM != NA ]; then
  if [ $SEX_SYSTEM == ZW ]; then
    HOMGAM=Male;
    HETGAM=Female;
  elif [ $SEX_SYSTEM == XY ]; then
    HOMGAM=Female;
    HETGAM=Male;
  fi;
fi;

cp $SAMPLE_INFO $OUTDIR0/INDS.txt;
INDS=$SAMPLE_INFO;
cat $INDS | cut -f5 | tail -n+2 > $OUTDIR0/BAMS.txt;
BAMS=$OUTDIR0/BAMS.txt;
if [ $SEX_SYSTEM != NA ]; then
  cat $INDS | tail -n+2 | awk -F'\t' '{if($3=="HETGAM") print $5}' \
  > $OUTDIR0/HETGAM_BAMS.txt;
  fi;
cat $REF.fai | cut -f1,2 | sort -k2,2 -nr | \
awk -F'\t' -v MIN_SCAFFOLD_SIZE=$MIN_SCAFFOLD_SIZE '{if($2 >= MIN_SCAFFOLD_SIZE) print}' \
> $OUTDIR0/SCAFFOLDS.txt;

if [ $CDS_SCAFFOLDS_ONLY == Yes ]; then
  for SCAFFOLD0 in $(cat $CDS | cut -f1 | sort -u); do
    cat $OUTDIR0/SCAFFOLDS.txt | awk -F'\t' -v SCAFFOLD0=$SCAFFOLD0 '$1==SCAFFOLD0 {print}';
  done > $OUTDIR0/SCAFFOLDS2.txt;
  mv $OUTDIR0/SCAFFOLDS2.txt $OUTDIR0/SCAFFOLDS.txt;
fi;

if [ $LIST_ONLY == Yes ]; then
  for SCAFFOLD0 in $(cat $PHASE_SCAFFOLD_LIST | cut -f1 | sort -u); do
    cat $OUTDIR0/SCAFFOLDS.txt | awk -F'\t' -v SCAFFOLD0=$SCAFFOLD0 '$1==SCAFFOLD0 {print}';
  done > $OUTDIR0/SCAFFOLDS2.txt;
  mv $OUTDIR0/SCAFFOLDS2.txt $OUTDIR0/SCAFFOLDS.txt;

fi;

#### Check if there are only haploids in the dataset
if [ $NONDIPLOIDS == Yes ]; then
  if [ $(cat $INDS | tail -n+2 | awk -F'\t' '{if($4==2) print}' | wc -l) == 0 ]; then
    HAPLOIDS_ONLY=TRUE;
    echo All individuals are already haploid. No need for phasing.
  else
    HAPLOIDS_ONLY=FALSE;
  fi;
else
  HAPLOIDS_ONLY=FALSE;
fi;
export HAPLOIDS_ONLY;

#### Calculate number of genomes in the dataset and set mac thresholds
GENOMES=$(cat $INDS | awk -F'\t' '{sum+=$4} END {print sum}')