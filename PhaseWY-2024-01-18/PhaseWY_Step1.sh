#!/bin/bash -l

# Version 2023-10-16
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

### Step 1. Get unmapped, low coverage, too high coverage, too high missingness, and repeat regions
mkdir $OUTDIR/1_target_region/$SCAFFOLD;
OUTDIR1=$OUTDIR/1_target_region/$SCAFFOLD;
mkdir $OUTDIR1/temp;

samtools depth -aa -r $SCAFFOLD -J -s $(cat $INDS | cut -f 5 | tail -n1) | cut -f2 \
> $OUTDIR1/temp/$SCAFF_NAME\_genome_coverage_base.txt;
NR=0;
for IND in $(cat $INDS | tail -n+2 | cut -f1); do \
  BAM=$(cat $INDS | tail -n+2 | awk -F'\t' -v IND=$IND '$1==IND {print $5}');
  base=$(find $OUTDIR1/temp/ -name "*_genome_coverage_*")
  NR=$(echo $NR + 1 | bc)
  samtools depth -aa -r $SCAFFOLD -J -s $BAM | cut -f 2- \
  > $OUTDIR1/temp/$SCAFF_NAME\_genome_coverage_$IND.txt;
  join -j 1 $base $OUTDIR1/temp/$SCAFF_NAME\_genome_coverage_$IND.txt \
  > $OUTDIR1/temp/$SCAFF_NAME\_genome_coverage_$NR.txt;
  rm $base $OUTDIR1/temp/$SCAFF_NAME\_genome_coverage_$IND.txt;  
done;

cat $(find $OUTDIR1/temp/ -name "*_genome_coverage_*") | tr ' ' '\t' | \
awk -F'\t' -v SCAFFOLD=$SCAFFOLD '{print SCAFFOLD"\t"$0}' | \
awk -F'\t' -v MIN_DP=$MIN_DP '{sum=0; missing=0; for(i=3; i<=NF; i++) \
{sum+=$i; if($i < MIN_DP) {missing++} } print $0"\t"sum/(NF-2)"\t"1-(missing/(NF-2))}'| \
awk -F'\t' -v MIN_MEAN=$MIN_MEAN -v MAX_MEAN=$MAX_MEAN -v MISSING=$MISSING \
'{if($(NF-1) >= MIN_MEAN && $(NF-1) <= MAX_MEAN && $NF >= MISSING) print $1"\t"($2-1)"\t"$2}' | \
bedtools merge > $OUTDIR1/$SCAFF_NAME\_target_depth.bed;

rm -r $OUTDIR1/temp;

if [ $REPEATS != NA ]; then
  #### Remove repeat regions from target depth
  cat $REPEATS | grep $SCAFFOLD \
  > $OUTDIR1/$SCAFF_NAME\_repeats.bed;
  bedtools subtract -a $OUTDIR1/$SCAFF_NAME\_target_depth.bed \
  -b $OUTDIR1/$SCAFF_NAME\_repeats.bed \
  > $OUTDIR1/$SCAFF_NAME\_target_region.bed;
  
  #### Remove temporary files
  rm $OUTDIR1/$SCAFF_NAME\_target_depth.bed \
  $OUTDIR1/$SCAFF_NAME\_repeats.bed;

  
else
  mv $OUTDIR1/$SCAFF_NAME\_target_depth.bed $OUTDIR1/$SCAFF_NAME\_target_region.bed;
fi;

#### Do not continue if there has been no alignments to the target region of the scaffold
ALIGNMENT=$(cat $OUTDIR1/$SCAFF_NAME\_target_region.bed | \
awk 'NR>0 {print; exit}' | wc -l);
if [ $ALIGNMENT -gt 0 ]; then
  TARGET_ALIGNMENT=TRUE;
else
  TARGET_ALIGNMENT=FALSE;
  echo No alignments to target regions in scaffold $SCAFFOLD. $SCAFFOLD will be dropped from the analysis.
fi;
