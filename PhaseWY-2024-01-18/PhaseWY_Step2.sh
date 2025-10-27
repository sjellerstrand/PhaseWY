#!/bin/bash -l

# Version 2023-07-03
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

### Step 2. Calculate normalized depth difference between sexes (heterogametic sex chromosome regions have not aligned to homogametic reference)

mkdir $OUTDIR/2_sex_depth_difference/$SCAFFOLD;
OUTDIR2=$OUTDIR/2_sex_depth_difference/$SCAFFOLD;

#### Get normalised individual depth based on mean vcf depth and bam depth
for IND in $(cat $INDS | tail -n+2 | cut -f1); do \
  IND_DP=$(cat $INDS | tail -n+2 | awk -F'\t' -v IND=$IND '$1==IND {print $2}');
  SEX=$(cat $INDS | tail -n+2 | awk -F'\t' -v IND=$IND '$1==IND {print $3}');
  BAM=$(cat $INDS | tail -n+2 | awk -F'\t' -v IND=$IND '$1==IND {print $5}');
  samtools depth -aa -b $OUTDIR1/$SCAFF_NAME\_target_region.bed \
  -r $SCAFFOLD -J -s $BAM | awk -F'\t' -v IND_DP=$IND_DP '{print $1"\t"$2"\t"($3/IND_DP)}' \
  > $OUTDIR2/$SCAFF_NAME\_$IND\_$SEX\_depth.txt;
done;

#### Calculate normalised depth ratio per site between sexes
Rscript $FUNCTIONS/sex_depth_diff_scaffold.r --no-save --args PROJECT=$PROJECT OUTDIR2=$OUTDIR2 \
SCAFFOLD=$SCAFFOLD;

#### Extract sex-linked regions due to sex depth differences
DEPTH_DIFF=0.75;
cat $OUTDIR2/$PROJECT\_$SCAFFOLD\_sex_diff.txt | \
awk -F'\t' -v DEPTH_DIFF=$DEPTH_DIFF '{if($3 < DEPTH_DIFF) print $1"\t"$2-1"\t"$2}' | \
bedtools sort -g $REF.fai | bedtools merge \
> $OUTDIR2/$PROJECT\_$SCAFFOLD\_hetgam_dropout.bed;

#### Remove temporary files
rm $OUTDIR2/$SCAFF_NAME\_*_*GAM_depth.txt \
$OUTDIR2/$PROJECT\_$SCAFFOLD\_sex_diff.txt;

#### Remove bed file if there are no heterogametic drop-out regions
if [ $(cat $OUTDIR2/$PROJECT\_$SCAFFOLD\_hetgam_dropout.bed | wc -l) -eq 0 ]; then
  rm $OUTDIR2/$PROJECT\_$SCAFFOLD\_hetgam_dropout.bed;
fi;
