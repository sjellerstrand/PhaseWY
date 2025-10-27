#!/bin/bash -l

# Version 2023-10-16
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

# Step 8. Summarize the whole genome with bed files
mkdir $OUTDIR/final_output/beds;
OUTDIRb=$OUTDIR/final_output/beds;

if [ $SEX_SYSTEM != NA ] && [ $PHASE == "Yes" ]; then

  ## Concatenate regions of the genome with hetrogametic drop-out due to sex depth differences
  cat $(find $OUTDIR/7_concat_vcfs/ -name "*_hetgam_dropout.bed") | \
  bedtools sort -g $REF.fai \
  > $OUTDIRb/$PROJECT\_hetgam_dropout.bed;
  rm -r $OUTDIR/2_sex_depth_difference/*;


  ## Concatenate sex-linked regions with unknown phase
  echo "###Header" > $OUTDIRb/$PROJECT\_unknown_phase.bed;
  cat $(find $OUTDIR/6_sex_linkage/ -name "*_phase_info.bed") | awk -F'\t' '{if($4=="Unknown") print $1"\t"$2"\t"$3}' | \
  bedtools sort -g $REF.fai | \
  bedtools subtract -a - -b $OUTDIRb/$PROJECT\_hetgam_dropout.bed \
  >> $OUTDIRb/$PROJECT\_unknown_phase.bed;

  ## Concatenate regions of the genome under analysis (alignment within target depth, target missingness, no repeat regions, no unknown phase)
  cat $(find $OUTDIR/1_target_region/ -name "*_target_region.bed") | bedtools sort -g $REF.fai | \
  bedtools subtract -a - -b  $OUTDIRb/$PROJECT\_unknown_phase.bed \
  > $OUTDIRb/$PROJECT\_target_region.bed;

else

  ## Concatenate regions of the genome under analysis (alignment within target depth, target missingness, no repeat regions)
  cat $(find $OUTDIR/1_target_region/ -name "*_target_region.bed") | bedtools sort -g $REF.fai \
  > $OUTDIRb/$PROJECT\_target_region.bed;
fi;
rm -r $OUTDIR/1_target_region/*;

## Concatenate regions with missing data
bedtools complement -i $OUTDIRb/$PROJECT\_target_region.bed -g $REF.fai \
> $OUTDIRb/$PROJECT\_missing_region.bed;

if [ $GENES != NA ]; then
  ## Features overlapping with target region
  bedtools intersect -a $GENES -b $OUTDIRb/$PROJECT\_target_region.bed -f 0.95 \
  > $OUTDIRb/$PROJECT\_target_features1.bed;
fi;


if [ $EXON != NA ]; then
  ## Exonic with 20 000 bp flanking region
  bedtools slop -i $EXON -g $REF.fai -b 20000 | bedtools merge \
  > $OUTDIRb/$PROJECT\_exon_with_flank.bed;
fi;

## Concatenate target regions of the genome that have been phased
if [ PHASE_SCAFFOLD_LIST != NA ] && [ LIST_ONLY == No ]; then
  for SCAFFOLD0 in $(cat $PHASE_SCAFFOLD_LIST | cut -f1 | sort -u); do
    cat $OUTDIRb/$PROJECT\_target_region.bed | grep $SCAFFOLD0;
  done > $OUTDIRb/$PROJECT\_target_region_phased.bed;

  for SCAFFOLD0 in $(find $OUTDIR/3_scaffold_input/ -name "*_biallelic.vcf.gz" | rev | cut -d'/' -f2 | rev); do
    sed -i "/$SCAFFOLD0/d" $OUTDIRb/$PROJECT\_target_region_phased.bed;
  done;

elif [ $PHASE == "Yes" ]; then
  cp $OUTDIRb/$PROJECT\_target_region.bed $OUTDIRb/$PROJECT\_target_region_phased.bed;
  for SCAFFOLD0 in $(find $OUTDIR/3_scaffold_input/ -name "*_biallelic.vcf.gz" | rev | cut -d'/' -f2 | rev); do
    sed -i "/$SCAFFOLD0/d" $OUTDIRb/$PROJECT\_target_region_phased.bed;
  done;
fi;

if [ $SEX_SYSTEM != NA  ] && [ $PHASE == "Yes" ]; then

  ### Extract sex-linked regions due to sex specfic clustering and phase info
  cat $(find $OUTDIR/6_sex_linkage/ -name "*_phase_windows.txt") | grep -v "chromStart" \
  > $OUTDIRb/$PROJECT\_phase_windows.txt;
  cat $(find $OUTDIR/6_sex_linkage/ -name "*_phase_info.bed") | grep -v "chromStart" \
  > $OUTDIRb/$PROJECT\_phase_info.bed;
  cat $OUTDIRb/$PROJECT\_phase_info.bed | grep "Sex-linked" | \
  cut -f 1,2,3 | bedtools sort -g $REF.fai | bedtools merge \
  > $OUTDIRb/$PROJECT\_phase_sex_linked.bed;
  bedtools intersect -a $OUTDIRb/$PROJECT\_phase_sex_linked.bed \
  -b $OUTDIRb/$PROJECT\_target_region_phased.bed \
  > $OUTDIRb/$PROJECT\_target_phase_sex_linked.bed;

  ### Make a bed file of sex linked regions within target regions
  bedtools intersect -a $OUTDIRb/$PROJECT\_hetgam_dropout.bed \
  -b $OUTDIRb/$PROJECT\_target_region_phased.bed \
  > $OUTDIRb/$PROJECT\_target_phase_hetgam_dropout.bed;
  cat $OUTDIRb/$PROJECT\_target_phase_hetgam_dropout.bed \
  $OUTDIRb/$PROJECT\_target_phase_sex_linked.bed | \
  bedtools sort -g $REF.fai | bedtools merge \
  > $OUTDIRb/$PROJECT\_sex_linked.bed;

  ### Make a bed file of homogametic (all sex-linked) regions
  cp $OUTDIRb/$PROJECT\_sex_linked.bed \
  $OUTDIRb/$PROJECT\_homogametic.bed;

  if [ $EXON != NA ]; then
    ### Make a bed file of homogametic (all sex-linked) neutral regions within target regions
    bedtools subtract -a $OUTDIRb/$PROJECT\_homogametic.bed \
    -b $OUTDIRb/$PROJECT\_exon_with_flank.bed \
    > $OUTDIRb/$PROJECT\_homogametic_neutral.bed;
  fi;

  if [ $CDS != NA ]; then
    ### Make a bed file of homogametic (all sex-linked) CDS regions within target regions
    bedtools intersect -a $OUTDIRb/$PROJECT\_homogametic.bed -b $CDS \
    > $OUTDIRb/$PROJECT\_homogametic_CDS.bed;
  fi;

  ### Make a bed file of target heterogametic regions (where there is no heterogametic drop-out, i.e. missing data)
  bedtools subtract -a $OUTDIRb/$PROJECT\_sex_linked.bed \
  -b  $OUTDIRb/$PROJECT\_hetgam_dropout.bed \
  > $OUTDIRb/$PROJECT\_heterogametic.bed;

  if [ $EXON != NA ]; then
    ### Make a bed file of sex linked heterogametic neutral regions
    bedtools subtract -a $OUTDIRb/$PROJECT\_heterogametic.bed \
    -b $OUTDIRb/$PROJECT\_exon_with_flank.bed \
    > $OUTDIRb/$PROJECT\_heterogametic_neutral.bed;
  fi;
  
  if [ $CDS != NA ]; then
    ### Make a bed file of sex linked heterogametic CDS regions
    bedtools intersect -a $OUTDIRb/$PROJECT\_heterogametic.bed -b $CDS \
    > $OUTDIRb/$PROJECT\_heterogametic_CDS.bed;
  fi;

  ### Make a bed file of target autosomal regions
  bedtools subtract -a $OUTDIRb/$PROJECT\_target_region.bed \
  -b $OUTDIRb/$PROJECT\_sex_linked.bed \
  > $OUTDIRb/$PROJECT\_autosomal.bed;

  if [ $EXON != NA ]; then
    ### Make a bed file of autosomal neutral regions within target regions
    bedtools subtract -a $OUTDIRb/$PROJECT\_autosomal.bed \
    -b $OUTDIRb/$PROJECT\_exon_with_flank.bed \
    > $OUTDIRb/$PROJECT\_autosomal_neutral.bed;
  fi;
  
  if [ $CDS != NA ]; then
    ### Make a bed file of autosomal CDS regions within target regions
    bedtools intersect -a $OUTDIRb/$PROJECT\_autosomal.bed -b $CDS \
    > $OUTDIRb/$PROJECT\_autosomal_CDS.bed;
  fi;
  
  if [ $GENES != NA ]; then

    ### Target features overlapping with homogametic (sex-linked) regions
    bedtools intersect -a $OUTDIRb/$PROJECT\_target_features1.bed \
    -b $OUTDIRb/$PROJECT\_homogametic.bed -f 0.95 -wa | \
    awk -F'\t' '!seen[$0]++' \
    > $OUTDIRb/$PROJECT\_homogametic_features.bed;

    ### Target features overlapping with heterogametic regions
    bedtools intersect -a $OUTDIRb/$PROJECT\_target_features1.bed \
    -b $OUTDIRb/$PROJECT\_heterogametic.bed -f 0.95 -wa | \
    awk -F'\t' '!seen[$0]++' \
    > $OUTDIRb/$PROJECT\_heterogametic_features.bed;

    ### Target features overlapping with autosomal regions
    bedtools intersect -a $OUTDIRb/$PROJECT\_target_features1.bed \
    -b $OUTDIRb/$PROJECT\_autosomal.bed -f 0.95 -wa | \
    awk -F'\t' '!seen[$0]++' \
    > $OUTDIRb/$PROJECT\_autosomal_features.bed;
    
fi;

  ### Problematic sex-linked sites due to signs of incomplete linage sorting (heterogametic dropout=ILS1 and unsuccessful phasing=ILS2)
  echo "###Header" > $OUTDIRb/$PROJECT\_sex_linked_ILS1.bed;
  cat $(find $OUTDIR/6_sex_linkage/ -name "*_sex_linked_ILS1.bed") \
  >> $OUTDIRb/$PROJECT\_sex_linked_ILS1.bed;
  echo "###Header" > $OUTDIRb/$PROJECT\_sex_linked_ILS1_not_hetgam_dropout.bed;
  bedtools subtract -a $OUTDIRb/$PROJECT\_sex_linked_ILS1.bed \
  -b $OUTDIRb/$PROJECT\_target_phase_hetgam_dropout.bed \
  >> $OUTDIRb/$PROJECT\_sex_linked_ILS1_not_hetgam_dropout.bed;
  rm -r $OUTDIR/6_sex_linkage/*;
  echo "###Header" > $OUTDIRb/$PROJECT\_sex_linked_ILS2.bed;
  cat $(find $OUTDIR/7_concat_vcfs/ -name "*_sex_linked_ILS2.bed") \
  >> $OUTDIRb/$PROJECT\_sex_linked_ILS2.bed;

fi;

if [ $EXON != NA ]; then
  ### Make a bed file of target neutral regions
  bedtools subtract -a $OUTDIRb/$PROJECT\_target_region.bed \
  -b $OUTDIRb/$PROJECT\_exon_with_flank.bed \
  > $OUTDIRb/$PROJECT\_neutral.bed;
fi;

if [ $CDS != NA ]; then
  ### Make a bed file of target CDS regions
  bedtools intersect -a $OUTDIRb/$PROJECT\_target_region.bed -b $CDS \
  > $OUTDIRb/$PROJECT\_CDS.bed;
fi;

if [ $GENES != NA ]; then
  ### Target features overlapping with target regions
  bedtools intersect -a $OUTDIRb/$PROJECT\_target_features1.bed \
  -b $OUTDIRb/$PROJECT\_target_region.bed -f 0.95 -wa | \
  awk -F'\t' '!seen[$0]++' \
  > $OUTDIRb/$PROJECT\_target_features.bed;
  rm $OUTDIRb/$PROJECT\_target_features1.bed;
fi;

