#!/bin/bash -l

# Version 2023-09-13
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand


### Step 6. Determine sex-linkage of each haplotype
mkdir $OUTDIR/6_sex_linkage/$SCAFFOLD;
OUTDIR6=$OUTDIR/6_sex_linkage/$SCAFFOLD;

bcftools view $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz | \
vcftools --vcf - \
--max-alleles 2 \
--recode --recode-INFO-all --stdout | \
vcffilter -f "AC > 1 & AC < ( $GENOMES - 1 )" | \
bgzip -c > $OUTDIR6/$SCAFF_NAME\_phased_all_variants_mac2.vcf.gz;
tabix $OUTDIR6/$SCAFF_NAME\_phased_all_variants_mac2.vcf.gz;

Rscript $FUNCTIONS/sex_linkage.r --no-save --args PROJECT=$PROJECT INDS=$INDS \
SCAFFOLD=$SCAFFOLD SCAFFOLD_LENGTH=$SCAFFOLD_LENGTH OUTDIR5=$OUTDIR5 \
OUTDIR6=$OUTDIR6 THRESHOLD=$THRESHOLD WINDOW=$WINDOW STEP=$STEP \
SEX_SYSTEM=$SEX_SYSTEM;

rm $OUTDIR6/$SCAFF_NAME\_phased_all_variants_mac2.vcf.gz*;

if [ $(cat $OUTDIR6/$SCAFF_NAME\_phase_info.bed | grep "Sex-linked" | wc -l) -gt 0 ]; then

  #### Identify problematic sex-linked sites. If heterogametes show homozygotic genotypes for both alleles it is likely a sign of heterogametic drop out or incomplete lineage sorting.
  cat $OUTDIR6/$SCAFF_NAME\_phase_info.bed | grep "#\|Sex-linked" > $OUTDIR6/$SCAFF_NAME\_sex_linked.bed;
  HETGAM=$(echo "-s "$(cat $INDS | awk -F'\t' '$3=="HETGAM" {print $1}' | tr '\n' ','));
  HETGAM=${HETGAM::-1};
  bcftools view $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz $HETGAM -R $OUTDIR6/$SCAFF_NAME\_sex_linked.bed | \
  vcffilter -f "AC > 0 & AF < 1" | \
  vcftools --vcf - --hardy --stdout | \
  tail -n+2 | cut -f1,2,3 | awk -F'\t|/' '{if($3 > 0 && $5 > 0) print $1"\t"$2-1"\t"$2}' \
  > $OUTDIR6/$SCAFF_NAME\_sex_linked_ILS1.bed;
fi;
