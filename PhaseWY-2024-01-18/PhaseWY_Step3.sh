#!/bin/bash -l

# Version 2023-09-22
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

### Step 3. Prepare scaffold for phasing.

mkdir $OUTDIR/3_scaffold_input/$SCAFFOLD;
OUTDIR3=$OUTDIR/3_scaffold_input/$SCAFFOLD;

#### Are there any variants to phase, and enough heterozygotes to phase the scaffold?
if [ $NONDIPLOIDS == Yes ]; then
  DIPLOIDS=$(cat $INDS | tail -n+2 | awk -F'\t' '{if($4==1) print "--remove-indv",$1}' | tr '\n' ' ');
  HETEROZYGOTES_INDS=$(bcftools view $VCF_IN.vcf.gz -r $SCAFFOLD -Ov  | vcftools --vcf - --hardy $DIPLOIDS \
  --stdout | tail -n+2 | cut -d'/' -f2 |  awk '$1>1 {print $1; exit}' | wc -l);
else
  HETEROZYGOTES_INDS=$(bcftools view $VCF_IN.vcf.gz -r $SCAFFOLD -Ov  | vcftools --vcf - --hardy --stdout | \
  tail -n+2 | cut -d'/' -f2 | awk '$1>1 {print $1; exit}' | wc -l);
fi;

if [ $HETEROZYGOTES_INDS -gt 0 ]; then
  HETEROZYGOTES=TRUE;
else
  HETEROZYGOTES=FALSE;
  echo No biallelic sites, or only one heterozygote per site at scaffold $SCAFFOLD. Cannot phase with shapeit4. $SCAFFOLD will not be phased.
fi;
export HETEROZYGOTES;

#### Extract scaffold and split multiallelic sites
bcftools annotate $VCF_IN.vcf.gz -r $SCAFFOLD -x FORMAT -Ou | \
bcftools norm -m - -Ov | \
bgzip -c > $OUTDIR3/$SCAFF_NAME\_biallelic.vcf.gz;
tabix $OUTDIR3/$SCAFF_NAME\_biallelic.vcf.gz;

if [ $HETEROZYGOTES == TRUE ]; then

  if [ $NONDIPLOIDS == Yes ]; then

    #### Seperate haploids out and convert diploid missing haplotypes to haploid missing genotypes
    HAPLOIDS=$(echo "-s "$(cat $INDS | tail -n+2 | awk -F'\t' '{if($4==1) print $1}' | tr '\n' ','));
    HAPLOIDS=${HAPLOIDS::-1};
    bcftools view $OUTDIR3/$SCAFF_NAME\_biallelic.vcf.gz $HAPLOIDS | \
    awk -F'\t' 'OFS="\t" {if($1 ~ "#") print; else for(i=10; i<=NF;i++) {if($i~/\.|\./ || $i~/\.\/\./) $i="."}; print}' | \
    bgzip -c > $OUTDIR3/$SCAFF_NAME\_biallelic_haploids.vcf.gz;
    tabix $OUTDIR3/$SCAFF_NAME\_biallelic_haploids.vcf.gz;

    #### Convert haploid genotypes to homozygote diploids
    python3 $FUNCTIONS/haploid2diploid.py -i $OUTDIR3/$SCAFF_NAME\_biallelic.vcf.gz \
    -o $OUTDIR3/$SCAFF_NAME\_biallelic.temp.vcf;
    bgzip $OUTDIR3/$SCAFF_NAME\_biallelic.temp.vcf;
    rm $OUTDIR3/$SCAFF_NAME\_biallelic.vcf.gz*;
    mv $OUTDIR3/$SCAFF_NAME\_biallelic.temp.vcf.gz $OUTDIR3/$SCAFF_NAME\_biallelic.vcf.gz;
    tabix $OUTDIR3/$SCAFF_NAME\_biallelic.vcf.gz;
  fi;

fi;
