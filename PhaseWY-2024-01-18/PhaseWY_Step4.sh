#!/bin/bash -l

# Version 2023-09-13
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

### Step 4. Perform read based phasing with whatshap
mkdir $OUTDIR/4_whatshap/$SCAFFOLD;
OUTDIR4=$OUTDIR/4_whatshap/$SCAFFOLD;
mkdir $OUTDIR4/temp1;
mkdir $OUTDIR4/temp2;

for IND in $(cat $INDS | tail -n+2 | cut -f1); do 
  bcftools view $OUTDIR3/$SCAFF_NAME\_biallelic.vcf.gz -s $IND | \
  bgzip -c > $OUTDIR4/temp1/$SCAFF_NAME\_biallelic_$IND.vcf.gz;
  tabix $OUTDIR4/temp1/$SCAFF_NAME\_biallelic_$IND.vcf.gz;

  whatshap phase -o $OUTDIR4/temp2/$SCAFF_NAME\_$IND\_whatshap.vcf.gz \
  --indels -r $REF $OUTDIR4/temp1/$SCAFF_NAME\_biallelic_$IND.vcf.gz \
  $(cat $BAMS | grep $IND) --sample $IND;
  tabix $OUTDIR4/temp2/$SCAFF_NAME\_$IND\_whatshap.vcf.gz;

  rm $OUTDIR4/temp1/$SCAFF_NAME\_biallelic_$IND.vcf.gz*;
done;

##### Activate bcftools 1.17 conda environment
source $conda/etc/profile.d/conda.sh;
conda activate PhaseWY_bcftools1.17;

##### Merge homogametic haplotypes
merge_files=$(ls $OUTDIR4/temp2/$SCAFF_NAME\_*_whatshap.vcf.gz);
bcftools merge $merge_files | \
bgzip -c > $OUTDIR4/$SCAFF_NAME\_whatshap.vcf.gz;
tabix $OUTDIR4/$SCAFF_NAME\_whatshap.vcf.gz;

###### Dectivate bcftools 1.17 conda environment and return to PhaseWY bioinfo environment again
conda deactivate;

rm -r $OUTDIR4/temp*;
