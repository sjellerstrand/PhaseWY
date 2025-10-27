#!/bin/bash -l

# Version 2023-09-13
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

### Step 5. Preform statistical phasing with shapeit4
mkdir $OUTDIR/5_shapeit4/$SCAFFOLD;
OUTDIR5=$OUTDIR/5_shapeit4/$SCAFFOLD;
INDS_ORDER=$(cat $INDS | tail -n+2 | cut -f1 | tr '\n' ',');
INDS_ORDER=${INDS_ORDER::-1};

if [ $NONDIPLOIDS == Yes ]; then

  ##### If there are haploids, use them as a haplotype refernce panel
  
  ###### Seperate haploids out for reference panel
  HAPLOIDS=$(echo "-s "$(cat $INDS | tail -n+2 | awk -F'\t' '{if($4==1) print $1}' | tr '\n' ','));
  HAPLOIDS=${HAPLOIDS::-1};
  bcftools view  $OUTDIR4/$SCAFF_NAME\_whatshap.vcf.gz $HAPLOIDS | \
  bgzip > $OUTDIR5/$SCAFF_NAME\_reference_panel.vcf.gz;
  tabix $OUTDIR5/$SCAFF_NAME\_reference_panel.vcf.gz;
  
  ###### Seperate diploids out for phasing
  DIPLOIDS=$(echo "-s "$(cat $INDS | tail -n+2 | awk -F'\t' '{if($4==2) print $1}' | tr '\n' ','));
  DIPLOIDS=${DIPLOIDS::-1};
  bcftools view  $OUTDIR4/$SCAFF_NAME\_whatshap.vcf.gz $DIPLOIDS | \
  bgzip > $OUTDIR5/$SCAFF_NAME\_shapeit4_temp1.vcf.gz;
  tabix $OUTDIR5/$SCAFF_NAME\_shapeit4_temp1.vcf.gz;
  
  ###### Activate shapeit4 conda environment
  conda activate PhaseWY_shapeit4;
  
  ###### Phase samples
  shapeit4 --input $OUTDIR5/$SCAFF_NAME\_shapeit4_temp1.vcf.gz \
  --thread 1 --region $SCAFFOLD --use-PS 0.0001 \
  --sequencing --output $OUTDIR5/$SCAFF_NAME\_shapeit4_temp2.vcf.gz \
  --reference $OUTDIR5/$SCAFF_NAME\_reference_panel.vcf.gz;

  ###### Dectivate shapeit4 conda environment and return to PhaseWY bioinfo environment again
  conda deactivate;
  rm $OUTDIR5/$SCAFF_NAME\_reference_panel.vcf.gz* \
  $OUTDIR5/$SCAFF_NAME\_shapeit4_temp1.vcf.gz*;
  tabix $OUTDIR5/$SCAFF_NAME\_shapeit4_temp2.vcf.gz;

  ##### Activate bcftools 1.17 conda environment
  conda activate PhaseWY_bcftools1.17;

  ##### Merge back haploids with haploid genotypes
  bcftools merge $OUTDIR5/$SCAFF_NAME\_shapeit4_temp2.vcf.gz \
  $OUTDIR3/$SCAFF_NAME\_biallelic_haploids.vcf.gz -Ou | \
  bcftools view -s $INDS_ORDER | \
  vcffixup - | vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR5/$SCAFF_NAME\_shapeit4_temp3.vcf.gz;

  ###### Dectivate bcftools 1.17 conda environment and return to PhaseWY bioinfo environment again
  conda deactivate;
  rm $OUTDIR5/$SCAFF_NAME\_shapeit4_temp2.vcf.gz*;
  tabix $OUTDIR5/$SCAFF_NAME\_shapeit4_temp3.vcf.gz;

  ##### Merge multiallelic sites and remove haploids with diploid genotypes
  bcftools norm $OUTDIR5/$SCAFF_NAME\_shapeit4_temp3.vcf.gz -m +any | \
  awk -F $'\t' 'BEGIN {OFS = FS} /^[#]/ {print; next} {for (i=10; i<=NF; i++) { gsub("/","|",$i)} print}' | \
  bgzip -c > $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz;
  tabix $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz;
  rm $OUTDIR5/$SCAFF_NAME\_shapeit4_temp3.vcf.gz*;

else

  ##### Activate shapeit4 conda environment
  conda activate PhaseWY_shapeit4;

  shapeit4 --input $OUTDIR4/$SCAFF_NAME\_whatshap.vcf.gz \
  --thread 1 --region $SCAFFOLD --use-PS 0.0001 \
  --sequencing --output $OUTDIR5/$SCAFF_NAME\_shapeit4.vcf.gz;

  ###### Dectivate shapeit4 conda environment and return to PhaseWY bioinfo environment again
  conda deactivate;
  tabix $OUTDIR5/$SCAFF_NAME\_shapeit4.vcf.gz;

  #### Modify phased vcf
  bcftools view $OUTDIR5/$SCAFF_NAME\_shapeit4.vcf.gz -s $INDS_ORDER -Ou | \
  bcftools norm -m +any | vcffixup - | vcffilter -f "AC > 0" | \
  awk -F $'\t' 'BEGIN {OFS = FS} /^[#]/ {print; next} {for (i=10; i<=NF; i++) { gsub("/","|",$i)} print}' | \
  bgzip -c > $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz;
  tabix $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz;

fi;

#### Remove input and whatshap files
rm -r $OUTDIR3 $OUTDIR4 $OUTDIR5/$SCAFF_NAME\_shapeit4*;
