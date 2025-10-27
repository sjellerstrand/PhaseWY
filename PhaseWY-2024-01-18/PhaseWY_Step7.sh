#!/bin/bash -l

# Version 2023-07-10
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

### Step 7. Split phased scaffold into autosomal, homogametic, and heterogametic genomic regions
source $conda/etc/profile.d/conda.sh;

#### Set up directory and files
mkdir $OUTDIR/7_concat_vcfs/$SCAFFOLD;
OUTDIR7=$OUTDIR/7_concat_vcfs/$SCAFFOLD;

#### If scaffold contains sex-linked regions, split pseudoautosmal and sex-linked regions, and heterogametes and homogametes

if [ $(cat $OUTDIR6/$SCAFF_NAME\_phase_info.bed | grep "Sex-linked" | wc -l) -gt 0 ]; then

  mkdir $OUTDIR7/$SCAFFOLD\_temp_sex_link;
  OUTDIR7t=$OUTDIR7/$SCAFFOLD\_temp_sex_link;
  INDS_ORDER_HETS=$(cat $INDS | awk -F'\t' '$3=="HETGAM" {print $1}' | tr '\n' ',');
  INDS_ORDER_HETS=${INDS_ORDER_HETS::-1};
  cat $OUTDIR6/$SCAFF_NAME\_phase_info.bed | grep "#\|Sex-linked" | cut -f1,2,3 \
  > $OUTDIR7t/sex_link.bed;
  HOMGAM=$(echo "-s "$(cat $INDS | awk -F'\t' '$3=="HOMGAM" {print $1}' | tr '\n' ','));
  HOMGAM=${HOMGAM::-1};
  HETGAM=$(echo "-s "$(cat $INDS | awk -F'\t' '$3=="HETGAM" {print $1}' | tr '\n' ','));
  HETGAM=${HETGAM::-1};

  ##### Split pseudoautosmal and sex-linked regions, and heterogametes and homogametes
  if [ $(cat $OUTDIR6/$SCAFF_NAME\_phase_info.bed | grep "Pseudoautosomal" | wc -l) -gt 0 ]; then
    cat $OUTDIR6/$SCAFF_NAME\_phase_info.bed | grep "#\|Pseudoautosomal" | cut -f1,2,3 \
    > $OUTDIR7t/pseudoautosomal.bed;

    if [ $(ls $OUTDIR2/$SCAFF_NAME\_hetgam_dropout.bed | wc -l) -gt 0 ]; then
      bedtools subtract -a $OUTDIR7t/pseudoautosomal.bed \
      -b $OUTDIR2/$SCAFF_NAME\_hetgam_dropout.bed \
      > $OUTDIR7t/autosomal.bed;
      cp $OUTDIR2/$SCAFF_NAME\_hetgam_dropout.bed $OUTDIR7/$SCAFF_NAME\_hetgam_dropout.bed;
      cat $OUTDIR7/$SCAFF_NAME\_hetgam_dropout.bed \
      $OUTDIR7t/sex_link.bed | \
      bedtools sort -g $REF.fai | bedtools merge \
      > $OUTDIR7t/sex_link_hetgam_dropout.bed;
      rm $OUTDIR7t/pseudoautosomal.bed;
    else
      mv $OUTDIR7t/pseudoautosomal.bed $OUTDIR7t/autosomal.bed;
    fi;

    if [ $(cat $OUTDIR7t/autosomal.bed | wc -l) -gt 0 ]; then
      bcftools view $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz -R $OUTDIR7t/autosomal.bed | \
      bgzip -c > $OUTDIR7/$SCAFF_NAME\_autosomal.vcf.gz;
      tabix $OUTDIR7/$SCAFF_NAME\_autosomal.vcf.gz;
    fi;

  elif [ $(ls $OUTDIR2/$SCAFF_NAME\_hetgam_dropout.bed | wc -l) -gt 0 ]; then

    cp $OUTDIR2/$SCAFF_NAME\_hetgam_dropout.bed $OUTDIR7/$SCAFF_NAME\_hetgam_dropout.bed;
    cat $OUTDIR7/$SCAFF_NAME\_hetgam_dropout.bed \
    $OUTDIR7t/sex_link.bed | \
    bedtools sort -g $REF.fai | bedtools merge \
    > $OUTDIR7t/sex_link_hetgam_dropout.bed;

  fi;

  if [ $(ls $OUTDIR7/$SCAFF_NAME\_hetgam_dropout.bed | wc -l) -gt 0 ]; then
    bcftools view $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz -R $OUTDIR7t/sex_link_hetgam_dropout.bed $HOMGAM | \
    bgzip -c > $OUTDIR7t/$SCAFF_NAME\_HOMGAM.vcf.gz;
    rm $OUTDIR7t/sex_link_hetgam_dropout.bed;
  else
    bcftools view $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz -R $OUTDIR7t/sex_link.bed $HOMGAM | \
    bgzip -c > $OUTDIR7t/$SCAFF_NAME\_HOMGAM.vcf.gz;
  fi;
  tabix $OUTDIR7t/$SCAFF_NAME\_HOMGAM.vcf.gz;
  bcftools view $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz -R $OUTDIR7t/sex_link.bed $HETGAM | \
  bgzip -c > $OUTDIR7t/$SCAFF_NAME\_HETGAM.vcf.gz;

  ##### Split heterogamets into heterogametic and homogametic genotypes
  python3 $FUNCTIONS/split_phase.py -i $OUTDIR7t/$SCAFF_NAME\_HETGAM.vcf.gz -p left \
  -o $OUTDIR7t/$SCAFF_NAME\_HETGAM_left.vcf;
  bgzip $OUTDIR7t/$SCAFF_NAME\_HETGAM_left.vcf;
  tabix $OUTDIR7t/$SCAFF_NAME\_HETGAM_left.vcf.gz;
  python3 $FUNCTIONS/split_phase.py -i $OUTDIR7t/$SCAFF_NAME\_HETGAM.vcf.gz -p right \
  -o $OUTDIR7t/$SCAFF_NAME\_HETGAM_right.vcf;
  bgzip $OUTDIR7t/$SCAFF_NAME\_HETGAM_right.vcf;
  tabix $OUTDIR7t/$SCAFF_NAME\_HETGAM_right.vcf.gz;
  rm $OUTDIR7t/$SCAFF_NAME\_HETGAM.vcf.gz;

  ##### Concatenate individual heterogametes into heterogametic and homogametic phase
  for IND in $(cat $INDS | awk -F'\t' '$3=="HETGAM" {print $1}'); do

    if [ $(cat $OUTDIR6/$SCAFF_NAME\_$IND\_het_left.bed | grep -v "#" | awk '{print; exit}' | wc -l) -gt 0 ]; then

      ###### Extract heterogametic sites with left phase
      bcftools view $OUTDIR7t/$SCAFF_NAME\_HETGAM_left.vcf.gz -s $IND \
      -R $OUTDIR6/$SCAFF_NAME\_$IND\_het_left.bed |
       bgzip -c > $OUTDIR7t/$SCAFF_NAME\_HETGAM_$IND\_hetsex1.vcf.gz;

      ###### Extract homogametic sites with right phase
      bcftools view $OUTDIR7t/$SCAFF_NAME\_HETGAM_right.vcf.gz -s $IND\
      -R $OUTDIR6/$SCAFF_NAME\_$IND\_het_left.bed |
      bgzip -c > $OUTDIR7t/$SCAFF_NAME\_HETGAM_$IND\_homsex1.vcf.gz;
    fi;

    if [ $(cat $OUTDIR6/$SCAFF_NAME\_$IND\_het_right.bed | grep -v "#" | awk '{print; exit}' | wc -l) -gt 0 ]; then

      ###### Extract homogametic sites with left phase
      bcftools view $OUTDIR7t/$SCAFF_NAME\_HETGAM_left.vcf.gz -s $IND \
      -R $OUTDIR6/$SCAFF_NAME\_$IND\_het_right.bed |
      bgzip -c > $OUTDIR7t/$SCAFF_NAME\_HETGAM_$IND\_homsex2.vcf.gz;

      ###### Extract heterogametic sites with right phase
      bcftools view $OUTDIR7t/$SCAFF_NAME\_HETGAM_right.vcf.gz -s $IND \
      -R $OUTDIR6/$SCAFF_NAME\_$IND\_het_right.bed  |
      bgzip -c > $OUTDIR7t/$SCAFF_NAME\_HETGAM_$IND\_hetsex2.vcf.gz;
     fi;

    ###### Concatenate heterogametic genotypes
    concat_files=$(ls $OUTDIR7t/$SCAFF_NAME\_HETGAM_$IND\_hetsex*.vcf.gz);
    bcftools concat $concat_files -Ou | bcftools sort -T $OUTDIR7/$SCAFF_NAME\_temp_merge | \
    bgzip -c > $OUTDIR7t/$SCAFF_NAME\_HETGAM_$IND\_hetsex.vcf.gz;
    tabix $OUTDIR7t/$SCAFF_NAME\_HETGAM_$IND\_hetsex.vcf.gz;
    rm $concat_files;

    ###### Concatenate homogametic genotypes
    concat_files=$(ls $OUTDIR7t/$SCAFF_NAME\_HETGAM_$IND\_homsex*.vcf.gz);
    bcftools concat $concat_files -Ou | bcftools sort -T $OUTDIR7/$SCAFF_NAME\_temp_merge | \
    bgzip -c > $OUTDIR7t/$SCAFF_NAME\_HETGAM_$IND\_homsex.vcf.gz;
    tabix $OUTDIR7t/$SCAFF_NAME\_HETGAM_$IND\_homsex.vcf.gz;
    rm $concat_files;
  done;

  ##### Remove temporary files
  rm $OUTDIR7t/$SCAFF_NAME\_HETGAM_left.vcf.gz* \
  $OUTDIR7t/$SCAFF_NAME\_HETGAM_right.vcf.gz*;

  ###### Activate bcftools 1.17 conda environment
  conda activate PhaseWY_bcftools1.17;

  ##### Merge heterogametic haplotypes
  merge_files=$(ls $OUTDIR7t/$SCAFF_NAME\_HETGAM_*_hetsex.vcf.gz);
  bcftools merge $merge_files -Ou | bcftools view -s $INDS_ORDER_HETS | \
  vcffixup - | vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR7/$SCAFF_NAME\_heterogametic_temp1.vcf.gz;
  rm $(ls $OUTDIR7t/$SCAFF_NAME\_HETGAM_*_hetsex.vcf.gz*);

  ###### Dectivate bcftools 1.17 conda environment and return to PhaseWY bioinfo environment again
  conda deactivate;

  ###### If there were heterogametic drop-out regions, remove these in the heterogametic file
  if [ $(ls $OUTDIR7/$SCAFF_NAME\_hetgam_dropout.bed | wc -l) -gt 0 ]; then

    vcftools --gzvcf $OUTDIR7/$SCAFF_NAME\_heterogametic_temp1.vcf.gz \
    --exclude-bed $OUTDIR7/$SCAFF_NAME\_hetgam_dropout.bed \
    --recode --recode-INFO-all --stdout | \
    bgzip -c > $OUTDIR7/$SCAFF_NAME\_heterogametic.vcf.gz;
    rm $OUTDIR7/$SCAFF_NAME\_heterogametic_temp1.vcf.gz;

  else
    mv $OUTDIR7/$SCAFF_NAME\_heterogametic_temp1.vcf.gz $OUTDIR7/$SCAFF_NAME\_heterogametic.vcf.gz;

  fi;

  tabix $OUTDIR7/$SCAFF_NAME\_heterogametic.vcf.gz;

  ###### Activate bcftools 1.17 conda environment
  conda activate PhaseWY_bcftools1.17;

  ##### Merge homogametic haplotypes
  merge_files=$(ls $OUTDIR7t/$SCAFF_NAME\_HETGAM_*_homsex.vcf.gz);
  bcftools merge $merge_files -Ou | bcftools view -s $INDS_ORDER_HETS | \
  bgzip -c > $OUTDIR7/$SCAFF_NAME\_homogametic_temp1.vcf.gz;
  tabix $OUTDIR7/$SCAFF_NAME\_homogametic_temp1.vcf.gz;
  rm $(ls $OUTDIR7t/$SCAFF_NAME\_HETGAM_*_homsex.vcf.gz*);

  ###### Dectivate bcftools 1.17 conda environment and return to PhaseWY bioinfo environment again
  conda deactivate;

  ###### If there were heterogametic drop-out variants, add these to the homogametic file for the heterogametic sex
  if [ $(ls $OUTDIR7/$SCAFF_NAME\_hetgam_dropout.bed | wc -l) -gt 0 ]; then

    bcftools view $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz \
    -R $OUTDIR7/$SCAFF_NAME\_hetgam_dropout.bed $HETGAM | \
    vcftools --vcf - \
    --exclude-bed $OUTDIR7t/sex_link.bed \
    --recode --recode-INFO-all --stdout | \
    bgzip -c > $OUTDIR7/$SCAFF_NAME\_homogametic_temp2.vcf.gz

    if [ $(bcftools view -H $OUTDIR7/$SCAFF_NAME\_homogametic_temp2.vcf.gz | awk '{print; exit}' | wc -l) -gt 0 ]; then

      ###### Re-code homozygose genotypes to haploid for the heterogametic sex
      python3 $FUNCTIONS/diploid2haploid.py -i $OUTDIR7/$SCAFF_NAME\_homogametic_temp2.vcf.gz \
      -o $OUTDIR7t/$SCAFF_NAME\_HETGAM_homogametic.vcf;
      bgzip $OUTDIR7t/$SCAFF_NAME\_HETGAM_homogametic.vcf;
      tabix $OUTDIR7t/$SCAFF_NAME\_HETGAM_homogametic.vcf.gz;

      ##### Concatenate homogametic sites for the heterogametic sex
      bcftools concat $OUTDIR7/$SCAFF_NAME\_homogametic_temp1.vcf.gz \
      $OUTDIR7t/$SCAFF_NAME\_HETGAM_homogametic.vcf.gz -Ou | \
      bcftools sort -T $OUTDIR7/$SCAFF_NAME\_temp_merge | \
      bgzip -c > $OUTDIR7/$SCAFF_NAME\_homogametic_temp3.vcf.gz;
      rm $OUTDIR7/$SCAFF_NAME\_homogametic_temp1.vcf.gz* \
      $OUTDIR7/$SCAFF_NAME\_homogametic_temp2.vcf.gz \
      $OUTDIR7t/$SCAFF_NAME\_HETGAM_homogametic.vcf.gz;

    else
      mv $OUTDIR7/$SCAFF_NAME\_homogametic_temp2.vcf.gz \
      $OUTDIR7/$SCAFF_NAME\_homogametic_temp3.vcf.gz;
      rm $OUTDIR7/$SCAFF_NAME\_homogametic_temp1.vcf.gz*;

    fi;

  else
    mv  $OUTDIR7/$SCAFF_NAME\_homogametic_temp1.vcf.gz \
    $OUTDIR7/$SCAFF_NAME\_homogametic_temp3.vcf.gz;

  fi;
  tabix $OUTDIR7/$SCAFF_NAME\_homogametic_temp3.vcf.gz;

  ###### Activate bcftools 1.17 conda environment
  conda activate PhaseWY_bcftools1.17;

  ##### Merge homogametic haplotypes from the homogametic and heterogametic sex
  bcftools merge $OUTDIR7t/$SCAFF_NAME\_HOMGAM.vcf.gz \
  $OUTDIR7/$SCAFF_NAME\_homogametic_temp3.vcf.gz | \
  bgzip -c > $OUTDIR7/$SCAFF_NAME\_homogametic_temp4.vcf.gz;
  tabix $OUTDIR7/$SCAFF_NAME\_homogametic_temp4.vcf.gz;

  ###### Dectivate bcftools 1.17 conda environment and return to PhaseWY bioinfo environment again
  conda deactivate;

  ###### Order individuals and filter sites
  bcftools view $OUTDIR7/$SCAFF_NAME\_homogametic_temp4.vcf.gz -s $INDS_ORDER | \
  vcffixup - | vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR7/$SCAFF_NAME\_homogametic.vcf.gz;
  tabix $OUTDIR7/$SCAFF_NAME\_homogametic.vcf.gz;
  rm -r $OUTDIR7/$SCAFF_NAME\_homogametic_temp3.vcf.gz* \
  $OUTDIR7/$SCAFF_NAME\_homogametic_temp4.vcf.gz*;

  #### Indentify problematic sex-linked sites. If both alleles occur on both sex chromosomes it is likely a sign of unsuccessful phasing or incomplete lineage sorting.
  bedtools intersect -a $OUTDIR7/$SCAFF_NAME\_heterogametic.vcf.gz \
  -b $OUTDIR7/$SCAFF_NAME\_homogametic.vcf.gz | awk -F'\t' '{print $1"\t"$2-1"\t"$2}' \
  > $OUTDIR7/$SCAFF_NAME\_overlaps.bed;

  if [ $(cat $OUTDIR7/$SCAFF_NAME\_overlaps.bed | awk 'NR>0 {print; exit}' | wc -l) -gt 0 ]; then

    bcftools query $OUTDIR7/$SCAFF_NAME\_heterogametic.vcf.gz \
    -R $OUTDIR7/$SCAFF_NAME\_overlaps.bed -f '%CHROM\t%POS\t%INFO/AN\t%INFO/AC\n' | \
    awk -F'\t|,' '{sum=0; for(i=4; i<=NF; ++i) sum+=$i} \
    {if(($4/$3==1) || sum==0) print $1"\t"$2"\tmonomorphic\t"; \
    else if(NF>4) print $1"\t"$2"\tmultiallelic"; \
    else print $1"\t"$2"\tpolymorphic"}' | \
    sort -k2,2 > $OUTDIR7/$SCAFF_NAME\_hetgam_var.txt;

    bcftools query $OUTDIR7/$SCAFF_NAME\_homogametic.vcf.gz \
    -R $OUTDIR7/$SCAFF_NAME\_overlaps.bed -f '%CHROM\t%POS\t%INFO/AN\t%INFO/AC\n' | \
    awk -F'\t|,' '{sum=0; for(i=4; i<=NF; ++i) sum+=$i} \
    {if(($4/$3==1) || sum==0) print $1"\t"$2"\tmonomorphic\t"; \
    else if(NF>4) print $1"\t"$2"\tmultiallelic"; \
    else print $1"\t"$2"\tpolymorphic"}' | \
    sort -k2,2 > $OUTDIR7/$SCAFF_NAME\_homgam_var.txt;

    join -j 2 -o 1.1,1.2,1.3,2.3 $OUTDIR7/$SCAFF_NAME\_hetgam_var.txt $OUTDIR7/$SCAFF_NAME\_homgam_var.txt | \
    awk '{if(($3 == "polymorphic" && $3 == $4) || \
    ($3 == "multiallelic" && $4 != "monomorphic") || \
    ($4 == "multiallelic" && $3 != "monomorphic")) \
    print $1"\t"$2-1"\t"$2"\t"$3"\t"$4}' \
    > $OUTDIR7/$SCAFF_NAME\_sex_linked_ILS2.bed;

    rm $OUTDIR7/$SCAFF_NAME\_overlaps.bed \
    $OUTDIR7/$SCAFF_NAME\_hetgam_var.txt \
    $OUTDIR7/$SCAFF_NAME\_homgam_var.txt;

  fi;

  ##### Remove temporary files
  rm -r $OUTDIR6/*_het_left.bed $OUTDIR6/*_het_right.bed $OUTDIR7t $OUTDIR2;

elif [ $(ls $OUTDIR2/$SCAFF_NAME\_hetgam_dropout.bed | wc -l) -gt 0 ]; then

  mkdir $OUTDIR7/$SCAFFOLD\_temp_sex_link;
  OUTDIR7t=$OUTDIR7/$SCAFFOLD\_temp_sex_link;
  HOMGAM=$(echo "-s "$(cat $INDS | awk -F'\t' '$3=="HOMGAM" {print $1}' | tr '\n' ','));
  HOMGAM=${HOMGAM::-1};
  HETGAM=$(echo "-s "$(cat $INDS | awk -F'\t' '$3=="HETGAM" {print $1}' | tr '\n' ','));
  HETGAM=${HETGAM::-1};

  ##### Split pseudoautosmal and sex-linked regions
  if [ $(cat $OUTDIR6/$SCAFF_NAME\_phase_info.bed | grep "Pseudoautosomal\|Autosomal" | wc -l) -gt 0 ]; then
    cat $OUTDIR6/$SCAFF_NAME\_phase_info.bed | grep "#\|Pseudoautosomal\|Autosomal" | cut -f1,2,3 \
   > $OUTDIR7t/pseudoautosomal.bed;

    bedtools subtract -a $OUTDIR7t/pseudoautosomal.bed \
    -b $OUTDIR2/$SCAFF_NAME\_hetgam_dropout.bed \
    > $OUTDIR7t/autosomal.bed;
    cp $OUTDIR2/$SCAFF_NAME\_hetgam_dropout.bed $OUTDIR7/$SCAFF_NAME\_hetgam_dropout.bed;
    rm $OUTDIR7t/pseudoautosomal.bed;

    if [ $(cat $OUTDIR7t/autosomal.bed | wc -l) -gt 0 ]; then
      bcftools view $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz -R $OUTDIR7t/autosomal.bed | \
      bgzip -c > $OUTDIR7/$SCAFF_NAME\_autosomal.vcf.gz;
      tabix $OUTDIR7/$SCAFF_NAME\_autosomal.vcf.gz;
    fi;

    ###### Split sex-linked regions into heterogametes and homogametes
    bcftools view $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz -R $OUTDIR7/$SCAFF_NAME\_hetgam_dropout.bed $HOMGAM | \
    bgzip -c > $OUTDIR7t/$SCAFF_NAME\_HOMGAM.vcf.gz;
    tabix $OUTDIR7t/$SCAFF_NAME\_HOMGAM.vcf.gz;
    bcftools view $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz -R $OUTDIR7/$SCAFF_NAME\_hetgam_dropout.bed $HETGAM | \
    bgzip -c > $OUTDIR7t/$SCAFF_NAME\_HETGAM.vcf.gz;

    if [ $(bcftools view -H $OUTDIR7t/$SCAFF_NAME\_HETGAM.vcf.gz | awk '{print; exit}' | wc -l) -gt 0 ]; then

      ###### Re-code homozygose genotypes to haploid for the heterogametic sex
      python3 $FUNCTIONS/diploid2haploid.py -i $OUTDIR7t/$SCAFF_NAME\_HETGAM.vcf.gz \
      -o $OUTDIR7t/$SCAFF_NAME\_HETGAM_homogametic.vcf;
      bgzip $OUTDIR7t/$SCAFF_NAME\_HETGAM_homogametic.vcf;
      tabix $OUTDIR7t/$SCAFF_NAME\_HETGAM_homogametic.vcf.gz;

      ###### Activate bcftools 1.17 conda environment
      conda activate PhaseWY_bcftools1.17;

      ###### Merge homogametic haplotypes from the homogametic and heterogametic sex
      bcftools merge $OUTDIR7t/$SCAFF_NAME\_HOMGAM.vcf.gz \
      $OUTDIR7t/$SCAFF_NAME\_HETGAM_homogametic.vcf.gz | \
      bgzip -c > $OUTDIR7t/$SCAFF_NAME\_homogametic_temp4.vcf.gz;
      tabix $OUTDIR7t/$SCAFF_NAME\_homogametic_temp4.vcf.gz;

      ###### Dectivate bcftools 1.17 conda environment and return to PhaseWY bioinfo environment again
      conda deactivate;

      ###### Order individuals and filter sites
      bcftools view $OUTDIR7t/$SCAFF_NAME\_homogametic_temp4.vcf.gz -s $INDS_ORDER | \
      vcffixup - | vcffilter -f "AC > 0" | \
      bgzip -c > $OUTDIR7/$SCAFF_NAME\_homogametic.vcf.gz;
      tabix $OUTDIR7/$SCAFF_NAME\_homogametic.vcf.gz;

    fi;

  fi;

  ##### Remove temporary files
  rm -r $OUTDIR7t $OUTDIR2;

else
  cp $OUTDIR5/$SCAFF_NAME\_phased_all_variants.vcf.gz \
  $OUTDIR7/$SCAFF_NAME\_autosomal.vcf.gz;
  tabix $OUTDIR7/$SCAFF_NAME\_autosomal.vcf.gz;
fi;
