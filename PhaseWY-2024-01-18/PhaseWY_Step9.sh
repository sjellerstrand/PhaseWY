#!/bin/bash -l

# Version 2024-01-18
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

# Step 9. Create new VCFs
mkdir $OUTDIR/final_output/vcfs;
OUTDIRv=$OUTDIR/final_output/vcfs;
mkdir $OUTDIRv/removed_variants;

if [ $PHASE == "Yes" ]; then

  if [ $HAPLOIDS_ONLY == TRUE ]; then

    ### Copy input vcf
    bcftools view $VCF_IN.vcf.gz | \
    bgzip -c > $OUTDIRv/$PROJECT\_phased_all_variants.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_phased_all_variants.vcf.gz;

  else

    ### Create list of vcf files with all phased variants
    find $OUTDIR/5_shapeit4 -name "*_phased_all_variants.vcf.gz" \
    > $OUTDIRv/$PROJECT\_phased_all_variants_vcf_list.txt;

    ## Merge all phased marker vcfs back together
    bcftools concat -f $OUTDIRv/$PROJECT\_phased_all_variants_vcf_list.txt -Ou | \
    bcftools sort -T $OUTDIRv/temp_merge | \
    bgzip -c > $OUTDIRv/$PROJECT\_phased_all_variants_temp.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_phased_all_variants_temp.vcf.gz;
    bcftools view $OUTDIRv/$PROJECT\_phased_all_variants_temp.vcf.gz \
    -R $OUTDIRb/$PROJECT\_target_region.bed | \
    bgzip -c > $OUTDIRv/$PROJECT\_phased_all_variants.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_phased_all_variants.vcf.gz;
    rm $OUTDIRv/$PROJECT\_phased_all_variants_temp.vcf.gz*;
    
    ## Save sites that could not be phased
    find $OUTDIR/3_scaffold_input/ -name "*_biallelic.vcf.gz" \
    > $OUTDIRv/removed_variants/$PROJECT\_nonphased_variants_vcf_list.txt;
    if [ $(cat $OUTDIRv/removed_variants/$PROJECT\_nonphased_variants_vcf_list.txt | wc -l) -gt 0 ]; then
      bcftools concat -f $OUTDIRv/removed_variants/$PROJECT\_nonphased_variants_vcf_list.txt -Ou | \
      bcftools sort -T $OUTDIRv/temp_merge -Ou | \
      bcftools norm -m +any | \
      bgzip -c > $OUTDIRv/removed_variants/$PROJECT\_nonphased_variants_temp.vcf.gz;
      tabix $OUTDIRv/removed_variants/$PROJECT\_nonphased_variants_temp.vcf.gz;
      bcftools view $OUTDIRv/removed_variants/$PROJECT\_nonphased_variants_temp.vcf.gz \
      -R $OUTDIRb/$PROJECT\_target_region.bed -Ou | \
      bgzip -c > $OUTDIRv/removed_variants/$PROJECT\_nonphased_variants.vcf.gz;
      tabix $OUTDIRv/removed_variants/$PROJECT\_nonphased_variants.vcf.gz;
      rm $OUTDIRv/removed_variants/$PROJECT\_nonphased_variants_temp.vcf.gz*;
    fi;
    rm -r $OUTDIR/3_scaffold_input/* \
    $OUTDIRv/removed_variants/$PROJECT\_nonphased_variants_vcf_list.txt;

    ## Remove temporary files
    rm -r $OUTDIR/5_shapeit4/* \
    $OUTDIRv/$PROJECT\_phased_all_variants_vcf_list.txt;

  fi;

  if [ $SEX_SYSTEM != NA  ]; then

    ## Filter problematic sites
    bcftools view $OUTDIRv/$PROJECT\_phased_all_variants.vcf.gz \
    -R $OUTDIRb/$PROJECT\_target_region.bed | \
    vcftools --vcf - \
    --exclude-bed $OUTDIRb/$PROJECT\_unknown_phase.bed \
    --recode --recode-INFO-all --stdout | \
    vcftools --vcf - \
    --exclude-bed $OUTDIRb/$PROJECT\_sex_linked_ILS2.bed \
    --recode --recode-INFO-all --stdout | \
    bgzip -c > $OUTDIRv/$PROJECT\_phased_all_variants_filtered.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_phased_all_variants_filtered.vcf.gz;

    ## Save all sites that could not be phased

    ### Save all problematic sex-linked sites removed due to unknown phasing
    bcftools view $OUTDIRv/$PROJECT\_phased_all_variants.vcf.gz \
    -R $OUTDIRb/$PROJECT\_unknown_phase.bed | \
    bgzip -c > $OUTDIRv/removed_variants/$PROJECT\_unknown_phase.vcf.gz;
    tabix $OUTDIRv/removed_variants/$PROJECT\_unknown_phase.vcf.gz;

    ## Create autosomal dataset

    ### Create list of autosomal vcf files
    find $OUTDIR/7_concat_vcfs -name "*_autosomal.vcf.gz" \
    > $OUTDIRv/$PROJECT\_autosomal_vcf_list.txt;

    ### Merge all autosomal vcfs back together
    bcftools concat -f $OUTDIRv/$PROJECT\_autosomal_vcf_list.txt -Ou | \
    bcftools sort -T $OUTDIRv/temp_merge | \
    bgzip -c > $OUTDIRv/$PROJECT\_autosomal_temp.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_autosomal_temp.vcf.gz; 
    bcftools view $OUTDIRv/$PROJECT\_autosomal_temp.vcf.gz \
    -R $OUTDIRb/$PROJECT\_autosomal.bed | \
    vcftools --vcf - \
    --exclude-bed $OUTDIRb/$PROJECT\_unknown_phase.bed \
    --recode --recode-INFO-all --stdout | \
    bgzip -c > $OUTDIRv/$PROJECT\_autosomal.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_autosomal.vcf.gz;
    rm $OUTDIRv/$PROJECT\_autosomal_temp.vcf.gz*;

    ### Create vcf with autosomal snps and mnps
    bcftools view $OUTDIRv/$PROJECT\_autosomal.vcf.gz | \
    vcfclassify - | \
    vcffilter -s -f "!( INS | DEL )" | \
    vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIRv/$PROJECT\_autosomal_snps.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_autosomal_snps.vcf.gz;

    if [ $EXON != NA ]; then
      ### Create vcf with autosomal, neutral and polymorphic snps
      bcftools view $OUTDIRv/$PROJECT\_autosomal_snps.vcf.gz \
      -R $OUTDIRb/$PROJECT\_autosomal_neutral.bed | \
      vcffilter -s -f "!MNP" | \
      vcffilter -f "AF < 1" | \
      bgzip -c > $OUTDIRv/$PROJECT\_autosomal_neutral.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_autosomal_neutral.vcf.gz;

      ### Create vcf with pruned neutral autosomal snps
      plink --vcf $OUTDIRv/$PROJECT\_autosomal_neutral.vcf.gz \
      --double-id --allow-extra-chr --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.1 --out $OUTDIRv/$PROJECT\_autosomal_neutral;
      cat $OUTDIRv/$PROJECT\_autosomal_neutral.prune.in | awk '{gsub(":", "\t"); print}' \
      > $OUTDIRv/$PROJECT\_autosomal_neutral.keep_sites;
      vcftools --gzvcf $OUTDIRv/$PROJECT\_autosomal_neutral.vcf.gz \
      --positions $OUTDIRv/$PROJECT\_autosomal_neutral.keep_sites \
      --recode --recode-INFO-all --stdout | \
      bgzip -c > $OUTDIRv/$PROJECT\_autosomal_neutral_pruned.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_autosomal_neutral_pruned.vcf.gz;
    fi;

    if [ $CDS != NA ]; then
      ### Create vcf with all autosomal CDS variants
      bcftools view $OUTDIRv/$PROJECT\_autosomal.vcf.gz \
      -R $OUTDIRb/$PROJECT\_autosomal_CDS.bed | \
      vcffilter -f "AC > 0" | \
      bgzip -c > $OUTDIRv/$PROJECT\_autosomal_CDS.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_autosomal_CDS.vcf.gz;
    fi;

    ## Remove temporary files
    rm $(cat $OUTDIRv/$PROJECT\_autosomal_vcf_list.txt) \
    $OUTDIRv/$PROJECT\_autosomal_vcf_list.txt;
    
    if [ $EXON != NA ]; then
      rm $OUTDIRv/$PROJECT\_*_neutral.keep_sites \
      $OUTDIRv/$PROJECT\_*_neutral.log \
      $OUTDIRv/$PROJECT\_*_neutral.nosex \
      $OUTDIRv/$PROJECT\_*_neutral.prune.in \
      $OUTDIRv/$PROJECT\_*_neutral.prune.out;
    fi;

    ## Create homogametic dataset

    ### Create list of homogametic vcf files
    find $OUTDIR/7_concat_vcfs -name "*_homogametic.vcf.gz" \
    > $OUTDIRv/$PROJECT\_homogametic_vcf_list.txt;

    ### Merge all homogametic vcfs back together
    bcftools concat -a -f $OUTDIRv/$PROJECT\_homogametic_vcf_list.txt -Ou | \
    bcftools sort -T $OUTDIRv/temp_merge | \
    bgzip -c > $OUTDIRv/$PROJECT\_homogametic_temp0.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_homogametic_temp0.vcf.gz;
    bcftools view $OUTDIRv/$PROJECT\_homogametic_temp0.vcf.gz \
    -R $OUTDIRb/$PROJECT\_homogametic.bed | \
    vcftools --vcf - \
    --exclude-bed $OUTDIRb/$PROJECT\_unknown_phase.bed \
    --recode --recode-INFO-all --stdout |
    bgzip -c > $OUTDIRv/$PROJECT\_homogametic_temp1.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_homogametic_temp1.vcf.gz;
    rm $OUTDIRv/$PROJECT\_homogametic_temp0.vcf.gz*;

    ### Save all problematic sites that are excluded from the final datasets

    #### Save all problematic sex-linked sites removed due to signs of incomplete lineage sorting (heterogametic dropout=ILS1)
    bcftools view $OUTDIRv/$PROJECT\_homogametic_temp1.vcf.gz \
    -R $OUTDIRb/$PROJECT\_sex_linked_ILS1.bed | \
    vcftools --vcf - \
    --exclude-bed $OUTDIRb/$PROJECT\_hetgam_dropout.bed \
    --recode --recode-INFO-all --stdout | \
    bcftools norm --rm-dup all | \
    bgzip -c > $OUTDIRv/removed_variants/$PROJECT\_homogametic_ILS1.vcf.gz;
    tabix $OUTDIRv/removed_variants/$PROJECT\_homogametic_ILS1.vcf.gz;

    #### Save all problematic sex-linked sites removed  due to signs of incomplete lineage sorting (unsuccessful phasing=ILS2)
    bcftools view $OUTDIRv/$PROJECT\_homogametic_temp1.vcf.gz \
    -R $OUTDIRb/$PROJECT\_sex_linked_ILS2.bed -Ou | \
    bcftools norm --rm-dup all -Ov | \
    bgzip -c > $OUTDIRv/removed_variants/$PROJECT\_homogametic_ILS2.vcf.gz;
    tabix $OUTDIRv/removed_variants/$PROJECT\_homogametic_ILS2.vcf.gz;

    #### Filter problematic sites
    vcftools --gzvcf $OUTDIRv/$PROJECT\_homogametic_temp1.vcf.gz \
    --exclude-bed $OUTDIRb/$PROJECT\_sex_linked_ILS2.bed \
    --recode --recode-INFO-all --stdout | \
    vcftools --vcf - \
    --exclude-bed $OUTDIRb/$PROJECT\_sex_linked_ILS1_not_hetgam_dropout.bed \
    --recode --recode-INFO-all --stdout | \
    vcfclassify - | \
    bgzip -c > $OUTDIRv/$PROJECT\_homogametic_temp2.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_homogametic_temp2.vcf.gz;
    rm $OUTDIRv/$PROJECT\_homogametic_temp1.vcf.gz*;

    ### Create homogametic datasets with homogametes only
    HOMGAM=$(echo "-s "$(cat $INDS | awk -F'\t' '$3=="HOMGAM" {print $1}' | tr '\n' ','));
    HOMGAM=${HOMGAM::-1};

    ### Create vcf with variants from homogametes
    bcftools view $OUTDIRv/$PROJECT\_homogametic_temp2.vcf.gz $HOMGAM | \
    vcffixup - | vcffilter -f "AC > 0" | \
    vcftools --vcf - \
    --max-missing $MISSING \
    --recode --recode-INFO-all --stdout | \
    bcftools norm --rm-dup all -Ov | \
    bgzip -c > $OUTDIRv/$PROJECT\_homogametic_homogametes.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_homogametic_homogametes.vcf.gz;

    ### Create vcf with homogametic snps and mnps from homogametes
    bcftools view $OUTDIRv/$PROJECT\_homogametic_homogametes.vcf.gz | \
    vcfclassify - | \
    vcffilter -s -f "!( INS | DEL )" | \
    vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIRv/$PROJECT\_homogametic_snps_homogametes.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_homogametic_snps_homogametes.vcf.gz;

    if [ $EXON != NA ]; then
      ### Create vcf with homogametic, neutral and polymorphic snps from homogametes
      bcftools view $OUTDIRv/$PROJECT\_homogametic_snps_homogametes.vcf.gz \
      -R $OUTDIRb/$PROJECT\_homogametic_neutral.bed | \
      vcffilter -s -f "!MNP" | \
      vcffilter -f "AF < 1" | \
      bgzip -c > $OUTDIRv/$PROJECT\_homogametic_neutral_homogametes.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_homogametic_neutral_homogametes.vcf.gz;

      ### Create vcf with pruned neutral homogametic snps from homogametes
      plink --vcf $OUTDIRv/$PROJECT\_homogametic_neutral_homogametes.vcf.gz \
      --double-id --allow-extra-chr --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.1 --out $OUTDIRv/$PROJECT\_homogametic_neutral_homogametes;
      cat $OUTDIRv/$PROJECT\_homogametic_neutral_homogametes.prune.in | awk '{gsub(":", "\t"); print}' \
      > $OUTDIRv/$PROJECT\_homogametic_neutral_homogametes.keep_sites;
      vcftools --gzvcf $OUTDIRv/$PROJECT\_homogametic_neutral_homogametes.vcf.gz \
      --positions $OUTDIRv/$PROJECT\_homogametic_neutral_homogametes.keep_sites \
      --recode --recode-INFO-all --stdout | \
      bgzip -c > $OUTDIRv/$PROJECT\_homogametic_neutral_pruned_homogametes.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_homogametic_neutral_pruned_homogametes.vcf.gz;
    fi;

    if [ $CDS != NA ]; then
      ### Create vcf with homogametic CDS variants from homogametes
      bcftools view $OUTDIRv/$PROJECT\_homogametic_homogametes.vcf.gz \
      -R $OUTDIRb/$PROJECT\_homogametic_CDS.bed | \
      bgzip -c > $OUTDIRv/$PROJECT\_homogametic_CDS_homogametes.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_homogametic_CDS_homogametes.vcf.gz;
    fi;

    if [ $EXON != NA ]; then
      ### Remove temporary files
      rm $OUTDIRv/$PROJECT\_*_neutral_homogametes.keep_sites \
      $OUTDIRv/$PROJECT\_*_neutral_homogametes.log \
      $OUTDIRv/$PROJECT\_*_neutral_homogametes.nosex \
      $OUTDIRv/$PROJECT\_*_neutral_homogametes.prune.in \
      $OUTDIRv/$PROJECT\_*_neutral_homogametes.prune.out;
    fi;

    ## Create homogametic datasets with heterogametes only
    HETGAM=$(echo "-s "$(cat $INDS | awk -F'\t' '$3=="HETGAM" {print $1}' | tr '\n' ','));
    HETGAM=${HETGAM::-1};

    ### Create vcf with variants from heterogametes
    bcftools view $OUTDIRv/$PROJECT\_homogametic_temp2.vcf.gz $HETGAM | \
    vcffixup - | vcffilter -f "AC > 0" | \
    vcftools --vcf - \
    --max-missing $MISSING \
    --recode --recode-INFO-all --stdout | \
    bcftools norm --rm-dup all -Ov | \
    bgzip -c > $OUTDIRv/$PROJECT\_homogametic_heterogametes.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_homogametic_heterogametes.vcf.gz;

    ### Create vcf with homogametic snps and mnps from heterogametes
    bcftools view $OUTDIRv/$PROJECT\_homogametic_heterogametes.vcf.gz | \
    vcfclassify - | \
    vcffilter -s -f "!( INS | DEL )" | \
    vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIRv/$PROJECT\_homogametic_snps_heterogametes.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_homogametic_snps_heterogametes.vcf.gz;

    if [ $EXON != NA ]; then
      ### Create vcf with homogametic, neutral and polymorphic snps from heterogametes
      bcftools view $OUTDIRv/$PROJECT\_homogametic_snps_heterogametes.vcf.gz \
      -R $OUTDIRb/$PROJECT\_homogametic_neutral.bed | \
      vcffilter -s -f "!MNP" | \
      vcffilter -f "AF < 1" | \
      bgzip -c > $OUTDIRv/$PROJECT\_homogametic_neutral_heterogametes.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_homogametic_neutral_heterogametes.vcf.gz;

      ### Create vcf with pruned neutral homogametic snps from heterogametes
      plink --vcf $OUTDIRv/$PROJECT\_homogametic_neutral_heterogametes.vcf.gz \
      --double-id --allow-extra-chr --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.1 --out $OUTDIRv/$PROJECT\_homogametic_neutral_heterogametes;
      cat $OUTDIRv/$PROJECT\_homogametic_neutral_heterogametes.prune.in | awk '{gsub(":", "\t"); print}' \
      > $OUTDIRv/$PROJECT\_homogametic_neutral_heterogametes.keep_sites;
      vcftools --gzvcf $OUTDIRv/$PROJECT\_homogametic_neutral_heterogametes.vcf.gz \
      --positions $OUTDIRv/$PROJECT\_homogametic_neutral_heterogametes.keep_sites \
      --recode --recode-INFO-all --stdout | \
      bgzip -c > $OUTDIRv/$PROJECT\_homogametic_neutral_pruned_heterogametes.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_homogametic_neutral_pruned_heterogametes.vcf.gz;
    fi;

    if [ $CDS != NA ]; then
      ### Create vcf with homogametic CDS variants from heterogametes
      bcftools view $OUTDIRv/$PROJECT\_homogametic_heterogametes.vcf.gz \
      -R $OUTDIRb/$PROJECT\_homogametic_CDS.bed | \
      bgzip -c > $OUTDIRv/$PROJECT\_homogametic_CDS_heterogametes.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_homogametic_CDS_heterogametes.vcf.gz;
    fi;

    if [ $EXON != NA ]; then
      ### Remove temporary files
      rm $OUTDIRv/$PROJECT\_*_neutral_heterogametes.keep_sites \
      $OUTDIRv/$PROJECT\_*_neutral_heterogametes.log \
      $OUTDIRv/$PROJECT\_*_neutral_heterogametes.nosex \
      $OUTDIRv/$PROJECT\_*_neutral_heterogametes.prune.in \
      $OUTDIRv/$PROJECT\_*_neutral_heterogametes.prune.out;
    fi;

    ## Create homogametic dataset with all individuals and missingness filters

    ### Filter missing variants
    bcftools view $OUTDIRv/$PROJECT\_homogametic_temp2.vcf.gz | \
    vcffixup - | vcffilter -f "AC > 0" | \
    vcftools --vcf - \
    --max-missing $MISSING \
    --recode --recode-INFO-all --stdout | \
    bcftools norm --rm-dup all -Ov | \
    bgzip -c > $OUTDIRv/$PROJECT\_homogametic.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_homogametic.vcf.gz;

    ### Create vcf with homogametic snps and mnps
    bcftools view $OUTDIRv/$PROJECT\_homogametic.vcf.gz \
    -R $OUTDIRb/$PROJECT\_homogametic.bed | \
    vcfclassify - | \
    vcffilter -s -f "!( INS | DEL )" | \
    vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIRv/$PROJECT\_homogametic_snps.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_homogametic_snps.vcf.gz;

    if [ $EXON != NA ]; then
      ### Create vcf with homogametic, neutral and polymorphic snps
      bcftools view $OUTDIRv/$PROJECT\_homogametic_snps.vcf.gz \
      -R $OUTDIRb/$PROJECT\_homogametic_neutral.bed | \
      vcffilter -s -f "!MNP" | \
      vcffilter -f "AF < 1" | \
      bgzip -c > $OUTDIRv/$PROJECT\_homogametic_neutral.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_homogametic_neutral.vcf.gz;

      ### Create vcf with pruned neutral homogametic snps
      plink --vcf $OUTDIRv/$PROJECT\_homogametic_neutral.vcf.gz \
      --double-id --allow-extra-chr --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.1 --out $OUTDIRv/$PROJECT\_homogametic_neutral;
      cat $OUTDIRv/$PROJECT\_homogametic_neutral.prune.in | awk '{gsub(":", "\t"); print}' \
      > $OUTDIRv/$PROJECT\_homogametic_neutral.keep_sites;
      vcftools --gzvcf $OUTDIRv/$PROJECT\_homogametic_neutral.vcf.gz \
      --positions $OUTDIRv/$PROJECT\_homogametic_neutral.keep_sites \
      --recode --recode-INFO-all --stdout | \
      bgzip -c > $OUTDIRv/$PROJECT\_homogametic_neutral_pruned.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_homogametic_neutral_pruned.vcf.gz;
    fi;

    if [ $CDS != NA ]; then
      ### Create vcf with homogametic CDS variants
      bcftools view $OUTDIRv/$PROJECT\_homogametic.vcf.gz \
      -R $OUTDIRb/$PROJECT\_homogametic_CDS.bed | \
      bgzip -c > $OUTDIRv/$PROJECT\_homogametic_CDS.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_homogametic_CDS.vcf.gz;
    fi;

    ### Remove temporary files
    rm $(cat $OUTDIRv/$PROJECT\_homogametic_vcf_list.txt) \
    $OUTDIRv/$PROJECT\_homogametic_vcf_list.txt \
    $OUTDIRv/$PROJECT\_homogametic_temp*;
      
    if [ $EXON != NA ]; then
      rm $OUTDIRv/$PROJECT\_*_neutral.keep_sites \
      $OUTDIRv/$PROJECT\_*_neutral.log \
      $OUTDIRv/$PROJECT\_*_neutral.nosex \
      $OUTDIRv/$PROJECT\_*_neutral.prune.in \
      $OUTDIRv/$PROJECT\_*_neutral.prune.out;
    fi;

    ## Create heterogametic dataset

    ### Create list of heterogametic vcf files
    find $OUTDIR/7_concat_vcfs -name "*_heterogametic.vcf.gz" \
    > $OUTDIRv/$PROJECT\_heterogametic_vcf_list.txt;

    ### Merge all heterogametic vcfs back together
    bcftools concat -f $OUTDIRv/$PROJECT\_heterogametic_vcf_list.txt -Ou | \
    bcftools sort -T $OUTDIRv/temp_merge | \
    bgzip -c > $OUTDIRv/$PROJECT\_heterogametic_temp0.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_heterogametic_temp0.vcf.gz;
    bcftools view $OUTDIRv/$PROJECT\_heterogametic_temp0.vcf.gz \
    -R $OUTDIRb/$PROJECT\_heterogametic.bed | \
    vcftools --vcf - \
    --exclude-bed $OUTDIRb/$PROJECT\_unknown_phase.bed \
    --recode --recode-INFO-all --stdout |
    bgzip -c > $OUTDIRv/$PROJECT\_heterogametic_temp1.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_heterogametic_temp1.vcf.gz;
    rm $OUTDIRv/$PROJECT\_heterogametic_temp0.vcf.gz*;

    ### Save all problematic sites that are excluded from the final datasets

    ### Save all problematic sex-linked sites removed  due to signs of incomplete linage sorting (heterogametic dropout=ILS1)
    bcftools view $OUTDIRv/$PROJECT\_heterogametic_temp1.vcf.gz \
    -R $OUTDIRb/$PROJECT\_sex_linked_ILS1.bed | \
    bgzip -c > $OUTDIRv/removed_variants/$PROJECT\_heterogametic_ILS1.vcf.gz;
    tabix $OUTDIRv/removed_variants/$PROJECT\_heterogametic_ILS1.vcf.gz;

    ### Save all problematic sex-linked sites removed  due to signs of incomplete linage sorting (unsuccessful phasing=ILS2)
    bcftools view $OUTDIRv/$PROJECT\_heterogametic_temp1.vcf.gz \
    -R $OUTDIRb/$PROJECT\_sex_linked_ILS2.bed | \
    bgzip -c > $OUTDIRv/removed_variants/$PROJECT\_heterogametic_ILS2.vcf.gz;
    tabix $OUTDIRv/removed_variants/$PROJECT\_heterogametic_ILS2.vcf.gz;

    #### Filter problematic sites
    vcftools --gzvcf $OUTDIRv/$PROJECT\_heterogametic_temp1.vcf.gz \
    --exclude-bed $OUTDIRb/$PROJECT\_sex_linked_ILS2.bed \
    --recode --recode-INFO-all --stdout | \
    vcftools --vcf - \
    --exclude-bed $OUTDIRb/$PROJECT\_sex_linked_ILS1.bed \
    --recode --recode-INFO-all --stdout | \
    vcfclassify - | \
    bgzip -c > $OUTDIRv/$PROJECT\_heterogametic.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_heterogametic.vcf.gz;
    rm $OUTDIRv/$PROJECT\_heterogametic_temp1.vcf.gz*;

    ### Create vcf with heterogametic snps and mnps
    bcftools view $OUTDIRv/$PROJECT\_heterogametic.vcf.gz \
    -R $OUTDIRb/$PROJECT\_heterogametic.bed | \
    vcfclassify - | \
    vcffilter -s -f "!( INS | DEL )" | \
    vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIRv/$PROJECT\_heterogametic_snps.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_heterogametic_snps.vcf.gz;

    if [ $EXON != NA ]; then
      ### Create vcf with heterogametic, neutral and polymorphic snps
      bcftools view $OUTDIRv/$PROJECT\_heterogametic_snps.vcf.gz \
      -R $OUTDIRb/$PROJECT\_heterogametic_neutral.bed | \
      vcffilter -s -f "!MNP" | \
      vcffilter -f "AF < 1" | \
      bgzip -c > $OUTDIRv/$PROJECT\_heterogametic_neutral.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_heterogametic_neutral.vcf.gz;

      ### Create vcf with pruned neutral heterogametic snps
      plink --vcf $OUTDIRv/$PROJECT\_heterogametic_neutral.vcf.gz \
      --double-id --allow-extra-chr --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.1 --out $OUTDIRv/$PROJECT\_heterogametic_neutral;
      cat $OUTDIRv/$PROJECT\_heterogametic_neutral.prune.in | awk '{gsub(":", "\t"); print}' \
      > $OUTDIRv/$PROJECT\_heterogametic_neutral.keep_sites;
      vcftools --gzvcf $OUTDIRv/$PROJECT\_heterogametic_neutral.vcf.gz \
      --positions $OUTDIRv/$PROJECT\_heterogametic_neutral.keep_sites \
      --recode --recode-INFO-all --stdout | \
      bgzip -c > $OUTDIRv/$PROJECT\_heterogametic_neutral_pruned.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_heterogametic_neutral_pruned.vcf.gz;
    fi;

    if [ $CDS != NA ]; then
      ### Create vcf with heterogametic CDS variants
      bcftools view $OUTDIRv/$PROJECT\_heterogametic.vcf.gz \
      -R $OUTDIRb/$PROJECT\_heterogametic_CDS.bed | \
      bgzip -c > $OUTDIRv/$PROJECT\_heterogametic_CDS.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_heterogametic_CDS.vcf.gz;
    fi;

    ### Remove temporary files
    rm -r $OUTDIRv/$PROJECT\_heterogametic_vcf_list.txt \
    $OUTDIR/7_concat_vcfs/*;
    
    if [ $EXON != NA ]; then
      $OUTDIRv/$PROJECT\_*_neutral.keep_sites \
      $OUTDIRv/$PROJECT\_*_neutral.log \
      $OUTDIRv/$PROJECT\_*_neutral.nosex \
      $OUTDIRv/$PROJECT\_*_neutral.prune.in \
      $OUTDIRv/$PROJECT\_*_neutral.prune.out;
    fi;

  else

    ## Create vcf with snps and mnps
    vcfclassify $OUTDIRv/$PROJECT\_phased_all_variants.vcf.gz | \
    vcffilter -s -f "!( INS | DEL )" | \
    vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIRv/$PROJECT\_snps.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_snps.vcf.gz;

    if [ $EXON != NA ]; then
      ## Create vcf with neutral and polymorphic snps
      bcftools view $OUTDIRv/$PROJECT\_snps.vcf.gz \
      -R $OUTDIRb/$PROJECT\_neutral.bed | \
      vcffilter -s -f "!MNP" | \
      vcffilter -f "AF < 1" | \
      bgzip -c > $OUTDIRv/$PROJECT\_neutral.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_neutral.vcf.gz;

      ## Create vcf with pruned neutral snps
      plink --vcf $OUTDIRv/$PROJECT\_neutral.vcf.gz \
      --double-id --allow-extra-chr --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.1 --out $OUTDIRv/$PROJECT\_neutral;
      cat $OUTDIRv/$PROJECT\_neutral.prune.in | awk '{gsub(":", "\t"); print}' \
      > $OUTDIRv/$PROJECT\_neutral.keep_sites;
      vcftools --gzvcf $OUTDIRv/$PROJECT\_neutral.vcf.gz \
      --positions $OUTDIRv/$PROJECT\_neutral.keep_sites \
      --recode --recode-INFO-all --stdout | \
      bgzip -c > $OUTDIRv/$PROJECT\_neutral_pruned.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_neutral_pruned.vcf.gz;
    fi;

    if [ $CDS != NA ]; then
      ## Create vcf with all CDS variants
      bcftools view $OUTDIRv/$PROJECT\_phased_all_variants.vcf.gz \
      -R $OUTDIRb/$PROJECT\_CDS.bed | \
      vcffilter -f "AC > 0" | \
      bgzip -c > $OUTDIRv/$PROJECT\_CDS.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_CDS.vcf.gz;
    fi;

    if [ $EXON != NA ]; then
      ## Remove temporary files
      rm $OUTDIRv/$PROJECT\_neutral.keep_sites \
      $OUTDIRv/$PROJECT\_neutral.log \
      $OUTDIRv/$PROJECT\_neutral.nosex \
      $OUTDIRv/$PROJECT\_neutral.prune.in \
      $OUTDIRv/$PROJECT\_neutral.prune.out;
    fi;

  fi;

  if [ PHASE_SCAFFOLD_LIST != NA ] && [ LIST_ONLY == No ]; then

    ### Create vcf with unphased variants
    bcftools filter -t ^$(cat $OUTDIR/0_file_info/SCAFFOLDS.txt | cut -f1 | tr '\n' ',') $VCF_IN.vcf.gz | \
    vcffixup - | vcffilter -f "AC > 0" | \
    vcftools --vcf - \
    --max-missing $MISSING \
    --recode --recode-INFO-all --stdout | \
    bgzip -c > $OUTDIRv/$PROJECT\_nonphased_variants.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_nonphased_variants.vcf.gz;

    ## Create vcf with snps and mnps
    vcfclassify $OUTDIRv/$PROJECT\_nonphased_variants.vcf.gz | \
    vcffilter -s -f "!( INS | DEL )" | \
    vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIRv/$PROJECT\_nonphased_snps.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_nonphased_snps.vcf.gz;

    if [ $EXON != NA ]; then
      ## Create vcf with neutral and polymorphic snps
      bcftools view $OUTDIRv/$PROJECT\_nonphased_snps.vcf.gz \
      -R $OUTDIRb/$PROJECT\_neutral.bed | \
      vcffilter -s -f "!MNP" | \
      vcffilter -f "AF < 1" | \
      bgzip -c > $OUTDIRv/$PROJECT\_nonphased_neutral.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_nonphased_neutral.vcf.gz;

      ## Create vcf with pruned neutral autosomal snps
      plink --vcf OUTDIRv/$PROJECT\_nonphased_neutral.vcf.gz \
      --double-id --allow-extra-chr --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.1 --out $OUTDIRv/$PROJECT\_nonphased_neutral;
      cat $OUTDIRv/$PROJECT\_nonphased_neutral.prune.in | awk '{gsub(":", "\t"); print}' \
      > $OUTDIRv/$PROJECT\_nonphased_neutral.keep_sites;
      vcftools --gzvcf $OUTDIRv/$PROJECT\_nonphased_neutral.vcf.gz \
      --positions $OUTDIRv/$PROJECT\_nonphased_neutral.keep_sites \
      --recode --recode-INFO-all --stdout | \
      bgzip -c > $OUTDIRv/$PROJECT\_nonphased_neutral_pruned.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_nonphased_neutral_pruned.vcf.gz;
    fi;

    if [ $CDS != NA ]; then
      ## Create vcf with all autosomal CDS variants
      bcftools view $OUTDIRv/$PROJECT\_nonphased_variants.vcf.gz \
      -R $OUTDIRb/$PROJECT\_CDS.bed | \
      vcffilter -f "AC > 0" | \
      bgzip -c > $OUTDIRv/$PROJECT\_nonphased_CDS.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_nonphased_CDS.vcf.gz;
    fi;

    if [ $EXON != NA ]; then
      ## Remove temporary files
      rm -r $OUTDIRv/$PROJECT\_nonphased_neutral.keep_sites \
      $OUTDIRv/$PROJECT\_nonphased_neutral.log \
      $OUTDIRv/$PROJECT\_nonphased_neutral.nosex \
      $OUTDIRv/$PROJECT\_nonphased_neutral.prune.in \
      $OUTDIRv/$PROJECT\_nonphased_neutral.prune.out;
    fi;

  fi;

else

    ## Copy input vcf
    cp $VCF_IN.vcf.gz $OUTDIRv/$PROJECT\_nonphased_all_variants.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_nonphased_all_variants.vcf.gz;

    ## Create vcf with snps and mnps
    vcfclassify $OUTDIRv/$PROJECT\_nonphased_all_variants.vcf.gz | \
    vcffilter -s -f "!( INS | DEL )" | \
    vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIRv/$PROJECT\_snps.vcf.gz;
    tabix $OUTDIRv/$PROJECT\_snps.vcf.gz;

    if [ $EXON != NA ]; then
      ## Create vcf with neutral and polymorphic snps
      bcftools view $OUTDIRv/$PROJECT\_snps.vcf.gz \
      -R $OUTDIRb/$PROJECT\_neutral.bed | \
      vcffilter -s -f "!MNP" | \
      vcffilter -f "AF < 1" | \
      bgzip -c > $OUTDIRv/$PROJECT\_neutral.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_neutral.vcf.gz;

      ## Create vcf with pruned neutral autosomal snps
      plink --vcf $OUTDIRv/$PROJECT\_neutral.vcf.gz \
      --double-id --allow-extra-chr --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.1 --out $OUTDIRv/$PROJECT\_neutral;
      cat $OUTDIRv/$PROJECT\_neutral.prune.in | awk '{gsub(":", "\t"); print}' \
      > $OUTDIRv/$PROJECT\_neutral.keep_sites;
      vcftools --gzvcf $OUTDIRv/$PROJECT\_neutral.vcf.gz \
      --positions $OUTDIRv/$PROJECT\_neutral.keep_sites \
      --recode --recode-INFO-all --stdout | \
      bgzip -c > $OUTDIRv/$PROJECT\_neutral_pruned.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_neutral_pruned.vcf.gz;
    fi;

    if [ $CDS != NA ]; then
      ## Create vcf with all autosomal CDS variants
      bcftools view $OUTDIRv/$PROJECT\_nonphased_all_variants.vcf.gz \
      -R $OUTDIRb/$PROJECT\_CDS.bed | \
      vcffilter -f "AC > 0" | \
      bgzip -c > $OUTDIRv/$PROJECT\_CDS.vcf.gz;
      tabix $OUTDIRv/$PROJECT\_CDS.vcf.gz;
    fi;

    if [ $EXON != NA ]; then
      ## Remove temporary files
      rm -r $OUTDIRv/$PROJECT\_neutral.keep_sites \
      $OUTDIRv/$PROJECT\_neutral.log \
      $OUTDIRv/$PROJECT\_neutral.nosex \
      $OUTDIRv/$PROJECT\_neutral.prune.in \
      $OUTDIRv/$PROJECT\_neutral.prune.out;
    fi;

fi;
