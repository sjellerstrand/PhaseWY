# Version 2023-07-10
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

# Step 10. Print some vcf stats
mkdir $OUTDIR/final_output/stats;
OUTDIRs=$OUTDIR/final_output/stats;

## Create list of vcf files
find $OUTDIRv -name "*.vcf.gz" | grep -v pruned | sort \
> $OUTDIRs/$PROJECT\_vcf_list.txt;

export OUTDIRs;

## Define function to run in parallell
print_stats() {

### Input parameter
FILE=$1;

### Set up analysis
FILE_PREFIX=$(echo $FILE | rev | cut -d'/' -f 1 | cut -d'.' -f 3- | rev);
mkdir $OUTDIRs/$FILE_PREFIX;
cp $FILE* $OUTDIRs/$FILE_PREFIX;
VCF_OUT=$OUTDIRs/$FILE_PREFIX/$FILE_PREFIX;

### Run scripts
source $FUNCTIONS/filter_stats.sh;
source $FUNCTIONS/pca.sh;

### Remove vcf
rm $VCF_OUT.vcf.gz*;

}

## Excecute function in parallell
export -f print_stats;
parallel 'print_stats {1}' :::: $OUTDIRs/$PROJECT\_vcf_list.txt;

## Remove temporary files
rm $OUTDIRs/$PROJECT\_vcf_list.txt;

if [ $SEX_SYSTEM != NA  ] && [ $PHASE == "Yes" ]; then

  ## Print genome summary
  mkdir $OUTDIRs/genome_summary;
  mkdir $OUTDIRs/genome_summary/focal;

  bedtools complement -i $OUTDIRb/$PROJECT\_target_region_phased.bed -g $REF.fai | \
  awk -F'\t' '{print $0"\tMissing data"}' \
  > $OUTDIRs/genome_summary/focal/$PROJECT\_Missing_data.bed;

  bedtools intersect -a $OUTDIRb/$PROJECT\_target_region_phased.bed \
  -b $OUTDIRb/$PROJECT\_autosomal.bed | \
  awk -F'\t' '{print $0"\tAutosomal"}' \
  > $OUTDIRs/genome_summary/focal/$PROJECT\_Autosomal.bed;

  bedtools intersect -a $OUTDIRb/$PROJECT\_target_phase_hetgam_dropout.bed \
  -b $OUTDIRb/$PROJECT\_target_phase_sex_linked.bed | \
  awk -F'\t' '{print $0"\tSex phase & depth difference"}' \
  > $OUTDIRs/genome_summary/focal/$PROJECT\_Sex_phase_depth_difference.bed;

  bedtools subtract -a $OUTDIRb/$PROJECT\_target_phase_hetgam_dropout.bed \
  -b $OUTDIRs/genome_summary/focal/$PROJECT\_Sex_phase_depth_difference.bed | \
  awk -F'\t' '{print $0"\tSex depth difference"}' \
  > $OUTDIRs/genome_summary/focal/$PROJECT\_Sex_depth_difference.bed;

  bedtools subtract -a $OUTDIRb/$PROJECT\_target_phase_sex_linked.bed \
  -b $OUTDIRs/genome_summary/focal/$PROJECT\_Sex_phase_depth_difference.bed | \
  awk -F'\t' '{print $0"\tSex phase difference"}' \
  > $OUTDIRs/genome_summary/focal/$PROJECT\_Sex_phase_difference.bed;

  cat $OUTDIRs/genome_summary/focal/$PROJECT\_Missing_data.bed \
  $OUTDIRs/genome_summary/focal/$PROJECT\_Autosomal.bed \
  $OUTDIRs/genome_summary/focal/$PROJECT\_Sex_phase_depth_difference.bed \
  $OUTDIRs/genome_summary/focal/$PROJECT\_Sex_depth_difference.bed \
  $OUTDIRs/genome_summary/focal/$PROJECT\_Sex_phase_difference.bed \
  > $OUTDIRs/genome_summary/focal/$PROJECT\_Genome_summary.bed;

  rm $OUTDIRs/genome_summary/focal/$PROJECT\_Missing_data.bed \
  $OUTDIRs/genome_summary/focal/$PROJECT\_Autosomal.bed \
  $OUTDIRs/genome_summary/focal/$PROJECT\_Sex_phase_depth_difference.bed \
  $OUTDIRs/genome_summary/focal/$PROJECT\_Sex_depth_difference.bed \
  $OUTDIRs/genome_summary/focal/$PROJECT\_Sex_phase_difference.bed;

 # GENOME=$REF.fai;

#  Rscript $FUNCTIONS/genome_summary.r --no-save --args PROJECT=$PROJECT OUTDIR=$OUTDIRs/genome_summary/focal GENOME=$GENOME;

#  if [ $(ls $SYNTENY | wc -l) -eq 1 ]; then

    #mkdir $OUTDIRs/genome_summary/synteny;
   # GENOME=
   # Rscript $FUNCTIONS/genome_summary.r --no-save --args PROJECT=$PROJECT OUTDIR=$OUTDIRs/genome_summary/synteny GENOME=$GENOME;
#  fi;

fi;
