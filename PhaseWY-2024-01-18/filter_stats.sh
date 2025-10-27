# Version 2023-03-07
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

### Retrive statistics
## Number of variants
bcftools view -H $VCF_OUT.vcf.gz | wc -l > $VCF_OUT.number_variants.txt;
## Site quality
vcftools --gzvcf $VCF_OUT.vcf.gz --site-quality --out $VCF_OUT;
## Mean depth per site
vcftools --gzvcf $VCF_OUT.vcf.gz --site-mean-depth --out $VCF_OUT;
## Proportion of missing data per site
vcftools --gzvcf $VCF_OUT.vcf.gz --missing-site --out $VCF_OUT;
## Allele frequency
vcftools --gzvcf $VCF_OUT.vcf.gz --freq2 --out $VCF_OUT --max-alleles 2;
## Mean depth per individual
vcftools --gzvcf $VCF_OUT.vcf.gz --depth --out $VCF_OUT;
## Heterozygosity and inbreeding coefficient per individual
if [ $NONDIPLOIDS == Yes ]; then
     DIPLOIDS=$(cat $METADATA/ploidy.txt | awk -F'\t' '{if($2==1) print "--remove-indv",$1}' | tr '\n' ' ');
     vcftools --gzvcf $VCF_OUT.vcf.gz --het $DIPLOIDS --out $VCF_OUT;
else
     vcftools --gzvcf $VCF_OUT.vcf.gz --het --out $VCF_OUT;
fi;
## Proportion of missing data per individual
vcftools --gzvcf $VCF_OUT.vcf.gz --missing-indv --out $VCF_OUT;
## Proportion of missing data per site
vcftools --gzvcf $VCF_OUT.vcf.gz --missing-site --out $VCF_OUT;

### Plot statistics
Rscript $FUNCTIONS/filter_stats.r \
--args VCF_OUT=$VCF_OUT;
