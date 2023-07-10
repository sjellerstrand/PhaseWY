# PhaseWY
A pipeline for phasing genomic data in general and homologous sex chromosome sequences in particular. The pipeline uses several females and males to identify and to seperate the sex chromosome sequences in the heterogametic sex (XY males or ZW females).

## Introduction


## Input data
- Variants: Variants in vcf-format from males and females that will be phased. Biallelic SNPs, as well as indels and MNPs are supported. It is up to the user to properly filter the vcf-file before running the pipeline. Disclaimer: The pipeline has so far only been tested for vcf-files produced by freebayes. However, other callers should be supported!
- Aligned sequences: The alignment of each sample in bam-format. The pipeline will use all reads available in the bam-files. It is therrefore up to the user to properly filter the bam files from duplicates and improper pairing (e.g. samtools -f2 -F260 -q20). The files are used to:
  * Identify callable regions
  * Sex-differences in sequence depth
  * Read based phasing<a/>
- Reference genome: The reference genome used to aligned the bam files and call variants.
- Variant statistics: The parameters used for filtering variants, i.e. minimum depth, minumum mean depth, maximum mean depth, and missingness should be provided to identify non-polymorphic regions of the genome as callable. Mean depth per individual should further be provided for normalisaion to identify sex-differences in sequence depth.
- Repetitve regions: Calling is assumed to be unreliable in repetitive regions. These are therefore excluded in Step1 and set as missing data (repetititve regions can be estimated with https://www.repeatmasker.org/).
