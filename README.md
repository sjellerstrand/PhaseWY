# PhaseWY
A pipeline for phasing genomic data in general and homologous sex chromosome sequences in particular. The pipeline uses several females and males to identify and to seperate the sex chromosome sequences in the heterogametic sex (XY males or ZW females).

## Introduction
In species with homologous sex chromosomes, retrieving sequence data from Y and W is challenging. This is unfortunate as such data would facilitate understanding interesting aspects of sex chromosome evolution such as, neo-sex chromosome formation, the degeneration of heterogametic Y and W chromosomes, dosage compensation on the X and Z chromosomes in response to Y and W loss of function, and the accumulation of sexually antagonistic mutations on X and Z. Here we present PhaseWY, a bioinformatic pipeline that aims to identify and extracts Y and W sequences by phasing and clustering haplotypes from multiple females and males.

![Plot](https://github.com/sjellerstrand/PhaseWY/blob/8d0f1a702563986c4c5611cc220fe53a14a27afe/Pipeline%20flowchart.jpg)

To run the pipeline, only knowledge of the sex determination system and the sex of each sequenced individual is required. The pipeline runs independently per scaffold and goes through several key steps. (1) It identifies callable regions of the genome, and (2) sex-linked regions based on sex differences in sequencing depth. Then, it (3) subsets variants per scaffold, it (4) performs read-based phasing with WhatsHap, and (5) statistical phasing with Shapeit4. It (6) identifies sex-linked regions based on clustering of phased haplotypes. Finally, it (7) re-genotypes the heterogametic sex to haploids for both the homogametic and the heterogametic sex chromosomes. The pipeline summarises the genomic regions in bed format as autosomal, homogametic, and heterogametic. It further outputs phased variants in vcf format for the corresponding regions. The pipeline is not limited to phasing sex chromosomes and performs phasing on strictly autosomal data. If the dataset contains haploid and diploid individuals (e.g., haplodiploid Hymenoptera), haploid individuals can be used as a haplotype reference panel for statistical phasing of the diploid individuals.

## Input data
- Variants: Variants in vcf-format from males and females that will be phased. Biallelic SNPs, as well as fixed alternate alleles, indels and MNPs are supported. It is up to the user to properly filter the vcf-file before running the pipeline. Disclaimer: The pipeline has so far only been tested for vcf-files produced by freebayes. However, other callers should be supported!
- Aligned sequences: The alignment of each sample in bam-format with corresponding sample IDs in SM-tag of the @RG header to the ID in the vcf-file. The pipeline will use all reads available in the bam-files. It is therrefore up to the user to properly filter the bam files from duplicates and improper pairing (e.g. samtools -f2 -F260 -q20). The files are used to:
  * Identify callable regions
  * Sex-differences in sequence depth
  * Read based phasing<a/>
- Masked regions: Calling is assumed to be unreliable in regions such as repetitive regions. The user may also want to mask some regions of the genome for other reasons. These are therefore excluded in Step 1 and set as missing data (repetititve regions can be estimated with https://www.repeatmasker.org/).
- Reference genome: The reference genome used to aligned the bam files and call variants.
- Variant filtering paramters: The parameters used for filtering variants, i.e. minimum depth, minumum mean depth, maximum mean depth, and missingness should be provided to identify non-polymorphic regions of the genome as callable. Mean depth per individual should further be provided for normalisaion to identify sex-differences in sequence depth.

## Current status and disclaimers
The core of the pipeline has been built and tested for several datasets with mostly optimistic result. The pipeline has since this version been improved and converted into a snakemake pipeline which will be published at the end of 2025 or during 2026.

The pipeline is developed for use with short-redas such as illumina. No other data types has so far been tested.

Current tests work very well at identifying sex linked region. Retrieving W-sequences from relatively recent neo-sex chromosomes so far works well for low diversity species, but seem to be a bit problematic for high diversity species.
