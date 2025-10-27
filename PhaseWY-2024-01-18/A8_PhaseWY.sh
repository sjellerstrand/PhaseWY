#!/bin/bash -l

#SBATCH -A snic2022-5-484
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 50:00:00
#SBATCH -J PhaseWY
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

# Version 2024-01-18
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

# Set parameters
PHASE=Yes                # Should phasing be performed? Alternatively, only output callable sites (step 1). "Yes" or "No" are accepted.
CDS_SCAFFOLDS_ONLY=Yes   # Analysis is only performed on scaffolds annotated with coding regions. Requres an annotation file. "Yes or "No" are accepted.
MIN_SCAFFOLD_SIZE=0     # Analysis is only performed on scaffolds longer than the specified nubmer of base pairs. Accepts any positive integer.
SEX_SYSTEM=ZW           # Sex-chromosome system. "ZW", "XY", or "NA" are accepted.
NONDIPLOIDS=No          # Are there any non-diploids? "Yes" or "No" are accepted.
THRESHOLD=1             # Ratio of heterogametes needed in smallest cluster for it to be called sex linked.
WINDOW=1000               # Number of variants used per window in the initial search.
STEP=250                # Number of variants for each sliding window step in the initial search.
PHASE_SCAFFOLD_LIST=NA    # List of scaffolds to be phased. Useful if sex chromosomes are already known, and there is no interest in phasing the autosomes. If "PHASE=Yes" and "LIST_ONLY=No", callable sites will still be identified and output for the autosomes. "NA" or path to list accepted. The list should contain one scaffold per row in the first column.
LIST_ONLY=No            # Should analysis only be performed on listed scaffolds given by "PHASE_SCAFFOLD_LIST"? "Yes" or "No" are accepted.

conda=/crex/proj/snic2020-2-25/nobackup/simon/conda;
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Rasolark_2021;
PROJECT2=Skylark_2021; # Annotation lift over
SYNTENY=Skylark_2021; # Synteny
WORKDIR=$MAINDIR/data/$PROJECT;
WORKDIR2=$MAINDIR/data/$PROJECT2;
VCF_IN=$WORKDIR/A7_filter_variants/filter3\
/$PROJECT\_3;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
REF=$WORKDIR/A4_make_new_reference/$PROJECT\_consensus_reference.fasta;
CDS=$WORKDIR2/B3_annotation_lift_over/Aarv_cds.bed;                              # If provided, bed and vcf files of callable CDS regions will be outputed. Also needed if CDS_SCAFFOLDS_ONLY=Yes is specified. Path or "NA" accepted.
EXON=$WORKDIR2/B3_annotation_lift_over/Aarv_exonic.bed;                          # If provided, bed and vcf files of callable regions assumed to be neutral (20 000 bp form closest exon will be outputed). Path or "NA" accepted.
GENES=NA;									 # If provided, annotated features within callable regions will be outputed, as well as their location within autosomal, homogametic and heterogametic regions. Path or "NA" accepted.
REPEATS=$WORKDIR2/B1_prepare_reference/$PROJECT2\_masked_repeats.bed;            # If provided, Repeat regions will be subtracted from callable regions . Path or "NA" accepted
SAMPLE_INFO=$WORKDIR/metadata/INDS.txt;
FUNCTIONS=/crex/proj/snic2020-2-25/nobackup/simon/bin/PhaseWY-2024-01-18;

# Import parameters from filters
MIN_MEAN=$(cat $WORKDIR/A7_filter_variants/filter2/filtering_parameters.txt | grep "MIN_MEAN" |cut -d'=' -f2); # The minimum accepted mean depth for a site to be kept in the analysis.
MAX_MEAN=$(cat $WORKDIR/A7_filter_variants/filter2/filtering_parameters.txt | grep "MAX_MEAN" |cut -d'=' -f2); # The maximum accepted mean depth for a site to be kept in the analysis.
MIN_DP=$(cat $WORKDIR/A7_filter_variants/filter2/filtering_parameters.txt | grep "MIN_DP" |cut -d'=' -f2); # The minimum accepted depth at an individual for a site to be kept in the analysis.
MISSING=$(cat $WORKDIR/A7_filter_variants/filter2/filtering_parameters.txt | grep "MISSING" |cut -d'=' -f2); # The minimum accepted missingness for a site to be kept in the analysis

# Load modules
module load bioinfo-tools R_packages/4.0.0 gnuparallel/20180822;

# Activate conda environment
conda activate PhaseWY_bioinfo;

# Load bcftools/1.17, commit b7b2a32
export PATH=$PATH:/crex/proj/snic2020-2-25/nobackup/simon/softwares/bcftools;

# Step 0. Setup analysis
source $FUNCTIONS/PhaseWY_Step0.sh;

# Preform haplotype phasing for each scaffold in parallel
export PHASE SEX_SYSTEM NONDIPLOIDS THRESHOLD WINDOW STEP PHASE_SCAFFOLD_LIST LIST_ONLY conda PROJECT \
SYNTENY VCF_IN REF OUTDIR REPEATS SAMPLE_INFO FUNCTIONS MIN_MEAN MAX_MEAN MIN_DP MISSING OUTDIR0 INDS BAMS GENOMES;

## Define function to run in parallell
phase_haplotypes() {

### Input scaffold parameters
SCAFFOLD=$1;
SCAFFOLD_LENGTH=$2;
SCAFF_NAME=$PROJECT\_$SCAFFOLD;
export SCAFFOLD SCAFFOLD_LENGTH SCAFF_NAME;

### Step 1. Get unmapped, low coverage, too high coverage, too high missingness, and repeat regions
source $FUNCTIONS/PhaseWY_Step1.sh;

if [ $TARGET_ALIGNMENT == TRUE ] && [ $PHASE == "Yes" ] && [ $HAPLOIDS_ONLY == FALSE ]; then

  if [ $PHASE_SCAFFOLD_LIST == NA ] || [ $(cat $PHASE_SCAFFOLD_LIST | grep $SCAFFOLD | wc -l) -gt 0 ]; then

    if [ $SEX_SYSTEM != NA ]; then

      ### Step 2. Calculate normalized depth difference between sexes (heterogametic sex chromosome regions have not aligned to homogametic reference)
      source $FUNCTIONS/PhaseWY_Step2.sh;
    fi;

    ### Step 3. Prepare scaffold for phasing.
    source $FUNCTIONS/PhaseWY_Step3.sh;

    if [ $HETEROZYGOTES == TRUE ]; then

      ### Step 4. Perform read based phasing with whatshap
      source $FUNCTIONS/PhaseWY_Step4.sh;

      ### Step 5. Preform statistical phasing with shapeit4
      source $FUNCTIONS/PhaseWY_Step5.sh;

      if [ $SEX_SYSTEM != NA  ]; then

        ### Step 6. Determine sex-linkage of each haplotype
        source $FUNCTIONS/PhaseWY_Step6.sh;

        ### Step 7. Split phased scaffold into autosomal, homogametic, and heterogametic genomic regions
        source $FUNCTIONS/PhaseWY_Step7.sh;

      fi;

    fi;

  fi;

fi;

### Update progress report

touch $OUTDIR/Scaffold_progress/$SCAFFOLD\_processed;

}

## Excecute function in parallell
export -f phase_haplotypes;
parallel --colsep '\t' 'phase_haplotypes {1} {2}' :::: $OUTDIR/0_file_info/SCAFFOLDS.txt;

# Step 8. Summarize the whole genome with bed files and plots
source $FUNCTIONS/PhaseWY_Step8.sh;

#### I steg 10 vill jag även skapa mer sammanfattande plottar som tex:
### Antal features some överlappar target och olika områden (gör även lista på festurs med % överlapp)
#### features och CDS i auto, och hetero vs homo regioner

# Step 9. Create new VCFs
source $FUNCTIONS/PhaseWY_Step9.sh;

# Step 10. Print some vcf stats
source $FUNCTIONS/PhaseWY_Step10.sh;

