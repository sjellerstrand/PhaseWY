#! /usr/bin/env python

# Version 2023-02-28
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

from sys import *
import os, time, argparse, re
import gzip

parser = argparse.ArgumentParser(description='Converts haploid genotypes into diploid homozygous genotypes')
parser.add_argument('-i', '--input', dest='i', help="input file in freebayes vcf format [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="output file [required]", required=True)
args = parser.parse_args()



# Check if the input vcf file is gzipped
if args.i.endswith(".gz"):
    inputF=gzip.open(args.i,'rt')

else:
    inputF=open(args.i,'r')

outputF=open(args.o, 'w')

for Line in inputF:
        # Check if the header section
        if re.match('^#',Line) is None:

                # Get the columns of that line
                columns=Line.strip("\n").split("\t")

                # Add the info to the site
                result=columns[0:9]

                # Get the genotypes
                genotypecolumns=range(9,len(columns))

                # Check each individual if it is haploid
                for ind in genotypecolumns:
                        genotype=columns[ind]
                        genotype=genotype.split(":")

                        # Check genotype field
                        if len(str(genotype[0]))==1:
                            genotype[0]=genotype[0]+"|"+genotype[0]
                            result.append(":".join(genotype))
                            
                        # If diploid genotype
                        else:
                                result.append(":".join(genotype))

                outputF.write('\t'.join(result)+"\n")

        # If it is a header line, just write it out
        else:
                outputF.write(Line)

inputF.close()
outputF.close()

