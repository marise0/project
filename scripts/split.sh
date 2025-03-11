#!/bin/bash

# Run Python script to export environment variables
source <(python -m scripts.export)

# Define VCF file path
vcf="$VCF_PATH/CAGI_exome_hg19.gatk.snps.vcf"

# Extract sample names and loop through them
gatk --java-options "-Xmx4g" VariantsToTable \
   -V data/CAGI_exome_hg19.gatk.snps.vcf \
   -F SAMPLE \
   -O sample_list.txt