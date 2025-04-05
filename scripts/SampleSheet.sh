#!/bin/bash

# Run Python script to export environment variables
source <(python -m scripts.export)
SAMPLESHEET="$VCF_PATH/samplesheet.csv"
echo "sampleset,path_prefix,chrom,format" > "$SAMPLESHEET"

# Process VCF files in the cases directory
for vcf in "$VCF_PATH"/cases/*.vcf; do
  echo "Processing VCF: $vcf"
  case_name=$(basename "$vcf" .vcf)
  output_vcf="$VCF_PATH/cases/$case_name.vcf.gz"
  bgzip -f "$vcf" -c > "$output_vcf"
  tabix -p vcf "$output_vcf"
  echo "$case_name,${VCF_PATH}/cases/${case_name},,vcf" >> "$SAMPLESHEET"
done

# Process VCF files in the control directory
for vcf in "$VCF_PATH"/control/*.vcf; do
  echo "Processing VCF: $vcf"
  control_name=$(basename "$vcf" .vcf)
  output_vcf="$VCF_PATH/control/$control_name.vcf.gz"
  bgzip -f "$vcf" -c > "$output_vcf"
  tabix -p vcf "$output_vcf"
  echo "$control_name,${VCF_PATH}/control/${control_name},,vcf" >> "$SAMPLESHEET"
done
echo "Samplesheet created at: $SAMPLESHEET"
