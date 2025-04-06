#!/bin/bash

# Run Python script to export environment variables
source <(python -m scripts.export)
# Define your input and output paths
vcf_file="$VCF_PATH/CAGI_exome_hg19.gatk.snps.vcf.gz"
output_vcf="$VCF_PATH/filteredMultisample.vcf.gz"
# List of rsIDs
rsIDs=(
  "rs2241880" "rs2066844" "rs2066845" "rs2476601" "rs3764147"
  "rs5743271" "rs41313262" "rs13107325" "rs4077515" "rs492602"
  "rs138629813" "rs16844401" "rs41267765" "rs34215892" "rs61759893"
  "rs56143179" "rs104895443" "rs2228015" "rs73166641"
)

# Convert the rsIDs into a filter expression: ID=="rs1" || ID=="rs2" ...
filter_expr=$(printf 'ID=="%s" || ' "${rsIDs[@]}")
filter_expr=${filter_expr% || }  # Remove trailing " || "
# Filter the VCF
bcftools view -i "$filter_expr" "$vcf_file" -Oz -o "$output_vcf"
# Index the result
tabix -p vcf "$output_vcf"
echo "Extracted variants written to: $output_vcf"
# Create the samplesheet
SAMPLESHEET="$VCF_PATH/samplesheet.csv"
echo "sampleset,path_prefix,chrom,format" > "$SAMPLESHEET"
echo "filteredMultisampleVCF,$VCF_PATH/filteredMultisample,,vcf" >> "$SAMPLESHEET"
echo "Samplesheet created at: $SAMPLESHEET"
