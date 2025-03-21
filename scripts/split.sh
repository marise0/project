#!/bin/bash

# Run Python script to export environment variables
source <(python -m scripts.export)

# Path to the input VCF file
vcf_file=$VCF_PATH"/CAGI_exome_hg19.gatk.snps.vcf.gz"

# Loop through the sample names
rsids=("rs2066844" "rs2066845" "rs5743271" "rs41313262" "rs76418789" "rs104895444" "rs104895467" "rs34536443" "rs2476601" "rs11209026" "rs3764147")

rsid_filter=$(printf 'ID=="%s" || ' "${rsids[@]}")
rsid_filter=${rsid_filter% || }  # Remove trailing ' || '

for sample in $(bcftools query -l $vcf_file); do
    if [[ $sample == CD* ]]; then
        # If the sample starts with "CD", it's a case
        bcftools view -s $sample $vcf_file | \
        bcftools filter -i 'FORMAT/GT!="0/0"' | \
        bcftools view -i "$rsid_filter" -Ov -o "$VCF_PATH/cases/${sample}.vcf"
    elif [[ $sample == H* ]]; then
        # If the sample starts with "H", it's a control
        bcftools view -s $sample $vcf_file | \
        bcftools filter -i 'FORMAT/GT!="0/0"' | \
        bcftools view -i "$rsid_filter" -Ov -o "$VCF_PATH/control/${sample}.vcf"
    fi
done
