#!/bin/bash

# Run Python script to export environment variables
source <(python -m scripts.export)

# Path to the input VCF file
vcf_file=$VCF_PATH"/CAGI_exome_hg19.gatk.snps.vcf.gz"

# Loop through the sample names
for sample in $(bcftools query -l $vcf_file); do
    if [[ $sample == CD* ]]; then
        # If the sample starts with "CD", it's a case
        bcftools view -s $sample -r chr16:50700000-50800000 $vcf_file | \
        bcftools filter -i 'FORMAT/GT!="0/0"' -Ov -o "$VCF_PATH/cases/${sample}.vcf"
    elif [[ $sample == H* ]]; then
        # If the sample starts with "H", it's a control
        bcftools view -s $sample -r chr16:50700000-50800000 $vcf_file | \
        bcftools filter -i 'FORMAT/GT!="0/0"' -Ov -o "$VCF_PATH/control/${sample}.vcf"
    fi
done