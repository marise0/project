#!/bin/bash

# Run Python script to export environment variables
source <(python -m scripts.export)
vcf_file=$VCF_PATH"/CAGI_exome_hg19.gatk.snps.vcf.gz"
rsIDs=("rs2241880" "rs2066844" "rs2066845" "rs2476601" "rs3764147" "rs5743271" "rs41313262" "rs13107325" "rs4077515" "rs492602" "rs138629813" "rs16844401" "rs41267765" "rs34215892" "rs61759893" "rs56143179" "rs104895443" "rs2228015" "rs73166641")
rsid_filter=$(printf 'ID=="%s" || ' "${rsIDs[@]}")
rsid_filter=${rsid_filter% || }  

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
