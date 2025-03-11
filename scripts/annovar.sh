#!/bin/bash

# Run Python script to export environment variables
source <(python -m scripts.export)

vcf="$VCF_PATH/CAGI_exome_hg19.gatk.snps.vcf"
echo "Processing VCF: $vcf"
"$DOC"/annovar/table_annovar.pl "$vcf" "$DOC"/annovar/humandb/ \
    -out "$OUTPUT_PATH"/$(basename "$vcf" .vcf) \
    -vcfinput -buildver hg19 -protocol refGene -operation g -nastring . -remove