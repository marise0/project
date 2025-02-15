#!/bin/bash

# Run Python script to export environment variables
source <(python -m scripts.export)

# Loop through all VCF files in the specified directory
for vcf in "$VCF_PATH"/cases/*.vcf; do
    echo "Processing VCF: $vcf"
    # Run the table_annovar.pl command for each VCF file
   "$DOC"/annovar/table_annovar.pl "$vcf" "$DOC"/annovar/humandb/ -out "$OUTPUT_PATH"/cases/$(basename "$vcf" .vcf) -vcfinput -buildver hg19 -protocol refGene -operation g -nastring . -remove
done

for vcf in "$VCF_PATH"/control/*.vcf; do
    echo "Processing VCF: $vcf"
    # Run the table_annovar.pl command for each VCF file
   "$DOC"/annovar/table_annovar.pl "$vcf" "$DOC"/annovar/humandb/ -out "$OUTPUT_PATH"/control/$(basename "$vcf" .vcf) -vcfinput -buildver hg19 -protocol refGene -operation g -nastring . -remove
done