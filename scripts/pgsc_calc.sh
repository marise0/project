#!/bin/bash

# Run Python script to export environment variables
source <(python -m scripts.export)

# Run the pgsc_calc Nextflow pipeline
nextflow run pgscatalog/pgsc_calc \
    -profile singularity \
    --input "$VCF_PATH/samplesheet.csv" \
    --scorefile "$VCF_PATH/PRS.txt" \
    --target_build GRCh37 \
    --outdir "$OUTPUT_PATH/results"  
    