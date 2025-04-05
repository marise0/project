#!/bin/bash

# Run Python script to export environment variables (like VCF_PATH)
source <(python -m scripts.export)

# Run the pgsc_calc Nextflow pipeline
nextflow run pgscatalog/pgsc_calc \
    -profile singularity \
    --input "$VCF_PATH/samplesheet.csv" \
    --scorefile "$VCF_PATH/PRS.txt" \
    --target_build GRCh37 \
    --run_ancestry pgsc_HGDP+1kGP_v1.tar.zst \
    --outdir "$OUTPUT_PATH/results"
