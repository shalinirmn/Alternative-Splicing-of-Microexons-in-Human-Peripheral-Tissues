#!/bin/bash

READ_DIR=$1
OUTDIR=${2:-"vast_align"}
VAST_TOOLS=${3:-"vast-tools"}
SPECIES="Hs2"
CORES=12

mkdir -p "$OUTDIR"

echo ">>> Starting VAST-TOOLS alignment for all _trim.fq.gz pairs..."
echo ">>> Using input directory: $READ_DIR"

# Loop through all *_1_trim.fq.gz files
for R1 in ${READ_DIR}/*_1_trim.fq.gz; do
    R2="${R1/_1_trim.fq.gz/_2_trim.fq.gz}"   
    SAMPLE=$(basename "$R1" _1_trim.fq.gz)

    # check if R2 exists
    if [[ ! -f "$R2" ]]; then
        echo "Warning: Missing pair for $R1 — skipping this sample!"
        continue
    fi

    echo ">>> Running VAST-TOOLS align for sample: $SAMPLE"
    echo "R1: $R1"
    echo "R2: $R2"

    # Run VAST-TOOLS alignment
    "$VAST_TOOLS" align -sp "$SPECIES" -o "$OUTDIR" "$R1" "$R2" --cores "$CORES"

    echo "Completed: $SAMPLE"
done

echo "All alignments completed. Output saved to: $OUTDIR/"
