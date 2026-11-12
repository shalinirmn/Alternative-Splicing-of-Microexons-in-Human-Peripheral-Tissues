#!/bin/bash

OUTDIR=${1:-"vast_align"}
SPECIES=${2:-"Hsa"}
VAST_TOOLS=${3:-"vast-tools"}

echo ">>> Combining VAST-TOOLS alignment outputs in: $OUTDIR (Species: $SPECIES)"

"$VAST_TOOLS" combine -o "$OUTDIR" -sp "$SPECIES"

echo "Combine completed. Inclusion table saved in: $OUTDIR/"
