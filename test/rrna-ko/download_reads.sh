#!/bin/bash

SRR_LIST="./SRR_Acc_List.txt"
OUTPUT_BASE="../data/delta-te-run"

mkdir -p "$OUTPUT_BASE"

while read -r SRR_ID; do
    [[ -z "$SRR_ID" || "$SRR_ID" =~ ^# ]] && continue

    echo "Processing $SRR_ID ..."

    prefetch --output-directory "$OUTPUT_BASE" "$SRR_ID"

    SAMPLE_DIR="${OUTPUT_BASE}/${SRR_ID}"
    mkdir -p "$SAMPLE_DIR"

    fasterq-dump --split-files --outdir "$SAMPLE_DIR" --progress "${OUTPUT_BASE}/${SRR_ID}"

    echo "Finished $SRR_ID"
    echo "-------------------------------------"
done < "$SRR_LIST"

echo "All downloads completed!"
