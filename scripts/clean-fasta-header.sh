#!/usr/bin/env bash
# Usage: bash clean_fasta_headers.sh input.fa output.fa

in_fasta="$1"
out_fasta="$2"

awk '
  /^>/ {
    # remove ">" and split by "|"
    split(substr($0, 2), a, "|")
    print ">" a[1]
    next
  }
  { print }
' "$in_fasta" > "$out_fasta"