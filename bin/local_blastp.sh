#!/bin/bash
set -e

QUERY=$1
DB=$2
OUT=$3

blastp -query "$QUERY" \
       -db "$DB" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -evalue 100000 \
       -max_target_seqs 10000 \
       -out "$OUT"
