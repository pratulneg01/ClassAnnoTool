#!/bin/bash

INPUT="$1"
OUTPUT="$2"

if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "Uso: fastq_to_fasta.sh <input.fastq> <output.fasta>"
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT")"

seqtk seq -A "$INPUT" > "$OUTPUT"

if [[ $? -eq 0 ]]; then
    echo "Conversión completada: $OUTPUT"
else
    echo "Error en la conversión con seqtk."
    exit 1
fi
