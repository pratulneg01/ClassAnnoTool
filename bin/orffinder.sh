#!/bin/bash

INPUT=$1

if [[ -z "$INPUT" || ! -f "$INPUT" ]]; then
    echo "Uso: $0 input.fasta"
    exit 1
fi

BASENAME=$(basename "$INPUT" .fasta)
OUTDIR="results/orfs"
mkdir -p "$OUTDIR"

# Ejecutar orfipy con salida en su carpeta por defecto
orfipy "$INPUT" --dna "${BASENAME}_orfs.fasta"

# Si ha funcionado y se ha generado el directorio de salida
DIR_OUTPUT="orfipy_${BASENAME}.fasta_out"

if [[ -d "$DIR_OUTPUT" ]]; then
    mv "$DIR_OUTPUT" "$OUTDIR/"
    echo "orfipy completado. Salida: $OUTDIR/$DIR_OUTPUT"
else
    echo "orfipy no generó salida. Revisa si el archivo contiene ORFs válidos."
    exit 1
fi