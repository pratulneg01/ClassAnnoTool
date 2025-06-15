#!/bin/bash

INPUT="$1"
ADAPTER="${2:-AGATCGGAAGAGC}"
OUTDIR="results/cutadapt"

if [[ -z "$INPUT" || ! -f "$INPUT" ]]; then
    echo "Error: archivo no especificado o no encontrado."
    exit 1
fi

mkdir -p "$OUTDIR"

# Detectar formato según la extensión del archivo o contenido
EXT=""
if [[ "$INPUT" == *.fastq || "$INPUT" == *.fq ]] || grep -q '^@' "$INPUT"; then
    EXT="fastq"
elif [[ "$INPUT" == *.fasta || "$INPUT" == *.fa ]] || grep -q '^>' "$INPUT"; then
    EXT="fasta"
else
    echo "Error: no se pudo determinar el formato del archivo (ni FASTQ ni FASTA)."
    exit 1
fi

BASENAME=$(basename "$INPUT" | sed 's/\.[^.]*$//')

# Ejecutar Cutadapt con formato de salida correspondiente
cutadapt -a "$ADAPTER" -o "$OUTDIR/${BASENAME}_cutadapt.${EXT}" "$INPUT"

if [[ $? -eq 0 ]]; then
    echo "Cutadapt finalizado con éxito. Archivo de salida: $OUTDIR/${BASENAME}_cutadapt.${EXT}"
else
    echo "Error al ejecutar Cutadapt."
    exit 1
fi
