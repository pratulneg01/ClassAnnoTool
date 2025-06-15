#!/bin/bash

# Uso: ./blastn.sh input.fasta output.txt

INPUT=$1
OUTPUT=$2

if [[ -z "$INPUT" || -z "$OUTPUT" || ! -f "$INPUT" ]]; then
    echo "Uso: $0 input.fasta output.txt"
    echo "Error: archivo de entrada no encontrado o parámetros incompletos."
    exit 1
fi

# Ejecutar BLASTn contra la base de datos remota 'nt'
blastn -query "$INPUT" -db nt -remote -out "$OUTPUT" -outfmt 6 -evalue 1e-5