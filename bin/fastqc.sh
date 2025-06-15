#!/bin/bash

INPUT=$1
OUTDIR="results/fastqc"

if [[ -z "$INPUT" || ! -f "$INPUT" ]]; then
    echo "Error: archivo no especificado o no encontrado."
    exit 1
fi

# Crear carpeta de salida si no existe
mkdir -p "$OUTDIR"

# Ejecutar FastQC
fastqc "$INPUT" --outdir="$OUTDIR"
if [[ $? -ne 0 ]]; then
    echo "Error en FastQC."
    exit 1
fi

# Extraer nombre base del archivo original sin extensión
ORIG_BASE=$(basename "$INPUT" | sed 's/\.[^.]*$//')
BASENAME="fastqc_${ORIG_BASE}"

# Renombrar archivos generados por FastQC
mv "$OUTDIR/${ORIG_BASE}_fastqc.html" "$OUTDIR/${BASENAME}.html"
mv "$OUTDIR/${ORIG_BASE}_fastqc.zip" "$OUTDIR/${BASENAME}.zip"

echo "FastQC completado correctamente. Resultados guardados como:"
echo " - $OUTDIR/${BASENAME}.html"
echo " - $OUTDIR/${BASENAME}.zip"
