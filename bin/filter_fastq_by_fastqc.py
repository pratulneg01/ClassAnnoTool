#!/usr/bin/env python3

import zipfile
import sys
import statistics
import os

def parse_quality_threshold(zip_path):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        txt_name = [name for name in zip_ref.namelist() if name.endswith('fastqc_data.txt')]
        if not txt_name:
            raise ValueError("No se encontró 'fastqc_data.txt' dentro del ZIP de FastQC.")

        with zip_ref.open(txt_name[0]) as f:
            lines = [line.decode('utf-8').strip() for line in f]

    # Buscar sección '>>Per base sequence quality'
    start, end = None, None
    for i, line in enumerate(lines):
        if line.startswith(">>Per base sequence quality"):
            start = i + 1
        elif start and line.startswith(">>END_MODULE"):
            end = i
            break

    if start is None or end is None:
        raise ValueError("No se encontró la sección 'Per base sequence quality' en el informe.")

    qual_scores = []
    for line in lines[start:end]:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split()
        try:
            qual = float(parts[0])
            count = int(float(parts[1]))
            qual_scores.extend([qual] * count)
        except (ValueError, IndexError):
            continue  # saltar líneas mal formateadas

    if not qual_scores:
        raise ValueError("No se pudo calcular el umbral: sin datos de calidad válidos.")

    return statistics.mean(qual_scores)

def filter_fastq(fastq_file, output_file, quality_threshold):
    def phred_score(char):
        return ord(char) - 33

    total = 0
    passed = 0

    with open(fastq_file, 'r') as f_in, open(output_file, 'w') as f_out:
        while True:
            header = f_in.readline()
            seq = f_in.readline()
            plus = f_in.readline()
            qual = f_in.readline()
            if not qual:
                break  # EOF

            total += 1

            if not seq or not qual.strip():
                continue  # saltar secuencias vacías

            scores = [phred_score(q) for q in qual.strip()]
            if not scores:
                continue  # evitar división por cero

            avg_score = sum(scores) / len(scores)

            if avg_score >= quality_threshold:
                passed += 1
                f_out.write(header + seq + plus + qual)

    print(f"Umbral de calidad promedio detectado: {quality_threshold:.2f}")
    print(f"Secuencias totales: {total}")
    print(f"Secuencias eliminadas: {total - passed}")
    print(f"Secuencias conservadas: {passed}")
    print(f"FASTQ filtrado guardado en: {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Uso: python3 filter_fastq_by_fastqc.py <archivo_fastqc.zip> <input.fastq> <output.fastq>")
        sys.exit(1)

    zip_file = sys.argv[1]
    fastq_file = sys.argv[2]
    output_file = sys.argv[3]

    try:
        threshold = parse_quality_threshold(zip_file)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        filter_fastq(fastq_file, output_file, threshold)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
