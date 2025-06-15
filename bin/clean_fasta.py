#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import os

def is_valid_sequence(seq, allowed="ACGTN"):
    return re.fullmatch(f"[{allowed}]+", seq.upper()) is not None

def clean_fasta(input_file, output_file, min_length=100):
    unique_seqs = set()
    total = 0
    written = 0

    cleaned_records = []

    # Crear carpeta log si no existe
    os.makedirs("log", exist_ok=True)

    # Archivos log
    too_short_log = open("log/removed_too_short.txt", "w")
    invalid_chars_log = open("log/removed_invalid_chars.txt", "w")
    duplicates_log = open("log/removed_duplicates.txt", "w")

    for record in SeqIO.parse(input_file, "fasta"):
        total += 1
        seq_str = str(record.seq).upper()

        if len(seq_str) < min_length:
            too_short_log.write(f"{record.id} | length: {len(seq_str)}\n")
            continue

        if not is_valid_sequence(seq_str):
            invalid_chars = set(re.sub(r'[ACGTN]', '', seq_str))
            invalid_chars_log.write(f"{record.id} | invalid characters: {', '.join(invalid_chars)}\n")
            continue

        if seq_str in unique_seqs:
            duplicates_log.write(f"{record.id} | duplicate sequence\n")
            continue

        unique_seqs.add(seq_str)
        cleaned_record = SeqRecord(record.seq.upper(), id=record.id, description="")
        cleaned_records.append(cleaned_record)
        written += 1

    # Cierre de logs
    too_short_log.close()
    invalid_chars_log.close()
    duplicates_log.close()

    SeqIO.write(cleaned_records, output_file, "fasta")

    print(f"Total sequences processed: {total}")
    print(f"Sequences written to output: {written}")
    print(f"Sequences removed: {total - written}")
    print("Logs written to:")
    print("  log/removed_too_short.txt")
    print("  log/removed_invalid_chars.txt")
    print("  log/removed_duplicates.txt")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean a FASTA file by filtering invalid or low-quality sequences.")
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("output", help="Output cleaned FASTA file")
    parser.add_argument("--min_length", type=int, default=100, help="Minimum sequence length (default: 100)")

    args = parser.parse_args()
    clean_fasta(args.input, args.output, args.min_length)
