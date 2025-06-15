import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Argumentos de línea de comandos
parser = argparse.ArgumentParser(description="Traducir secuencias al frame más largo")
parser.add_argument("input_fasta", help="Archivo FASTA con secuencias de ADN")
parser.add_argument("output_fasta", help="Archivo FASTA de salida con proteínas traducidas")
args = parser.parse_args()

translated_records = []

for record in SeqIO.parse(args.input_fasta, "fasta"):
    seq = record.seq
    best_translation = ""
    
    for frame in range(3):
        translated = seq[frame:].translate(to_stop=False)
        if len(translated) > len(best_translation):
            best_translation = translated

    new_record = SeqRecord(
        Seq(str(best_translation)),
        id=record.id,
        description="translated_best_frame"
    )
    translated_records.append(new_record)

SeqIO.write(translated_records, args.output_fasta, "fasta")
