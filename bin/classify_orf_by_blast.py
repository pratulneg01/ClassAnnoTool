import argparse
from collections import defaultdict
from Bio import SeqIO

# Argumentos de línea de comandos
parser = argparse.ArgumentParser(description="Clasifica ORFs como anotados o no anotados según resultados BLAST")
parser.add_argument("blast_file", help="Archivo BLAST en formato outfmt 6")
parser.add_argument("query_fasta", help="Archivo FASTA original con los ORFs")
parser.add_argument("--identity", type=float, default=80.0, help="Umbral de identidad mínima (porcentaje)")
parser.add_argument("--coverage", type=float, default=70.0, help="Umbral de cobertura mínima (longitud de alineamiento en nt)")
parser.add_argument("--evalue", type=float, default=1e-5, help="Valor E máximo permitido")

args = parser.parse_args()

# Guardar hits que cumplen criterios
blast_hits = defaultdict(list)

with open(args.blast_file) as bf:
    for line in bf:
        cols = line.strip().split("\t")
        if len(cols) < 11:
            continue  # línea malformada
        qid = cols[0]
        identity = float(cols[2])
        aln_len = int(cols[3])
        evalue = float(cols[10])

        if identity >= args.identity and aln_len >= args.coverage and evalue <= args.evalue:
            blast_hits[qid].append((identity, aln_len, evalue))

# Clasificar IDs
all_records = list(SeqIO.parse(args.query_fasta, "fasta"))
all_ids = set(rec.id for rec in all_records)
annotated_ids = set(blast_hits.keys())
unannotated_ids = all_ids - annotated_ids

# Guardar listas de IDs
with open("annotated_ids.txt", "w") as out1:
    for aid in sorted(annotated_ids):
        out1.write(f"{aid}\n")

with open("unannotated_ids.txt", "w") as out2:
    for uid in sorted(unannotated_ids):
        out2.write(f"{uid}\n")

# Guardar FASTA con secuencias separadas
with open("annotated.fasta", "w") as aout, open("unannotated.fasta", "w") as uout:
    for rec in all_records:
        if rec.id in annotated_ids:
            SeqIO.write(rec, aout, "fasta")
        else:
            SeqIO.write(rec, uout, "fasta")

# Informe final
print("Clasificación completada.")
print(f"ORFs anotados:     {len(annotated_ids)} → annotated.fasta")
print(f"ORFs no anotados:  {len(unannotated_ids)} → unannotated.fasta")
