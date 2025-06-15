import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Classify proteins as annotated or unannotated based on BLASTp results.")
    parser.add_argument("fasta", help="FASTA file with translated protein sequences")
    parser.add_argument("blast", help="BLASTp results in outfmt 6")
    parser.add_argument("--identity", type=float, default=40.0, help="Minimum identity percentage")
    parser.add_argument("--coverage", type=float, default=70.0, help="Minimum query coverage percentage")
    parser.add_argument("--evalue", type=float, default=1e-3, help="Maximum e-value allowed")
    return parser.parse_args()

def read_hits(blast_file, identity_thresh, coverage_thresh, evalue_thresh):
    hits = set()
    with open(blast_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue
            qseqid = cols[0]
            pident = float(cols[2])
            align_len = int(cols[3])
            qstart = int(cols[6])
            qend = int(cols[7])
            evalue = float(cols[10])

            qlen = abs(qend - qstart) + 1
            coverage = (qlen / align_len) * 100  # Aproximación de la cobertura

            if pident >= identity_thresh and coverage >= coverage_thresh and evalue <= evalue_thresh:
                hits.add(qseqid)
    return hits

def split_fasta_by_hits(fasta_file, hits):
    annotated = []
    unannotated = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in hits:
            annotated.append(record)
        else:
            unannotated.append(record)
    return annotated, unannotated

def main():
    args = parse_args()
    hits = read_hits(args.blast, args.identity, args.coverage, args.evalue)
    annotated, unannotated = split_fasta_by_hits(args.fasta, hits)

    SeqIO.write(annotated, "annotated_proteins.fasta", "fasta")
    SeqIO.write(unannotated, "unannotated_proteins.fasta", "fasta")
    print(f"Total sequences: {len(annotated) + len(unannotated)}")
    print(f"Annotated: {len(annotated)}")
    print(f"Unannotated: {len(unannotated)}")

if __name__ == "__main__":
    main()