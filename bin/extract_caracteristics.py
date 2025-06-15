#!/usr/bin/env python3
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import sys
import os

def extract_features(fasta_file, output_file):
    features = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        analysed_seq = ProteinAnalysis(seq)
        comp = analysed_seq.get_amino_acids_percent()
        mw = analysed_seq.molecular_weight()
        aromaticity = analysed_seq.aromaticity()
        instability_index = analysed_seq.instability_index()
        isoelectric_point = analysed_seq.isoelectric_point()

        aa_order = sorted(comp.keys())
        comp_list = [comp[aa] for aa in aa_order]

        feature_dict = {
            "id": record.id,
            "mw": mw,
            "aromaticity": aromaticity,
            "instability_index": instability_index,
            "isoelectric_point": isoelectric_point
        }
        feature_dict.update({f"aa_percent_{aa}": comp[aa] for aa in aa_order})
        features.append(feature_dict)

    df = pd.DataFrame(features)
    df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    output_file = sys.argv[2]
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    extract_features(fasta_file, output_file)