import matplotlib.pyplot as plt
from Bio import SeqIO
import math
import os

FASTA_FILE = "output/orfs.fasta"
BLAST_FILE = "output/blast_results.tsv"

def get_validated_orfs(blast_path):
    validated_ids = set()
    with open(blast_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            qseqid = line.strip().split("\t")[0]
            validated_ids.add(qseqid)
    return validated_ids

def load_orf_lengths(fasta_path, validated_ids):
    true_lengths = []
    false_lengths = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        orf_id = record.id
        orf_len = len(record.seq)
        if orf_id in validated_ids:
            true_lengths.append(orf_len)
        else:
            false_lengths.append(orf_len)
    return true_lengths, false_lengths

def plot_histograms(true_lengths, false_lengths):
    plt.figure(figsize=(10, 6))
    bins = range(0, max(true_lengths + false_lengths) + 100, 60)

    plt.hist(false_lengths, bins=bins, color='red', alpha=0.6, label='False Positives')
    plt.hist(true_lengths, bins=bins, color='green', alpha=0.7, label='True Positives')

    plt.title("ORF length distribution")
    plt.xlabel("Length (nucleotides)")
    plt.ylabel("ORF Amount")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("orf_length_distribution.png")
    plt.show()

def estimate_geometric_expectation():
    p = 3 / 64
    expected_codons = 1 / p
    expected_length_nt = expected_codons * 3
    print(f"  - p (codon stop) = {p:.4f}")
    print(f"  - E[X] = {expected_codons:.2f} codons")
    print(f"  - Expected mean length ≈ {expected_length_nt:.0f} nucleotides")

def main():
    print("Loading validated ORFs")
    validated_ids = get_validated_orfs(BLAST_FILE)

    print("Loading ORF sequences")
    true_lengths, false_lengths = load_orf_lengths(FASTA_FILE, validated_ids)

    print("Building histogram (i hope it looks nice (it doesn't))")
    plot_histograms(true_lengths, false_lengths)

    estimate_geometric_expectation()
    
    print(f"Valided ORFs (true positives) : {len(true_lengths)}")
    print(f"Non-Validated ORFs (false positives) : {len(false_lengths)}")
    print(f"Total amount of analyzed ORFs : {len(true_lengths) + len(false_lengths)}")


if __name__ == "__main__":
    main()

