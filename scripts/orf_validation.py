import os
from collections import defaultdict

BLAST_FILE = "output/blast_results.tsv"
GFF_FILE = "output/orfs.gff"
UPDATED_GFF_FILE = "output/orfs_validated.gff"
FASTA_FILE = "output/orfs.fasta"

EVALUE_THRESHOLD = 1e-5


def parse_blast_results(blast_path):
    validated_orfs = {}
    with open(blast_path) as blast:
        for line in blast:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = cols[:12]

            evalue = float(evalue)
            if evalue <= EVALUE_THRESHOLD:
                qstart, qend = int(qstart), int(qend)
                validated_orfs[qseqid] = {
                    "hit": sseqid,
                    "evalue": evalue,
                    "qstart": min(qstart, qend),
                    "qend": max(qstart, qend)
                }
    return validated_orfs


def annotate_gff(gff_path, validated_orfs, output_path):
    new_entries = []
    all_orf_ids = set()

    with open(gff_path) as gff:
        original_lines = gff.readlines()

    for line in original_lines:
        if line.startswith("#") or not line.strip():
            continue

        cols = line.strip().split("\t")
        if len(cols) < 9:
            continue

        source = cols[0]
        feature_type = cols[2]
        start = int(cols[3])
        end = int(cols[4])
        attributes = cols[8]

        if "ID=" in attributes:
            orf_id = attributes.split("ID=")[-1].split(";")[0]
            all_orf_ids.add(orf_id)

            if orf_id in validated_orfs:
                hit_info = validated_orfs[orf_id]
                cds_start = start + (hit_info["qstart"] - 1) * 3
                cds_end = start + (hit_info["qend"] - 1) * 3
                cds_entry = [source, "BLAST", "CDS", str(cds_start), str(cds_end), ".", "+", "0",
                             f"ID=CDS_{orf_id};Parent={orf_id};hit={hit_info['hit']}"]
                new_entries.append("\t".join(cds_entry))

                if cds_start > start:
                    utr5_entry = [source, "BLAST", "five_prime_UTR", str(start), str(cds_start - 1), ".", "+", ".", f"Parent={orf_id}"]
                    new_entries.append("\t".join(utr5_entry))
                if cds_end < end:
                    utr3_entry = [source, "BLAST", "three_prime_UTR", str(cds_end + 1), str(end), ".", "+", ".", f"Parent={orf_id}"]
                    new_entries.append("\t".join(utr3_entry))

    with open(output_path, "w") as out:
        out.writelines(original_lines)
        out.write("\n".join(new_entries))
        out.write("\n")

    return all_orf_ids


def assess_false_positives(validated_orfs, all_orf_ids):
    validated = set(validated_orfs.keys())
    total = len(all_orf_ids)
    true_positives = len(validated)
    false_positives = total - true_positives
    fpr = false_positives / total if total else 0
    print(f"\nORF Validation Summary:")
    print(f"Total ORFs        : {total}")
    print(f"Validated (CDS)   : {true_positives}")
    print(f"False Positives   : {false_positives}")
    print(f"False Positive Rate: {fpr:.2%}")
    return fpr


def main():
    print("Parsing BLAST results")
    validated_orfs = parse_blast_results(BLAST_FILE)

    print("Annotating GFF file")
    all_orfs = annotate_gff(GFF_FILE, validated_orfs, UPDATED_GFF_FILE)

    print("Assessing false positive rate")
    assess_false_positives(validated_orfs, all_orfs)

    print(f"\nGFF saved to: {UPDATED_GFF_FILE}")


if __name__ == "__main__":
    main()
