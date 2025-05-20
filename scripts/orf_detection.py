from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data.CodonTable import standard_dna_table
import os

class ORF:
    def __init__(self, seq, src, start, end, frame):
        self.seq = seq
        self.src = src
        self.start = start
        self.end = end
        self.frame = frame
        self.length = end - start + 1
        self.id = None
        self.alt_starts = []
        self.nested_orfs = []

    def add_alt_starts(self, positions):
        if isinstance(positions, int):
            positions = [positions]
        self.alt_starts.extend(positions)

    def add_nested(self, orf):
        self.nested_orfs.append(orf)

    def to_gff(self):
        return f"{self.src}\tcustom_ORF_finder\tORF\t{self.start}\t{self.end}\t.\t{self.frame[0]}\t{int(self.frame[-1]) - 1}\tID={self.id}"


class ORFExtractor:
    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.orfs = []

    def _find_codons(self, seq, codons, strand):
        positions = {}
        if isinstance(codons, str): codons = [codons]

        for i in range(3):
            frame_key = f"{strand}{i + 1}"
            frame_seq = seq[i:]
            pos_list = [i + pos * 3 for pos in range(len(frame_seq) // 3) if frame_seq[pos * 3:pos * 3 + 3] in codons]
            if pos_list:
                positions[frame_key] = pos_list
        return positions

    def _extract_orfs_from_frame(self, starts, stops, frame, seq, rev_seq):
        while starts and stops:
            alt_starts = []
            start = starts.pop(0)

            while stops and stops[0] < start:
                stops.pop(0)

            if stops:
                end = stops[0] + 2
                while starts and starts[0] < end:
                    alt_starts.append(starts.pop(0))

                if frame.startswith("+"):
                    orf_seq = seq.seq[start:end + 1]
                    start_adj, end_adj = start + 1, end + 1
                else:
                    orf_seq = rev_seq.seq[start:end + 3]
                    start_adj, end_adj = len(seq) - end, len(seq) - start

                orf = ORF(orf_seq, seq.id, start_adj, end_adj, frame)
                orf.add_alt_starts([start_adj] + [s + 1 for s in alt_starts])
                self.orfs.append(orf)

    def extract(self):
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            rev_record = record.reverse_complement()
            starts = self._find_codons(record.seq, "ATG", "+")
            starts.update(self._find_codons(rev_record.seq, "ATG", "-"))

            stops = self._find_codons(record.seq, standard_dna_table.stop_codons, "+")
            stops.update(self._find_codons(rev_record.seq, standard_dna_table.stop_codons, "-"))

            common_frames = set(starts).intersection(stops)

            for frame in common_frames:
                self._extract_orfs_from_frame(starts[frame], stops[frame], frame, record, rev_record)

        for i, orf in enumerate(self.orfs, 1):
            orf.id = f"ORF_{i}"

    def export(self, output_prefix):
        with open(output_prefix + ".gff", "w") as gff:
            for orf in self.orfs:
                gff.write(orf.to_gff() + "\n")

        with open(output_prefix + ".fasta", "w") as fasta:
            for orf in self.orfs:
                record = SeqRecord(orf.seq, id=orf.id,
                                   description=f"length={orf.length} source={orf.src} start={orf.start} end={orf.end} frame={orf.frame}")
                SeqIO.write(record, fasta, "fasta")


fasta_input = "data/Homo_sapiens_cdna_assembed.fasta"
output_name = "output/orfs"

extractor = ORFExtractor(fasta_input)
extractor.extract()
extractor.export(output_name)

print(f"ORFs found : {len(extractor.orfs)}")
print(f"Saved files as : {output_name}.fasta and {output_name}.gff")
