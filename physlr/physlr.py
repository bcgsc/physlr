#!/usr/bin/env python3
"""
Physlr: Physical Mapping of Linked Reads
"""

import argparse
import itertools
import sys

from physlr.minimerize import minimerize
from physlr.read_fasta import read_fasta

class Physlr:
    """
    Physlr: Physical Mapping of Linked Reads
    """

    def physlr_indexfa(self):
        "Index a set of sequences. The output file format is TSV."
        for filename in self.args.FASTA:
            with open(filename) as fin:
                for name, seq, _ in read_fasta(fin):
                    print(name, "\t", sep="", end="")
                    print(*minimerize(self.args.k, self.args.w, seq.upper()))

    def physlr_indexlr(self):
        "Index a set of linked reads. The output file format is TSV."
        for filename in self.args.FASTA:
            with open(filename) as fin:
                for _, seq, bx in read_fasta(fin):
                    print(bx, "\t", sep="", end="")
                    print(*minimerize(self.args.k, self.args.w, seq.upper()))

    def physlr_overlap(self):
        "Read a sketch of linked reads and find overlapping barcodes."

        # Read a dictionary of barcodes to minimizers.
        bxtomin = {}
        for filename in self.args.FASTA:
            with open(filename) as fin:
                for line in fin:
                    fields = line.split(None, 1)
                    if len(fields) < 2:
                        continue
                    bx = fields[0]
                    minimizers = fields[1].split()
                    if bx not in bxtomin:
                        bxtomin[bx] = set()
                    bxtomin[bx].update(minimizers)

        # Construct a dictionary of minimizers to barcodes.
        mintobx = {}
        for bx, minimizers in bxtomin.items():
            for x in minimizers:
                if x not in mintobx:
                    mintobx[x] = set()
                mintobx[x].add(bx)

        # Construct a set of barcode pairs that share a minimizer.
        bxpairs = set()
        for bxs in mintobx.values():
            for u, v in itertools.combinations(bxs, 2):
                bxpairs.add((min(u, v), max(u, v)))

        # Output overlapping barcodes.
        print("U\tV\tUn\tVn\tn")
        for u, v in bxpairs:
            print(u, v, len(bxtomin[u]), len(bxtomin[v]), len(bxtomin[u] & bxtomin[v]), sep="\t")

    @staticmethod
    def parse_arguments():
        "Parse the command line arguments."
        argparser = argparse.ArgumentParser()
        argparser.add_argument(
            "-k", "--k", action="store", type=int, required=True,
            help="size of a k-mer (bp)")
        argparser.add_argument(
            "-w", "--window", action="store", dest="w", type=int, required=True,
            help="number of k-mers in a window of size k + w - 1 bp")
        argparser.add_argument(
            "command",
            help="A command: indexfa, indexlr, overlap")
        argparser.add_argument(
            "FASTA", nargs="+",
            help="FASTA/FASTQ file of linked reads")
        return argparser.parse_args()

    def __init__(self):
        "Create a new instance of Physlr."
        self.args = self.parse_arguments()
        self.args.FASTA = ["/dev/stdin" if s == "-" else s for s in self.args.FASTA]

    def main(self):
        "Run Physlr."
        if self.args.command == "indexfa":
            self.physlr_indexfa()
        elif self.args.command == "indexlr":
            self.physlr_indexlr()
        elif self.args.command == "overlap":
            self.physlr_overlap()
        else:
            print("Unrecognized command:", self.args.command, file=sys.stderr)
            exit(1)

def main():
    "Run Physlr."
    Physlr().main()

if __name__ == "__main__":
    main()
