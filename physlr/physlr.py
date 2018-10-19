#!/usr/bin/env python3
"""
Physlr: Physical Mapping of Linked Reads
"""

import argparse
import sys

from physlr.minimerize import minimerize
from physlr.read_fasta import read_fasta

class Physlr:
    """
    Physlr: Physical Mapping of Linked Reads
    """

    def physlr_indexfa(self):
        "Index a set of sequences. The output file format is JSON."
        print("{")
        for filename in self.args.FASTA:
            with open(filename) as fin:
                for name, seq, _ in read_fasta(fin):
                    print(
                        '"', name, '": ',
                        minimerize(self.args.k, self.args.w, seq.upper()),
                        sep="")
        print("}")

    def physlr_indexlr(self):
        "Index a set of linked reads. The output file format is JSON."
        print("{")
        for filename in self.args.FASTA:
            with open(filename) as fin:
                for name, seq, bx in read_fasta(fin):
                    print(
                        '"', bx, '": ',
                        minimerize(self.args.k, self.args.w, seq.upper()),
                        sep="")
        print("}")

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
            help="A command: indexfa, indexlr")
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
        else:
            print("Unrecognized command:", self.args.command, file=sys.stderr)
            exit(1)

def main():
    "Run Physlr."
    Physlr().main()

if __name__ == "__main__":
    main()
