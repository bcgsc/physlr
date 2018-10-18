#!/usr/bin/env pypy3

import argparse
import sys

from physlr.minimerize import minimerize
from physlr.read_fasta import read_fasta

class Physlr:
    """
    Physlr: Physical Mapping of Linked Reads
    """

    def physlr_index(self):
        "Index a set of sequences. The output file format is JSON."
        print("{")
        for filename in self.args.FASTA:
            with open(filename) as fin:
                for name, seq in read_fasta(fin):
                    print('"', name, '": ',
                        minimerize(self.args.k, self.args.w, seq.upper()),
                        sep ="")
        print("}")

    def parse_arguments(self):
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
            help="A command: assemble, graph, index")
        argparser.add_argument(
            "FASTA", nargs="+",
            help="FASTA/FASTQ file of linked reads")
        return argparser.parse_args()

    def main(self):
        "Process each file specified on the command line"
        self.args = self.parse_arguments()
        self.args.FASTA = ["/dev/stdin" if s == "-" else s for s in self.args.FASTA]
        if self.args.command == "index":
            self.physlr_index()
        else:
            print("Unrecognized command: ", self.args.command, file=sys.stderr)
            exit(1)

def main():
    Physlr().main()

if __name__ == "__main__":
    main()
