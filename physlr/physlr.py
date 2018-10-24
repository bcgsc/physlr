#!/usr/bin/env python3
"""
Physlr: Physical Mapping of Linked Reads
"""

import argparse
import csv
import itertools
import sys
import networkx as nx

from networkx.algorithms.shortest_paths.weighted import single_source_dijkstra_path_length
from physlr.minimerize import minimerize
from physlr.benv.graph import Graph
from physlr.read_fasta import read_fasta

class Physlr:
    """
    Physlr: Physical Mapping of Linked Reads
    """

    @staticmethod
    def read_tsv(filenames):
        "Read a graph in TSV format."
        g = nx.Graph()
        for filename in filenames:
            with open(filename) as fin:
                header = fin.readline()
                if header != "U\tV\tUn\tVn\tn\n":
                    print("Unexpected header:", header, file=sys.stderr)
                    exit(1)
                tsvin = csv.reader(fin, delimiter="\t")
                for u, v, u_minimizers, v_minimizers, shared_minimizers in tsvin:
                    g.add_edge(
                        u, v,
                        Un=int(u_minimizers), Vn=int(v_minimizers),
                        n=int(shared_minimizers))
        return g

    @staticmethod
    def read_graphviz(filename):
        "Read a GraphViz file."
        g = nx.drawing.nx_agraph.read_dot(filename)
        for _, eprop in g.edges().items():
            eprop["Un"] = int(eprop["Un"])
            eprop["Vn"] = int(eprop["Vn"])
            eprop["n"] = int(eprop["n"])
        return g

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

    def physlr_tsvtogv(self):
        "Convert a graph from TSV to GraphViz."
        g = nx.Graph()
        for filename in self.args.FASTA:
            with open(filename) as fin:
                header = fin.readline()
                if header != "U\tV\tUn\tVn\tn\n":
                    print("Unexpected header:", header, file=sys.stderr)
                    exit(1)
                tsvin = csv.reader(fin, delimiter="\t")
                for u, v, u_minimizers, v_minimizers, shared_minimizers in tsvin:
                    g.add_edge(u, v, Un=u_minimizers, Vn=v_minimizers, n=shared_minimizers)
        nx.drawing.nx_agraph.write_dot(g, sys.stdout)

    def physlr_mst(self):
        "Determine the maximum spanning tree."
        g = self.read_tsv(self.args.FASTA)
        gmst = nx.algorithms.tree.mst.maximum_spanning_tree(g, weight="n")
        nx.drawing.nx_agraph.write_dot(gmst, sys.stdout)

    def physlr_backbone(self):
        "Determine the backbone path of the graph."
        g = self.read_graphviz(self.args.FASTA[0])
        backbone = self.determine_backbone(g)
        print(*backbone)

    def physlr_graph(self, fmt):
        "Generate a graph from the minimizer index."
        graph = Graph()
        for filename in self.args.FASTA:
            with open(filename) as fin:
                graph.read_index(fin)
                graph.output_graph(pmin=0, fmt=fmt)

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
            help="A command: indexfa, indexlr, graphtsv, graphgv, overlap, tsvtogv, backbone")
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
        if self.args.command == "backbone":
            self.physlr_backbone()
        elif self.args.command == "indexfa":
            self.physlr_indexfa()
        elif self.args.command == "indexlr":
            self.physlr_indexlr()
        elif self.args.command == "graphtsv":
            self.physlr_graph("tsv")
        elif self.args.command == "graphgv":
            self.physlr_graph("graphviz")
        elif self.args.command == "mst":
            self.physlr_mst()
        elif self.args.command == "overlap":
            self.physlr_overlap()
        elif self.args.command == "tsvtogv":
            self.physlr_tsvtogv()
        else:
            print("Unrecognized command:", self.args.command, file=sys.stderr)
            exit(1)

def main():
    "Run Physlr."
    Physlr().main()

if __name__ == "__main__":
    main()
