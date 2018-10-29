#!/usr/bin/env python3
"""
Physlr: Physical Mapping of Linked Reads
"""

import argparse
import itertools
import sys
import networkx as nx

from networkx.algorithms.shortest_paths.weighted import single_source_dijkstra_path_length
from physlr.minimerize import minimerize
from physlr.benv.graph import Graph
from physlr.read_fasta import read_fasta

def quantile(quantiles, xs):
    "Return the specified quantiles p of xs."
    sorted_xs = sorted(xs)
    return [sorted_xs[round(p * (len(sorted_xs)-1))] for p in quantiles]

class Physlr:
    """
    Physlr: Physical Mapping of Linked Reads
    """

    @staticmethod
    def write_tsv(g, fout):
        "Write a graph in TSV format."
        if "m" in next(iter(g.nodes.values())):
            print("U\tn\tm", file=fout)
        else:
            print("U\tn", file=fout)
        for u, prop in g.nodes.items():
            if "m" in prop:
                print(u, prop["n"], prop["m"], sep="\t", file=fout)
            else:
                print(u, prop["n"], sep="\t", file=fout)
        print("\nU\tV\tn", file=fout)
        for e, prop in g.edges.items():
            print(e[0], e[1], prop["n"], sep="\t", file=fout)

    @staticmethod
    def write_graph(g, fout, graph_format):
        "Write a graph."
        if graph_format == "gv":
            nx.drawing.nx_agraph.write_dot(g, sys.stdout)
        elif graph_format == "tsv":
            Physlr.write_tsv(g, fout)
        else:
            print("Unknown graph format:", graph_format, file=sys.stderr)
            exit(1)

    @staticmethod
    def read_tsv(g, filename):
        "Read a graph in TSV format."
        with open(filename) as fin:
            line = fin.readline()
            if line not in ["U\tn\n", "U\tn\tm\n"]:
                print("Unexpected header:", line, file=sys.stderr)
                exit(1)
            reading_vertices = True
            for line in fin:
                if line == "\n":
                    line = fin.readline()
                    if line == "U\tV\tn\n":
                        reading_vertices = False
                    else:
                        print("Unexpected header:", line, file=sys.stderr)
                        exit(1)
                    line = fin.readline()
                xs = line.split()
                if reading_vertices:
                    if len(xs) == 2:
                        g.add_node(xs[0], n=int(xs[1]))
                    elif len(xs) == 3:
                        g.add_node(xs[0], n=int(xs[1]), m=int(xs[2]))
                    else:
                        print("Unexpected row:", line, file=sys.stderr)
                        exit(1)
                else:
                    if len(xs) == 3:
                        g.add_edge(xs[0], xs[1], n=int(xs[2]))
                    else:
                        print("Unexpected row:", line, file=sys.stderr)
                        exit(1)
        return g

    @staticmethod
    def read_graphviz(g, filename):
        "Read a GraphViz file."
        graph = nx.drawing.nx_agraph.read_dot(filename)
        for vprop in graph.nodes().values():
            vprop["n"] = int(vprop["n"])
        for eprop in graph.edges().values():
            eprop["n"] = int(eprop["n"])
        return nx.algorithms.operators.binary.compose(g, graph)

    @staticmethod
    def sort_vertices(g):
        """
        Sort the vertices of a graph by name.
        There may be more than one tree with the same minimum or maximum weight.
        Which spanning tree is chosen depends on the order of the vertices.
        Sort the vertices of a graph by name to ensure consistent results.
        """
        gsorted = nx.Graph()
        gsorted.add_nodes_from(sorted(g.nodes().items()))
        gsorted.add_edges_from((e[0], e[1], eprops) for e, eprops in g.edges().items())
        return gsorted

    @staticmethod
    def read_graph(filenames):
        "Read a graph in either GraphViz or TSV format."
        g = nx.Graph()
        for filename in filenames:
            with open(filename) as fin:
                c = fin.read(1)
                if c == "s":
                    g = Physlr.read_graphviz(g, filename)
                elif c == "U":
                    g = Physlr.read_tsv(g, filename)
                else:
                    print("Unexpected graph format", c + fin.readline(), file=sys.stderr)
                    sys.exit(1)
        return Physlr.sort_vertices(g)

    @staticmethod
    def read_minimizers(filenames):
        "Read minimizers in TSV format."
        bxtomin = {}
        for filename in filenames:
            with open(filename) as fin:
                for line in fin:
                    fields = line.split(None, 1)
                    if len(fields) < 2:
                        continue
                    bx = fields[0]
                    if bx not in bxtomin:
                        bxtomin[bx] = set()
                    bxtomin[bx].update(int(x) for x in fields[1].split())
        return bxtomin

    @staticmethod
    def construct_minimizers_to_barcodes(bxtomin):
        "Construct a dictionary of minimizers to barcodes."
        mintobx = {}
        for bx, minimizers in bxtomin.items():
            for x in minimizers:
                if x not in mintobx:
                    mintobx[x] = set()
                mintobx[x].add(bx)
        return mintobx

    @staticmethod
    def determine_backbone(g):
        "Determine the backbone of the maximum spanning tree."
        g.remove_nodes_from([u for u, deg in g.degree if deg == 0])
        ecc = nx.algorithms.distance_measures.eccentricity(g)
        diameter = nx.algorithms.distance_measures.diameter(g, e=ecc)
        sources = [u for u, d in ecc.items() if d == diameter]
        u, v, _ = max(
            (
                (u, *max(
                    single_source_dijkstra_path_length(g, u, weight="n").items(),
                    key=lambda x: (x[1], x)))
                for u in sources),
            key=lambda x: (x[2], x))
        return nx.algorithms.shortest_paths.generic.shortest_path(g, u, v)

    @staticmethod
    def count_adajcent_components(g):
        "Estimate the nubmer of adjacent components of each vertex."
        adjcomps = {}
        for u in g:
            vs = g.neighbors(u)
            adjcomps[u] = nx.algorithms.components.number_connected_components(g.subgraph(vs))
        return adjcomps

    def physlr_filter(self):
        "Filter a graph."
        g = self.read_graph(self.args.FILES)
        if self.args.n > 0:
            edges = [e for e, eprop in g.edges().items() if eprop["n"] < self.args.n]
            g.remove_edges_from(edges)
            print("Removed", len(edges), "edges with fewer than ", self.args.n, "common markers.",
                  file=sys.stderr)
        if self.args.M > 0:
            vertices = [u for u, prop in g.nodes().items() if prop["M"] >= self.args.M]
            g.remove_nodes_from(vertices)
            print("Removed", len(vertices), "vertices with", self.args.M, "or more molecules.",
                  file=sys.stderr)
        self.write_graph(g, sys.stdout, self.args.graph_format)

    def physlr_indexfa(self):
        "Index a set of sequences. The output file format is TSV."
        for filename in self.args.FILES:
            with open(filename) as fin:
                for name, seq, _ in read_fasta(fin):
                    print(name, "\t", sep="", end="")
                    print(*minimerize(self.args.k, self.args.w, seq.upper()))

    def physlr_indexlr(self):
        "Index a set of linked reads. The output file format is TSV."
        for filename in self.args.FILES:
            with open(filename) as fin:
                for _, seq, bx in read_fasta(fin):
                    print(bx, "\t", sep="", end="")
                    print(*minimerize(self.args.k, self.args.w, seq.upper()))

    def physlr_count_markers(self):
        "Count the frequency of each minimizer."
        if self.args.n == 0:
            self.args.n = 2
        bxtomin = self.read_minimizers(self.args.FILES)
        freq = {}
        for xs in bxtomin.values():
            for x in xs:
                if x in freq:
                    freq[x] += 1
                else:
                    freq[x] = 1
        print("Marker\tCount")
        for x, c in sorted(freq.items(), key=lambda x: x[1]):
            if c >= self.args.n:
                print(x, c, sep="\t")

    def physlr_intersect(self):
        "Print the minimizers in the intersection of each pair of barcodes."
        bxtomin = self.read_minimizers(self.args.FILES)
        mintobx = self.construct_minimizers_to_barcodes(bxtomin)
        seen = set()
        for bxs in mintobx.values():
            for u, v in itertools.combinations(bxs, 2):
                if not (u, v) in seen:
                    seen.add((u, v))
                    common = bxtomin[u] & bxtomin[v]
                    if len(common) >= self.args.n:
                        print(u, v, "", sep="\t", end="")
                        print(*common)

    def physlr_overlap(self):
        "Read a sketch of linked reads and find overlapping barcodes."
        bxtomin = self.read_minimizers(self.args.FILES)
        mintobx = self.construct_minimizers_to_barcodes(bxtomin)

        # Ignore markers that were seen only once.
        freqs = [len(bxs) for bxs in mintobx.values() if len(bxs) >= 2]
        q1, q2, q3 = quantile([0.25, 0.5, 0.75], freqs)
        whisker = int(q3 + self.args.coef * (q3 - q1))
        if self.args.C is None:
            self.args.C = whisker
        print(
            "Removing markers that are frequent and likely repetitive.\n",
            "Q1=", q1, " Q2=", q2, " Q3=", q3,
            " Q3+", self.args.coef, "*(Q3-Q1)=", whisker,
            " C=", self.args.C, sep="", file=sys.stderr)

        # Remove frequent (likely repetitive) minimizers.
        unique = set(x for x, bxs in mintobx.items() if len(bxs) < self.args.C)
        for xs in bxtomin.values():
            xs &= unique

        # Add the vertices.
        g = nx.Graph()
        for u, minimizers in sorted(bxtomin.items()):
            g.add_node(u, n=len(minimizers))

        # Add the overlap edges.
        for x, bxs in mintobx.items():
            if x not in unique:
                continue
            for u, v in itertools.combinations(bxs, 2):
                if not g.has_edge(u, v):
                    g.add_edge(u, v, n=len(bxtomin[u] & bxtomin[v]))

        # Write the graph.
        self.write_tsv(g, sys.stdout)

    def physlr_mst(self):
        "Determine the maximum spanning tree."
        g = self.read_graph(self.args.FILES)
        gmst = nx.algorithms.tree.mst.maximum_spanning_tree(g, weight="n")
        self.write_graph(gmst, sys.stdout, self.args.graph_format)

    def physlr_backbone(self):
        "Determine the backbone path of the graph."
        g = self.read_graph(self.args.FILES)
        gmst = nx.algorithms.tree.mst.maximum_spanning_tree(g, weight="n")
        backbone = self.determine_backbone(gmst)
        print(*backbone)

    def physlr_backbone_graph(self):
        "Determine the backbone-induced subgraph."
        g = self.read_graph(self.args.FILES)
        gmst = nx.algorithms.tree.mst.maximum_spanning_tree(g, weight="n")
        backbone = self.determine_backbone(gmst)
        subgraph = g.subgraph(backbone)
        self.write_graph(subgraph, sys.stdout, self.args.graph_format)

    def physlr_tiling_graph(self):
        "Determine the minimum-tiling-path-induced subgraph."
        g = self.read_graph(self.args.FILES)
        gmst = nx.algorithms.tree.mst.maximum_spanning_tree(g, weight="n")
        backbone = self.determine_backbone(gmst)
        if self.args.n == 0:
            self.args.n = min(g[u][v]["n"] for u, v in zip(backbone, backbone[1:]))
            print("Using n=", self.args.n, sep="", file=sys.stderr)
        g.remove_edges_from([e for e, eprop in g.edges().items() if eprop["n"] < self.args.n])
        u, v = backbone[0], backbone[-1]
        tiling_path = nx.algorithms.shortest_paths.generic.shortest_path(g, u, v)
        subgraph = g.subgraph(tiling_path)
        self.write_graph(subgraph, sys.stdout, self.args.graph_format)

    def physlr_molecules(self):
        "Estimate the nubmer of molecules per barcode."
        g = self.read_graph(self.args.FILES)
        g.remove_edges_from([e for e, eprop in g.edges().items() if eprop["n"] < self.args.n])
        for u, nmolecules in self.count_adajcent_components(g).items():
            g.nodes[u]["m"] = nmolecules
        self.write_graph(g, sys.stdout, self.args.graph_format)

    def physlr_graph(self, fmt):
        "Generate a graph from the minimizer index."
        graph = Graph()
        for filename in self.args.FILES:
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
            "-c", "--coef", action="store", dest="coef", type=float, default=1.5,
            help="ignore markers that occur in Q3+c*(Q3-Q1) or more barcodes [0]")
        argparser.add_argument(
            "-C", "--max-count", action="store", dest="C", type=int,
            help="ignore markers that occur in C or more barcodes [NA]")
        argparser.add_argument(
            "-M", "--max-molecules", action="store", dest="M", type=int,
            help="remove barcodes with M or more molecules [NA]")
        argparser.add_argument(
            "-n", "--min-n", action="store", dest="n", type=int,
            help="remove edges with fewer than n shared markers [0]")
        argparser.add_argument(
            "-O", "--output-format", action="store", dest="graph_format", default="tsv",
            help="the output graph file format [tsv]")
        argparser.add_argument(
            "command",
            help="A command")
        argparser.add_argument(
            "FILES", nargs="+",
            help="FASTA/FASTQ, TSV, or GraphViz format")
        return argparser.parse_args()

    def __init__(self):
        "Create a new instance of Physlr."
        self.args = self.parse_arguments()
        self.args.FILES = ["/dev/stdin" if s == "-" else s for s in self.args.FILES]

    def main(self):
        "Run Physlr."
        if self.args.command == "backbone":
            self.physlr_backbone()
        elif self.args.command == "backbone-graph":
            self.physlr_backbone_graph()
        elif self.args.command == "count-markers":
            self.physlr_count_markers()
        elif self.args.command == "filter":
            self.physlr_filter()
        elif self.args.command == "indexfa":
            self.physlr_indexfa()
        elif self.args.command == "indexlr":
            self.physlr_indexlr()
        elif self.args.command == "graphtsv":
            self.physlr_graph("tsv")
        elif self.args.command == "graphgv":
            self.physlr_graph("graphviz")
        elif self.args.command == "intersect":
            self.physlr_intersect()
        elif self.args.command == "molecules":
            self.physlr_molecules()
        elif self.args.command == "mst":
            self.physlr_mst()
        elif self.args.command == "overlap":
            self.physlr_overlap()
        elif self.args.command == "tiling-graph":
            self.physlr_tiling_graph()
        else:
            print("Unrecognized command:", self.args.command, file=sys.stderr)
            exit(1)

def main():
    "Run Physlr."
    Physlr().main()

if __name__ == "__main__":
    main()
