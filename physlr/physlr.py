#!/usr/bin/env python3
"""
Physlr: Physical Mapping of Linked Reads
"""

import argparse
import itertools
import multiprocessing
import os
import random
import re
import statistics
import sys
import timeit
from collections import Counter


import networkx as nx
import tqdm

from physlr.minimerize import minimerize
from physlr.read_fasta import read_fasta

t0 = timeit.default_timer()

def quantile(quantiles, xs):
    "Return the specified quantiles p of xs."
    sorted_xs = sorted(xs)
    return [sorted_xs[round(p * (len(sorted_xs)-1))] for p in quantiles]

def progress_bar_for_file(fin):
    "Return a progress bar for a file."
    return tqdm.tqdm(
        total=os.fstat(fin.fileno()).st_size,
        mininterval=1, smoothing=0.1,
        bar_format="{percentage:4.1f}% {elapsed} ETA {remaining} {bar}")

def progress(iterator):
    "Return an iterator that displays a progress bar."
    return tqdm.tqdm(
        iterator, mininterval=1, smoothing=0.1,
        bar_format="{percentage:4.1f}% {elapsed} ETA {remaining} {bar}")

class Physlr:
    """
    Physlr: Physical Mapping of Linked Reads
    """

    @staticmethod
    def print_graph_stats(g, fout=sys.stderr):
        "Print graph stats."
        v = g.number_of_nodes()
        e = g.number_of_edges()
        print(int(timeit.default_timer() - t0), f"V={v} E={e} E/V={round(e/v, 2)}", file=fout)

    @staticmethod
    def read_bed(filenames):
        "Read BED files."
        bed = []
        for filename in filenames:
            print(int(timeit.default_timer() - t0), "Reading", filename, file=sys.stderr)
            with open(filename) as fin:
                progressbar = progress_bar_for_file(fin)
                for line in fin:
                    progressbar.update(len(line))
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 5:
                        print("physlr: expected five or more BED fields:", line, file=sys.stderr)
                        exit(1)
                    tname, tstart, tend, qname, score = fields[0:5]
                    orientation = fields[5] if len(fields) >= 6 else "."
                    bed.append((tname, int(tstart), int(tend), qname, int(score), orientation))
                progressbar.close()
            print(int(timeit.default_timer() - t0), "Read", filename, file=sys.stderr)
        return bed

    @staticmethod
    def read_fastas(filenames):
        "Read FASTA files. Return a dictionary of names to sequences."
        seqs = {}
        for filename in filenames:
            print(int(timeit.default_timer() - t0), "Reading", filename, file=sys.stderr)
            with open(filename) as fin:
                for name, seq, _, _ in read_fasta(fin):
                    seqs[name] = seq
            print(
                int(timeit.default_timer() - t0),
                "Read", len(seqs), "sequences", file=sys.stderr)
        return seqs

    @staticmethod
    def read_path(filenames):
        "Read path files."
        paths = []
        for filename in filenames:
            print(int(timeit.default_timer() - t0), "Reading", filename, file=sys.stderr)
            with open(filename) as fin:
                progressbar = progress_bar_for_file(fin)
                for line in fin:
                    progressbar.update(len(line))
                    paths.append(line.split())
                progressbar.close()
            print(int(timeit.default_timer() - t0), "Read", filename, file=sys.stderr)
        return paths

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
            u, v = sorted(e)
            print(u, v, prop["n"], sep="\t", file=fout)

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
            progressbar = progress_bar_for_file(fin)
            line = fin.readline()
            progressbar.update(len(line))
            if line not in ["U\tn\n", "U\tn\tm\n"]:
                print("Unexpected header:", line, file=sys.stderr)
                exit(1)
            reading_vertices = True
            for line in fin:
                progressbar.update(len(line))
                if line == "\n":
                    line = fin.readline()
                    progressbar.update(len(line))
                    if line == "U\tV\tn\n":
                        reading_vertices = False
                    else:
                        print("Unexpected header:", line, file=sys.stderr)
                        exit(1)
                    line = fin.readline()
                    progressbar.update(len(line))
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
        progressbar.close()
        return g

    @staticmethod
    def read_graphviz(g, filename):
        "Read a GraphViz file."
        graph = nx.drawing.nx_agraph.read_dot(filename)
        for vprop in graph.nodes().values():
            vprop["n"] = int(vprop["n"])
        for _, _, eprop in graph.edges.data():
            eprop["n"] = int(eprop["n"])
        return nx.algorithms.operators.binary.compose(g, graph)

    # Complement nucleotides.
    TRANSLATE_COMPLEMENT = str.maketrans(
        "ACGTUNMRWSYKVHDBacgtunmrwsykvhdb",
        "TGCAANKYWSRMBDHVtgcaankywsrmbdhv")

    @staticmethod
    def reverse_complement(seq):
        "Return the reverse complement of this sequence."
        return seq[::-1].translate(Physlr.TRANSLATE_COMPLEMENT)

    @staticmethod
    def get_oriented_sequence(sequences, name_orientation):
        "Fetch and orient the specified sequence."
        name = name_orientation[0:-1]
        orientation = name_orientation[-1]
        if orientation == "+":
            return sequences[name]
        if orientation == "-":
            return Physlr.reverse_complement(sequences[name])
        print("physlr: Unexpected orientation:", orientation, file=sys.stderr)
        sys.exit(1)

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
        print(int(timeit.default_timer() - t0), "Reading", *filenames, file=sys.stderr)
        read_gv = False
        g = nx.Graph()
        for filename in filenames:
            with open(filename) as fin:
                c = fin.read(1)
                if c == "s":
                    g = Physlr.read_graphviz(g, filename)
                    read_gv = True
                elif c == "U":
                    g = Physlr.read_tsv(g, filename)
                else:
                    print("Unexpected graph format", c + fin.readline(), file=sys.stderr)
                    sys.exit(1)
        print(int(timeit.default_timer() - t0), "Read", *filenames, file=sys.stderr)
        if read_gv:
            print(int(timeit.default_timer() - t0), "Sorting the vertices", file=sys.stderr)
            g = Physlr.sort_vertices(g)
            print(int(timeit.default_timer() - t0), "Sorted the vertices", file=sys.stderr)
        return g

    @staticmethod
    def remove_singletons(g):
        "Remove singletons (isolated vertices) and return the number removed."
        singletons = [u for u, deg in g.degree if deg == 0]
        g.remove_nodes_from(singletons)
        return len(singletons)

    @staticmethod
    def filter_edges(g, arg_n):
        "Remove edges with n < arg_n."
        if arg_n == 0:
            return
        edges = [(u, v) for u, v, n in progress(g.edges(data="n")) if n < arg_n]
        print(
            int(timeit.default_timer() - t0),
            "Removed", len(edges), "edges with fewer than", arg_n,
            "common markers of", g.number_of_edges(),
            f"({round(100 * len(edges) / g.number_of_edges(), 2)}%)", file=sys.stderr)
        g.remove_edges_from(edges)

        num_singletons = Physlr.remove_singletons(g)
        print(
            int(timeit.default_timer() - t0),
            "Removed", num_singletons, "isolated vertices.", file=sys.stderr)

    @staticmethod
    def read_minimizers(filenames):
        "Read minimizers in TSV format. Returns unordered set."
        bxtomin = {}
        for filename in filenames:
            print(int(timeit.default_timer() - t0), "Reading", filename, file=sys.stderr)
            with open(filename) as fin:
                progressbar = progress_bar_for_file(fin)
                for line in fin:
                    progressbar.update(len(line))
                    fields = line.split(None, 1)
                    if len(fields) < 2:
                        continue
                    bx = fields[0]
                    if bx not in bxtomin:
                        bxtomin[bx] = set()
                    bxtomin[bx].update(int(x) for x in fields[1].split())
                progressbar.close()
            print(int(timeit.default_timer() - t0), "Read", filename, file=sys.stderr)
        return bxtomin

    @staticmethod
    def read_minimizers_list(filenames):
        "Read minimizers in TSV format. Returns ordered list."
        bxtomin = {}
        for filename in filenames:
            print(int(timeit.default_timer() - t0), "Reading", filename, file=sys.stderr)
            with open(filename) as fin:
                progressbar = progress_bar_for_file(fin)
                for line in fin:
                    progressbar.update(len(line))
                    fields = line.split(None, 1)
                    if len(fields) < 2:
                        continue
                    bx = fields[0]
                    if bx in bxtomin:
                        print("Error: Expected single id per in file", file=sys.stderr)
                        exit(1)
                    bxtomin[bx] = [int(x) for x in fields[1].split()]
                progressbar.close()
            print(int(timeit.default_timer() - t0), "Read", filename, file=sys.stderr)
        return bxtomin

    @staticmethod
    def count_molecules_per_bx(moltomin):
        "Iterate over minimizers dictionary, track # molecules per barcode"
        mol_counts = Counter()
        bx_match = re.compile(r'^(\S+)_(\d+)$')
        for bx_mol in moltomin:
            bx_mol_match = re.search(bx_match, bx_mol)
            if bx_mol_match:
                mol_counts[bx_mol_match.group(1)] = max(mol_counts[bx_mol_match.group(1)],
                                                        int(bx_mol_match.group(2)) + 1)
        return mol_counts

    @staticmethod
    def construct_minimizers_to_barcodes(bxtomin):
        "Construct a dictionary of minimizers to barcodes."
        mintobx = {}
        for bx, minimizers in progress(bxtomin.items()):
            for x in minimizers:
                if x not in mintobx:
                    mintobx[x] = set()
                mintobx[x].add(bx)
        print(
            int(timeit.default_timer() - t0),
            "Indexed", len(mintobx), "minimizers", file=sys.stderr)
        return mintobx

    @staticmethod
    def triconnected_components(g):
        "Return the triconnected components of the graph."
        components = []
        for component in nx.biconnected_components(g):
            if len(component) < 3:
                components.append(component)
                continue
            try:
                cuts = next(nx.all_node_cuts(g.subgraph(component), k=2))
                if len(cuts) > 2:
                    components.append(component)
                    continue
                assert len(cuts) == 2
                subcomponents = list(nx.connected_components(g.subgraph(component - cuts)))
                if len(subcomponents) == 1:
                    components.append(component)
                    continue
                components += subcomponents
                components.append(cuts)
            except StopIteration:
                components.append(component)
        return components

    @staticmethod
    def diameter_of_tree(g, weight=None):
        """
        Compute the diameter of a tree.
        The diameter of an arbitrary component is returned if there are multiple components."
        """
        u = next(iter(g.nodes))
        paths = nx.shortest_path_length(g, u, weight=weight)
        u, _ = max(paths.items(), key=lambda x: x[1])
        paths = nx.shortest_path_length(g, u, weight=weight)
        v, diameter = max(paths.items(), key=lambda x: x[1])
        return (u, v, diameter)

    @staticmethod
    def determine_backbones_of_trees(g):
        "Determine the backbones of the maximum spanning trees."
        paths = []
        for component in nx.connected_components(g):
            gcomponent = g.subgraph(component)
            u, v, _ = Physlr.diameter_of_tree(gcomponent, weight="n")
            path = nx.shortest_path(gcomponent, u, v, weight="n")
            paths.append(path)
        paths.sort(key=len, reverse=True)
        return paths

    @staticmethod
    def determine_backbones(g):
        "Determine the backbones of the graph."
        g = g.copy()
        backbones = []
        while not nx.is_empty(g):
            gmst = nx.maximum_spanning_tree(g, weight="n")
            paths = Physlr.determine_backbones_of_trees(gmst)
            backbones.extend(paths)
            vertices = [u for path in paths for u in path]
            neighbors = [v for u in vertices for v in g.neighbors(u)]
            g.remove_nodes_from(vertices)
            g.remove_nodes_from(neighbors)
            Physlr.remove_singletons(g)
        backbones.sort(key=len, reverse=True)
        print(int(timeit.default_timer() - t0), "Determined the backbone paths", file=sys.stderr)
        return backbones

    @staticmethod
    def print_flesh_path(backbone, backbone_insertions):
        "Print out the backbone path with 'flesh' barcodes added"
        for i, mol in enumerate(backbone):
            print(mol, file=sys.stdout, end=' ')
            if i in backbone_insertions:
                insertions = ",".join(backbone_insertions[i])
                print("(" + insertions + ")", file=sys.stdout, end=' ')
        print(backbone[len(backbone)-1], file=sys.stdout)

    # RE match for read header
    header_prefix_re = re.compile(r'^(\S+?)(\/[1-2])?')

    @staticmethod
    def is_valid_pair(bx1, bx2, name1, name2):
        "Checks if read pair's headers match (sanity check), and they have associated barcodes"
        header_prefix_match_r1 = re.search(Physlr.header_prefix_re, name1)
        header_prefix_match_r2 = re.search(Physlr.header_prefix_re, name2)
        return bx1 is not None and bx2 is not None and bx1 == bx2 and \
            header_prefix_match_r1.group(1) == header_prefix_match_r2.group(1)

    def physlr_filter(self):
        "Filter a graph."
        g = self.read_graph(self.args.FILES)
        Physlr.filter_edges(g, self.args.n)
        if self.args.M is not None:
            vertices = [u for u, prop in g.nodes().items() if prop["m"] >= self.args.M]
            g.remove_nodes_from(vertices)
            print(
                int(timeit.default_timer() - t0),
                "Removed", len(vertices), "vertices with", self.args.M, "or more molecules.",
                file=sys.stderr)
        if self.args.min_component_size > 0:
            ncomponents, nvertices = 0, 0
            vertices = set()
            for component in nx.connected_components(g):
                if len(component) < self.args.min_component_size:
                    vertices.update(component)
                    ncomponents += 1
                    nvertices += len(component)
            g.remove_nodes_from(vertices)
            print(
                int(timeit.default_timer() - t0),
                "Removed", nvertices, "vertices in", ncomponents, "components",
                "with fewer than", self.args.min_component_size, "vertices in a component.",
                file=sys.stderr)
        self.write_graph(g, sys.stdout, self.args.graph_format)

    def physlr_flesh_backbone(self):
        "Flesh out the barcodes in the backbone paths"
        g = self.read_graph([self.args.FILES[0]])
        backbones_raw = self.read_path([self.args.FILES[1]])
        backbones = [b for b in backbones_raw if len(b) >= self.args.min_component_size]
        print(
            int(timeit.default_timer() - t0),
            "Removed", len(backbones_raw) - len(backbones),
            "backbones of", len(backbones_raw),
            "with fewer than", self.args.min_component_size,
            "vertices", file=sys.stderr)
        for backbone in backbones:
            backbone_insertions = {}
            neighbours = {v for u in backbone for v in g.neighbors(u)}
            # Find where the neighbours should go in the backbone path
            for neighbour in neighbours:
                if neighbour in backbone:
                    continue
                (max_n, max_index) = (float("-inf"), float("-inf"))
                for i, k in enumerate(backbone):
                    if not g.has_edge(neighbour, k):
                        continue
                    n = g[neighbour][k]["n"]
                    if n > max_n:
                        (max_n, max_index) = (n, i)
                # Put R or L of node with most shared markers?
                if max_index == 0:
                    insert_index = max_index
                elif max_index == len(backbone) - 1:
                    insert_index = max_index - 1
                else:
                    r_n = g[neighbour][backbone[max_index+1]]['n'] \
                        if g.has_edge(neighbour, backbone[max_index+1]) else 0
                    l_n = g[neighbour][backbone[max_index-1]]['n'] \
                        if g.has_edge(neighbour, backbone[max_index-1]) else 0
                    if l_n > r_n:
                        insert_index = max_index - 1
                    else:
                        insert_index = max_index
                if insert_index not in backbone_insertions:
                    backbone_insertions[insert_index] = []
                backbone_insertions[insert_index].append(neighbour)
            self.print_flesh_path(backbone, backbone_insertions)

    def physlr_subgraph(self):
        "Extract a vertex-induced subgraph."
        if self.args.d not in (0, 1):
            exit("physlr subgraph: error: Only -d0 and -d1 are currently supported.")
        vertices = set(self.args.v.split())
        exclude_vertices = set(self.args.exclude_vertices.split())
        g = self.read_graph(self.args.FILES)
        if self.args.d == 1:
            vertices.update(v for u in vertices for v in g.neighbors(u))
        subgraph = g.subgraph(vertices - exclude_vertices)
        print(int(timeit.default_timer() - t0), "Extracted subgraph", file=sys.stderr)
        self.write_graph(subgraph, sys.stdout, self.args.graph_format)
        print(int(timeit.default_timer() - t0), "Wrote graph", file=sys.stderr)

    def physlr_subgraphs(self):
        "Extract multiple vertex-induced subgraphs."
        if self.args.output is None:
            exit("physlr subgraphs: missing parameter: --output is need but not provided.")
        if self.args.d not in (0, 1):
            exit("physlr subgraphs: error: Only -d0 and -d1 are currently supported.")
        vertices = set(self.args.v.split(","))
        exclude_vertices = set(self.args.exclude_vertices.split(","))
        g = self.read_graph(self.args.FILES)
        num_empty_subgraphs = 0
        print(int(timeit.default_timer() - t0),
              "Extracting and writing (sub)graphs",
              file=sys.stderr)
        if not os.path.exists(self.args.output):
            os.makedirs(self.args.output)
        for u in progress(vertices):
            vertices_u = set()
            if self.args.exclude_source == 0:
                vertices_u.add(u)
            if self.args.d == 1:
                vertices_u.update(v for v in g.neighbors(u))
            subgraph = g.subgraph(vertices_u - exclude_vertices)
            if subgraph.number_of_nodes() == 0:
                num_empty_subgraphs += 1
            else:
                fout = open(self.args.output+"/"+u+"."+self.args.graph_format, "w+")
                self.write_graph(subgraph, fout, self.args.graph_format)
                fout.close()
        print(int(timeit.default_timer() - t0),
              "Number of empty subgraphs (not written):", num_empty_subgraphs,
              file=sys.stderr)
        print(int(timeit.default_timer() - t0), "Wrote graphs", file=sys.stderr)

    def physlr_indexfa(self):
        "Index a set of sequences. The output file format is TSV."
        for filename in self.args.FILES:
            with open(filename) as fin:
                for name, seq, _, _ in read_fasta(fin):
                    print(name, "\t", sep="", end="")
                    print(*minimerize(self.args.k, self.args.w, seq.upper()))

    def physlr_indexlr(self):
        "Index a set of linked reads. The output file format is TSV."
        for filename in self.args.FILES:
            with open(filename) as fin:
                for _, seq, bx, _ in read_fasta(fin):
                    print(bx, "\t", sep="", end="")
                    print(*minimerize(self.args.k, self.args.w, seq.upper()))

    def physlr_count_markers(self):
        "Count the frequency of each minimizer."
        bxtomin = self.read_minimizers(self.args.FILES)
        marker_counts = Counter(x for markers in progress(bxtomin.values()) for x in markers)
        print(
            int(timeit.default_timer() - t0),
            "Counted", len(marker_counts), "minimizers", file=sys.stderr)

        print("Marker\tCount")
        for x, n in sorted(marker_counts.items(), key=lambda x: x[1]):
            if n >= 2:
                print(x, n, sep="\t")

    def physlr_intersect(self):
        "Print the minimizers in the intersection of each pair of barcodes."
        if self.args.n == 0:
            self.args.n = 1
        bxtomin = self.read_minimizers(self.args.FILES)
        mintobx = self.construct_minimizers_to_barcodes(bxtomin)
        if self.args.v:
            pairs = itertools.combinations(self.args.v.split(), 2)
        else:
            pairs = {(u, v) for bxs in mintobx.values() for u, v in itertools.combinations(bxs, 2)}
        for u, v in pairs:
            common = bxtomin[u] & bxtomin[v]
            if len(common) >= self.args.n:
                print(u, v, "", sep="\t", end="")
                print(*common)

    @staticmethod
    def remove_singleton_markers(bxtomin):
        """
        Remove markers that occur only once.
        Return the counts of markers after removing singletons.
        """
        marker_counts = Counter(x for markers in progress(bxtomin.values()) for x in markers)
        print(
            int(timeit.default_timer() - t0),
            "Counted", len(marker_counts), "minimizers", file=sys.stderr)

        singletons = {x for x, n in progress(marker_counts.items()) if n < 2}
        for markers in progress(bxtomin.values()):
            markers -= singletons
        print(
            int(timeit.default_timer() - t0),
            "Removed", len(singletons), "minimizers that occur only once of", len(marker_counts),
            f"({round(100 * len(singletons) / len(marker_counts), 2)}%)", file=sys.stderr)
        for marker in singletons:
            del marker_counts[marker]
        return marker_counts

    def physlr_filter_barcodes(self):
        """
        Filter barcodes by number of markers.
        Read a TSV file of barcodes to minimizers.
        Remove markers that occur only once.
        Remove barkers with too few or too many markers.
        Write a TSV file of barcodes to minimizers.
        """
        bxtomin = self.read_minimizers(self.args.FILES)
        Physlr.remove_singleton_markers(bxtomin)

        q0, q1, q2, q3, q4 = quantile(
            [0, 0.25, 0.5, 0.75, 1], (len(markers) for markers in bxtomin.values()))
        low_whisker = int(q1 - self.args.coef * (q3 - q1))
        high_whisker = int(q3 + self.args.coef * (q3 - q1))
        if self.args.n == 0:
            self.args.n = max(q0, low_whisker)
        if self.args.N is None:
            self.args.N = min(1 + q4, high_whisker)

        print(
            int(timeit.default_timer() - t0), " Counted markers per barcode\n",
            f"    Markers per barcode: Q0={q0} Q1={q1} Q2={q2} Q3={q3} Q4={q4} Q3-Q1={q3 - q1}\n",
            f"    Q3-{self.args.coef}*(Q3-Q1)={low_whisker} n={self.args.n}\n",
            f"    Q3+{self.args.coef}*(Q3-Q1)={high_whisker} N={self.args.N}",
            sep="", file=sys.stderr)

        too_few, too_many = 0, 0
        for bx, markers in progress(bxtomin.items()):
            if len(markers) < self.args.n:
                too_few += 1
            elif len(markers) >= self.args.N:
                too_many += 1
            else:
                print(bx, "\t", sep="", end="")
                print(*markers)
        print(
            "    Discarded", too_few, "barcodes with too few markers of", len(bxtomin),
            f"({round(100 * too_few / len(bxtomin), 2)}%)", file=sys.stderr)
        print(
            "    Discarded", too_many, "barcodes with too many markers of", len(bxtomin),
            f"({round(100 * too_many / len(bxtomin), 2)}%)", file=sys.stderr)
        print(
            int(timeit.default_timer() - t0),
            "Wrote", len(bxtomin) - too_few - too_many, "barcodes", file=sys.stderr)

    def physlr_filter_minimizers(self):
        "Filter minimizers by depth of coverage. Remove repetitive minimizers."
        bxtomin = self.read_minimizers(self.args.FILES)
        marker_counts = Physlr.remove_singleton_markers(bxtomin)

        # Identify frequent markers.
        q1, q2, q3 = quantile([0.25, 0.5, 0.75], marker_counts.values())
        low_whisker = int(q1 - self.args.coef * (q3 - q1))
        high_whisker = int(q3 + self.args.coef * (q3 - q1))
        if self.args.C is None:
            self.args.C = high_whisker
        print(
            int(timeit.default_timer() - t0),
            " Minimizer frequency: Q1=", q1, " Q2=", q2, " Q3=", q3,
            " Q1-", self.args.coef, "*(Q3-Q1)=", low_whisker,
            " Q3+", self.args.coef, "*(Q3-Q1)=", high_whisker,
            " C=", self.args.C, sep="", file=sys.stderr)

        # Remove frequent markers.
        frequent_markers = {x for x, count in marker_counts.items() if count >= self.args.C}
        num_empty_barcodes = 0
        for bx, markers in progress(bxtomin.items()):
            markers -= frequent_markers
            if not markers:
                num_empty_barcodes += 1
                continue
            print(bx, "\t", sep="", end="")
            print(*markers)

        print(
            int(timeit.default_timer() - t0),
            "Removed", len(frequent_markers), "most frequent minimizers of", len(marker_counts),
            f"({round(100 * len(frequent_markers) / len(marker_counts), 2)}%)", file=sys.stderr)
        print(
            int(timeit.default_timer() - t0),
            "Removed", num_empty_barcodes, "empty barcodes of", len(bxtomin),
            f"({round(100 * num_empty_barcodes / len(bxtomin), 2)}%)", file=sys.stderr)
        print(
            int(timeit.default_timer() - t0),
            "Wrote", len(bxtomin) - num_empty_barcodes, "barcodes", file=sys.stderr)

    def remove_repetitive_minimizers(self, bxtomin, mintobx):
        "Remove repetitive minimizers."

        # Remove markers that occur only once.
        num_markers = len(mintobx)
        singletons = {x for x, bxs in progress(mintobx.items()) if len(bxs) < 2}
        for x in singletons:
            del mintobx[x]
        for markers in progress(bxtomin.values()):
            markers -= singletons
        print(
            int(timeit.default_timer() - t0),
            "Removed", len(singletons), "minimizers that occur only once of", num_markers,
            f"({round(100 * len(singletons) / num_markers, 2)}%)", file=sys.stderr)

        # Identify repetitive markers.
        q1, q2, q3 = quantile([0.25, 0.5, 0.75], (len(bxs) for bxs in mintobx.values()))
        whisker = int(q3 + self.args.coef * (q3 - q1))
        if self.args.C is None:
            self.args.C = whisker
        print(
            int(timeit.default_timer() - t0),
            " Minimizer frequency: Q1=", q1, " Q2=", q2, " Q3=", q3,
            " Q3+", self.args.coef, "*(Q3-Q1)=", whisker,
            " C=", self.args.C, sep="", file=sys.stderr)

        # Remove frequent (likely repetitive) minimizers.
        num_markers = len(mintobx)
        repetitive = {x for x, bxs in mintobx.items() if len(bxs) >= self.args.C}
        for x in repetitive:
            del mintobx[x]
        for xs in progress(bxtomin.values()):
            xs -= repetitive
        print(
            int(timeit.default_timer() - t0),
            "Removed", len(repetitive), "most frequent minimizers of", num_markers,
            f"({round(100 * len(repetitive) / num_markers, 2)}%)", file=sys.stderr)

    def physlr_overlap(self):
        "Read a sketch of linked reads and find overlapping barcodes."

        bxtomin = self.read_minimizers(self.args.FILES)
        mintobx = self.construct_minimizers_to_barcodes(bxtomin)

        # Add the vertices.
        g = nx.Graph()
        for u, minimizers in sorted(progress(bxtomin.items())):
            if len(minimizers) >= self.args.n:
                g.add_node(u, n=len(minimizers))
        print(
            int(timeit.default_timer() - t0),
            "Added", g.number_of_nodes(), "barcodes to the graph", file=sys.stderr)

        # Add the overlap edges.
        edges = Counter(
            (u, v) for bxs in progress(mintobx.values()) for u, v in itertools.combinations(bxs, 2))
        print(int(timeit.default_timer() - t0), "Loaded", len(edges), "edges", file=sys.stderr)

        for (u, v), n in progress(edges.items()):
            if n >= self.args.n:
                g.add_edge(u, v, n=n)
        num_removed = len(edges) - g.number_of_edges()
        print(
            int(timeit.default_timer() - t0),
            "Removed", num_removed, "edges with fewer than", self.args.n,
            "common markers of", len(edges),
            f"({round(100 * num_removed / len(edges), 2)}%)", file=sys.stderr)

        num_singletons = Physlr.remove_singletons(g)
        print(
            int(timeit.default_timer() - t0),
            "Removed", num_singletons, "isolated vertices.", file=sys.stderr)

        Physlr.print_graph_stats(g)

        # Write the graph.
        self.write_tsv(g, sys.stdout)
        print(int(timeit.default_timer() - t0), "Wrote the graph", file=sys.stderr)

    def physlr_mst(self):
        "Determine the maximum spanning tree."
        g = self.read_graph(self.args.FILES)
        gmst = nx.algorithms.tree.mst.maximum_spanning_tree(g, weight="n")
        self.write_graph(gmst, sys.stdout, self.args.graph_format)

    def physlr_backbone(self):
        "Determine the backbone path of the graph."
        g = self.read_graph(self.args.FILES)
        backbones = self.determine_backbones(g)
        for backbone in backbones:
            print(*backbone)

    def physlr_backbone_graph(self):
        "Determine the backbone-induced subgraph."
        g = self.read_graph(self.args.FILES)
        Physlr.remove_singletons(g)
        backbones = self.determine_backbones(g)
        backbone = (u for path in backbones for u in path)
        subgraph = self.sort_vertices(g.subgraph(backbone))
        self.write_graph(subgraph, sys.stdout, self.args.graph_format)
        print(int(timeit.default_timer() - t0), "Output the backbone subgraph", file=sys.stderr)

    def physlr_biconnected_components(self):
        "Separate a graph into its biconnected components by removing its cut vertices."
        g = self.read_graph(self.args.FILES)
        cut_vertices = list(nx.articulation_points(g))
        g.remove_nodes_from(cut_vertices)
        print(
            int(timeit.default_timer() - t0),
            "Removed", len(cut_vertices), "cut vertices.", file=sys.stderr)
        self.write_graph(g, sys.stdout, self.args.graph_format)
        print(int(timeit.default_timer() - t0), "Wrote graph", file=sys.stderr)

    def physlr_tiling_graph(self):
        "Determine the minimum-tiling-path-induced subgraph."
        g = self.read_graph(self.args.FILES)
        Physlr.filter_edges(g, self.args.n)
        backbones = self.determine_backbones(g)
        tiling = {u for path in backbones for u in nx.shortest_path(g, path[0], path[-1])}
        subgraph = g.subgraph(tiling)
        self.write_graph(subgraph, sys.stdout, self.args.graph_format)

    def physlr_count_molecules(self):
        "Estimate the nubmer of molecules per barcode."
        g = self.read_graph(self.args.FILES)
        Physlr.filter_edges(g, self.args.n)
        print(
            int(timeit.default_timer() - t0),
            "Separating barcodes into molecules", file=sys.stderr)

        for u, prop in progress(g.nodes.items()):
            subgraph = g.subgraph(g.neighbors(u))
            # Ignore K3 (triangle) components.
            prop["m"] = sum(
                1 for component in nx.biconnected_components(subgraph)
                if len(component) >= 4)
        self.write_graph(g, sys.stdout, self.args.graph_format)

    @staticmethod
    def split_minimizers_bx(bx, g, bxtomin):
        "Partition the minimizers of the given barcode"
        bx_match = re.compile(r'^(\S+)_\d+$')
        bx_min = bxtomin[bx]
        mol = 0
        mol_list = []
        while g.has_node(bx + "_" + str(mol)):
            bxmol = bx + "_" + str(mol)
            neighbour_minimizers_list = [bxtomin[re.search(bx_match, v).group(1)] \
                                         for v in g.neighbors(bxmol) \
                                         if re.search(bx_match, v).group(1) in bxtomin]
            if not neighbour_minimizers_list:
                neighbour_minimizers_list = [set()]
            neighbour_minimizers_set = set.union(*neighbour_minimizers_list)
            molec_minimizers = set.intersection(bx_min, neighbour_minimizers_set)
            mol_list.append((bxmol, molec_minimizers))
            mol += 1
        return mol_list

    @staticmethod
    def split_minimizers_bx_process(bx):
        """
        Partition the minimizers of this barcode.
        The Graph and bx->min dictionary are passed as class variables.
        """
        return Physlr.split_minimizers_bx(bx, Physlr.graph, Physlr.bxtomin)

    def physlr_split_minimizers(self):
        "Given the molecule overlap graph, split the minimizers into molecules"
        if len(self.args.FILES) < 2:
            exit("physlr split-minimizers: error: graph file and bx to minimizer inputs required")
        g = self.read_graph([self.args.FILES[0]])
        bxtomin = self.read_minimizers([self.args.FILES[1]])

        if self.args.threads == 1:
            moltomin = [self.split_minimizers_bx(bx, g, bxtomin) for bx in progress(bxtomin)]
            moltomin = dict(x for l in moltomin for x in l)

        else:
            Physlr.graph = g
            Physlr.bxtomin = bxtomin
            with multiprocessing.Pool(self.args.threads) as pool:
                moltomin = dict(x for l in pool.map(self.split_minimizers_bx_process,
                                                    progress(bxtomin), chunksize=100) for x in l)
            Physlr.graph = None
            Physlr.bxtomin = None

        empty_ct = 0
        for mol in moltomin:
            if not moltomin[mol]:
                empty_ct += 1
                print("%s\t%s" % (mol, ""), file=sys.stdout)
                print("Warning:", mol, "has no associated minimizers", file=sys.stderr)
            else:
                print("%s\t%s" % (mol, " ".join(map(str, moltomin[mol]))), file=sys.stdout)

    def physlr_split_reads_molecules(self):
        "Given the molecule -> minimizers table and the reads, partition reads into molecules"
        if len(self.args.FILES) < 3:
            exit("physlr split-reads-molecules: error: molecule minimizers,\
            bx minimizers, reads inputs required")
        moltomin = self.read_minimizers([self.args.FILES[0]])
        mol_counts = self.count_molecules_per_bx(moltomin)
        bxtomin = self.args.FILES[1]
        num_pairs, num_valid_pairs, num_no_min, num_equal_min, num_no_int_min = 0, 0, 0, 0, 0

        readmin = open(bxtomin, 'r')

        read_count = 0
        for readfile in self.args.FILES[2:]:
            with open(readfile, 'r') as fin:
                for name, seq, bx, qual in read_fasta(fin):
                    if read_count == 0:
                        (name1, seq1, bx1, qual1) = (name, seq, bx, qual)
                        read_count += 1
                    else:
                        num_pairs += 1
                        (min_info1, min_info2) = (readmin.readline().strip(),
                                                  readmin.readline().strip())
                        if self.is_valid_pair(bx1, bx, name1, name) and bx in mol_counts:
                            (bx1_mol, minimizers1) = self.parse_minimizer_line(min_info1)
                            (bx2_mol, minimizers2) = self.parse_minimizer_line(min_info2)
                            if bx1_mol != bx2_mol or bx1_mol != bx1:
                                print("Should match: ", bx1_mol, bx2_mol, bx1, file=sys.stderr)
                                exit("Error: Minimizer TSV order doesn't match reads fq file")

                            (mol, inc_no_int_min, inc_equal_min) = self.assign_read_molecule(
                                set.union(minimizers1, minimizers2), moltomin, mol_counts, bx1)
                            bx_mol = bx1 + mol
                            num_no_int_min += inc_no_int_min
                            num_equal_min += inc_equal_min
                            self.print_read(name1 + " BX:Z:" + bx_mol, seq1, qual1)
                            self.print_read(name + " BX:Z:" + bx_mol, seq, qual)
                            num_valid_pairs += 1
                        elif self.is_valid_pair(bx1, bx, name1, name) and \
                                not self.args.molecules_bx_only:
                            self.print_read(name1 + " BX:Z:" + bx1, seq1, qual1)
                            self.print_read(name + " BX:Z:" + bx, seq, qual)
                            num_no_min += 1
                        elif not self.args.molecules_bx_only:
                            self.print_read(name1, seq1, qual1)
                            self.print_read(name, seq, qual)
                        read_count = 0

        print("Saw", num_pairs, "read pairs, saw",
              num_valid_pairs, "valid read pairs with associated barcodes.",
              num_no_min, "read pairs' barcodes had no split minimizers.",
              num_no_int_min, "read pairs' barcodes had no intersecting minimizers",
              num_equal_min, "read pairs had multiple molecule assignments",
              file=sys.stderr)
        readmin.close()

    @staticmethod
    def print_read(name, seq, qual):
        "Prints a read to stdout"
        print("@", name, "\n", seq, "\n+\n", qual, sep="")

    @staticmethod
    def parse_minimizer_line(min_line):
        "Given a minimizer line from a read, parse out the barcode_mol and set of minimizers"
        min_line = min_line.split("\t", 1)
        if len(min_line) == 2:
            return (min_line[0], set(map(int, min_line[1].split())))
        return (min_line[0], set())

    @staticmethod
    def assign_read_molecule(minimizers, moltomin, mol_counts, bx):
        "Given the minimizers of a read and the barcode, assign to a molecule"
        intersections = {mol: len(set.intersection(minimizers, moltomin[bx + "_" + str(mol)]))
                         for mol in range(0, mol_counts[bx]) if bx + "_" + str(mol) in moltomin}
        if not intersections.values():
            return "", 1, 0
        max_intersection = max(intersections.values())
        if max_intersection == 0:
            return "", 1, 0
        max_hits = [mol for mol in intersections if intersections[mol] == max_intersection]
        if len(max_hits) == 1:
            return "_" + str(max_hits[0]), 0, 0
        return "_" + str(random.choice(max_hits)), 0, 1

    @staticmethod
    def determine_molecules(g, u):
        "Assign the neighbours of this vertex to molecules."
        cut_vertices = set(nx.articulation_points(g.subgraph(g.neighbors(u))))
        components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)) - cut_vertices)))
        components.sort(key=len, reverse=True)
        return u, {v: i for i, vs in enumerate(components) if len(vs) > 1 for v in vs}

    @staticmethod
    def determine_molecules_process(u):
        """
        Assign the neighbours of this vertex to molecules.
        The graph is passed in the class variable Physlr.graph.
        """
        return Physlr.determine_molecules(Physlr.graph, u)

    def physlr_molecules(self):
        "Separate barcodes into molecules."
        gin = self.read_graph(self.args.FILES)
        Physlr.filter_edges(gin, self.args.n)
        print(
            int(timeit.default_timer() - t0),
            "Separating barcodes into molecules", file=sys.stderr)

        # Parition the neighbouring vertices of each barcode into molecules.
        if self.args.threads == 1:
            molecules = dict(self.determine_molecules(gin, u) for u in progress(gin))
        else:
            Physlr.graph = gin
            with multiprocessing.Pool(self.args.threads) as pool:
                molecules = dict(pool.map(
                    self.determine_molecules_process, progress(gin), chunksize=100))
            Physlr.graph = None
        print(int(timeit.default_timer() - t0), "Identified molecules", file=sys.stderr)

        # Add vertices.
        gout = nx.Graph()
        for u, vs in sorted(molecules.items()):
            n = gin.nodes[u]["n"]
            nmolecules = 1 + max(vs.values()) if vs else 0
            for i in range(nmolecules):
                gout.add_node(f"{u}_{i}", n=n)

        print(
            int(timeit.default_timer() - t0),
            "Identified", gout.number_of_nodes(), "molecules in",
            gin.number_of_nodes(), "barcodes.",
            round(gout.number_of_nodes() / gin.number_of_nodes(), 2), "mean molecules per barcode",
            file=sys.stderr)

        # Add edges.
        for (u, v), prop in gin.edges.items():
            # Skip singleton and cut vertices, which are excluded from the partition.
            if v not in molecules[u] or u not in molecules[v]:
                continue
            u_molecule = molecules[u][v]
            v_molecule = molecules[v][u]
            gout.add_edge(f"{u}_{u_molecule}", f"{v}_{v_molecule}", n=prop["n"])
        print(int(timeit.default_timer() - t0), "Separated molecules", file=sys.stderr)

        self.write_graph(gout, sys.stdout, self.args.graph_format)
        print(int(timeit.default_timer() - t0), "Wrote graph", file=sys.stderr)

    @staticmethod
    def index_markers_in_backbones(backbones, bxtomin):
        "Index the positions of the markers in the backbones."
        markertopos = {}
        for tid, path in enumerate(progress(backbones)):
            for pos, (u, v) in enumerate(zip(path, path[1:])):
                u = u.split("_", 1)[0]
                v = u.split("_", 1)[0]
                for marker in bxtomin[u] & bxtomin[v]:
                    markertopos.setdefault(marker, set()).add((tid, pos))
        print(
            int(timeit.default_timer() - t0),
            "Indexed", len(markertopos), "markers", file=sys.stderr)
        return markertopos

    @staticmethod
    def determine_orientation(x, y, z):
        "Determine the orientation of an alignment."
        if x is not None and z is not None:
            return "." if x == y == z else "+" if x <= y <= z else "-" if x >= y >= z else "."
        if x is not None:
            return "+" if x < y else "-" if x > y else "."
        if z is not None:
            return "+" if y < z else "-" if y > z else "."
        return "."

    def physlr_filter_markers(self):
        """
        Removes repeats from minimizers and outputs result to stdout
        """
        target_filenames = [self.args.FILES[0]]
        bxtomin = self.read_minimizers(target_filenames)
        mintobx = self.construct_minimizers_to_barcodes(bxtomin)
        self.remove_repetitive_minimizers(bxtomin, mintobx)
        for bx, markers in progress(bxtomin.items()):
            print(bx, "\t", sep="", end="")
            print(*markers)

    def physlr_map_mkt(self):
        """
        Map sequences to a physical map.
        Usage: physlr map TGRAPH.tsv TMARKERS.tsv QMARKERS.tsv... >MAP.bed
        """
        import physlr.mkt
        import numpy

        if len(self.args.FILES) < 3:
            exit("physlr map: error: at least three file arguments are required")
        graph_filenames = [self.args.FILES[0]]
        target_filenames = [self.args.FILES[1]]
        query_filenames = self.args.FILES[2:]

        g = self.read_graph(graph_filenames)
        bxtomin = self.read_minimizers(target_filenames)
        query_markers = self.read_minimizers_list(query_filenames)

        # Index the positions of the markers in the backbone.
        backbones = Physlr.determine_backbones(g)
        markertopos = Physlr.index_markers_in_backbones(backbones, bxtomin)

        # Map the query sequences to the physical map.
        num_mapped = 0
        for qid, markers in progress(query_markers.items()):
            # Count the number of markers mapped to each target position.
            tidpos_to_n = Counter(pos for marker in markers for pos in markertopos.get(marker, ()))
            # Map each target position to a query position.
            #tid->tpos->qpos_list
            tid_to_qpos = {}
            for qpos, marker in enumerate(markers):
                for (tid, tpos) in markertopos.get(marker, ()):
                    if not tid in tid_to_qpos:
                        tid_to_qpos[tid] = {}
                    if not tpos in tid_to_qpos[tid]:
                        tid_to_qpos[tid][tpos] = []
                    tid_to_qpos[tid][tpos].append(qpos)

            tid_to_mkt = {}
            for (tid, tpos_to_qpos) in tid_to_qpos.items():
                #build array of the time points of measurements
                #build array containing the measurements corresponding to entries of time
                timepoints = []
                measurements = []
                for (tpos, qpos_list) in tpos_to_qpos.items():
                    #do not use islands (noise?)
                    if tpos + 1 in tpos_to_qpos or tpos - 1 in tpos_to_qpos:
                        for qpos in qpos_list:
                            timepoints.append(tpos)
                            measurements.append(qpos)
                if timepoints:
                    tid_to_mkt[tid] = physlr.mkt.test(numpy.array(timepoints), \
                                                      numpy.array(measurements), \
                                                      1, self.args.p, "upordown")
            mapped = False
            for (tid, tpos), score in tidpos_to_n.items():
                if score >= self.args.n:
                    orientation = "."
                    if tid in tid_to_mkt:
                        #mk: string of test result
                        #m: slope
                        #c: intercept
                        #p: significance
                        result = tid_to_mkt[tid]
                        mapped = True
                        if result[3] < self.args.p and result[1] != 0:
                            orientation = "+" if result[1] > 0 else "-"
                    print(tid, tpos, tpos + 1, qid, score, orientation, sep="\t")
            if mapped:
                num_mapped += 1
        print(
            int(timeit.default_timer() - t0),
            "Mapped", num_mapped, "sequences of", len(query_markers),
            f"({round(100 * num_mapped / len(query_markers), 2)}%)", file=sys.stderr)

    def physlr_map(self):
        """
        Map sequences to a physical map.
        Usage: physlr map TGRAPH.tsv TMARKERS.tsv QMARKERS.tsv... >MAP.bed
        """

        if len(self.args.FILES) < 3:
            exit("physlr map: error: at least three file arguments are required")
        graph_filenames = [self.args.FILES[0]]
        target_filenames = [self.args.FILES[1]]
        query_filenames = self.args.FILES[2:]

        g = self.read_graph(graph_filenames)
        bxtomin = self.read_minimizers(target_filenames)
        query_markers = bxtomin if target_filenames == query_filenames else \
            self.read_minimizers(query_filenames)

        # Index the positions of the markers in the backbone.
        backbones = Physlr.determine_backbones(g)
        markertopos = Physlr.index_markers_in_backbones(backbones, bxtomin)

        # Map the query sequences to the physical map.
        num_mapped = 0
        for qid, markers in progress(query_markers.items()):
            # Map each target position to a query position.
            tidpos_to_qpos = {}
            for qpos, marker in enumerate(markers):
                for tidpos in markertopos.get(marker, ()):
                    tidpos_to_qpos.setdefault(tidpos, []).append(qpos)
            for tidpos, qpos in tidpos_to_qpos.items():
                tidpos_to_qpos[tidpos] = statistics.median_low(qpos)

            # Count the number of markers mapped to each target position.
            tidpos_to_n = Counter(pos for marker in markers for pos in markertopos.get(marker, ()))

            mapped = False
            for (tid, tpos), score in tidpos_to_n.items():
                if score >= self.args.n:
                    mapped = True
                    orientation = Physlr.determine_orientation(
                        tidpos_to_qpos.get((tid, tpos - 1), None),
                        tidpos_to_qpos.get((tid, tpos + 0), None),
                        tidpos_to_qpos.get((tid, tpos + 1), None))
                    print(tid, tpos, tpos + 1, qid, score, orientation, sep="\t")
            if mapped:
                num_mapped += 1
        print(
            int(timeit.default_timer() - t0),
            "Mapped", num_mapped, "sequences of", len(query_markers),
            f"({round(100 * num_mapped / len(query_markers), 2)}%)", file=sys.stderr)

    def physlr_bed_to_path(self):
        """
        Convert a BED file of mappings to scaffolds paths.
        Usage: physlr bed-to-path BED... >PATH
        """

        # Construct a dictionary of query names to their best position on the physical map.
        qnames = {}
        for tname, tstart, _, qname, score, orientation in \
                progress(Physlr.read_bed(self.args.FILES)):
            if score < self.args.n:
                continue
            qnames.setdefault(qname, {}).setdefault(tname, []).append((tstart, orientation))

        # Order and orient the queries on the physical map.
        scaffolds = []
        for qname, targets in progress(qnames.items()):
            tname, positions = max(targets.items(), key=lambda target: len(target[1]))
            tstart = statistics.median_low(x[0] for x in positions)
            orientations = Counter(x[1] for x in positions).most_common(2)
            if len(orientations) == 1 or orientations[0][1] > orientations[1][1]:
                orientation = orientations[0][0]
            else:
                orientation = "."
            scaffolds.append((tname, tstart, orientation, qname))
        scaffolds.sort()
        print(int(timeit.default_timer() - t0), f"Ordered and oriented queries.", file=sys.stderr)

        num_scaffolds = 0
        num_contigs = 0
        prev_tname = None
        for tname, _, orientation, qname in progress(scaffolds):
            num_contigs += 1
            if prev_tname:
                if tname != prev_tname:
                    num_scaffolds += 1
                    print()
                else:
                    print(" ", end="")
            print(qname, orientation, sep="", end="")
            prev_tname = tname
        num_scaffolds += 1
        print()
        print(
            int(timeit.default_timer() - t0),
            f"Wrote {num_contigs} contigs in {num_scaffolds} scaffolds.",
            file=sys.stderr)

    def physlr_path_to_fasta(self):
        """
        Produce sequences in FASTA format from paths.
        Usage: physlr path-to-fasta FASTA PATH... >FASTA
        """
        if len(self.args.FILES) < 2:
            exit("physlr path-to-fasta: error: at least two file arguments are required")
        fasta_filenames = self.args.FILES[0:1]
        path_filenames = self.args.FILES[1:]
        seqs = Physlr.read_fastas(fasta_filenames)
        paths = Physlr.read_path(path_filenames)

        num_scaffolds = 0
        num_contigs = 0
        num_bases = 0
        for path in progress(paths):
            # Remove unoriented sequences.
            path = [name for name in path if name[-1] != "."]
            if not path:
                continue

            seq = "NNNNNNNNNN".join(Physlr.get_oriented_sequence(seqs, name) for name in path)
            if len(seq) < self.args.min_length:
                continue
            print(f">{str(num_scaffolds).zfill(7)} LN:i:{len(seq)} xn:i:{len(path)}\n{seq}")
            num_scaffolds += 1
            num_contigs += len(path)
            num_bases += len(seq)
        print(
            int(timeit.default_timer() - t0),
            f"Wrote {num_bases} bases in {num_contigs} contigs in {num_scaffolds} scaffolds.",
            file=sys.stderr)

    def physlr_filter_bed(self):
        """
        Select records from a BED file in the order given by a .path file.
        Usage: physlr filter-bed BED PATH... >BED
        """
        if len(self.args.FILES) < 2:
            exit("physlr filter-bed: error: at least two file arguments are required")
        bed_filenames = [self.args.FILES[0]]
        path_filenames = self.args.FILES[1:]

        num_beds = 0
        bxtobeds = {}
        for bed in progress(Physlr.read_bed(bed_filenames)):
            num_beds += 1
            qname = bed[3]
            bxtobeds.setdefault(qname, []).append(bed)
        print(
            int(timeit.default_timer() - t0),
            "Indexed", len(bxtobeds), "barcodes in", num_beds, "BED records", file=sys.stderr)

        num_paths = 0
        num_barcodes = 0
        num_beds = 0
        num_too_small = 0
        num_missing = 0
        for i, path in enumerate(progress(Physlr.read_path(path_filenames))):
            path = [x for subpath in path for x in subpath.strip("()").split(",")]
            if len(path) < self.args.min_component_size:
                num_too_small += 1
                continue
            num_paths += 1
            if i != 0:
                print("NA\tNA\tNA\tNA\tNA\tNA")
            for u in path:
                num_barcodes += 1
                bx = u
                if not self.args.molecule_bed:
                    bx = u.split("_", 1)[0]
                if bx not in bxtobeds:
                    num_missing += 1
                    if self.args.verbose >= 1:
                        print("warning:", bx, "not found in BED", *bed_filenames, file=sys.stderr)
                    continue
                beds = bxtobeds[bx]
                num_beds += len(beds)
                for bed in beds:
                    print(*bed, sep="\t")
        print(
            int(timeit.default_timer() - t0),
            "Skipped", num_too_small, "paths shorter than",
            self.args.min_component_size, "vertices", file=sys.stderr)
        if num_missing > 0:
            print(
                int(timeit.default_timer() - t0),
                "Missing", num_missing, "barcodes of", num_barcodes,
                f"({round(100 * num_missing / num_barcodes, 2)}%)", file=sys.stderr)
        print(
            int(timeit.default_timer() - t0),
            "Wrote", num_beds, "BED records in", num_barcodes, "barcodes in", num_paths, "paths",
            file=sys.stderr)

    @staticmethod
    def parse_arguments():
        "Parse the command line arguments."
        argparser = argparse.ArgumentParser()
        argparser.add_argument(
            "-t", "--threads", action="store", dest="threads", type=int,
            default=min(16, os.cpu_count()),
            help="number of threads [16 or number of CPU]")
        argparser.add_argument(
            "-k", "--k", action="store", type=int,
            help="size of a k-mer (bp)")
        argparser.add_argument(
            "-w", "--window", action="store", dest="w", type=int,
            help="number of k-mers in a window of size k + w - 1 bp")
        argparser.add_argument(
            "-c", "--coef", action="store", dest="coef", type=float, default=1.5,
            help="ignore markers that occur in Q3+c*(Q3-Q1) or more barcodes [0]")
        argparser.add_argument(
            "-C", "--max-count", action="store", dest="C", type=int,
            help="ignore markers that occur in C or more barcodes [None]")
        argparser.add_argument(
            "-M", "--max-molecules", action="store", dest="M", type=int,
            help="remove barcodes with M or more molecules [None]")
        argparser.add_argument(
            "-n", "--min-n", action="store", dest="n", type=int, default=0,
            help="remove edges with fewer than n shared markers [0]")
        argparser.add_argument(
            "-N", "--max-n", action="store", dest="N", type=int, default=None,
            help="remove edges with at least N shared markers [None]")
        argparser.add_argument(
            "--min-length", action="store", dest="min_length", type=int, default=0,
            help="remove sequences with length less than N bp [0]")
        argparser.add_argument(
            "--min-component-size", action="store", dest="min_component_size", type=int, default=0,
            help="remove components with fewer than N vertices [0]")
        argparser.add_argument(
            "--molecule-bed", action="store", dest="molecule_bed", type=int, default=0,
            help="Retain molecule splits in filtered BED (0 or 1) [0]")
        argparser.add_argument(
            "--molecules-bx-only", action="store", dest="molecules_bx_only", type=int, default=1,
            help="Only print reads with barcodes that have been split to molecules (0 or 1) [1]")
        argparser.add_argument(
            "-v", "--vertices", action="store", dest="v",
            help="list of vertices [None]")
        argparser.add_argument(
            "-V", "--exclude-vertices", action="store", dest="exclude_vertices", default="",
            help="list of vertices to exclude [None]")
        argparser.add_argument(
            "--exclude-source", action="store", dest="exclude_source", type=int, default=1,
            help="exclude the barcode itself from the subgraph (0 or 1) [1]")
        argparser.add_argument(
            "-d", "--distance", action="store", dest="d", type=int, default=0,
            help="include vertices within d edges away [0]")
        argparser.add_argument(
            "-o", "--output", action="store", dest="output",
            help="the output file or directory")
        argparser.add_argument(
            "-O", "--output-format", action="store", dest="graph_format", default="tsv",
            help="the output graph file format [tsv]")
        argparser.add_argument(
            "-p", "--min_p_val", action="store", dest="p", type=float, default=0.01,
            help="Minimum significance threshold (FPR) for Mann-Kendall Test")
        argparser.add_argument(
            "--verbose", action="store", dest="verbose", type=int, default="0",
            help="the level of verbosity [0]")
        argparser.add_argument(
            "--version", action="version", version="physlr 0.1.0")
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
        method_name = "physlr_" + self.args.command.replace("-", "_")
        if not hasattr(Physlr, method_name):
            print("physlr: error: unrecognized command:", self.args.command, file=sys.stderr)
            exit(1)
        getattr(Physlr, method_name)(self)

def main():
    "Run Physlr."
    Physlr().main()

if __name__ == "__main__":
    main()
