#!/usr/bin/env python3
"""
Physlr: Physical Mapping of Linked Reads
"""

import argparse
import itertools
import math
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

# The time at which execution started.
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
    if Physlr.args.verbose < 2:
        return iterator
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
        "Read BED files. Also able to read PAF files (columns 1, 5, 6, 8, 9, 10)."
        bed = []
        for filename in filenames:
            print(int(timeit.default_timer() - t0), "Reading", filename, file=sys.stderr)
            with open(filename) as fin:
                if Physlr.args.verbose >= 2:
                    progressbar = progress_bar_for_file(fin)
                for line in fin:
                    if Physlr.args.verbose >= 2:
                        progressbar.update(len(line))
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) == 5:
                        # BED format, five columns
                        tname, tstart, tend, qname, score = fields[0:5]
                        orientation = "."
                    elif len(fields) == 6:
                        # BED format, six columns
                        tname, tstart, tend, qname, score, orientation = fields[0:6]
                    elif len(fields) >= 12:
                        # PAF format
                        qname, _, _, _, orientation, tname, _, tstart, tend, score = fields[0:10]
                    else:
                        print("physlr: expected 5 or 6 BED fields, or 12 or more PAF fields:",
                              line, file=sys.stderr)
                        sys.exit(1)
                    bed.append((tname, int(tstart), int(tend), qname, int(score), orientation))
                if Physlr.args.verbose >= 2:
                    progressbar.close()
            print(int(timeit.default_timer() - t0), "Read", filename, file=sys.stderr)
        return bed

    @staticmethod
    def read_paf(filenames):
        """Read PAF files."""
        paf = []
        for filename in filenames:
            print(int(timeit.default_timer() - t0), "Reading", filename, file=sys.stderr)
            with open(filename) as fin:
                if Physlr.args.verbose >= 2:
                    progressbar = progress_bar_for_file(fin)
                for line in fin:
                    if Physlr.args.verbose >= 2:
                        progressbar.update(len(line))
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 12:
                        print("physlr: expected 12 or more PAF fields:", line, file=sys.stderr)
                        sys.exit(1)
                    qname, qlength, qstart, qend, orientation, \
                        tname, tlength, tstart, tend, score, length, mapq = fields[0:12]
                    paf.append((
                        qname, int(qlength), int(qstart), int(qend), orientation,
                        tname, int(tlength), int(tstart), int(tend),
                        int(score), int(length), int(mapq)))
                if Physlr.args.verbose >= 2:
                    progressbar.close()
            print(int(timeit.default_timer() - t0), "Read", filename, file=sys.stderr)
        return paf

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
    def read_paths(filenames):
        "Read path files."
        paths = []
        for filename in filenames:
            print(int(timeit.default_timer() - t0), "Reading", filename, file=sys.stderr)
            with open(filename) as fin:
                if Physlr.args.verbose >= 2:
                    progressbar = progress_bar_for_file(fin)
                for line in fin:
                    if Physlr.args.verbose >= 2:
                        progressbar.update(len(line))
                    paths.append(line.split())
                if Physlr.args.verbose >= 2:
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
            sys.exit(1)

    @staticmethod
    def read_tsv(g, filename):
        "Read a graph in TSV format."
        with open(filename) as fin:
            if Physlr.args.verbose >= 2:
                progressbar = progress_bar_for_file(fin)
            line = fin.readline()
            if Physlr.args.verbose >= 2:
                progressbar.update(len(line))
            if line not in ["U\tn\n", "U\tn\tm\n"]:
                print("Unexpected header:", line, file=sys.stderr)
                sys.exit(1)
            reading_vertices = True
            for line in fin:
                if Physlr.args.verbose >= 2:
                    progressbar.update(len(line))
                if line == "\n":
                    line = fin.readline()
                    if Physlr.args.verbose >= 2:
                        progressbar.update(len(line))
                    if line == "U\tV\tn\n":
                        reading_vertices = False
                    else:
                        print("Unexpected header:", line, file=sys.stderr)
                        sys.exit(1)
                    line = fin.readline()
                    if not line: # a graph with no edges
                        print("Warning: input graph has no edges, input file:\n",
                              filename, file=sys.stderr)
                        break
                    if Physlr.args.verbose >= 2:
                        progressbar.update(len(line))
                xs = line.split()
                if reading_vertices:
                    if len(xs) == 2:
                        g.add_node(xs[0], n=int(xs[1]))
                    elif len(xs) == 3:
                        g.add_node(xs[0], n=int(xs[1]), m=int(xs[2]))
                    else:
                        print("Unexpected row:", line, file=sys.stderr)
                        sys.exit(1)
                else:
                    if len(xs) == 3:
                        g.add_edge(xs[0], xs[1], n=int(xs[2]))
                    else:
                        print("Unexpected row:", line, file=sys.stderr)
                        sys.exit(1)
        if Physlr.args.verbose >= 2:
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
            "common minimizers of", g.number_of_edges(),
            f"({round(100 * len(edges) / g.number_of_edges(), 2)}%)", file=sys.stderr)
        g.remove_edges_from(edges)

        num_singletons = Physlr.remove_singletons(g)
        print(
            int(timeit.default_timer() - t0),
            "Removed", num_singletons, "isolated vertices.", file=sys.stderr)

    @staticmethod
    def keep_best_edges(g, bestn):
        """Keep the best edges of each vertex."""
        if bestn is None:
            return
        num_edges = g.number_of_edges()
        num_removed = 0
        us = list(g.nodes)
        us.sort(key=g.degree, reverse=True)
        for u in progress(us):
            vs = list(g[u])
            vs.sort(key=lambda v, u=u: g[u][v]["n"], reverse=True)
            for v in vs[bestn:]:
                g.remove_edge(u, v)
                num_removed += 1
        print(
            int(timeit.default_timer() - t0),
            "Removed", num_removed, "non-best edges of", g.number_of_edges(),
            f"({round(100 * num_removed / num_edges, 2)}%)", file=sys.stderr)

        num_singletons = Physlr.remove_singletons(g)
        print(
            int(timeit.default_timer() - t0),
            "Removed", num_singletons, "isolated vertices.", file=sys.stderr)

    @staticmethod
    def remove_unsupported_edges(g):
        "Remove edges with no common neighbours."
        unsupported = [(u, v) for u, v in progress(g.edges())
                       if next(nx.common_neighbors(g, u, v), None) is None]
        print(
            int(timeit.default_timer() - t0),
            "Removed", len(unsupported), "unsupported edges of", g.number_of_edges(),
            f"({round(100 * len(unsupported) / g.number_of_edges(), 2)}%)", file=sys.stderr)
        g.remove_edges_from(unsupported)

        num_singletons = Physlr.remove_singletons(g)
        print(
            int(timeit.default_timer() - t0),
            "Removed", num_singletons, "isolated vertices.", file=sys.stderr)

    @staticmethod
    def read_minimizers(filenames):
        "Read minimizers in TSV format. Returns unordered set."
        bxtomxs = {}
        for filename in filenames:
            print(int(timeit.default_timer() - t0), "Reading", filename, file=sys.stderr)
            with open(filename) as fin:
                if Physlr.args.verbose >= 2:
                    progressbar = progress_bar_for_file(fin)
                for line in fin:
                    if Physlr.args.verbose >= 2:
                        progressbar.update(len(line))
                    fields = line.split(None, 1)
                    if len(fields) < 2:
                        continue
                    bx = fields[0]
                    if bx not in bxtomxs:
                        bxtomxs[bx] = set()
                    bxtomxs[bx].update(int(mx.split(":", 1)[0]) for mx in fields[1].split())
                if Physlr.args.verbose >= 2:
                    progressbar.close()
            print(int(timeit.default_timer() - t0), "Read", filename, file=sys.stderr)
        return bxtomxs

    @staticmethod
    def read_minimizers_list(filenames):
        "Read minimizers in TSV format. Returns ordered list."
        bxtomxs = {}
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
                    if bx in bxtomxs:
                        print("Error: Expected single id per in file", file=sys.stderr)
                        sys.exit(1)
                    bxtomxs[bx] = [int(mx.split(":", 1)[0]) for mx in fields[1].split()]
                progressbar.close()
            print(int(timeit.default_timer() - t0), "Read", filename, file=sys.stderr)
        return bxtomxs

    @staticmethod
    def read_minimizers_pos(filenames):
        "Read minimizers with positions. Returns an ordered list of (position, hash) tuples."
        nametomxs = {}
        for filename in filenames:
            print(int(timeit.default_timer() - t0), "Reading", filename, file=sys.stderr)
            with open(filename) as fin:
                progressbar = progress_bar_for_file(fin)
                for line in fin:
                    progressbar.update(len(line))
                    fields = line.split("\t", 1)
                    if len(fields) < 2:
                        continue
                    name = fields[0]
                    if name in nametomxs:
                        print("Error: Duplicate sequence name:", name, "in", filename, \
                            file=sys.stderr)
                        sys.exit(1)
                    posmxs = []
                    for mx_pos in fields[1].split():
                        if ":" not in mx_pos:
                            print("Error: Minimizers do not include positions:", filename, \
                                file=sys.stderr)
                            sys.exit(1)
                        mx, pos = mx_pos.split(":", 1)
                        posmxs.append((int(pos), int(mx)))
                    nametomxs[name] = posmxs
                progressbar.close()
            print(int(timeit.default_timer() - t0), "Read", filename, file=sys.stderr)
        return nametomxs

    @staticmethod
    def count_molecules_per_bx(moltomxs):
        "Iterate over minimizers dictionary, track # molecules per barcode"
        mol_counts = Counter()
        bx_match = re.compile(r'^(\S+)_(\d+)$')
        for bx_mol in moltomxs:
            bx_mol_match = re.search(bx_match, bx_mol)
            if bx_mol_match:
                mol_counts[bx_mol_match.group(1)] = max(mol_counts[bx_mol_match.group(1)],
                                                        int(bx_mol_match.group(2)) + 1)
        return mol_counts

    @staticmethod
    def construct_minimizers_to_barcodes(bxtomxs):
        "Construct a dictionary of minimizers to barcodes."
        mxtobxs = {}
        for bx, mxs in progress(bxtomxs.items()):
            for mx in mxs:
                if mx not in mxtobxs:
                    mxtobxs[mx] = set()
                mxtobxs[mx].add(bx)
        print(
            int(timeit.default_timer() - t0),
            "Indexed", len(mxtobxs), "minimizers", file=sys.stderr)
        return mxtobxs

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
    def detect_junctions_of_tree(gcomponent, minor_branch_size):
        """"
        Detect junctions in the tree, with at least 3 branches larger than minor_branch_size.
        """
        branch_lengths = Physlr.measure_branch_length(gcomponent)
        candidate_junctions = [node
                               for node in list(gcomponent.nodes)
                               if gcomponent.degree(node) >= 3]
        junctions = []
        for candidate_junction in candidate_junctions:
            lengths = [branch_lengths[(candidate_junction, neighbor)]
                       for neighbor in gcomponent.neighbors(candidate_junction)]
            lengths.sort()
            if lengths[-3] >= minor_branch_size:
                # If the 3rd largest branch is considerably large, this node is a junction.
                junctions.append(candidate_junction)
        return junctions

    @staticmethod
    def split_junctions_of_tree(prune_junctions, gin, keep_largest=0):
        """"
        Detect and split junctions of trees, with at least 3 branches larger than prune_junctions.
        For each junction, keep the two incident edges with the largest weight, and remove the rest.
        """
        if prune_junctions == 0:
            return gin, 0
        junctions = Physlr.detect_junctions_of_tree(gin, prune_junctions)
        g = gin.copy()
        for junction in junctions:
            edges = list(g.edges(junction, data="n"))
            edges.sort(key=lambda e: e[2], reverse=True)
            if Physlr.args.verbose >= 3:
                print("Junction:", junction, "Edges:", *edges, file=sys.stderr)
            # Keep the two incident edges with the largest weight, and remove the rest.
            if keep_largest:
                for u, v, _ in edges[2:]:
                    g.remove_edge(u, v)
            else:
                for u, v, _ in edges:
                    g.remove_edge(u, v)
        return (g, len(junctions))

    @staticmethod
    def determine_backbones_of_trees(g, prune_junctions):
        """"
        Determine backbones of the MSTs.
        Resolve junctions of >=3 branches of size >= prune_junctions.
        """
        paths = []
        removed_count = 0
        for component in nx.connected_components(g):
            gcomponents, rem_count =\
                Physlr.split_junctions_of_tree(prune_junctions, g.subgraph(component))
            removed_count += rem_count
            for subcomponent in nx.connected_components(gcomponents):
                gsubcomponent = g.subgraph(subcomponent)
                u, v, _ = Physlr.diameter_of_tree(gsubcomponent, weight="n")
                path = nx.shortest_path(gsubcomponent, u, v, weight="n")
                paths.append(path)
        paths.sort(key=len, reverse=True)
        print(
            int(timeit.default_timer() - t0),
            "Pruned", removed_count, "junctions in the backbone",
            file=sys.stderr)
        return paths

    @staticmethod
    def identify_chimera(g, backbones, distance, min_support):
        "Identify poorly supported (possibly chimeric) barcodes."
        print(
            int(timeit.default_timer() - t0),
            "Identifying poorly supported barcodes.", file=sys.stderr, flush=True)
        fout = open(Physlr.args.output, "w") if Physlr.args.output else None
        if fout:
            print("Tname", "Pos", "U", "V", "W", "Overlap", "Depth", "Support", sep="\t", file=fout)
        chimera = []
        for tname, backbone in enumerate(backbones):
            # Skip small components. Removing poorly supported barcodes can take many iterations.
            if len(backbone) < Physlr.args.min_path_size:
                continue

            for i, v in enumerate(backbone):
                us = set(backbone[max(0, i - distance) : i])
                ws = set(backbone[i + 1 : i + 1 + distance])
                if not us or not ws:
                    continue

                # Count the number molecules that directly overlap.
                overlappers = set(u for e in nx.edge_boundary(g, us, ws) for u in e)
                overlapping = len(overlappers)
                if overlapping > 0 and not fout:
                    continue

                # Count the number molecules that span the position.
                uneighbors = set()
                for u in us:
                    uneighbors.update(g.neighbors(u))
                wneighbors = set()
                for w in ws:
                    wneighbors.update(g.neighbors(w))
                uneighbors.difference_update(us, [v], ws)
                wneighbors.difference_update(us, [v], ws)
                spanners = uneighbors & wneighbors
                depth = overlapping + len(spanners)
                if depth >= min_support and not fout:
                    continue

                # Count the number of pairs of molecules that span the position.
                uneighbors -= spanners
                wneighbors -= spanners
                uboundary = set()
                wboundary = set()
                for u, w in nx.edge_boundary(g, uneighbors, wneighbors):
                    uboundary.add(u)
                    wboundary.add(w)
                support = depth + min(len(uboundary), len(wboundary))

                if fout:
                    print(tname, i, len(us), v, len(ws), overlapping, depth, support,
                          sep="\t", file=fout)
                if overlapping == 0 and depth < min_support and support < min_support:
                    chimera.append(v)
                    chimera += spanners
        if fout:
            fout.close()
        if Physlr.args.verbose >= 1:
            print("Chimera:", *chimera, file=sys.stderr)
        print(
            int(timeit.default_timer() - t0),
            "Identified", len(chimera), "poorly supported barcodes.", file=sys.stderr, flush=True)
        return chimera

    @staticmethod
    def determine_backbones_and_remove_chimera(g):
        """Remove chimeric vertices iteratively until none remain."""
        print(int(timeit.default_timer() - t0),
              "Determining the backbone-induced subgraph.", file=sys.stderr)
        if Physlr.args.d == 0:
            Physlr.args.d = 2
        if Physlr.args.s == 0:
            return Physlr.determine_backbones(g)

        iterations = 0
        nchimera = 0
        chimera = True
        while chimera:
            backbones = Physlr.determine_backbones(g)
            chimera = Physlr.identify_chimera(
                g, backbones, distance=Physlr.args.d, min_support=Physlr.args.s)
            g.remove_nodes_from(chimera)
            if chimera:
                iterations += 1
                nchimera += len(chimera)
        print(
            int(timeit.default_timer() - t0),
            "Removed", nchimera, "chimeric barcodes in", iterations, "iterations.",
            file=sys.stderr, flush=True)
        return backbones

    @staticmethod
    def physlr_cut_chimera(g):
        "Cut chimeric paths to correct misassemblies."
        if Physlr.args.d == 0:
            Physlr.args.d = 2
        g = Physlr.read_graph([Physlr.args.FILES[0]])
        backbones = Physlr.read_paths([Physlr.args.FILES[1]])
        chimera = Physlr.identify_chimera(
            g, backbones, distance=Physlr.args.d, min_support=Physlr.args.s)
        sep = ""
        for backbone in backbones:
            for u in backbone:
                if u in chimera:
                    sep = "\n"
                else:
                    print(sep, u, sep="", end="")
                    sep = " "
            sep = "\n"
        if sep:
            print()

    @staticmethod
    def determine_pruned_mst(g):
        """Return the pruned maximum spanning tree of the graph."""
        # having tested kruskal and prim, we found the former is faster in our case
        gmst = nx.maximum_spanning_tree(g, algorithm="kruskal", weight="n")
        print(
            int(timeit.default_timer() - t0),
            "Determined the maximum spanning tree.", file=sys.stderr, flush=True)
        Physlr.prune_mst(gmst, Physlr.args.prune_branches)
        return gmst

    @staticmethod
    def identify_contiguous_paths(g):
        """Return the contiguous paths and junctions of the graph."""
        g = g.copy()
        junctions = [u for u, deg in g.degree() if deg >= 3]
        g.remove_nodes_from(junctions)
        paths = list(nx.connected_components(g))
        paths.sort(key=len, reverse=True)
        return paths, junctions

    @staticmethod
    def identify_bridges(g, bridge_length):
        """Return the bridges of the graph"""
        paths, junctions = Physlr.identify_contiguous_paths(g)
        bridges = [path for path in paths if len(path) < bridge_length
                   and all(g.degree(u) == 2 for u in path)]
        bridges += g.subgraph(junctions).edges
        return bridges

    @staticmethod
    def remove_bridges_once(g, bridge_length):
        """
        Remove bridges from this graph.
        Bridges are non-blunt contigs shorter than a specified length.
        """
        gmst = Physlr.determine_pruned_mst(g)
        bridges = Physlr.identify_bridges(gmst, bridge_length)
        for bridge in bridges:
            g.remove_nodes_from(bridge)
        if Physlr.args.verbose >= 3:
            print("Bridges:", *bridges, file=sys.stderr, flush=True)
        print(
            int(timeit.default_timer() - t0),
            "Removed", sum(map(len, bridges)), "vertices in", len(bridges), "bridges",
            "shorter than", bridge_length, "vertices.",
            file=sys.stderr, flush=True)
        return len(bridges)

    @staticmethod
    def remove_bridges(g, bridge_length):
        """Remove bridges iteratively until none remain."""
        iterations = 0
        total_nbridges = 0
        while True:
            nbridges = Physlr.remove_bridges_once(g, bridge_length)
            if nbridges == 0:
                break
            iterations += 1
            total_nbridges += nbridges
        print(
            int(timeit.default_timer() - t0),
            "Removed", total_nbridges, "bridges in", iterations, "iterations.",
            file=sys.stderr, flush=True)

    @staticmethod
    def determine_backbones(g):
        "Determine the backbones of the graph."
        g = g.copy()
        if Physlr.args.prune_bridges > 0:
            Physlr.remove_bridges(g, Physlr.args.prune_bridges)
        backbones = []
        gmst = Physlr.determine_pruned_mst(g)
        while not nx.is_empty(gmst):
            paths = Physlr.determine_backbones_of_trees(gmst, Physlr.args.prune_junctions)
            backbones += (path for path in paths if len(path) >= Physlr.args.prune_branches)
            for path in paths:
                gmst.remove_nodes_from(path)
        backbones.sort(key=len, reverse=True)
        print(
            int(timeit.default_timer() - t0),
            "Assembled", sum(map(len, backbones)), "molecules in", len(backbones), "paths.",
            file=sys.stderr, flush=True)
        return backbones

    @staticmethod
    def wrap_up_messages_and_pass(mst, messages, sender, receiver):
        """
        Wrap up all incoming messages to this node (sender) except the one from receiver;
        and set (pass) the message from sender to receiver.
        """
        if mst.degree(sender) == 1:
            length = 0
        else:
            length = max(messages[(sender, u)] for u in mst.neighbors(sender) if u != receiver)
        messages[(receiver, sender)] = 1 + length

    @staticmethod
    def measure_branch_length(mst):
        """
        Measure the lengths of branches.
        The branch length of an edge (u,v) is the longest path from (u,v) to a leaf vertex.
        """
        dfs = list(nx.dfs_edges(mst))
        if not dfs:
            return {}
        stack = [dfs[0][0]]
        branch_lengths = {}
        # Gather
        for edge in dfs:
            while stack[-1] != edge[0]:
                Physlr.wrap_up_messages_and_pass(mst, branch_lengths, stack.pop(), stack[-1])
            stack.append(edge[1])
        while len(stack) != 1:
            Physlr.wrap_up_messages_and_pass(mst, branch_lengths, stack.pop(), stack[-1])
        # Distribute
        for edge in dfs:
            Physlr.wrap_up_messages_and_pass(mst, branch_lengths, edge[0], edge[1])
        return branch_lengths

    @staticmethod
    def prune_mst_once(g, branch_size):
        """"
        Prune branches smaller than branch_size.
        Return the number of pruned vertices.
        """
        g0 = g.copy()
        for component in nx.connected_components(g0):
            branch_lengths = Physlr.measure_branch_length(g0.subgraph(component))
            edges = [(u, v) for (u, v), length in branch_lengths.items()
                     if g0.degree(u) >= 3 and length < branch_size]
            g.remove_edges_from(edges)
            for _, v in edges:
                if g.has_node(v):
                    g.remove_nodes_from(nx.node_connected_component(g, v))
        n = g0.number_of_nodes()
        pruned = n - g.number_of_nodes()
        if Physlr.args.verbose >= 3:
            print(
                int(timeit.default_timer() - t0),
                "Pruned", pruned, "vertices of", n, f"({round(100 * pruned / n, 2)}%)",
                file=sys.stderr, flush=True)
        return pruned

    @staticmethod
    def prune_mst(g, branch_size):
        """"
        Iteratively prune branches smaller than branch_size.
        First prune branches shorter than branch_size, then iteratively prune very small branches.
        """
        n = g.number_of_nodes()
        iterations = 0
        total_pruned = 0
        pruned = True
        while pruned and not nx.is_empty(g):
            pruned = Physlr.prune_mst_once(g, branch_size)
            total_pruned += pruned
            if pruned:
                iterations += 1
        print(
            int(timeit.default_timer() - t0),
            "Pruned", total_pruned, "vertices of", n, f"({round(100 * total_pruned / n, 2)}%)",
            "in", iterations, "iterations.",
            file=sys.stderr, flush=True)

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

    @staticmethod
    def remove_small_components(g, min_component_size):
        "Remove comonents smaller than min_component_size"
        if min_component_size < 2:
            return
        ncomponents, nvertices = 0, 0
        vertices = set()
        for component in nx.connected_components(g):
            if len(component) < min_component_size:
                vertices.update(component)
                ncomponents += 1
                nvertices += len(component)
        g.remove_nodes_from(vertices)
        print(
            int(timeit.default_timer() - t0),
            "Removed", nvertices, "vertices in", ncomponents, "components",
            "with fewer than", min_component_size, "vertices in a component.",
            file=sys.stderr)

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
        Physlr.remove_small_components(g, self.args.min_component_size)
        self.write_graph(g, sys.stdout, self.args.graph_format)

    def physlr_best_edges(self):
        """Keep the best edges of each vertex."""
        g = self.read_graph(self.args.FILES)
        Physlr.filter_edges(g, self.args.n)
        Physlr.keep_best_edges(g, self.args.bestn)
        self.write_graph(g, sys.stdout, self.args.graph_format)
        print(int(timeit.default_timer() - t0), "Wrote graph", file=sys.stderr)

    def physlr_flesh_backbone(self):
        "Flesh out the barcodes in the backbone paths"
        g = self.read_graph([self.args.FILES[0]])
        backbones_raw = self.read_paths([self.args.FILES[1]])
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
                # Put R or L of node with most shared minimizers?
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
            sys.exit("physlr subgraph: error: Only -d0 and -d1 are currently supported.")
        vertices = set(self.args.v.split(","))
        if not self.args.exclude_vertices and len(vertices) == 1 and self.args.d == 1:
            self.args.exclude_vertices = self.args.v
        exclude_vertices = set(self.args.exclude_vertices.split(","))
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
            sys.exit("physlr subgraphs: missing parameter: --output is need but not provided.")
        if self.args.d not in (0, 1):
            sys.exit("physlr subgraphs: error: Only -d0 and -d1 are currently supported.")
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

    def physlr_count_minimizers(self):
        "Count the frequency of each minimizer."
        bxtomxs = self.read_minimizers(self.args.FILES)
        mx_counts = Counter(mx for mxs in progress(bxtomxs.values()) for mx in mxs)
        print(
            int(timeit.default_timer() - t0),
            "Counted", len(mx_counts), "minimizers", file=sys.stderr)

        count = 0
        print("Depth")
        for n in sorted(progress(mx_counts.values())):
            if n >= self.args.c:
                count += 1
                print(n)
        print(int(timeit.default_timer() - t0), "Wrote", count, "minimizers", file=sys.stderr)

    def physlr_intersect(self):
        "Print the minimizers in the intersection of each pair of barcodes."
        if self.args.n == 0:
            self.args.n = 1
        bxtomxs = self.read_minimizers(self.args.FILES)
        mxtobxs = self.construct_minimizers_to_barcodes(bxtomxs)
        if self.args.v:
            pairs = itertools.combinations(self.args.v.split(","))
        else:
            pairs = {(u, v) for bxs in mxtobxs.values() for u, v in itertools.combinations(bxs, 2)}
        for u, v in pairs:
            common = bxtomxs[u] & bxtomxs[v]
            if len(common) >= self.args.n:
                print(u, v, "", sep="\t", end="")
                print(*common)

    @staticmethod
    def remove_singleton_minimizers(bxtomxs):
        """
        Remove minimizers that occur only once.
        Return the counts of minimizers after removing singletons.
        """
        mx_counts = Counter(mx for mxs in progress(bxtomxs.values()) for mx in mxs)
        print(
            int(timeit.default_timer() - t0),
            "Counted", len(mx_counts), "minimizers", file=sys.stderr)

        singletons = {mx for mx, n in progress(mx_counts.items()) if n < 2}
        for mxs in progress(bxtomxs.values()):
            mxs -= singletons
        print(
            int(timeit.default_timer() - t0),
            "Removed", len(singletons), "minimizers that occur only once of", len(mx_counts),
            f"({round(100 * len(singletons) / len(mx_counts), 2)}%)", file=sys.stderr)
        for mx in singletons:
            del mx_counts[mx]
        return mx_counts

    def physlr_filter_barcodes(self):
        """
        Filter barcodes by number of minimizers.
        Read a TSV file of barcodes to minimizers.
        Remove minimizers that occur only once.
        Remove barkers with too few or too many minimizers.
        Write a TSV file of barcodes to minimizers.
        """
        bxtomxs = self.read_minimizers(self.args.FILES)
        Physlr.remove_singleton_minimizers(bxtomxs)

        q0, q1, q2, q3, q4 = quantile(
            [0, 0.25, 0.5, 0.75, 1], (len(mxs) for mxs in bxtomxs.values()))
        low_whisker = int(q1 - self.args.coef * (q3 - q1))
        high_whisker = int(q3 + self.args.coef * (q3 - q1))
        if self.args.n == 0:
            self.args.n = max(q0, low_whisker)
        if self.args.N is None:
            self.args.N = min(1 + q4, high_whisker)

        print(
            int(timeit.default_timer() - t0), " Counted minimizers per barcode\n",
            f"    Markers per barcode: Q0={q0} Q1={q1} Q2={q2} Q3={q3} Q4={q4} Q3-Q1={q3 - q1}\n",
            f"    Q3-{self.args.coef}*(Q3-Q1)={low_whisker} n={self.args.n}\n",
            f"    Q3+{self.args.coef}*(Q3-Q1)={high_whisker} N={self.args.N}",
            sep="", file=sys.stderr)

        too_few, too_many = 0, 0
        for bx, mxs in progress(bxtomxs.items()):
            if len(mxs) < self.args.n:
                too_few += 1
            elif len(mxs) >= self.args.N:
                too_many += 1
            else:
                print(bx, "\t", sep="", end="")
                print(*mxs)
        print(
            "    Discarded", too_few, "barcodes with too few minimizers of", len(bxtomxs),
            f"({round(100 * too_few / len(bxtomxs), 2)}%)", file=sys.stderr)
        print(
            "    Discarded", too_many, "barcodes with too many minimizers of", len(bxtomxs),
            f"({round(100 * too_many / len(bxtomxs), 2)}%)", file=sys.stderr)
        print(
            int(timeit.default_timer() - t0),
            "Wrote", len(bxtomxs) - too_few - too_many, "barcodes", file=sys.stderr)

    def physlr_filter_minimizers(self):
        "Filter minimizers by depth of coverage. Remove repetitive minimizers."
        bxtomxs = self.read_minimizers(self.args.FILES)
        mx_counts = Physlr.remove_singleton_minimizers(bxtomxs)

        # Identify frequent minimizers.
        q1, q2, q3 = quantile([0.25, 0.5, 0.75], mx_counts.values())
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

        # Remove frequent minimizers.
        frequent_mxs = {mx for mx, count in mx_counts.items() if count >= self.args.C}
        num_empty_barcodes = 0
        for bx, mxs in progress(bxtomxs.items()):
            mxs -= frequent_mxs
            if not mxs:
                num_empty_barcodes += 1
                continue
            print(bx, "\t", sep="", end="")
            print(*mxs)

        print(
            int(timeit.default_timer() - t0),
            "Removed", len(frequent_mxs), "most frequent minimizers of", len(mx_counts),
            f"({round(100 * len(frequent_mxs) / len(mx_counts), 2)}%)", file=sys.stderr)
        print(
            int(timeit.default_timer() - t0),
            "Removed", num_empty_barcodes, "empty barcodes of", len(bxtomxs),
            f"({round(100 * num_empty_barcodes / len(bxtomxs), 2)}%)", file=sys.stderr)
        print(
            int(timeit.default_timer() - t0),
            "Wrote", len(bxtomxs) - num_empty_barcodes, "barcodes", file=sys.stderr)

    def remove_repetitive_minimizers(self, bxtomxs, mxtobxs):
        "Remove repetitive minimizers."

        # Remove minimizers that occur only once.
        num_mxs = len(mxtobxs)
        singletons = {mx for mx, bxs in progress(mxtobxs.items()) if len(bxs) < 2}
        for mx in singletons:
            del mxtobxs[mx]
        for mxs in progress(bxtomxs.values()):
            mxs -= singletons
        print(
            int(timeit.default_timer() - t0),
            "Removed", len(singletons), "minimizers that occur only once of", num_mxs,
            f"({round(100 * len(singletons) / num_mxs, 2)}%)", file=sys.stderr)

        # Identify repetitive minimizers.
        q1, q2, q3 = quantile([0.25, 0.5, 0.75], (len(bxs) for bxs in mxtobxs.values()))
        whisker = int(q3 + self.args.coef * (q3 - q1))
        if self.args.C is None:
            self.args.C = whisker
        print(
            int(timeit.default_timer() - t0),
            " Minimizer frequency: Q1=", q1, " Q2=", q2, " Q3=", q3,
            " Q3+", self.args.coef, "*(Q3-Q1)=", whisker,
            " C=", self.args.C, sep="", file=sys.stderr)

        # Remove frequent (likely repetitive) minimizers.
        num_mxs = len(mxtobxs)
        repetitive = {mx for mx, bxs in mxtobxs.items() if len(bxs) >= self.args.C}
        for mx in repetitive:
            del mxtobxs[mx]
        for xs in progress(bxtomxs.values()):
            xs -= repetitive
        print(
            int(timeit.default_timer() - t0),
            "Removed", len(repetitive), "most frequent minimizers of", num_mxs,
            f"({round(100 * len(repetitive) / num_mxs, 2)}%)", file=sys.stderr)

    def physlr_overlap(self):
        "Read a sketch of linked reads and find overlapping barcodes."

        bxtomxs = self.read_minimizers(self.args.FILES)
        mxtobxs = self.construct_minimizers_to_barcodes(bxtomxs)

        # Add the vertices.
        g = nx.Graph()
        for u, mxs in sorted(progress(bxtomxs.items())):
            if len(mxs) >= self.args.n:
                g.add_node(u, n=len(mxs))
        print(
            int(timeit.default_timer() - t0),
            "Added", g.number_of_nodes(), "barcodes to the graph", file=sys.stderr)

        # Add the overlap edges.
        edges = Counter(
            (u, v) for bxs in progress(mxtobxs.values()) for u, v in itertools.combinations(bxs, 2))
        print(int(timeit.default_timer() - t0), "Loaded", len(edges), "edges", file=sys.stderr)

        for (u, v), n in progress(edges.items()):
            if n >= self.args.n:
                g.add_edge(u, v, n=n)
        num_removed = len(edges) - g.number_of_edges()
        print(
            int(timeit.default_timer() - t0),
            "Removed", num_removed, "edges with fewer than", self.args.n,
            "common minimizers of", len(edges),
            f"({round(100 * num_removed / len(edges), 2)}%)", file=sys.stderr)

        num_singletons = Physlr.remove_singletons(g)
        print(
            int(timeit.default_timer() - t0),
            "Removed", num_singletons, "isolated vertices.", file=sys.stderr)

        Physlr.print_graph_stats(g)

        # Write the graph.
        self.write_tsv(g, sys.stdout)
        print(int(timeit.default_timer() - t0), "Wrote the graph", file=sys.stderr)

    def physlr_filter_overlap(self):
        "Read a Physlr overlap graph and filter edges."

        edge_weight = []

        at_edges = False
        print(int(timeit.default_timer() - t0), "Processing Nodes", file=sys.stderr)
        with open(self.args.FILES[0], "r") as overlap_input:
            for line in overlap_input:
                if at_edges and self.args.minimizer_overlap != 0:
                    columns = line.rstrip().split("\t")
                    edge_weight.append(int(columns[2]))
                elif at_edges and self.args.minimizer_overlap == 0:
                    print(line.rstrip())
                else:
                    print(line.rstrip())
                    if not line.strip():
                        print(next(overlap_input).rstrip())
                        print(int(timeit.default_timer() - t0), "Processing Edges", file=sys.stderr)
                        at_edges = True

        if self.args.minimizer_overlap == 0:
            return

        print(int(timeit.default_timer() - t0), "Sorting Edges", file=sys.stderr)
        edge_weight.sort()
        lower_threshold = edge_weight[int(len(edge_weight) * self.args.minimizer_overlap / 100) - 1]

        # Faster to read tsv again than to store edges as a dictionary
        print(int(timeit.default_timer() - t0), "Filtering Edges", file=sys.stderr)
        print(int(timeit.default_timer() - t0), "Lower Threshold", lower_threshold, file=sys.stderr)
        at_edges = False
        with open(self.args.FILES[0], "r") as overlap_input:
            for line in overlap_input:
                if at_edges:
                    columns = line.rstrip().split("\t")
                    if int(columns[2]) > lower_threshold:
                        print(line.rstrip())
                else:
                    if not line.strip():
                        at_edges = True
                        next(overlap_input)

    def physlr_degree(self):
        "Print the degree of each vertex."
        g = self.read_graph(self.args.FILES)
        Physlr.filter_edges(g, self.args.n)
        print("U\tn\tDegree")
        for u, prop in progress(g.nodes.items()):
            print(u, prop["n"], g.degree(u), sep="\t")
        print(int(timeit.default_timer() - t0), "Wrote degrees of vertices", file=sys.stderr)

    def physlr_common_neighbours(self):
        """Count the number of common neighbours for each edge."""
        g = self.read_graph(self.args.FILES)
        for u, v in g.edges:
            g[u][v]["n"] = len(list(nx.common_neighbors(g, u, v)))
        print(int(timeit.default_timer() - t0), "Counted common neighbours.",
              file=sys.stderr, flush=True)
        self.write_graph(g, sys.stdout, self.args.graph_format)
        print(int(timeit.default_timer() - t0), "Wrote the graph.", file=sys.stderr)

    def physlr_mst(self):
        """Determine the maximum spanning tree pruned for small branches."""
        g = self.read_graph(self.args.FILES)
        Physlr.remove_singletons(g)
        if Physlr.args.prune_bridges > 0:
            Physlr.remove_bridges(g, Physlr.args.prune_bridges)
        gmst = Physlr.determine_pruned_mst(g)

        self.write_graph(gmst, sys.stdout, self.args.graph_format)
        print(int(timeit.default_timer() - t0), "Wrote the MST.", file=sys.stderr)

    def physlr_mst_troubleshoot(self):
        """Determine the maximum spanning tree pruned for small branches."""
        g = self.read_graph(self.args.FILES)
        Physlr.remove_singletons(g)
        if Physlr.args.prune_bridges > 0:
            Physlr.remove_bridges(g, Physlr.args.prune_bridges)
        gmst = Physlr.determine_pruned_mst(g)

        print(int(timeit.default_timer() - t0), "Measuring branches.", file=sys.stderr, flush=True)
        for component in nx.connected_components(gmst):
            gcomponent = gmst.subgraph(component)
            branch_lengths = Physlr.measure_branch_length(gcomponent)
            for u, v in branch_lengths:
                gmst[u][v]["l"] = min(branch_lengths[(u, v)], branch_lengths[(v, u)])
        print(int(timeit.default_timer() - t0), "Measured branches.", file=sys.stderr, flush=True)

        self.write_graph(gmst, sys.stdout, self.args.graph_format)
        print(int(timeit.default_timer() - t0), "Wrote the MST.", file=sys.stderr)

    def physlr_report_junctions_graph(self):
        """
        Report junctions in the MST of the graph and output a list of junc barcodes.
        """
        g = self.read_graph(self.args.FILES)
        gmst = Physlr.determine_pruned_mst(g)
        print(int(timeit.default_timer() - t0), "Searching for junctions...", file=sys.stderr)
        tree_junctions = []
        for component in nx.connected_components(gmst):
            tree_junctions += Physlr.detect_junctions_of_tree(
                gmst.subgraph(component), Physlr.args.prune_junctions)
        print(int(timeit.default_timer() - t0),
              "Found", len(tree_junctions), "junctions.", file=sys.stderr)
        junctions = []
        if self.args.junction_depth > 0:
            print(int(timeit.default_timer() - t0),
                  "Exapnding junctions, depth:", self.args.junction_depth, file=sys.stderr)
            if self.args.junction_depth > 1:
                tree_junctions_expanded = set()
                for tree_junction in tree_junctions:
                    tree_junctions_expanded.update(
                        nx.bfs_tree(gmst, source=tree_junction,
                                    depth_limit=self.args.junction_depth-1))
            else:
                tree_junctions_expanded = set(tree_junctions)
            junctions = {m for n in tree_junctions_expanded for m in g.neighbors(n)}
            junctions.update(tree_junctions_expanded)
            print(int(timeit.default_timer() - t0),
                  "Exapnded to", len(junctions), "junctions.", file=sys.stderr)
        else:
            junctions = tree_junctions
        for junction in junctions:
            print(junction, file=sys.stdout)
        print(int(timeit.default_timer() - t0), "Wrote junctions.", file=sys.stderr)

    def physlr_remove_bridges_graph(self):
        """
        Iteratively remove bridges in the MST of the graph and output the graph.
        """
        g = self.read_graph(self.args.FILES)
        if Physlr.args.prune_bridges > 0:
            Physlr.remove_bridges(g, Physlr.args.prune_bridges)
        self.write_graph(g, sys.stdout, self.args.graph_format)

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
        backbones = Physlr.determine_backbones_and_remove_chimera(g)
        subgraph = nx.Graph()
        for backbone in backbones:
            gbackbone = g.subgraph(backbone)
            subgraph.add_nodes_from(gbackbone.nodes.data())
            subgraph.add_edges_from(gbackbone.edges.data())
        subgraph = self.sort_vertices(subgraph)
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
    def split_minimizers_bx(bx, g, bxtomxs):
        "Partition the minimizers of the given barcode"
        bx_match = re.compile(r'^(\S+)_\d+$')
        mxs = bxtomxs[bx]
        mol = 0
        mol_list = []
        while g.has_node(bx + "_" + str(mol)):
            bxmol = bx + "_" + str(mol)
            neighbour_mxs_list = [bxtomxs[re.search(bx_match, v).group(1)] \
                                  for v in g.neighbors(bxmol) \
                                  if re.search(bx_match, v).group(1) in bxtomxs]
            if not neighbour_mxs_list:
                neighbour_mxs_list = [set()]
            neighbour_mxs_set = set.union(*neighbour_mxs_list)
            molec_mxs = set.intersection(mxs, neighbour_mxs_set)
            mol_list.append((bxmol, molec_mxs))
            mol += 1
        return mol_list

    @staticmethod
    def split_minimizers_bx_process(bx):
        """
        Partition the minimizers of this barcode.
        The Graph and bx->min dictionary are passed as class variables.
        """
        return Physlr.split_minimizers_bx(bx, Physlr.graph, Physlr.bxtomxs)

    def physlr_split_minimizers(self):
        "Given the molecule overlap graph, split the minimizers into molecules"
        if len(self.args.FILES) < 2:
            msg = "physlr split-minimizers: error: graph file and bx to minimizer inputs required"
            sys.exit(msg)
        g = self.read_graph([self.args.FILES[0]])
        bxtomxs = self.read_minimizers([self.args.FILES[1]])

        if self.args.threads == 1:
            moltomxs = [self.split_minimizers_bx(bx, g, bxtomxs) for bx in progress(bxtomxs)]
            moltomxs = dict(x for l in moltomxs for x in l)

        else:
            Physlr.graph = g
            Physlr.bxtomxs = bxtomxs
            with multiprocessing.Pool(self.args.threads) as pool:
                moltomxs = dict(x for l in pool.map(self.split_minimizers_bx_process,
                                                    progress(bxtomxs), chunksize=100) for x in l)
            Physlr.graph = None
            Physlr.bxtomxs = None

        empty_ct = 0
        for mol in moltomxs:
            if not moltomxs[mol]:
                empty_ct += 1
                print("%s\t%s" % (mol, ""), file=sys.stdout)
                print("Warning:", mol, "has no associated minimizers", file=sys.stderr)
            else:
                print("%s\t%s" % (mol, " ".join(map(str, moltomxs[mol]))), file=sys.stdout)

    def physlr_split_reads_molecules(self):
        "Given the molecule -> minimizers table and the reads, partition reads into molecules"
        if len(self.args.FILES) < 3:
            sys.exit("physlr split-reads-molecules: error: molecule minimizers,\
            bx minimizers, reads inputs required")
        moltomxs = self.read_minimizers([self.args.FILES[0]])
        mol_counts = self.count_molecules_per_bx(moltomxs)
        bxtomxs_filename = self.args.FILES[1]
        num_pairs, num_valid_pairs, num_no_mx, num_equal_mx, num_no_int_mx = 0, 0, 0, 0, 0

        bxtomxs_file = open(bxtomxs_filename, 'r')

        read_count = 0
        for readfile in self.args.FILES[2:]:
            with open(readfile, 'r') as fin:
                for name, seq, bx, qual in read_fasta(fin):
                    if read_count == 0:
                        (name1, seq1, bx1, qual1) = (name, seq, bx, qual)
                        read_count += 1
                    else:
                        num_pairs += 1
                        mx_info1 = bxtomxs_file.readline().strip()
                        mx_info2 = bxtomxs_file.readline().strip()
                        if self.is_valid_pair(bx1, bx, name1, name) and bx in mol_counts:
                            (bx1_mol, mxs1) = self.parse_minimizer_line(mx_info1)
                            (bx2_mol, mxs2) = self.parse_minimizer_line(mx_info2)
                            if bx1_mol != bx2_mol or bx1_mol != bx1:
                                print("Should match: ", bx1_mol, bx2_mol, bx1, file=sys.stderr)
                                sys.exit("Error: Minimizer TSV order doesn't match reads fq file")

                            (mol, inc_no_int_mx, inc_equal_mx) = self.assign_read_molecule(
                                set.union(mxs1, mxs2), moltomxs, mol_counts, bx1)
                            bx_mol = bx1 + mol
                            num_no_int_mx += inc_no_int_mx
                            num_equal_mx += inc_equal_mx
                            self.print_read(name1 + " BX:Z:" + bx_mol, seq1, qual1)
                            self.print_read(name + " BX:Z:" + bx_mol, seq, qual)
                            num_valid_pairs += 1
                        elif self.is_valid_pair(bx1, bx, name1, name) and \
                                not self.args.molecules_bx_only:
                            self.print_read(name1 + " BX:Z:" + bx1, seq1, qual1)
                            self.print_read(name + " BX:Z:" + bx, seq, qual)
                            num_no_mx += 1
                        elif not self.args.molecules_bx_only:
                            self.print_read(name1, seq1, qual1)
                            self.print_read(name, seq, qual)
                        read_count = 0

        print("Saw", num_pairs, "read pairs, saw",
              num_valid_pairs, "valid read pairs with associated barcodes.",
              num_no_mx, "read pairs' barcodes had no split minimizers.",
              num_no_int_mx, "read pairs' barcodes had no intersecting minimizers",
              num_equal_mx, "read pairs had multiple molecule assignments",
              file=sys.stderr)
        bxtomxs_file.close()

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
    def assign_read_molecule(mxs, moltomxs, mol_counts, bx):
        "Given the minimizers of a read and the barcode, assign to a molecule"
        intersections = {mol: len(set.intersection(mxs, moltomxs[bx + "_" + str(mol)]))
                         for mol in range(0, mol_counts[bx]) if bx + "_" + str(mol) in moltomxs}
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
    def detect_communities_biconnected_components(g, node_set):
        """Separate bi-connected components. Return components."""
        cut_vertices = set(nx.articulation_points(g.subgraph(node_set)))
        components = list(nx.connected_components(g.subgraph(node_set - cut_vertices)))
        return components

    @staticmethod
    def detect_communities_common_neighbours(g, node_set, cn_threshold):
        """
        Filter edges with few common neighbours (local bridges).
        Separate bi-connected components. Return components.
        """
        subgraph = g.subgraph(node_set).copy()
        weak_edges = [(u, v) for u, v in progress(subgraph.edges())
                      if len(list(nx.common_neighbors(subgraph, u, v))) < cn_threshold]
        subgraph.remove_edges_from(weak_edges)
        cut_vertices = list(nx.articulation_points(subgraph))
        subgraph.remove_nodes_from(cut_vertices)
        return list(nx.connected_components(subgraph))

    @staticmethod
    def detect_communities_k_clique(g, node_set, k):
        """Apply k-clique community detection. Return communities."""
        return list(nx.algorithms.community.k_clique_communities(g.subgraph(node_set), k))

    @staticmethod
    def detect_communities_louvain(g, node_set, init_communities=None):
        """Apply Louvain community detection on a single component. Return communities."""
        if len(node_set) < 2:
            return []
        import community as louvain
        partition = louvain.best_partition(g.subgraph(node_set), init_communities)
        return [{node for node in partition.keys() if partition[node] == com}
                for com in set(partition.values())]

    @staticmethod
    def detect_communities_cosine_of_squared(g, node_set, squaring=True, threshold=0.7):
        """
        Square the adjacency matrix and then use cosine similarity to detect communities.
        Return communities.
        """
        import scipy as sp
        import numpy as np
        from sklearn.metrics.pairwise import cosine_similarity

        communities = []
        if len(node_set) > 1:
            adj_array = nx.adjacency_matrix(g.subgraph(node_set)).toarray()
            if squaring:
                new_adj = np.multiply(
                    cosine_similarity(
                        sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) >= threshold, adj_array)
            else:
                new_adj = np.multiply(cosine_similarity(adj_array) >= threshold, adj_array)
            edges_to_remove = np.argwhere(new_adj != adj_array)
            barcode_dict = dict(zip(range(len(node_set)), list(node_set)))
            edges_to_remove_barcode = [(barcode_dict[i], barcode_dict[j])
                                       for i, j in edges_to_remove]
            sub_graph_copy = nx.Graph(g.subgraph(node_set))
            sub_graph_copy.remove_edges_from(edges_to_remove_barcode)
            cos_components = list(nx.connected_components(sub_graph_copy))
            for com in cos_components:
                communities.append(com)
        return communities

    @staticmethod
    def partition_subgraph_into_bins_randomly(node_set, max_size=50):
        """
        Partition the subgraph into bins randomly for faster processing. Return bins.
        Warning: This function is not deterministic.
        """
        bins_count = 1 + len(node_set) // max_size
        node_list = list(node_set)
        random.shuffle(node_list)
        size, leftover = divmod(len(node_set), bins_count)
        bins = [node_list[0 + size * i: size * (i + 1)] for i in range(bins_count)]
        edge = size * bins_count
        for i in range(leftover):
            bins[i % bins_count].append(node_list[edge + i])
        return [set(x) for x in bins]

    @staticmethod
    def merge_communities(g, communities, node_set=0, strategy=0, cutoff=20):
        """Merge communities if appropriate."""
        mode = 1
        if cutoff == -1:  # no merging
            return communities
        if len(communities) == 1 and (node_set == 0 or strategy != 1):
            return communities
        if strategy == 1:  # Merge by Initializing Louvain with the communities
            return Physlr.detect_communities_louvain(g, node_set, communities)
        # Ad-hoc Merge (default - strategy = 0)
        merge_network = nx.Graph()
        for i in range(len(communities)):
            merge_network.add_node(i)
        for i, com1 in enumerate(communities):
            for j, com2 in enumerate(communities):
                if i >= j:
                    continue
                if mode == 1:  # disjoint input communities.
                    if nx.number_of_edges(
                            g.subgraph(com1.union(com2))) - \
                            nx.number_of_edges(g.subgraph(com1)) - \
                            nx.number_of_edges(g.subgraph(com2)) > cutoff:
                        merge_network.add_edge(i, j)
                else:  # overlapping input communities.
                    if nx.number_of_edges(
                            g.subgraph(com1.union(com2))) - \
                            len(set(g.subgraph(com1).edges()).union(
                                set(g.subgraph(com2).edges()))) \
                            > cutoff:
                        merge_network.add_edge(i, j)
        return [{barcode for j in i for barcode in communities[j]}
                for i in nx.connected_components(merge_network)]

    @staticmethod
    def determine_molecules_partition_split_merge(g, component):
        """
        Assign the neighbours of this vertex to molecules using fast heuristic community detection.
        Pipeline: bi-connected (bc) + partition + bc + k-cliques communities + merge.
        """
        return [merged
                for bi_connected_component in
                Physlr.detect_communities_biconnected_components(g, component)
                for merged in
                Physlr.merge_communities(
                    g, [cluster
                        for bin_set in
                        Physlr.partition_subgraph_into_bins_randomly(
                            bi_connected_component)
                        for bi_con2 in
                        Physlr.detect_communities_biconnected_components(g, bin_set)
                        for cluster in
                        Physlr.detect_communities_k_clique(g, bi_con2, k=3)],
                    bi_connected_component, strategy=0)
                ]

    @staticmethod
    def determine_molecules(g, u, junctions, strategy):
        """Assign the neighbours of this vertex to molecules."""
        communities = [g[u].keys()]
        if junctions:
            if u not in junctions:
                strategy = "bc"
        alg_list = strategy.split("+")
        for algorithm in alg_list:
            communities_temp = []
            if algorithm == "bc":
                for component in communities:
                    communities_temp.extend(
                        Physlr.detect_communities_biconnected_components(g, component))
            elif algorithm == "cn2":
                for component in communities:
                    communities_temp.extend(
                        Physlr.detect_communities_common_neighbours(g, component, cn_threshold=2))
            elif algorithm == "cn3":
                for component in communities:
                    communities_temp.extend(
                        Physlr.detect_communities_common_neighbours(g, component, cn_threshold=3))
            elif algorithm == "k3":
                for component in communities:
                    communities_temp.extend(
                        Physlr.detect_communities_k_clique(g, component, k=3))
            elif algorithm == "k4":
                for component in communities:
                    communities_temp.extend(
                        Physlr.detect_communities_k_clique(g, component, k=4))
            elif algorithm == "cos":
                for component in communities:
                    communities_temp.extend(
                        Physlr.detect_communities_cosine_of_squared(
                            g, component, squaring=False, threshold=0.4))
            elif algorithm == "sqcos":
                for component in communities:
                    communities_temp.extend(
                        Physlr.detect_communities_cosine_of_squared(g, component))
            elif algorithm == "louvain":
                for component in communities:
                    communities_temp.extend(
                        Physlr.detect_communities_louvain(g, component))
            if algorithm == "distributed":
                for component in communities:
                    communities_temp.extend(
                        Physlr.determine_molecules_partition_split_merge(g, component))
            communities = communities_temp

        communities.sort(key=len, reverse=True)
        return u, {v: i for i, vs in enumerate(communities) if len(vs) > 1 for v in vs}

    @staticmethod
    def determine_molecules_process(u):
        """
        Assign the neighbours of this vertex to molecules.
        The graph is passed in the class variable Physlr.graph.
        """
        return Physlr.determine_molecules(Physlr.graph, u, Physlr.junctions, Physlr.args.strategy)

    def physlr_molecules(self):
        "Separate barcodes into molecules."
        alg_white_list = {"bc", "cn2", "cn3", "k3", "k4", "cos", "sqcos", "louvain", "distributed"}
        alg_list = self.args.strategy.split("+")
        if not alg_list:
            sys.exit("Error: physlr molecule: missing parameter --separation-strategy")
        if not set(alg_list).issubset(alg_white_list):
            exit_message = "Error: physlr molecule: wrong input parameter(s) " + \
                      "--separation-strategy: " + str(set(alg_list) - alg_white_list)
            sys.exit(exit_message)
        junctions = []
        if len(self.args.FILES) > 1:
            gin = self.read_graph([self.args.FILES[0]])
            with open(self.args.FILES[1]) as fin:
                for line in fin:
                    junctions.append(line.split()[0])
            print(
                int(timeit.default_timer() - t0),
                "Separating junction-causing barcodes into molecules "
                "using the following algorithm(s):\n\t",
                self.args.strategy.replace("+", " + "),
                "\n\tand other barcodes with bc.",
                file=sys.stderr)
        else:
            gin = self.read_graph(self.args.FILES)
            print(
                int(timeit.default_timer() - t0),
                "Separating barcodes into molecules using the following algorithm(s):\n\t",
                self.args.strategy.replace("+", " + "),
                file=sys.stderr)

        Physlr.filter_edges(gin, self.args.n)

        # Partition the neighbouring vertices of each barcode into molecules.
        if self.args.threads == 1:
            molecules = dict(
                self.determine_molecules(
                    gin, u, junctions, self.args.strategy) for u in progress(gin))
        else:
            Physlr.graph = gin
            Physlr.junctions = junctions
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

        num_singletons = Physlr.remove_singletons(gout)
        print(
            int(timeit.default_timer() - t0),
            "Removed", num_singletons, "isolated vertices.", file=sys.stderr)
        self.write_graph(gout, sys.stdout, self.args.graph_format)
        print(int(timeit.default_timer() - t0), "Wrote graph", file=sys.stderr)

    @staticmethod
    def write_subgraphs_stats(g, fout):
        "Write statistics of the subgraphs."
        print("Barcode\tNodes\tEdges\tDensity", file=fout)
        for i in progress(g):
            print(i, g[i][0], g[i][1], g[i][2], sep="\t", file=fout)

    @staticmethod
    def subgraph_stats(g, u):
        """Extract the statistics of the vertex-induced subgraph with the vertex being u."""
        sub_graph = g.subgraph(g.neighbors(u))
        nodes_count = sub_graph.number_of_nodes()
        edges_count = sub_graph.number_of_edges()
        if nodes_count < 2:
            return u, [nodes_count, edges_count, 0.0]
        return u, (nodes_count, edges_count,
                   (edges_count*2.0/(nodes_count*(nodes_count-1))))

    @staticmethod
    def subgraph_stats_process(u):
        """
        Extract the statistics of the subgraph of neighbours of this vertex.
        The graph is passed in the class variable Physlr.graph.
        """
        return Physlr.subgraph_stats(Physlr.graph, u)

    def physlr_subgraphs_stats(self):
        "Retrieve subgraphs' stats."
        gin = self.read_graph(self.args.FILES)
        Physlr.filter_edges(gin, self.args.n)
        print(
            int(timeit.default_timer() - t0),
            "Computing statistics of the subgraphs...", file=sys.stderr)
        if self.args.threads == 1:
            stats = dict(self.subgraph_stats(gin, u) for u in progress(gin))
        else:
            Physlr.graph = gin
            with multiprocessing.Pool(self.args.threads) as pool:
                stats = dict(pool.map(
                    self.subgraph_stats_process, progress(gin), chunksize=100))
            Physlr.graph = None
        print(int(timeit.default_timer() - t0), "Extracted subgraphs' statistics.", file=sys.stderr)
        self.write_subgraphs_stats(stats, sys.stdout)

    @staticmethod
    def index_minimizers_in_backbones(backbones, bxtomxs):
        "Index the positions of the minimizers in the backbones."
        mxtopos = {}
        for tid, path in enumerate(progress(backbones)):
            for pos, u in enumerate(path):
                if u not in bxtomxs:
                    u = u.rsplit("_", 1)[0]
                if u not in bxtomxs:
                    u = u.rsplit("_", 1)[0]
                for mx in bxtomxs[u]:
                    mxtopos.setdefault(mx, set()).add((tid, pos))
        print(
            int(timeit.default_timer() - t0),
            "Indexed", len(mxtopos), "minimizers", file=sys.stderr)
        return mxtopos

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

    def map_indexing(self):
        "Load data structures and indexes required for mapping."
        if len(self.args.FILES) < 3:
            sys.exit("physlr map: error: at least three file arguments are required")
        path_filenames = [self.args.FILES[0]]
        target_filenames = [self.args.FILES[1]]
        query_filenames = self.args.FILES[2:]

        # Index the positions of the minimizers in the backbone.
        moltomxs = Physlr.read_minimizers(target_filenames)
        query_mxs = moltomxs if target_filenames == query_filenames else \
            Physlr.read_minimizers_list(query_filenames)

        # Index the positions of the markers in the backbone.
        backbones = Physlr.read_paths(path_filenames)
        backbones = [backbone for backbone in backbones
                     if len(backbone) >= self.args.min_component_size]
        mxtopos = Physlr.index_minimizers_in_backbones(backbones, moltomxs)

        return query_mxs, mxtopos, backbones

    def physlr_map_mkt(self):
        """
        Map sequences to a physical map.
        Usage: physlr map TGRAPH.path TMARKERS.tsv QMARKERS.tsv... >MAP.bed
        """
        import physlr.mkt
        import numpy

        query_mxs, mxtopos, _backbones = self.map_indexing()

        # Map the query sequences to the physical map.
        num_mapped = 0
        for qid, mxs in progress(query_mxs.items()):
            # Count the number of minimizers mapped to each target position.
            tidpos_to_n = Counter(pos for mx in mxs for pos in mxtopos.get(mx, ()))
            # Map each target position to a query position.
            # tid->tpos->qpos_list
            tid_to_qpos = {}
            for qpos, mx in enumerate(mxs):
                for (tid, tpos) in mxtopos.get(mx, ()):
                    if not tid in tid_to_qpos:
                        tid_to_qpos[tid] = {}
                    if not tpos in tid_to_qpos[tid]:
                        tid_to_qpos[tid][tpos] = []
                    tid_to_qpos[tid][tpos].append(qpos)

            tid_to_mkt = {}
            for (tid, tpos_to_qpos) in tid_to_qpos.items():
                # build array of the time points of measurements
                # build array containing the measurements corresponding to entries of time
                timepoints = []
                measurements = []
                num_tpos = 0
                for (tpos, qpos_list) in tpos_to_qpos.items():
                    # do not use islands (noise?)
                    # determine count of non-island sequences
                    if tpos + 1 in tpos_to_qpos or tpos - 1 in tpos_to_qpos:
                        for qpos in qpos_list:
                            timepoints.append(tpos)
                            measurements.append(qpos)
                        num_tpos += 1
                if num_tpos > self.args.mkt_median_threshold or len(timepoints) > 50000:
                    # print("Warning ", len(timepoints), " minimizers positions in ", \
                    #       num_tpos, " backbone positions seen for scaffold ", qid, \
                    #       " to backbone ", tid, file=sys.stderr)
                    timepoints = []
                    measurements = []
                    for (tpos, qpos_list) in tpos_to_qpos.items():
                        if tpos + 1 in tpos_to_qpos or tpos - 1 in tpos_to_qpos:
                            timepoints.append(tpos)
                            measurements.append(statistics.median_low(qpos_list))
                tid_to_mkt[tid] = physlr.mkt.test(numpy.array(timepoints), \
                                                  numpy.array(measurements), \
                                                  1, self.args.p, "upordown")
            mapped = False
            for (tid, tpos), score in tidpos_to_n.items():
                if score >= self.args.n:
                    orientation = "."
                    if tid in tid_to_mkt:
                        # mk: string of test result
                        # m: slope
                        # c: intercept
                        # p: significance
                        result = tid_to_mkt[tid]
                        mapped = True
                        if result[3] < self.args.p and result[1] != 0:
                            orientation = "+" if result[1] > 0 else "-"
                    print(tid, tpos, tpos + 1, qid, score, orientation, sep="\t")
            if mapped:
                num_mapped += 1
        print(
            int(timeit.default_timer() - t0),
            "Mapped", num_mapped, "sequences of", len(query_mxs),
            f"({round(100 * num_mapped / len(query_mxs), 2)}%)", file=sys.stderr)

    def physlr_map(self):
        """
        Map sequences to a physical map.
        Usage: physlr map TPATHS.path TMARKERS.tsv QMARKERS.tsv... >MAP.bed
        """

        query_mxs, mxtopos, _backbones = self.map_indexing()

        # Map the query sequences to the physical map.
        num_mapped = 0
        for qid, mxs in progress(query_mxs.items()):
            # Map each target position to a query position.
            tidpos_to_qpos = {}
            for qpos, mx in enumerate(mxs):
                for tidpos in mxtopos.get(mx, ()):
                    tidpos_to_qpos.setdefault(tidpos, []).append(qpos)
            for tidpos, qpos in tidpos_to_qpos.items():
                tidpos_to_qpos[tidpos] = statistics.median_low(qpos)

            # Count the number of minimizers mapped to each target position.
            tidpos_to_n = Counter(pos for mx in mxs for pos in mxtopos.get(mx, ()))

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
            "Mapped", num_mapped, "sequences of", len(query_mxs),
            f"({round(100 * num_mapped / len(query_mxs), 2)}%)", file=sys.stderr)

    def physlr_map_paf(self):
        """
        Map sequences to a physical map and output a PAF file.
        Usage: physlr map TGRAPH.path TMARKERS.tsv QMARKERS.tsv... >MAP.paf
        """

        query_mxs, mxtopos, backbones = self.map_indexing()

        # Map the query sequences to the physical map.
        num_mapped = 0
        for qid, mxs in progress(query_mxs.items()):
            # Map each target position to a query position.
            tidpos_to_qpos = {}
            for qpos, mx in enumerate(mxs):
                for tidpos in mxtopos.get(mx, ()):
                    tidpos_to_qpos.setdefault(tidpos, []).append(qpos)
            for tidpos, qpos in tidpos_to_qpos.items():
                q0, q1, q2, q3, q4 = quantile([0, 0.25, 0.5, 0.75, 1], qpos)
                low_whisker = max(q0, int(q1 - self.args.coef * (q3 - q1)))
                high_whisker = min(q4, int(q3 + self.args.coef * (q3 - q1)))
                tidpos_to_qpos[tidpos] = (low_whisker, q2, high_whisker)

            # Count the number of minimizers mapped to each target position.
            tidpos_to_n = Counter(pos for mx in mxs for pos in mxtopos.get(mx, ()))

            mapped = False
            for (tid, tpos), score in tidpos_to_n.items():
                if score >= self.args.n:
                    mapped = True
                    # The qpos tuple is (low_whisker, median, high_whisker).
                    qmedian_before = tidpos_to_qpos.get((tid, tpos - 1), (None, None, None))[1]
                    qstart, qmedian, qend = tidpos_to_qpos[tid, tpos]
                    qmedian_after = tidpos_to_qpos.get((tid, tpos + 1), (None, None, None))[1]
                    orientation = Physlr.determine_orientation(
                        qmedian_before, qmedian, qmedian_after)
                    qlength = len(mxs)
                    tlength = len(backbones[tid])
                    mapq = int(100 * score / (qend - qstart))
                    print(
                        qid, qlength, qstart, qend,
                        orientation,
                        tid, tlength, tpos, tpos + 1,
                        score, qend - qstart, mapq, sep="\t")
            if mapped:
                num_mapped += 1
        print(
            int(timeit.default_timer() - t0),
            "Mapped", num_mapped, "sequences of", len(query_mxs),
            f"({round(100 * num_mapped / len(query_mxs), 2)}%)", file=sys.stderr)

    def physlr_liftover_paf(self):
        """Lift over query coordinates of a PAF file from minimzer index to nucleotide position."""
        if len(self.args.FILES) < 2:
            sys.exit("physlr liftover-paf: error: At least two file arguments are required")
        query_filenames = [self.args.FILES[0]]
        paf_filenames = self.args.FILES[1:]

        # Read the minimizer positions of the query sequence.
        liftover = {}
        qlengths = {}
        for qname, posmxs in Physlr.read_minimizers_pos(query_filenames).items():
            assert qname not in liftover
            liftover[qname] = [pos for pos, _ in posmxs]
            qlengths[qname] = max(liftover[qname])

        # Lift over the query coordinates of the PAF file.
        for qname, _qlength, qstart, qend, orientation, \
                tname, tlength, tstart, tend, score, length, mapq in \
                progress(Physlr.read_paf(paf_filenames)):
            print(
                qname, qlengths[qname], liftover[qname][qstart], liftover[qname][qend],
                orientation, tname, tlength, tstart, tend, score, length, mapq, sep="\t")

    @staticmethod
    def chr_isdecimal(x):
        "Return true if this string contains only decimal digits, ignoring a prefix of chr."
        return str.isdecimal(x[3:]) if x.startswith("chr") else str.isdecimal(x)

    @staticmethod
    def chr_int(x):
        "Convert this string to an integer, ignoring a prefix of chr."
        if not Physlr.chr_isdecimal(x):
            return None
        return int(x[3:]) if x.startswith("chr") else int(x)

    def physlr_annotate_graph(self):
        """
        Annotate a graph with a BED or PAF file of mappings.
        Usage: physlr annotate-graph GRAPH PATH BED_OR_PAF... >ANNOTATED-GRAPHVIZ
        """

        if len(self.args.FILES) < 2:
            sys.exit("physlr annotate-graph: error: at least two file arguments are required")
        graph_filenames = [self.args.FILES[0]]
        path_filenames = [self.args.FILES[1]]
        bed_filenames = self.args.FILES[2:]

        g = self.read_graph(graph_filenames)
        Physlr.remove_small_components(g, self.args.min_component_size)
        backbones = Physlr.read_paths(path_filenames)

        # Associate vertices with mapped query names.
        utomapping = {}
        for tname, tstart, _, qname, score, _ in \
                progress(Physlr.read_bed(bed_filenames)):
            if score < self.args.n:
                continue
            tname = int(tname)
            tstart = int(tstart)
            u = backbones[tname][tstart]
            if u not in utomapping or score > utomapping[u][3]:
                utomapping[u] = (tname, tstart, qname, score)

        # Map reference names to integers.
        qnames = list({qname for _, _, qname, _ in utomapping.values()})
        qnames.sort(key=lambda s: (
            not Physlr.chr_isdecimal(s),
            Physlr.chr_int(s) if Physlr.chr_isdecimal(s) else s))
        qnametoindex = {s: i for i, s in enumerate(qnames)}
        qnames = None

        # Output the annotated graph in GraphViz format.
        print("strict graph {\nnode [style=filled width=1 height=1]\nedge [color=lightgrey]")
        for u, prop in g.nodes.items():
            if u in utomapping:
                tname, tstart, qname, _ = utomapping[u]
                hue = round(qnametoindex[qname] / len(qnametoindex), 3)
                print(f'"{u}" [label="{tname}_{tstart}" color="{hue},1,1"]')
            else:
                print(f'"{u}"')
        for e, prop in g.edges.items():
            u, v = sorted(e)
            print(f'"{u}" -- "{v}" [label={prop["n"]}]')
        print("}")
        print(int(timeit.default_timer() - t0), "Wrote graph", file=sys.stderr)

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

    @staticmethod
    def read_arcs_pair(filename):
        """
        Read ARCS pair tsv. Return a dictionary of pairs to orientation evidence.
        """
        pairs = {}
        with open(filename, "r") as arks_pair_file:
            for line in arks_pair_file:
                columns = line.rstrip().split("\t")
                (u, v, hh, ht, th, tt) = [int(column) if idx > 1 else column
                                          for idx, column in enumerate(columns)]
                pairs[(u, v)] = [hh, ht, th, tt]
                pairs[(v, u)] = [hh, th, ht, tt]
        print("Read ARCS pairs", file=sys.stderr)
        return pairs

    @staticmethod
    def read_dist_est(filename, dist_type):
        """
        Read ARCS distance estimates. Return a dictionary of pairs to distance estimates.
        """
        dist = {}
        dist_type_to_idx = {"min":2, "avg":3, "max":4}
        if dist_type not in dist_type_to_idx:
            print("invalid --dist-type parameters. Acceptable values are: min, avg, max",
                  file=sys.stderr)
            sys.exit(1)
        idx = dist_type_to_idx[dist_type]
        with open(filename, "r") as dist_est_file:
            for line in dist_est_file:
                columns = line.rstrip().split("\t")
                if columns[0] == "contig1":
                    continue

                contig1 = columns[0].rstrip("-+")
                contig2 = columns[1].rstrip("-+")
                dist[(contig1, contig2)] = int(columns[idx])
                dist[(contig2, contig1)] = int(columns[idx])

        print("Read Distance Estimates", file=sys.stderr)
        return dist

    @staticmethod
    def normal_estimation(x, probability, n):
        """
        Normal approximation to the binomial distribution
        """
        mean_val = n * probability
        std_dev = float(math.sqrt(mean_val * (1 - probability)))
        return 0.5 * (1 + math.erf((x - mean_val) / (std_dev * math.sqrt(2))))

    @staticmethod
    def check_link_significance(ori_list):
        """
        Check if barcode link is strong enough
        """
        max_idx = ori_list.index(max(ori_list))
        ori_list.sort()
        max_val = ori_list[-1]
        sum_with_second_max = ori_list[-1] + ori_list[-2]
        normal_cdf = Physlr.normal_estimation(max_val, 0.5, sum_with_second_max)
        if 1 - normal_cdf < 0.05:
            return max_idx
        return -1

    @staticmethod
    def orient_part_of_path_backward(pairs, path, unoriented, curr_pos, name):
        """
        Orient small part of path based on ARCS scaffold pairing information going backwards
        """
        idxtojoin = {0:"-+", 1:"--", 2:"++", 3:"+-"}
        while unoriented:
            prev_pos = curr_pos - 1
            pair = (path[prev_pos][:-1], name[:-1])
            oriented = False
            if pair in pairs:
                join_ori = pairs[pair]
                max_idx = Physlr.check_link_significance(join_ori)
                if max_idx != -1:
                    curr_ori = idxtojoin[max_idx][1]
                    if curr_ori == name[-1]:
                        prev_ori = idxtojoin[max_idx][0]
                        path[prev_pos] = path[prev_pos][0:-1] + prev_ori
                        unoriented.pop()
                        curr_pos = prev_pos
                        name = path[prev_pos]
                        oriented = True
            if not oriented:
                unoriented.clear()

        return path, unoriented

    @staticmethod
    def orient_part_of_path_forward(pairs, path, unoriented, curr_pos, name):
        """
        Orient small part of path based on ARCS scaffold pairing information going forwards
        """
        idxtojoin = {0:"-+", 1:"--", 2:"++", 3:"+-"}
        prev_pos = curr_pos - 1
        prev_name = path[prev_pos]
        pair = (prev_name[:-1], name[:-1])
        oriented = False
        if pair in pairs:
            join_orientation = pairs[pair]
            max_idx = Physlr.check_link_significance(join_orientation)
            if max_idx != -1:
                if idxtojoin[max_idx][0] == prev_name[-1]:
                    path[curr_pos] = name[:-1] + idxtojoin[max_idx][1]
                    oriented = True
        if not oriented:
            unoriented.clear()

        return path, unoriented

    @staticmethod
    def orient_path(path, pairs):
        """
        Orient path based on ARCS scaffold pairing information
        """
        unoriented = []

        if len(path) == 1 and path[0][-1] == ".":
            path[0] = path[0][0:-1] + "+"
            return path

        for curr_pos, name in enumerate(path):
            if curr_pos == 0:
                if name[-1] == ".":
                    unoriented.append(curr_pos)
            else:
                if name[-1] != ".":
                    if unoriented:
                        temp_curr_pos = curr_pos
                        temp_name = name

                        path, unoriented = Physlr.orient_part_of_path_backward(pairs, path,
                                                                               unoriented, curr_pos,
                                                                               name)

                        curr_pos = temp_curr_pos
                        name = temp_name
                else:
                    if not unoriented:
                        path, unoriented = Physlr.orient_part_of_path_forward(pairs, path,
                                                                              unoriented, curr_pos,
                                                                              name)

                    else:
                        unoriented.append(curr_pos)
        return path

    @staticmethod
    def orient_paths(paths, pairs):
        """
        Orient paths based on ARCS scaffold pairing information
        """
        if not pairs:
            return paths
        for idx, path in enumerate(paths):
            paths[idx] = Physlr.orient_path(path, pairs)
        return paths

    @staticmethod
    def generate_seq_with_dist(seqs, dist, path, gaps):
        """
        Return scaffold using distance estimates based on ARCS distance estimation
        """
        seq = ""
        for idx, name in enumerate(path):
            if seq == "":
                if name[-1] == ".":
                    seq += ("N" * len(seqs[name[0:-1]]))
                else:
                    seq += Physlr.get_oriented_sequence(seqs, name)
            else:
                pair = (name[:-1], path[idx - 1][:-1])
                if name[-1] != "." and path[idx - 1][-1] != ".":
                    if pair in dist:
                        seq += (dist[pair] * "N")
                    else:
                        seq += gaps
                    seq += Physlr.get_oriented_sequence(seqs, name)
                elif name[-1] != "." and path[idx - 1][-1] == ".":
                    seq += gaps
                    seq += Physlr.get_oriented_sequence(seqs, name)
                else:
                    seq += gaps
                    seq += ("N" * len(seqs[name[0:-1]]))
        return seq

    @staticmethod
    def path_to_fasta_no_arcs(seqs, paths, gaps):
        """
        Convert path to fasta when --arcs-pair isn't used
        """

        num_scaffolds = 0
        num_contigs = 0
        num_bases = 0

        for path in progress(paths):
            if not path:
                continue

            all_unoriented = all([name[-1] == "." for name in path])
            if all_unoriented:
                continue

            seq = gaps.join(Physlr.get_oriented_sequence(seqs, name)
                            if name[-1] != "." else ("N" * len(seqs[name[0:-1]]))
                            for name in path)

            if len(seq) < Physlr.args.min_length:
                continue
            num_scaffolds += 1
            print(f">{str(num_scaffolds).zfill(7)} LN:i:{len(seq)} xn:i:{len(path)}\n{seq}")
            num_contigs += len(path)
            num_bases += len(seq)

        return num_scaffolds, num_contigs, num_bases

    @staticmethod
    def path_to_fasta_with_arcs(seqs, paths, gaps):
        """
        Convert path to fasta when --arcs-pair is used
        """

        pairs = Physlr.read_arcs_pair(Physlr.args.arcs_pair)
        dist = Physlr.read_dist_est(Physlr.args.dist_est, Physlr.args.dist_type)

        num_scaffolds = 0
        num_contigs = 0
        num_bases = 0

        num_unoriented = sum([1 for path in paths for name in path if name[-1] == "."])
        print(num_unoriented, "unoriented pieces before using ARCS pair information",
              file=sys.stderr)

        paths = Physlr.orient_paths(paths, pairs)

        num_unoriented = sum([1 for path in paths for name in path if name[-1] == "."])
        print(num_unoriented, "unoriented pieces after using ARCS pair information",
              file=sys.stderr)

        for path in paths:
            if not dist:
                seq = gaps.join(Physlr.get_oriented_sequence(seqs, name)
                                if name[-1] != "." else ("N" * len(seqs[name[0:-1]]))
                                for name in path)
            else:
                seq = Physlr.generate_seq_with_dist(seqs, dist, path, gaps)

            if len(seq) < Physlr.args.min_length:
                continue
            num_scaffolds += 1
            print(f">{str(num_scaffolds).zfill(7)} LN:i:{len(seq)} xn:i:{len(path)}\n{seq}")
            num_contigs += len(path)
            num_bases += len(seq)

        return num_scaffolds, num_contigs, num_bases, paths

    def physlr_path_to_fasta(self):
        """
        Produce sequences in FASTA format from paths.
        Usage: physlr path-to-fasta FASTA PATH... >FASTA
        """
        if len(self.args.FILES) < 2:
            sys.exit("physlr path-to-fasta: error: at least two file arguments are required")
        fasta_filenames = self.args.FILES[0:1]
        path_filenames = self.args.FILES[1:]
        seqs = Physlr.read_fastas(fasta_filenames)
        paths = Physlr.read_paths(path_filenames)

        gaps = "N" * self.args.gap_size

        if self.args.arcs_pair == "":
            num_scaffolds, num_contigs, num_bases = Physlr.path_to_fasta_no_arcs(seqs, paths, gaps)
        else:
            num_scaffolds, num_contigs, num_bases, paths = \
                Physlr.path_to_fasta_with_arcs(seqs, paths, gaps)

        used_seqs = {name[0:-1] for path in paths for name in path if name[-1] != "."}

        for name in seqs:
            if name not in used_seqs:
                seq = seqs[name]
                if len(seq) < self.args.min_length:
                    continue
                num_scaffolds += 1
                print(f">{str(num_scaffolds).zfill(7)} LN:i:{len(seq)} xn:i:{1}\n{seq}")
                num_contigs += 1
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
            sys.exit("physlr filter-bed: error: at least two file arguments are required")
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
        for i, path in enumerate(progress(Physlr.read_paths(path_filenames))):
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
                    if self.args.verbose >= 3:
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
    def compute_ngxx(xs, g, proportion):
        "Compute the NGxx metric. xs must be sorted from largest to smallest."
        target = g * proportion
        running_sum = 0
        for x in xs:
            running_sum += x
            if running_sum >= target:
                return x
        return 0

    def physlr_path_metrics(self):
        "Report assembly metrics of a path file."
        if self.args.g is None:
            sys.exit("physlr path-metrics: error: You must specify -g, --expected-molecules")
        print("Max\tNG25\tNG50\tNG75\tMin\tPaths\tNodes\tFile")
        for filename in self.args.FILES:
            paths = self.read_paths([filename])
            xs = [len(path) for path in paths if len(path) >= self.args.min_component_size]
            if not xs:
                print(f"0\t0\t0\t0\t0\t0\t0\t{filename}")
                continue
            xs.sort(reverse=True)
            ng25 = Physlr.compute_ngxx(xs, self.args.g, 0.25)
            ng50 = Physlr.compute_ngxx(xs, self.args.g, 0.50)
            ng75 = Physlr.compute_ngxx(xs, self.args.g, 0.75)
            print(max(xs), ng25, ng50, ng75, min(xs), len(xs), sum(xs), filename, sep="\t")

    def physlr_paf_metrics(self):
        "Report assembly metrics of a PAF file."
        if self.args.g is None:
            sys.exit("physlr paf-metrics: error: You must specify -g, --genome-size")
        print("Max\tNG25\tNG50\tNG75\tMin\tCount\tSum\tFile")
        for filename in self.args.FILES:
            xs = []
            for _qname, _qlength, qstart, qend, _orientation, \
                    _tname, _tlength, _tstart, _tend, _score, _length, _mapq in \
                    progress(Physlr.read_paf([filename])):
                xs.append(qend - qstart)
            if not xs:
                print(f"0\t0\t0\t0\t0\t0\t0\t{filename}")
                continue
            xs.sort(reverse=True)
            ng25 = Physlr.compute_ngxx(xs, self.args.g, 0.25)
            ng50 = Physlr.compute_ngxx(xs, self.args.g, 0.50)
            ng75 = Physlr.compute_ngxx(xs, self.args.g, 0.75)
            print(max(xs), ng25, ng50, ng75, min(xs), len(xs), sum(xs), filename, sep="\t")

    def physlr_fasta_gaps(self):
        """Print the coordinates in BED format of gaps in a FASTA file."""
        seqs = Physlr.read_fastas(self.args.FILES)
        gap_regex = re.compile(r'NN*')
        gap = 0
        for name, seq in seqs.items():
            for match in re.finditer(gap_regex, seq):
                gap += 1
                if match.end() - match.start() >= self.args.n:
                    print(name, match.start(), match.end(), f"gap{gap}", sep="\t")

    def physlr_find_ntcard_mode(self):
        """Find the first mode after minimum in ntCard histogram output."""
        # Assumption: There is no negative slope to the right of the first local minimum
        freq_count = [int(line.rstrip().split("\t")[2]) for line in open(self.args.FILES[0])
                      if line[0] != "k"]
        min_val = freq_count[0]
        for idx, freq in enumerate(freq_count):
            if freq > min_val:
                min_idx = idx - 1
                break
            min_val = freq
        freq_count = freq_count[min_idx:]
        print(freq_count.index(max(freq_count)) + 1 + min_idx)

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
            "--separation-strategy", action="store", dest="strategy", default="bc+k3",
            help="strategy for barcode to molecule separation [bc+k3]. Use a combination"
                 " of bc, k3, cos, sqcos, louvain, and distributed"
                 " concatenated plus sign (example:bc+k3+bc)")
        argparser.add_argument(
            "--coef", action="store", dest="coef", type=float, default=1.5,
            help="ignore minimizers that occur in Q3+c*(Q3-Q1) or more barcodes [0]")
        argparser.add_argument(
            "-c", "--min-count", action="store", dest="c", type=int, default=2,
            help="ignore minimizers that occur in less than c barcodes [2]")
        argparser.add_argument(
            "-C", "--max-count", action="store", dest="C", type=int,
            help="ignore minimizers that occur in C or more barcodes [None]")
        argparser.add_argument(
            "-g", "--expected-molecules", action="store", dest="g", type=int,
            help="the expected number of molecules in the assembly, for computing NGxx [None]")
        argparser.add_argument(
            "-M", "--max-molecules", action="store", dest="M", type=int,
            help="remove barcodes with M or more molecules [None]")
        argparser.add_argument(
            "-n", "--min-n", action="store", dest="n", type=int, default=0,
            help="remove edges with fewer than n shared minimizers [0]")
        argparser.add_argument(
            "-N", "--max-n", action="store", dest="N", type=int, default=None,
            help="remove edges with at least N shared minimizers [None]")
        argparser.add_argument(
            "--bestn", action="store", dest="bestn", type=int, default=None,
            help="Keep the best n edges of each vertex [None]")
        argparser.add_argument(
            "--min-length", action="store", dest="min_length", type=int, default=0,
            help="remove sequences with length less than N bp [0]")
        argparser.add_argument(
            "--min-component-size", action="store", dest="min_component_size", type=int, default=0,
            help="remove components with fewer than N vertices [0]")
        argparser.add_argument(
            "--min-path-size", action="store", dest="min_path_size", type=int, default=200,
            help="skip paths with fewer than N barcodes [200]")
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
            "--exclude-vertices", action="store", dest="exclude_vertices", default="",
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
            "-s", "--support", action="store", dest="s", type=int, default=3,
            help="Minimum number of barcodes required to support a position")
        argparser.add_argument(
            "--mkt-median-threshold", action="store", dest="mkt_median_threshold",
            type=int, default=50,
            help="Max number of backbones before using only medians in Mann-Kendall Test")
        argparser.add_argument(
            "-V", "--verbose", action="store", dest="verbose", type=int, default="2",
            help="the level of verbosity: 0:silent, 1:periodic, 2:progress, 3:verbose [2]")
        argparser.add_argument(
            "--version", action="version", version="physlr 0.1.0")
        argparser.add_argument(
            "command",
            help="A command")
        argparser.add_argument(
            "FILES", nargs="+",
            help="FASTA/FASTQ, TSV, or GraphViz format")
        argparser.add_argument(
            "--prune-branches", action="store", dest="prune_branches", type=int, default=10,
            help="size of the branches to be pruned [10]. set to 0 to skip prunning.")
        argparser.add_argument(
            "--prune-bridges", action="store", dest="prune_bridges", type=int, default=0,
            help="size of the bridges to be pruned [0]. set to 0 to skip bridge prunning.")
        argparser.add_argument(
            "--prune-junctions", action="store", dest="prune_junctions", type=int, default=0,
            help="split a backbone path when the alternative branch is longer than"
                 "prune-junctions [0]. set to 0 to skip.")
        argparser.add_argument(
            "--junction-depth", action="store", dest="junction_depth", type=int, default=0,
            help="depth for expanding the junctions by collecting all neighbors of this depth [0].")
        argparser.add_argument(
            "--gap-size", action="store", dest="gap_size", type=int, default=100,
            help="gap size used in scaffolding [100].")
        argparser.add_argument(
            "--minimizer-overlap", action="store", dest="minimizer_overlap", type=float, default=0,
            help="Percent of edges to remove [0].")
        argparser.add_argument(
            "--arcs-pair", action="store", dest="arcs_pair", type=str, default="",
            help="ARCS scaffold pairing file.")
        argparser.add_argument(
            "--dist-est", action="store", dest="dist_est", type=str, default="",
            help="ARCS scaffold pairing distance estimation file.")
        argparser.add_argument(
            "--dist-type", action="store", dest="dist_type", type=str, default="avg",
            help="ARCS scaffold pairing distance type."
                 "Acceptable values are: min, avg, max")
        return argparser.parse_args()

    def __init__(self):
        "Create a new instance of Physlr."
        self.args = self.parse_arguments()
        Physlr.args = self.args
        self.args.FILES = ["/dev/stdin" if s == "-" else s for s in self.args.FILES]

    def main(self):
        "Run Physlr."
        method_name = "physlr_" + self.args.command.replace("-", "_")
        if not hasattr(Physlr, method_name):
            print("physlr: error: unrecognized command:", self.args.command, file=sys.stderr)
            sys.exit(1)
        getattr(Physlr, method_name)(self)

def main():
    "Run Physlr."
    Physlr().main()

if __name__ == "__main__":
    main()
