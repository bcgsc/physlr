#!/usr/bin/env python3
import networkx as nx
import argparse
from intervaltree import *
from networkx.drawing.nx_agraph import write_dot

'''
Given a molecule TSV file from Tigmint, produce an overlap graph of barcodes (Nodes: barcodes; Edges: molecule overlap based on alignments)
Also output TSV file summarizing the molecule extents per barcode
'''

#Given the molecule TSV file, read through and store all the molecule extents in an Interval  Tree
def createIntervalTree(TSVfilename):
	interval_trees = {} #chr -> interval_tree

	with open(TSVfilename, 'r') as TSV:
		header = TSV.readline()
		for extent in TSV:
			extent = extent.strip().split("\t")
			(chr, start, end, BX) = (extent[0], int(extent[1]), int(extent[2]), extent[4])
			if chr not in interval_trees:
				interval_trees[chr] = IntervalTree()
			interval_trees[chr][start:end] = BX

	return interval_trees

#Go through the molecule extents, and produce a graph of overlaps (Nodes: barcodes, Edges: molecule(s) of that barcode overlap)
def createOverlapGraph(barcode_trees, TSVfilename, prefix, min_overlap):
	# Go over each molecule, and determine its overlaps
	overlap_graph = nx.Graph()
	out_overlapCoords = open(prefix + ".overlaps.tsv", 'w')
	out_overlapCoords.write("B1\tB2\tchr_overlap\tStart_overlap\tEnd_overlap\n")
	# Node attributes: BX (name)
	# Edge attributes: Length of overlap, jaccard score
	with open(TSVfilename, 'r') as TSV:
		header = TSV.readline()
		for extent in TSV:
			extent = extent.strip().split("\t")
			(chr, start, end, BX) = (extent[0], int(extent[1]), int(extent[2]), extent[4])
			interval_tree = barcode_trees[chr]
			overlaps = interval_tree[start:end]
			for overlap in overlaps:
				# Compute the overlap between the extents
				len_overlap = min(end, overlap.end) - max(start, overlap.begin)
				if len_overlap < min_overlap:
					continue
				# Decide if this is an overlap we want to store
				# Store overlap if start of overlap is less than others. If equal starts, store lexicographically smaller barcode
				if start < overlap.begin or (start == overlap.begin and BX < overlap.data):
					if overlap_graph.has_edge(BX, overlap.data):
						print("WARNING: Already an edge between %s and %s" % (BX, overlap.data))
						overlap_graph[BX][overlap.data]['l'] += len_overlap
					else:
						overlap_graph.add_edge(BX, overlap.data, l=len_overlap)
				out_overlapCoords.write("%s\t%s\t%s\t%d\t%d\n" % (BX, overlap.data, chr, max(start, overlap.begin), min(end, overlap.end)))
	write_dot(overlap_graph, prefix + "overlap_graph.dot")

def main():
	parser = argparse.ArgumentParser(description="Produce a barcode overlap graph based on read alignments")
	parser.add_argument("TSV", type=str, help="Molecule TSV file from Tigmint, sorted by position")
	parser.add_argument("-m", type=int, help='Minimum overlap between molecules to create edge', default=500, required=False)
	parser.add_argument("-p", type=str, help="Prefix for output files", default="barcode_overlap", required=False)
	args = parser.parse_args()

	barcode_trees = createIntervalTree(args.TSV)
	createOverlapGraph(barcode_trees, args.TSV, args.p, args.m)

if __name__ == "__main__":
	main()