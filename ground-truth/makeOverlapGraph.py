#!/usr/bin/env python3
import networkx as nx
import argparse

'''
Given a molecule TSV file from Tigmint, produce an overlap graph of barcodes (Nodes: barcodes; Edges: molecule overlap based on alignments)
Also output TSV file summarizing the molecule extents per barcode
'''

def readBED(BEDfilename, min_overlap, prefix):
	overlap_graph = nx.Graph()
	with open(BEDfilename, 'r') as BED:
		for bed_entry in BED:
			bed_entry = bed_entry.strip().split("\t")
			(chr1, start1, end1, bx1, mi1) = (bed_entry[0], int(bed_entry[1]), int(bed_entry[2]), bed_entry[3], int(bed_entry[4]))
			(chr2, start2, end2, bx2, mi2) = (bed_entry[5], int(bed_entry[6]), int(bed_entry[7]), bed_entry[8], int(bed_entry[9]))
			overlap = int(bed_entry[10])
			if overlap < min_overlap:
				continue
			if start1 < start2 or (start1 == start2 and bx1 < bx2):
				if overlap_graph.has_edge(bx1, bx2):
					# print("WARNING: Already an edge between %s and %s" % (BX, overlap.data))
					overlap_graph[bx1][bx2]['l'] += overlap
				else:
					overlap_graph.add_edge(bx1, bx2, l=overlap)
	print("DONE building barcode graph")
	out_graph = open(prefix + ".overlap_graph.dot", 'w')
	out_tsv = open(prefix + ".overlap_graph.tsv", 'w')
	out_graph.write("strict graph  {\n")
	out_tsv.write("U\tV\tl\n")
	for edge in overlap_graph.edges():
		out_str = "\t\"%s\" -- \"%s\"\t[l=%d];\n" % (edge[0], edge[1], overlap_graph[edge[0]][edge[1]]["l"])
		out_graph.write(out_str)
		out_str = "%s\t%s\t%d\n" % (edge[0], edge[1], overlap_graph[edge[0]][edge[1]]["l"])
		out_tsv.write(out_str)
	out_graph.write("}\n")
	print("DONE writing barcode graph")
	out_graph.close()
	out_tsv.close()

def readBED_molec(BEDfilename, min_overlap, prefix):
	out_graph = open(prefix + ".overlap_graph_molec.dot", 'w')
	out_graph.write("strict graph  {\n")

	out_tsv = open(prefix + ".overlap_graph_molec.tsv", 'w')
	out_tsv.write("U\tV\tl\n")

	with open(BEDfilename, 'r') as BED:
		for bed_entry in BED:
			bed_entry = bed_entry.strip().split("\t")
			(chr1, start1, end1, bx1, mi1) = (bed_entry[0], int(bed_entry[1]), int(bed_entry[2]), bed_entry[3], int(bed_entry[4]))
			(chr2, start2, end2, bx2, mi2) = (bed_entry[5], int(bed_entry[6]), int(bed_entry[7]), bed_entry[8], int(bed_entry[9]))
			overlap = int(bed_entry[10])
			if overlap < min_overlap:
				continue
			if start1 < start2 or (start1 == start2 and bx1 < bx2):
				out_graph.write("\t\"%s_%s\" -- \"%s_%s\"\t[l=%d];\n" % (bx1, mi1, bx2, mi2, overlap))
				out_tsv.write("%s_%s\t%s_%s\t%d\n" % (bx1, mi1, bx2, mi2, overlap))
	out_graph.write("}\n")
	print("DONE building molecule graph")
	print("DONE writing molecule graph")
	out_graph.close()
	out_tsv.close()


def analyze_successors(G):
	out_file = open("degree.tsv", 'w')
	for node in G.nodes():
		deg = G.degree(node)
		out_file.write("%s\t%d\n" % (node, deg))
	out_file.close()


def main():
	parser = argparse.ArgumentParser(description="Produce a barcode overlap graph based on read alignments")
	parser.add_argument("TSV", type=str, help="Molecule TSV file from Tigmint, sorted by position")
	parser.add_argument("-m", type=int, help='Minimum overlap between molecules to create edge', default=500, required=False)
	parser.add_argument("-p", type=str, help="Prefix for output files", default="barcode_overlap", required=False)
	args = parser.parse_args()

	#barcode_trees = createIntervalTree(args.TSV)
	#createOverlapGraph(barcode_trees, args.TSV, args.p, args.m)
	G = readBED(args.TSV, args.m, args.p)
	readBED_molec(args.TSV, args.m, args.p)
	analyze_successors(G)

if __name__ == "__main__":
	main()