import sys, os, time
import networkx as nx
from networkx.drawing.nx_pydot import write_dot

"""
Build overlap graph.

INPUT: Input is a .tsv file of reads with their start and end coordinates as one
of many pieces of information. The file is assumed sorted in ascending order by
the start coordinates. The program goes through the file line by line computing
the overlap of a read on a line with reads in the subsequent lines that have
start coordinates lower than or equal to the read's end coordinate.

Overlap is the count of the number of units of a read's length that are common
between two reads.

OUTPUT: A .dot file for input to Graphivz.

WARNING: Never tested for performace or correctness.
"""

def main(argv):

    G = nx.Graph()

    # Read in all the lines.
    lines = [line.rstrip().split() for line in open('sample.tsv')]
    n = len(lines)

    # Go over line by line, till the penultimate line.
    for i in range(1, n-1):

        # Log the line number in progress.
        print(str(i+1) + ":", end="")
        
        v1 = i
        l1 = int(lines[i][1])
        u1 = int(lines[i][2])

        j = i+1; c = 0;

        v2 = j
        l2 = int(lines[j][1])
        u2 = int(lines[j][2])

        # Check subsequent reads for overlap and compute weight.
        while l2 < u1:
            d1 = u1 - l2
            d2 = u1 - u2
            w = d1 + 1
            if d2 > 0:
                w = w - d2
            G.add_edge(v1, v2, weight=w)
            j += 1; c += 1
            if j == n:
                break
            v2 = j
            l2 = int(lines[j][1])
            u2 = int(lines[j][2])
        # Log the number of overlapping reads.
        print(c)

    # Write out the graph to disk in Graphviz format.
    pos = nx.nx_agraph.graphviz_layout(G)
    nx.draw(G, pos=pos)
    write_dot(G, 'overlap.dot')

 
if __name__ == "__main__":
    main(sys.argv)
