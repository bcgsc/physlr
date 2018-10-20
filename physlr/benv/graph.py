#!/usr/bin/env pypy3
"""
Build a barcode overlap graph from the minimizer index.
"""

import math

class Graph:

    barcodeToIndex = {}
    nextBarcodeIndex = 0
    barcodeToHashes = []
    barcodes = []
    hashToBarcodes = {}

    def read_index(self, fin):
        "Build graph data structures from the minimizer index."
        for line in fin:
            barcode, hashlist = line.split('\t')
            hashes = [int(x) for x in hashlist.split()]
            barcodeIndex = self.barcodeToIndex.get(barcode)
            if barcodeIndex == None:
                self.barcodeToIndex[barcode] = self.nextBarcodeIndex
                self.barcodes.append(barcode)
                self.barcodeToHashes.append([])
                barcodeIndex = self.nextBarcodeIndex
                self.nextBarcodeIndex = self.nextBarcodeIndex + 1
            self.barcodeToHashes[-1].extend(hashes)
            for h in hashes:
                self.hashToBarcodes.setdefault(h, []).append(barcodeIndex)
        # collapse duplicate entries
        for i, hashes in enumerate(self.barcodeToHashes):
            self.barcodeToHashes[i] = list(set(hashes))

    def output_graph(self, pmin, fmt="graphviz"):
        "Output the graph, with p-values on the edges."
        if fmt == "tsv":
            print("barcode1", "barcode2", "n1", "n2", "intersection",
                  "p", sep="\t")
        else:
            print("graph {")
        # universe size of possible minimizers
        H = len(self.hashToBarcodes)
        for i1, hashes1 in enumerate(self.barcodeToHashes):
            if i1 > len(self.barcodeToHashes)/2:
                break
            candidates = []
            for h in hashes1:
                candidates.extend(self.hashToBarcodes[h])
            candidates = list(set(candidates))
            for i2 in candidates:
                if i2 == i1:
                    continue
                barcode1 = self.barcodes[i1]
                barcode2 = self.barcodes[i2]
                hashes2 = self.barcodeToHashes[i2]
                # c1/c2: number of minimizers for barcodes i1/i2
                c1 = len(hashes1)
                c2 = len(hashes2)
                cmin = min(c1, c2)
                # space of possible hash values
                # E(c): expected intersection by chance
                # (hypergeometric distribution)
                E = c1 * c2 / H
                # c: intersection of hashes
                c = len(set(hashes1).intersection(set(hashes2)))
                # r: difference between expected and true intersection
                r = (c - E) / cmin
                # p: p-value to see this intersection size by chance
                p = math.exp(-2 * r**2 * cmin)
                if fmt == "tsv":
                    print(barcode1, barcode2, c1, c2, c,
                          "{:.3e}".format(p), sep="\t")
                else:
                    print('{} -- {} [p={:.3e}]'.format(self.barcodes[i1],
                                                       self.barcodes[i2], p))
        if fmt != "tsv":
            print("}")
