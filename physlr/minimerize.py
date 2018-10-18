#!/usr/bin/env pypy3
"""
Hash k-mers and integers
Written by Shaun Jackman @sjackman
"""

from physlr.hash import hash_kmer
import re

ACGT = re.compile("^[ACGT]+$")

def kmerize(k, seq):
    "Iterator over the kmers of a string."
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i : i + k]
        if ACGT.match(kmer):
            yield kmer

def minimerize(k, w, seq):
    "Return the minimizers of a string."
    hashes = [hash_kmer(kmer) for kmer in kmerize(k, seq)]
    minimizers = []
    previous_minimizer = -1
    for i in range(0, len(hashes) - w + 1):
        minimizer, minimizer_i = min((x, j) for (j, x) in enumerate(hashes[i : i + w]))
        minimizer_i += i
        if minimizer_i > previous_minimizer:
            previous_minimizer = minimizer_i
            minimizers.append(minimizer)
    return minimizers
