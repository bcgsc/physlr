# Physlr: Construct a Physical Map from Linked Reads

Physlr constructs a *de novo* physical map using linked reads from 10X Genomics or stLFR. This physical map can then be used to scaffold an existing assembly to yield chromosomal level contiguity.

## Build Status
[![Build Status](https://dev.azure.com/bcgsc/btl/_apis/build/status/bcgsc.physlr?branchName=master)](https://dev.azure.com/bcgsc/btl/_build/latest?definitionId=1&branchName=master)


# Dependencies

* [ntCard](https://github.com/bcgsc/ntCard)
* [ntHits](https://github.com/bcgsc/ntHits)
* GCC 5 or newer with openmp
* Python 3.5 or newer and the following packages
    * [community](https://python-louvain.readthedocs.io/en/latest/api.html)
    * [networkx](https://networkx.github.io/)
    * [numpy](https://numpy.org/)
    * [scipy](https://www.scipy.org/)
    * [sklearn](https://scikit-learn.org/stable/)
    * [tqdm](https://tqdm.github.io/)

Additionally, we recommend using pypy3 over regular python3 for speed.


## Optional dependencies

- [pigz](https://zlib.net/pigz/) for parallel gzip
- [zsh](https://sourceforge.net/projects/zsh/) for reporting time and memory usage


# Compiling Physlr from source

```
cd src
make
```

# Running Physlr

## Generating Physlr Physical Map with stLFR reads

To construct a physical map, you would need linked reads from 10X Genomics or stLFR. To visualize the correctness and contiguity of the physical map, you wouls also need a reference genome. In this example, the linked reads and reference genome are called `linkedReads.fq.gz` and `reference.fa`, respectively. The linked reads are from stLFR so we specify `minimizer_overlap=stLFR` to use the default value for stLFR reads.

```
cd experiment
ln -s data/Makefile
ls
Makefile    linkedReads.fq.gz #stLFR   reference.fa
make physical-map lr=linkedReads ref=ref minimizer_overlap=stLFR
```

## Scaffolding a draft assembly with Physlr Physical Map

To scaffold a draft assembly, you would need linked reads from 10X Genomics or stLFR, and an existing assembly. To obtain Quast metrics for the Physlr scaffolded assembly, you wouls also need a reference genome. In this example, the linked reads, draft assembly, and reference genome are called `linkedReads.fq.gz`, `draft.fa`, `reference.fa`, respectively. The linked reads are from 10X so we specify `minimizer_overlap=10X` to use the default value for 10X reads.

```
cd experiment
ln -s data/Makefile
ls
Makefile    linkedReads.fq.gz #10X   reference.fa    draft.fa
make scaffolds lr=linkedReads ref=reference draft=draft minimizer_overlap=10X
```

See the help page for further information.

# Output files

* `lr.physlr.physical-map.path`: Paths of barcodes (backbones) describing chromosomal length sequences.
* `lr.physlr.physical-map.ref.n10.paf.gz.*.pdf`: Various graphs showing the contiguity and correctness of the backbones with respect to the reference.
* `draft.physlr.fa`: Physlr scaffolded assembly using the physical map.
* `draft.physlr.quast.tsv`: Quast metrics comparing the Physlr scaffolded assembly against the reference.

# Testing Compiled Physlr executables

```
cd src
make check
```

# Acknowledgements

This projects uses:
* [btl_bloomfilter](https://github.com/bcgsc/btl_bloomfilter) BTL C/C++ Common bloom filters for bioinformatics projects implemented by Justin Chu
* [klib](https://github.com/attractivechaos/klib) A standalone and lightweight C library implemented by Attractive Chaos
* [nthash](https://github.com/bcgsc/ntHash) rolling hash implementation by Hamid Mohamadi
* [readfq](https://github.com/Tessil/robin-map) Fast multi-line FASTA/Q reader API implemented by Heng Li
* [robin-map](https://github.com/Tessil/robin-map) C++ implementation of a fast hash map and hash set using robin hood hashing by Thibaut G.
