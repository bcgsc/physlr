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
cd ~/physlr/src
make
```

# Running Physlr

## Generating Physlr Physical Map with stLFR reads

```
cd ~/experiment
ln -s ~/physlr/data/Makefile
ls
Makefile    lr.fq.gz #stLFR   ref.fa
make physical-map lr=lr ref=ref minimizer_overlap=stLFR
```

## Scaffolding a draft assembly with Physlr Physical Map

```
cd ~/experiment
ln -s ~/physlr/data/Makefile
ls
Makefile    lr.fq.gz #10X   ref.fa    draft.fa
make scaffolds lr=lr ref=ref draft=draft minimizer_overlap=10X
```

# Output files

* `lr.physlr.physical-map.path`: Paths of barcodes (backbones) describing chromosomal length sequences.
* `lr.physlr.physical-map.ref.n10.paf.gz.*.pdf`: Various visualizations of backbones against the reference.
* `draft.physlr.fa`: Physlr scaffolded assembly with physical-map.
* `draft.physlr.quast.tsv`: Quast metrics using the Physlr scaffolded assembly as query against the reference.

# Testing Compiled Physlr executables

```
cd ~/physlr/src
make check
```

# Acknowledgements

This projects uses:
* [btl_bloomfilter](https://github.com/bcgsc/btl_bloomfilter) BTL C/C++ Common bloom filters for bioinformatics projects implemented by Justin Chu
* [klib](https://github.com/attractivechaos/klib) A standalone and lightweight C library implemented by Attractive Chaos
* [nthash](https://github.com/bcgsc/ntHash) rolling hash implementation by Hamid Mohamadi
* [readfq](https://github.com/Tessil/robin-map) Fast multi-line FASTA/Q reader API implemented by Heng Li
* [robin-map](https://github.com/Tessil/robin-map) C++ implementation of a fast hash map and hash set using robin hood hashing by Thibaut G.
