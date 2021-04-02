# Physlr: Construct a Physical Map from Linked Reads

Physlr constructs a *de novo* physical map using linked reads from 10X Genomics or stLFR. This physical map can then be used to scaffold an existing assembly to yield chromosomal level contiguity.

# Dependencies

* [ntCard](https://github.com/bcgsc/ntCard)
* [ntHits](https://github.com/bcgsc/ntHits)
* GCC 5 or newer with [OpenMP](https://www.openmp.org) and [boost](https://www.boost.org)
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
pip3 install --user git+https://github.com/bcgsc/physlr
git clone https://github.com/bcgsc/physlr
cd physlr/src && make install
```

To install Physlr in a specified directory:

```
pip3 install --user git+https://github.com/bcgsc/physlr
git clone https://github.com/bcgsc/physlr
cd physlr/src && make install PREFIX=/opt/physlr
```

# Running Physlr

## Generating Physlr Physical Map with stLFR reads

To construct a physical map, you need linked reads from 10X Genomics or stLFR. In addition, to visualize the correctness and contiguity of the physical map, you need a reference genome.
In this example, the linked reads and reference genome are called `linkedreads.fq.gz` and `reference.fa`, respectively. The linked reads are from stLFR so we specify `protocol=stlfr` to use the default value for stLFR reads.

```
cd experiment
bin/physlr-make physical-map lr=linkedreads ref=reference protocol=stlfr
```

## Scaffolding a draft assembly with Physlr Physical Map

To scaffold a draft assembly, you need linked reads from 10X Genomics or stLFR, and an existing assembly. In addition, to calculate Quast summary metrics for the Physlr scaffolded assembly, you need a reference genome.
In this example, the linked reads, draft assembly, and reference genome are called `linkedreads.fq.gz`, `draft.fa`, `reference.fa`, respectively. The linked reads are from 10X Genomics so we specify `protocol=10x` to use the default value for 10X Genomics reads.

```
cd experiment
bin/physlr-make scaffolds lr=linkedreads ref=reference draft=draft protocol=10x
```

See the help page for further information.
`bin/physlr-make help`

# Output files

* `lr.physlr.physical-map.path`: Paths of barcodes (backbones).
* `lr.physlr.physical-map.ref.n10.paf.gz.*.pdf`: Various graphs showing the contiguity and correctness of the backbones with respect to the reference.
* `draft.physlr.fa`: Physlr scaffolded assembly using the physical map.
* `draft.physlr.quast.tsv`: Quast metrics comparing the Physlr scaffolded assembly against the reference.

# Acknowledgements

This projects uses:
* [btl_bloomfilter](https://github.com/bcgsc/btl_bloomfilter) BTL C/C++ Common bloom filters for bioinformatics projects implemented by Justin Chu
* [nthash](https://github.com/bcgsc/ntHash) rolling hash implementation by Hamid Mohamadi
* [readfq](https://github.com/lh3/readfq) Fast multi-line FASTA/Q reader API implemented by Heng Li
* [robin-map](https://github.com/Tessil/robin-map) C++ implementation of a fast hash map and hash set using robin hood hashing by Thibaut G.
