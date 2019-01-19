
# physlr-indexlr

### Compile
```sh
make
```
### Usage
```
physlr-indexlr  -k K -w W [-v] [-o file] file...

  -k K       use K as k-mer size
  -w W       use W as sliding-window size
  -v         enable verbose output
  -o file    write output to file, default is stdout
  --help     display this help and exit
  file       space separated list of FASTQ files
  ```
### Example
```sh
./physlr-indexlr -k100 -w5 data/tiny.fq >data/tiny.physlr.tsv
```
### Output
Each line of output is a barcode followed by a tab, and then a list of space-separated minimizers.
```
barcode	minimizers...
```
### Test
```sh
make check
```
