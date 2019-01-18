# physlr-indexlr

### Compile
```sh
make
```
### Usage
```sh
./physlr-indexlr -k K -w W [-v] file...
```
### Run
```sh
./physlr-indexlr -k100 -w5 data/tiny.fq >out
```
### Output
Each line of output is a barcode followed by a tab, and then a list of space-separated minimizers.
```
barcode   minimizers...
```
