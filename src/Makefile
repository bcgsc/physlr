CXXFLAGS=-g -O3 -Wall -Werror -std=c++11 -pthread -I./include -I./vendor
LDLIBS=-lz

all: physlr-filter-barcodes physlr-indexlr

# Check the C++ source code for errors with clang-tidy.
lint:
	clang-tidy -warnings-as-errors='*' include/*.h *.cc -- -std=c++11 -x c++ -I./include -I./vendor

check: check-physlr-index-fasta check-physlr-index-fastq check-physlr-filter-barcodes

check-physlr-index-fasta: all
	./physlr-indexlr -t16 -k100 -w5 data/mt.fa | diff -q - data/mt.physlr.tsv.good

check-physlr-index-fastq: all
	./physlr-indexlr -t16 -k100 -w5 data/tiny.fq | diff -q - data/tiny.physlr.tsv.good

check-physlr-filter-barcodes: all
	./physlr-filter-barcodes head50.tsv -o cpp.head50.n150-900.tsv -n150 -N900
	./checkMx.sh cpp.head50.n150-900.tsv py.head50.n150-900.tsv 

clean:
	rm -f physlr-indexlr physlr-filter-barcodes
