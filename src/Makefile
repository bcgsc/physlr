CXXFLAGS=-g -O3 -Wall -Werror -std=c++11 -pthread -I./include -I./vendor -fopenmp
LDLIBS=-lz

all: physlr-filter-barcodes physlr-indexlr overlap/physlr-overlap

# Check the C++ source code for errors with clang-tidy.
lint:
	clang-tidy -warnings-as-errors='*' include/*.h *.cc overlap/*.cc -- -std=c++11 -x c++ -I./include -I./vendor

check: check-physlr-index-fasta check-physlr-index-fastq check-physlr-filter-barcodes

check-physlr-index-fasta: all
	./physlr-indexlr -t16 -k100 -w5 data/mt.fa | diff -q - data/mt.physlr.tsv.good

check-physlr-index-fastq: all
	./physlr-indexlr -t16 -k100 -w5 data/tiny.fq | diff -q - data/tiny.physlr.tsv.good

check-physlr-filter-barcodes: all
	./physlr-filter-barcodes data/f1chr4-head50.tsv -o data/f1chr4-head50.n100-1000.tsv -n100 -N1000
	./data/checkMx.sh data/f1chr4-head50.n100-1000.good.tsv data/f1chr4-head50.n100-1000.tsv 
clean:
	rm -f physlr-indexlr physlr-filter-barcodes physlr-overlap
