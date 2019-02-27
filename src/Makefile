CXXFLAGS=-g -O3 -Wall -Werror -std=c++11

all: physlr-indexlr

# Check the C++ source code for errors with clang-tidy.
lint:
	clang-tidy -warnings-as-errors='*' *.cc -- -std=c++11

check: check-physlr-index-fasta check-physlr-index-fastq

check-physlr-index-fasta: all
	./physlr-indexlr -k100 -w5 data/mt.fa | diff -q - data/mt.physlr.tsv.good

check-physlr-index-fastq: all
	./physlr-indexlr -k100 -w5 data/tiny.fq | diff -q - data/tiny.physlr.tsv.good

clean:
	rm -f physlr-indexlr
