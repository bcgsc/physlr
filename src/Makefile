CXXFLAGS=-g -O3 -Wall -Werror -Wno-error=unknown-pragmas -std=c++11 -pthread -I./include -I./vendor -fopenmp
LDLIBS=-lz

all: physlr-indexlr physlr-filter-barcodes physlr-overlap physlr-filter-bxmx

clean:
	rm -f physlr-indexlr physlr-filter-barcodes physlr-filter-bxmx physlr-overlap

# Check the C++ source code for errors.
lint: clang-format clang-tidy

# Check the C++ source code for errors with clang-tidy.
clang-tidy:
	clang-tidy -warnings-as-errors='*' include/*.h *.cc -- -std=c++11 -x c++ -I./include -I./vendor

# Check the C++ source code for white-space errors with clang-format.
clang-format:
	for i in *.cc include/*.h; do clang-format -style=file $$i >$$i.fixed; done
	for i in *.cc include/*.h; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi

# Fix white-space errors with clang-format.
fix:
	for i in *.cc; do clang-format -i -style=file $$i; done

check: check-physlr-index-fasta check-physlr-index-pos-fasta check-physlr-index-fastq check-physlr-filter-barcodes check-physlr-filter-bxmx-first check-physlr-filter-bxmx-second

check-physlr-index-fasta: all
	./physlr-indexlr -t16 -k100 -w5 data/mt.fa | diff -q - data/mt.physlr.tsv.good

check-physlr-index-pos-fasta: all
	./physlr-indexlr -t16 -k100 -w5 --pos data/mt.fa | diff -q - data/mt.pos.physlr.tsv.good

check-physlr-index-fastq: all
	./physlr-indexlr -t16 -k100 -w5 data/tiny.fq | diff -q - data/tiny.physlr.tsv.good

check-physlr-filter-barcodes: all
	./physlr-filter-barcodes data/f1chr4-head50.tsv -o data/f1chr4-head50.n100-1000.tsv -n100 -N1000
	./data/checkMx.sh data/f1chr4-head50.n100-1000.good.tsv data/f1chr4-head50.n100-1000.tsv

check-physlr-filter-bxmx-first: all
	./physlr-filter-barcodes data/f1.tiny.indexlr.physlr.tsv -n2 -N10 -o fb-cpp
	env PYTHONPATH=.. ../bin/physlr filter-minimizers -V0 -C4 fb-cpp >fm-py
	./physlr-filter-bxmx data/f1.tiny.indexlr.physlr.tsv -o ff-cpp -n2 -N10 -C4
	./data/checkMx.sh fm-py ff-cpp
	rm fb-cpp fm-py ff-cpp

check-physlr-filter-bxmx-second: all
	./physlr-filter-barcodes data/f1.tiny.indexlr.physlr.tsv -n2 -N10 -o fb-cpp
	env PYTHONPATH=.. ../bin/physlr filter-minimizers -V0 fb-cpp >fm-py
	./physlr-filter-bxmx data/f1.tiny.indexlr.physlr.tsv -o ff-cpp -n2 -N10
	./data/checkMx.sh fm-py ff-cpp
	rm fb-cpp fm-py ff-cpp
