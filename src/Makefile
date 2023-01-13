CXXFLAGS=-g -O3 -Wall -Werror -Wno-error=unknown-pragmas -std=c++11 -pthread -I./include -I./vendor -fopenmp
LDLIBS=-lz -l btllib
PREFIX=$(HOME)/.local

.PHONY: all clean check install lint clang-format clang-tidy fix check-physlr-index-fasta check-physlr-index-pos-fasta check-physlr-index-fastq check-physlr-filter-barcodes check-physlr-filter-bxmx-first check-physlr-filter-bxmx-second check-physlr-index-stLFR check-physlr-index-stLFR-10X check-physlr-makebf check-physlr-overlap check-physlr-split-minimizers

all: physlr-filter-barcodes physlr-overlap physlr-filter-bxmx physlr-makebf physlr-split-minimizers # physlr-molecules

clean:
	rm -f physlr-filter-barcodes physlr-filter-bxmx physlr-overlap physlr-makebf physlr-split-minimizers # physlr-molecules

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

check: check-physlr-index-fasta check-physlr-index-pos-fasta check-physlr-index-fastq check-physlr-filter-barcodes check-physlr-filter-bxmx-first check-physlr-filter-bxmx-second check-physlr-index-stLFR check-physlr-index-stLFR-10X check-physlr-makebf check-physlr-overlap # check-physlr-molecules

check-physlr-filter-barcodes: all
	./physlr-filter-barcodes data/f1chr4-head50.tsv -o data/f1chr4-head50.n100-1000.tsv -n100 -N1000
	./data/checkMx.sh data/f1chr4-head50.n100-1000.good.tsv data/f1chr4-head50.n100-1000.tsv

check-physlr-filter-bxmx-first: all
	./physlr-filter-bxmx data/f1.tiny.indexlr.n2-10.physlr.tsv -o f1.tiny.indexlr.n2-10.c2-4.physlr.tsv -n2 -N10 -C4
	./data/checkMx.sh f1.tiny.indexlr.n2-10.c2-4.physlr.tsv data/f1.tiny.indexlr.n2-10.c2-4.physlr.tsv.good

check-physlr-filter-bxmx-second: all
	./physlr-filter-bxmx data/f1.tiny.indexlr.n2-10.physlr.tsv -o f1.tiny.indexlr.n2-10.c2-x.physlr.tsv -n2 -N10
	./data/checkMx.sh f1.tiny.indexlr.n2-10.c2-x.physlr.tsv data/f1.tiny.indexlr.n2-10.c2-x.physlr.tsv.good

#check-physlr-index-stLFR: all
#	indexlr --bx -t16 -k8 -w1 data/stLFR.tiny.fq | diff -q - data/stLFR.tiny.physlr.tsv.good

#check-physlr-index-stLFR-10X: all
#	indexlr --bx -t16 -k8 -w1 data/stLFR.10X.combined.tiny.fq | diff -q - data/stLFR.10X.combined.tiny.physlr.tsv.good

check-physlr-makebf: all
	./physlr-makebf -t16 -k8 -b1000 -o tiny.bf data/tiny_k8.rep
	indexlr --bx -t16 -k8 -w1 -r tiny.bf data/stLFR.tiny.fq | diff -q - data/stLFR.tiny.filtered.physlr.tsv.good

check-physlr-overlap: all
	./physlr-overlap -t4 -m1 data/tiny.overlap.input.tsv | sort --ignore-case | diff -q - data/tiny.physlr.overlap.n1.sorted.tsv
	./physlr-overlap -t4 -m2 data/tiny.overlap.input.tsv | sort --ignore-case | diff -q - data/tiny.physlr.overlap.n2.sorted.tsv

# check-physlr-molecules: all
# 	./physlr-molecules -s bc data/tiny.mol.input.tsv | diff -q - data/tiny.mol.tsv.good

check-physlr-split-minimizers: all
	./physlr-split-minimizers -t4 data/tiny.split-minimizers.ext.mol.tsv data/tiny.split-minimizers.physlr.tsv | sort --ignore-case |diff -q - data/tiny.split-minimizers.ext.mol.split.tsv.good

install: physlr-filter-barcodes physlr-overlap physlr-filter-bxmx physlr-makebf physlr-split-minimizers # physlr-molecules
	install -d $(DESTDIR)$(PREFIX)/bin
	install $^ $(DESTDIR)$(PREFIX)/bin
