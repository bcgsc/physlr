CXXFLAGS=-g -O3 -Wall -Werror -std=c++11

all: physlr-indexlr

check:
	./physlr-indexlr -k100 -w5 data/tiny.fq | diff -q - data/tiny.physlr.tsv.good

clean:
	rm -f physlr-indexlr
