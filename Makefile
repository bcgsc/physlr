# Size of a k-mer.
k=32

# Number of k-mers in a minimizer window.
w=32

.DELETE_ON_ERROR:
.SECONDARY:

all: check

check: mt/mt.fa.physlr.json

clean:
	rm -f mt/mt.fa.physlr.json

# Download the human mitochondrial genome.
mt/mt.fa:
	mkdir -p $(@D)
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=NC_012920.1&db=nucleotide&rettype=fasta' | seqtk seq >$@

# Physlr

# Index a FASTA file using Physlr.
%.fa.physlr.json: %.fa
	PYTHONPATH=. bin/physlr index -k$k -w$w $< >$@
