# Size of a k-mer.
k=32

# Number of k-mers in a minimizer window.
w=32

# Number of threads.
t=16

.DELETE_ON_ERROR:
.SECONDARY:

all: lint

# Check the source code for errors with Pylint.
lint:
	pylint physlr

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

################################################################################
# Fly
# See https://support.10xgenomics.com/de-novo-assembly/datasets/2.1.0/fly
# and https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/performance

# Download the fly linked reads from 10x Genomics.
fly/fly.tar:
	mkdir -p $(@D)
	curl -o $@ http://s3-us-west-2.amazonaws.com/10x.files/samples/assembly/2.1.0/fly/fly_fastqs.tar

################################################################################
# Picea sitchensis

# Download the Picea sitchensis plastid FASTA.
psitchensiscp/psitchensiscp.fa:
	mkdir -p $(@D)
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=KU215903.2&db=nucleotide&rettype=fasta' | seqtk seq >$@

################################################################################
# BWA

# Index the target genome.
%.fa.bwt: %.fa
	bwa index $<

# Align linked reads to the draft genome and do not sort.
%.$(lr).sortn.bam: %.fa.bwt $(lr).fq.gz
	bwa mem -t$t -pC $*.fa $(lr).fq.gz | samtools view -@$t -h -F4 -o $@

################################################################################
# samtools

# Convert a BAM file to FASTQ.
%.sortn.bam.fq.gz: %.sortn.bam
	samtools fastq -@$t -TBX $< | $(gzip) -p$t >$@
