# Size of a k-mer.
k=32

# Number of k-mers in a minimizer window.
w=32

# Number of threads.
t=16

# Compress in parallel.
gzip=pigz -p$t

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: fly humanmt psitchensiscp

all: lint

# Check the source code for errors with Pylint.
lint:
	pylint physlr

check: humanmt

clean:
	rm -f humanmt/mt.physlr.json

################################################################################
# Human mitochondrion

# Test Phsylr using the human mitochondrion.
humanmt: \
	humanmt/mt.physlr.json

# Download the human mitochondrial genome.
humanmt/mt.fa:
	mkdir -p $(@D)
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=NC_012920.1&db=nucleotide&rettype=fasta' | seqtk seq >$@

################################################################################
# Fly
# See https://support.10xgenomics.com/de-novo-assembly/datasets/2.1.0/fly
# and https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/performance

# Test Phsylr using the fly data.
fly: \
	fly/fly.f1.sortn.bam \
	fly/fly.f1.sort.bam.bai \
	fly/fly.f1.sortbx.bam \
	fly/fly.f1.sortbxn.bam

# Download the fly genome from NCBI.
fly/fly.fa:
	curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz | seqtk seq >$@

# Download the fly linked reads from 10x Genomics.
fly/f1.tar:
	mkdir -p $(@D)
	curl -o $@ http://s3-us-west-2.amazonaws.com/10x.files/samples/assembly/2.1.0/fly/fly_fastqs.tar

################################################################################
# Picea sitchensis plastid

# Test Phsylr using the Picea sitchensis plastid data.
psitchensiscp: \
	psitchensiscp/psitchensiscp.HYN5VCCXX_4.sortn.bam \
	psitchensiscp/psitchensiscp.HYN5VCCXX_4.sort.bam.bai \
	psitchensiscp/psitchensiscp.HYN5VCCXX_4.sortbx.bam \
	psitchensiscp/psitchensiscp.HYN5VCCXX_4.sortbxn.bam

# Download the Picea sitchensis plastid genome.
psitchensiscp/psitchensiscp.fa:
	mkdir -p $(@D)
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=KU215903.2&db=nucleotide&rettype=fasta' | seqtk seq >$@

# Symlink the plastid reads.
psitchensiscp/HYN5VCCXX_4cp.fq.gz: psitchensiscp/psitchensiscp.HYN5VCCXX_4.sortbxn.dropse.fq.gz
	ln -sf $< $@

################################################################################
# Trimadap

# Trim adapter sequences using trimadap.
%.trimadap.fq.gz: %.fq.gz
	trimadap-mt -p$t -t1 $< | sed 's/^X$$/N/' | $(gzip) >$@

################################################################################
# BWA

# Index the target genome.
%.fa.bwt: %.fa
	bwa index $<

# Align linked reads to a target genome.
%.$(lr).sortn.bam: %.fa.bwt $(lr).fq.gz
	bwa mem -t$t -pC $*.fa $(lr).fq.gz | samtools view -@$t -F4 -o $@

################################################################################
# minimap2

# Align linked reads to a target genome.
%.$(lr).minimap2.sortn.bam: %.fa.bwt $(lr).fq.gz
	minimap2 -t$t -a -xsr -y $*.fa $(lr).fq.gz | samtools view -@$t -F4 -o $@

################################################################################
# samtools

# Sort a BAM file by position.
%.sort.bam: %.sortn.bam
	samtools sort -@$t -T$$(mktemp -u -t $(@F).XXXXXX) -o $@ $<

# Sort a BAM file by BX tag and position.
%.sortbx.bam: %.sortn.bam
	samtools sort -@$t -tBX -T$$(mktemp -u -t $(@F).XXXXXX) -o $@ $<

# Sort a BAM file by BX tag and query name.
%.sortbxn.bam: %.sortn.bam
	samtools sort -@$t -tBX -n -T$$(mktemp -u -t $(@F).XXXXXX) -o $@ $<

# Index a BAM file.
%.bam.bai: %.bam
	samtools index -@$t $<

# Convert a BAM file to FASTQ.
%.sortbxn.fq.gz: %.sortbxn.bam
	samtools fastq -@$t -TBX $< | $(gzip) >$@

################################################################################
# seqtk

# Drop single-end reads.
%.dropse.fq.gz: %.fq.gz
	seqtk dropse $< | $(gzip) >$@

# Merge paired-end reads.
%.fq.gz: %.1.fq.gz %.2.fq.gz
	seqtk mergepe $^ | $(gzip) >$@

################################################################################
# EMA

# Download the barcode white list.
4M-with-alts-february-2016.txt:
	curl -o $@ https://raw.githubusercontent.com/10XGenomics/supernova/master/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt

# Count barcodes.
%.ema-ncnt: %.fq.gz 4M-with-alts-february-2016.txt
	ema count -w 4M-with-alts-february-2016.txt -o $* $<

# Extract the barcode to BX:Z tag using ema preproc.
%.bx.fq.gz: %.fq.gz %.ema-ncnt
	gunzip -c $< | ema preproc -t$t -b -n1 -w 4M-with-alts-february-2016.txt -o $*.ema $*.ema-ncnt
	$(gzip) <$*.ema/ema-bin-000 >$@
	rm -rf $*.ema

# Align linked reads to the draft genome using EMA and sort by position.
%.$(lr).bx.ema.sort.bam: $(lr).bx.fq.gz %.fa.bwt
	$(time) ema align -t$t -r $*.fa -1 $< \
	| samtools view -@$t -u -F4 \
	| samtools sort -@$t -T$$(mktemp -u -t $@.XXXXXX) -o $@

################################################################################
# ntHash

# Count k-mers using ntCard.
%.ntcard_k32.hist: %.fq.gz
	ntcard -t$t -c1000 -k 32,64,96,128 -p $*.ntcard $<

# Convert a .hist to a .histo file for GenomeScope.
%.histo: %.hist
	sed -n 's/^f//p' $< | tr '\t' ' ' >$@

################################################################################
# Physlr

# Index a FASTA file using Physlr.
%.physlr.json: %.fa
	PYTHONPATH=. bin/physlr indexfa -k$k -w$w $< >$@

# Index linked reads using Physlr.
%.physlr.json: %.fq.gz
	gunzip -c $< | PYTHONPATH=. bin/physlr indexlr -k$k -w$w - >$@
