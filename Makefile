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

check: humanmt psitchensiscp fly

clean:
	rm -f humanmt/mt.physlr.tsv

################################################################################
# Human mitochondrion

# Test Phsylr using the human mitochondrion.
humanmt: \
	humanmt/mt.physlr.tsv

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
	f1chr4.physlr.overlap.backbone.label.n50.gv.pdf \
	f1chr4.physlr.overlap.backbone.path.fly.molecule.bed.png \
	f1chr4.physlr.overlap.molecules.M2.backbone.label.n50.gv.pdf \
	f1chr4.physlr.overlap.molecules.M2.backbone.path.fly.molecule.bed.png \
	f1chr2R.physlr.overlap.backbone.label.n50.gv.pdf \
	f1chr2R.physlr.overlap.backbone.path.fly.molecule.bed.png \
	f1chr2R.physlr.overlap.molecules.M2.backbone.label.n50.gv.pdf \
	f1chr2R.physlr.overlap.molecules.M2.backbone.path.fly.molecule.bed.png

# Download the fly genome from NCBI.
fly/fly.fa:
	mkdir -p $(@D)
	curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz | seqtk seq >$@

# Download the fly linked reads from 10x Genomics.
fly/f1.tar:
	mkdir -p $(@D)
	curl -o $@ http://s3-us-west-2.amazonaws.com/10x.files/samples/assembly/2.1.0/fly/fly_fastqs.tar

# Extract the reads that map to chromosome 4.
%.chr4.sortbxn.bam: %.sortbxn.bam
	samtools view -h $< | awk '/^@/ || $$3 == "NC_004353.4"' | samtools view -@$t -o $@

# Symlink the chromosome 4 reads.
f1chr4.fq.gz: fly/fly.f1.chr4.sortbxn.dropse.bx100-200.fq.gz
	ln -sf $< $@

# Extract the reads that map to chromosome 2R.
%.chr2R.sortbxn.bam: %.sortbxn.bam
	samtools view -h $< | awk '/^@/ || $$3 == "NT_033778.4"' | samtools view -@$t -o $@

# Symlink the chromosome 2R reads.
f1chr2R.fq.gz: fly/fly.f1.chr2R.sortbxn.dropse.bx100-200.fq.gz
	ln -sf $< $@

################################################################################
# Picea sitchensis plastid

# Test Phsylr using the Picea sitchensis plastid data.
psitchensiscp: \
	HYN5VCCXX_4cp.physlr.overlap.backbone.label.n50.gv.pdf \
	HYN5VCCXX_4cp.physlr.overlap.backbone.path.psitchensiscp.molecule.bed.png

# Download the Picea sitchensis plastid genome.
psitchensiscp/psitchensiscp.fa:
	mkdir -p $(@D)
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=KU215903.2&db=nucleotide&rettype=fasta' | seqtk seq >$@

# Symlink the plastid reads.
HYN5VCCXX_4cp.fq.gz: psitchensiscp/psitchensiscp.HYN5VCCXX_4.sortbxn.dropse.bx100-200.fq.gz
	ln -sf $< $@

################################################################################
# Filter reads

# Count the number of reads per barcode.
%.bx.tsv: %.fq.gz
	gunzip -c $< | sed -n 's/^.*BX:Z://p' \
	| mlr -p --otsvlite rename 1,BX then count-distinct -f BX -o Reads then sort -nr Reads >$@

# Determine a set of barcodes containing on-target molecules.
%.bx100-200.txt: %.fq.gz
	gunzip -c $< | sed -n 's/^.*BX:Z://p' \
	| mlr -p rename 1,BX then count-distinct -f BX then filter '$$count >= 100 && $$count < 200' then cut -f BX then sort -f BX >$@

# Separate a set of on-target reads.
%.bx100-200.fq.gz: %.fq.gz %.bx100-200.txt
	gunzip -c $< | paste - - - - | grep -Ff $*.bx100-200.txt | tr '\t' '\n' | $(gzip) >$@

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
# EMA

# Map linked reads to the draft genome using EMA.
# Filter out reads without barcodes.
%.$(lr).ema.sortn.bam: $(lr).fq.gz %.fa.bwt
	gunzip -c $< | paste - - - - - - - - | grep "BX:Z:" | tr '\t' '\n' \
	| $(time) ema align -t$t -r $*.fa -1 /dev/stdin | samtools view -@$t -h -F4 -o $@

################################################################################
# minimap2

# Align linked reads to a target genome.
%.$(lr).minimap2.sortn.bam: %.fa.bwt $(lr).fq.gz
	minimap2 -t$t -a -xsr -y $*.fa $(lr).fq.gz | samtools view -@$t -F4 -o $@

################################################################################
# samtools

# Index a FASTA file.
%.fa.fai: %.fa
	samtools faidx $<

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
%.pe.fq.gz: %.1.fq.gz %.2.fq.gz
	seqtk mergepe $^ | $(gzip) >$@

# Select the first read of the read pair.
%.1.fq.gz: %.fq.gz
	seqtk dropse $< | seqtk seq -1 | $(gzip) >$@

# Select the second read of the read pair.
%.2.fq.gz: %.fq.gz
	seqtk dropse $< | seqtk seq -2 | $(gzip) >$@

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
# Tigmint
as=0.65
dist=10000
nm=5
mapq=1
minsize=2000

# Create molecule extents BED.
%.a$(as).d$(dist).n$(nm).q$(mapq).s$(minsize).molecule.bed: %.sortbx.bam
	tigmint-molecule -a$(as) -n$(nm) -q$(mapq) -d$(dist) -s$(minsize) -o $@ $<

# Create molecule extents TSV.
%.a$(as).d$(dist).n$(nm).q$(mapq).s$(minsize).molecule.tsv: %.sortbx.bam
	tigmint-molecule -a$(as) -n$(nm) -q$(mapq) -d$(dist) -s$(minsize) --tsv -o $@ $<

################################################################################
# ntHash

# Count k-mers using ntCard.
%.ntcard_k32.hist: %.fq.gz
	ntcard -t$t -c1000 -k 32,64,96,128 -p $*.ntcard $<

# Convert a .hist to a .histo file for GenomeScope.
%.histo: %.hist
	sed -n 's/^f//p' $< | tr '\t' ' ' >$@

################################################################################
# Unicycler

# Assembled paired-end reads.
%.unicycler.gfa: %.1.fq.gz %.2.fq.gz
	unicycler -t$t --mode bold -o $*.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz
	ln -s $*.unicycler/assembly.gfa $@

################################################################################
# Bandage

# Plot the assembly graph using Bandage.
%.gfa.png: %.gfa
	Bandage image $< $@

# Plot the assembly graph using Bandage.
%.gfa.svg: %.gfa
	Bandage image $< $@

################################################################################
# Physlr

# Index a FASTA file.
%.physlr.tsv: %.fa
	PYTHONPATH=. bin/physlr indexfa -k$k -w$w $< >$@

# Index linked reads.
%.physlr.tsv: %.fq.gz
	gunzip -c $< | PYTHONPATH=. bin/physlr indexlr -k$k -w$w - >$@

# Count the frequency of the markers.
%.physlr.markers.tsv: %.physlr.tsv
	PYTHONPATH=. bin/physlr count-markers $< >$@
	
# Identify the overlapping markers of each pair of barcodes.
%.physlr.intersect.tsv: %.physlr.tsv
	PYTHONPATH=. bin/physlr intersect $< >$@

# Determine overlaps and output the graph in TSV.
%.physlr.overlap.tsv: %.physlr.tsv
	PYTHONPATH=. bin/physlr overlap $< >$@

# Determine the maximum spanning tree.
%.physlr.overlap.mst.tsv: %.physlr.overlap.tsv
	PYTHONPATH=. bin/physlr mst $< >$@

# Determine the backbone graph from the overlap TSV.
%.backbone.tsv: %.tsv
	PYTHONPATH=. bin/physlr backbone-graph $< >$@

# Determine the backbone path of the backbone graph.
%.backbone.path: %.backbone.tsv
	PYTHONPATH=. bin/physlr backbone $< >$@

# Determine the minimum tiling graph of the backbone graph.
%.physlr.overlap.backbone.tiling.tsv: %.physlr.overlap.backbone.tsv
	PYTHONPATH=. bin/physlr tiling-graph $< >$@

# Estimate the nubmer of molecules per barcode.
%.physlr.overlap.molecules.tsv: %.physlr.overlap.tsv
	PYTHONPATH=. bin/physlr molecules -n10 $< >$@

# Remove barcodes with more than one molecule.
%.physlr.overlap.molecules.M2.tsv: %.physlr.overlap.molecules.tsv
	PYTHONPATH=. bin/physlr filter -M2 $< >$@

# Convert a graph from TSV to GraphViz.
%.gv: %.tsv
	PYTHONPATH=. bin/physlr filter -Ogv $< >$@

# Extract a BED file of the backbone barcodes.
%.backbone.path.$(ref).molecule.bed: %.backbone.path $(ref)/$(ref).$(lr).a0.65.d10000.n5.q1.s2000.molecule.bed
	for i in $$(<$<); do grep $$i $(ref)/$(ref).$(lr).a0.65.d10000.n5.q1.s2000.molecule.bed || true; done >$@

# Plot a BED file.
%.bed.png: %.bed
	Rscript -e 'rmarkdown::render("plotbed.rmd", "html_document", "$*.plotbed.html", params = list(input_bed="$<"))'
	cp -a $*.plotbed_files/figure-html/plot-bed-1.png $@

################################################################################
# Bedtools

# Compute genome coverage.
%.bed.$(ref).cov: %.bed $(ref)/$(ref).fa.fai
	bedtools genomecov -max 1 -g $(ref)/$(ref).fa.fai -i $< \
	| awk '$$2 != 0 || $$5 != 1' >$@

################################################################################
# GraphViz

# Label the edges with edge weight.
%.label.gv: %.gv
	gvpr -c 'E { label = n }' $< >$@

# Filter a graph by edge weight.
%.n5.gv: %.gv
	gvpr 'E[n >= 5]' $< >$@

# Filter a graph by edge weight.
%.n10.gv: %.gv
	gvpr 'E[n >= 10]' $< >$@

# Filter a graph by edge weight.
%.n20.gv: %.gv
	gvpr 'E[n >= 20]' $< >$@

# Filter a graph by edge weight.
%.n50.gv: %.gv
	gvpr 'E[n >= 50]' $< >$@

# Layout and render an undirected graph to PDF.
%.gv.pdf: %.gv
	neato -Goverlap=scale -Gsize=100,100 -Tpdf -o $@ $<

# Layout and render an undirected graph to PNG.
%.gv.png: %.gv
	neato -Goverlap=scale -Tpng -o $@ $<
