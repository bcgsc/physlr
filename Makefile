# Size of a k-mer.
k=32

# Number of k-mers in a minimizer window.
w=32

# Number of threads.
t=16

# Compress in parallel.
gzip=pigz -p$t

# Python interpreter.
python=PYTHONPATH=. pypy3

# Reference genome.
ref=fly

# Linked reads.
lr=f1

# Draft genome assembly.
draft=f1.supernova.scaftigs

SHELL=bash -e -o pipefail
ifeq ($(shell zsh -e -o pipefail -c 'true' 2>/dev/null; echo $$?), 0)
# Set pipefail to ensure that all commands of a pipe succeed.
SHELL=zsh -e -o pipefail
# Report run time and memory usage with zsh.
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J
endif

.DELETE_ON_ERROR:
.SECONDARY:

all: lint

# Check the source code for errors with Pylint.
lint:
	pylint physlr

check: f1

################################################################################
# Fly
# See https://support.10xgenomics.com/de-novo-assembly/datasets/2.1.0/fly
# and https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/performance

# Assemble a physical map of fly chromosome 4.
f1chr4: f1chr4.physlr.stamp

# Assemble a physical map of fly chromosome 2R.
f1chr2R: f1chr2R.physlr.stamp

# Assemble a physical map of the fly genome.
f1: \
	f1.n100-2000.physlr.overlap.n118.mol.backbone.path.fly.molecule.bed.fly.cov.tsv \
	f1.n100-2000.physlr.overlap.n118.mol.backbone.path.fly.molecule.bed.pdf \
	f1.n100-2000.physlr.overlap.n118.mol.backbone.label.gv.pdf \
	f1.n100-2000.physlr.overlap.n118.mol.backbone.map.f1.supernova.scaftigs.n10.sort.best.bed.pdf \
	f1.n100-2000.physlr.overlap.n118.mol.backbone.map.f1.supernova.scaftigs.n10.sort.best.bed.path.fly.paf.pdf \
	f1.n100-2000.physlr.overlap.n118.mol.backbone.map.f1.supernova.scaftigs.n10.sort.best.bed.path.quast.tsv

# Download the fly genome from NCBI.
fly/fly.fa:
	mkdir -p $(@D)
	curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz | seqtk seq >$@

# Download the fly annotation from NCBI.
fly/fly.gff:
	mkdir -p $(@D)
	curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz | gunzip -c >$@

# Download the Supernova assembly of the linked reads from 10x Genomics.
f1.supernova.fa:
	curl http://cf.10xgenomics.com/samples/assembly/2.1.0/fly/fly_pseudohap.fasta.gz | gunzip -c >$@

# Download the fly linked reads from 10x Genomics.
fly/f1.tar:
	mkdir -p $(@D)
	curl -o $@ http://s3-us-west-2.amazonaws.com/10x.files/samples/assembly/2.1.0/fly/fly_fastqs.tar

# Extract the tar file of fly FASTQ reads.
fly/f1.fq.gz: fly/f1.tar
	tar --wildcards -Oxf fly/f1.tar 'fly/H3C7LDMXX/read-RA*.fastq.gz' >$@

# Symlink the fly reads.
f1.fq.gz: fly/f1.fq.gz
	ln -s $< $@

# Download the fly linked reads from 10x Genomics.
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

# Symlink the subsampled chromosome 2R reads.
f1chr2R-bx%.fq.gz: fly/fly.f1.chr2R.sortbxn.dropse.bx%.fq.gz
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

# Determine a set of barcodes containing on-target molecules.
%.bx200-300.txt: %.fq.gz
	gunzip -c $< | sed -n 's/^.*BX:Z://p' \
	| mlr -p rename 1,BX then count-distinct -f BX then filter '$$count >= 200 && $$count < 300' then cut -f BX then sort -f BX >$@

# Separate a set of on-target reads.
%.bx200-300.fq.gz: %.fq.gz %.bx200-300.txt
	gunzip -c $< | paste - - - - | grep -Ff $*.bx200-300.txt | tr '\t' '\n' | $(gzip) >$@

# Determine a set of barcodes containing on-target molecules.
%.bx300-400.txt: %.fq.gz
	gunzip -c $< | sed -n 's/^.*BX:Z://p' \
	| mlr -p rename 1,BX then count-distinct -f BX then filter '$$count >= 300 && $$count < 400' then cut -f BX then sort -f BX >$@

# Separate a set of on-target reads.
%.bx300-400.fq.gz: %.fq.gz %.bx300-400.txt
	gunzip -c $< | paste - - - - | grep -Ff $*.bx300-400.txt | tr '\t' '\n' | $(gzip) >$@

# Determine a set of barcodes containing on-target molecules.
%.bx400-500.txt: %.fq.gz
	gunzip -c $< | sed -n 's/^.*BX:Z://p' \
	| mlr -p rename 1,BX then count-distinct -f BX then filter '$$count >= 400 && $$count < 500' then cut -f BX then sort -f BX >$@

# Separate a set of on-target reads.
%.bx400-500.fq.gz: %.fq.gz %.bx400-500.txt
	gunzip -c $< | paste - - - - | grep -Ff $*.bx400-500.txt | tr '\t' '\n' | $(gzip) >$@

# Determine a set of barcodes containing on-target molecules.
%.bx100-500.txt: %.fq.gz
	gunzip -c $< | sed -n 's/^.*BX:Z://p' \
	| mlr -p rename 1,BX then count-distinct -f BX then filter '$$count >= 100 && $$count < 500' then cut -f BX then sort -f BX >$@

# Separate a set of on-target reads.
%.bx100-500.fq.gz: %.fq.gz %.bx100-500.txt
	gunzip -c $< | paste - - - - | grep -Ff $*.bx100-500.txt | tr '\t' '\n' | $(gzip) >$@

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

# Align a FASTA file to the reference genome and produce a PAF file.
%.$(ref).paf.gz: $(ref)/$(ref).fa %.fa
	$(time) minimap2 -t$t -xasm10 $^ | $(gzip) >$@

################################################################################
# miniasm

# Draw a dot plot of a PAF file.
# Skip alignments to non-chromosomal sequences.
%.paf.ps: %.paf.gz
	gunzip -c $< | grep -v NW_ | minidot /dev/stdin >$@

# Convert Postscript to PDF
%.pdf: %.ps
	ps2pdf $< $@

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

# Cut scaffolds at Ns to create scaftigs.
%.scaftigs.fa: %.fa
	seqtk cutN -n1 $< | tr :- _ | seqtk seq >$@

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
	$(python) bin/physlr indexfa -k$k -w$w $< >$@

# Index linked reads.
%.physlr.tsv: %.fq.gz
	gunzip -c $< | $(python) bin/physlr indexlr -k$k -w$w - >$@

# Count the frequency of the markers.
%.physlr.markers.tsv: %.physlr.tsv
	$(python) bin/physlr count-markers $< >$@
	
# Identify the overlapping markers of each pair of barcodes.
%.physlr.intersect.tsv: %.physlr.tsv
	$(python) bin/physlr intersect $< >$@

# Filter barcodes by number of markers.
%.n100-2000.physlr.tsv: %.physlr.tsv
	$(python) bin/physlr filter-barcodes -n100 -N2000 $< >$@

# Determine overlaps and output the graph in TSV.
%.physlr.overlap.tsv: %.physlr.tsv
	$(python) bin/physlr overlap -n10 $< >$@

# Determine the maximum spanning tree.
%.physlr.overlap.mst.tsv: %.physlr.overlap.tsv
	$(python) bin/physlr mst $< >$@

# Separate a graph into its biconnected components by removing its cut vertices.
%.bic.tsv: %.tsv
	$(python) bin/physlr biconnected-components $< >$@

# Determine the backbone graph from the overlap TSV.
%.backbone.tsv: %.tsv
	$(python) bin/physlr backbone-graph $< >$@

# Determine the backbone path of the backbone graph.
%.path: %.tsv
	$(python) bin/physlr backbone $< >$@

# Flesh out the backbone path
%.backbone.fleshed.path: %.tsv %.backbone.path
	$(python) bin/physlr flesh-backbone --min-component-size=50 $< $*.backbone.path >$@

# Determine the minimum tiling graph of the backbone graph.
%.backbone.tiling.tsv: %.backbone.tsv
	$(python) bin/physlr tiling-graph $< >$@

# Map the barcodes to the backbone graph.
%.backbone.map.$(lr).n100-2000.n10.bed: %.backbone.tsv $(lr).n100-2000.physlr.tsv
	$(python) bin/physlr map -n10 $^ $(lr).n100-2000.physlr.tsv >$@

# Map the draft assembly to the backbone graph.
%.overlap.n50.mol.backbone.map.$(draft).n10.bed: %.overlap.n50.mol.backbone.tsv %.tsv $(draft).physlr.tsv
	$(python) bin/physlr map -n10 $^ >$@

# Map the draft assembly to the backbone graph.
%.overlap.n118.mol.backbone.map.$(draft).n10.bed: %.overlap.n118.mol.backbone.tsv %.tsv $(draft).physlr.tsv
	$(python) bin/physlr map -n10 $^ >$@

# Filter a BED file by score.
%.n100.bed: %.n10.bed
	awk '$$5 >= 100' $< >$@

# Estimate the number of molecules per barcode.
%.physlr.overlap.n20.countmol.tsv: %.physlr.overlap.tsv
	$(python) bin/physlr count-molecules -n20 $< >$@

# Remove barcodes with more than one molecule.
%.physlr.overlap.molecules.M2.tsv: %.physlr.overlap.molecules.tsv
	$(python) bin/physlr filter -M2 $< >$@

# Filter edges n >= 10 using Miller.
%.n10.tsv: %.tsv
	mlr --tsvlite filter '$$n >= 10' $< >$@

# Filter edges n >= 20 using Miller.
%.n20.tsv: %.tsv
	mlr --tsvlite filter '$$n >= 20' $< >$@

# Filter edges n >= 50 using Miller.
%.n50.tsv: %.tsv
	mlr --tsvlite filter '$$n >= 50' $< >$@

# Filter edges n >= 100 using Miller.
%.n100.tsv: %.tsv
	mlr --tsvlite filter '$$n >= 100' $< >$@

# Filter edges n >= 118 using Miller.
%.n118.tsv: %.tsv
	mlr --tsvlite filter '$$n >= 118' $< >$@

# Separate barcodes into molecules.
%.mol.tsv: %.tsv
	$(python) bin/physlr molecules -t$t $< >$@

# Convert a graph from TSV to GraphViz.
# Filter out small components.
%.gv: %.tsv
	$(python) bin/physlr filter --min-component-size=50 -Ogv $< >$@

# Extract a BED file of the backbone barcodes.
# Filter out small components.
%.path.$(ref).molecule.bed: $(ref)/$(ref).$(lr).a0.65.d10000.n5.q1.s2000.molecule.bed %.path
	$(python) bin/physlr filter-bed --min-component-size=50 $^ >$@

# Extract a BED file of the fleshed-out backbone barcodes.
# Filter out small components.
%.backbone.fleshed.path.$(ref).molecule.bed: $(ref)/$(ref).$(lr).a0.65.d10000.n5.q1.s2000.molecule.bed %.backbone.fleshed.path
	$(python) bin/physlr filter-bed --min-component-size=50 $^ >$@

# Sort a BED file.
%.sort.bed: %.bed
	sort -k1,1n -k1,1 -k2,2n -k3,3n -k5,5nr -k4,4 $< >$@

# Keep the best record at each position.
%.sort.best.bed: %.sort.bed
	awk '{ keep = $$1 " " $$2 " " $$3 != x; x = $$1 " " $$2 " " $$3 } keep' $< >$@

# Extract scaffolds paths from a BED file.
%.bed.path: %.bed
	$(python) bin/physlr bed-to-path $^ >$@

# Produce sequences in FASTA format from paths.
%.$(draft).n10.sort.best.bed.path.fa: $(draft).fa %.$(draft).n10.sort.best.bed.path
	$(python) bin/physlr path-to-fasta --min-length 100000 $^ >$@

# Plot a BED file.
%.bed.pdf: %.bed
	Rscript -e 'rmarkdown::render("plotbed.rmd", "html_document", "$*.plotbed.html", params = list(input_bed="$<"))'

# Assemble a physical map.
%.physlr.stamp: \
		%.n100-2000.physlr.overlap.n50.mol.backbone.path.$(ref).molecule.bed.$(ref).cov.tsv \
		%.n100-2000.physlr.overlap.n50.mol.backbone.path.$(ref).molecule.bed.pdf \
		%.n100-2000.physlr.overlap.n50.mol.backbone.label.gv.pdf \
		%.n100-2000.physlr.overlap.n50.mol.backbone.fleshed.path.$(ref).molecule.bed.pdf \
		%.n100-2000.physlr.overlap.n50.mol.backbone.map.f1.supernova.scaftigs.n10.sort.best.bed.pdf \
		%.n100-2000.physlr.overlap.n50.mol.backbone.map.f1.supernova.scaftigs.n10.sort.best.bed.path.fly.paf.pdf
	touch $@

################################################################################
# Bedtools

# Compute genome coverage.
%.bed.$(ref).cov.tsv: %.bed $(ref)/$(ref).fa.fai
	grep -v NA $< | sort -k1,1 -k2,2n -k3,3n | bedtools genomecov -max 1 -g $(ref)/$(ref).fa.fai -i - | awk '$$2 != 0 || $$5 != 1' >$@

################################################################################
# QUAST

# Calculate assembly contiguity and correctness metrics using QUAST.
%.quast.tsv: %.fa $(ref)/$(ref).fa $(ref)/$(ref).gff
	quast-lg -t$t -es --fast --large --scaffold-gap-max-size 100000 --min-identity 95 -R $(ref)/$(ref).fa -o $*.quast $<
	cp $*.quast/transposed_report.tsv $@

# Aggregate QUAST metrics.
f1.quast.tsv: \
		f1.abyss.scaftigs.quast.tsv \
		f1.abyss.quast.tsv \
		f1.supernova.scaftigs.quast.tsv \
		f1.supernova.quast.tsv \
		f1.n100-2000.physlr.overlap.n118.mol.backbone.map.f1.abyss.n10.sort.best.bed.path.quast.tsv \
		f1.n100-2000.physlr.overlap.n118.mol.backbone.map.f1.supernova.scaftigs.n10.sort.best.bed.path.quast.tsv
	mlr --tsvlite cut -x -f NG75,NGA75,LG75,LGA75 $^ >$@

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

################################################################################
# RMarkdown reports

# Compare assembly metrics.
%.quast.html: %.quast.tsv
	Rscript -e 'rmarkdown::render("quast.rmd", "html_document", "$*.quast.html", params = list(input_tsv="$<", output_tsv="$*.quast.table.tsv"))'
