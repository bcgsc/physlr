# Size of a k-mer.
k=32

# Number of k-mers in a minimizer window.
w=32

# Number of threads.
t=16

# Path to the Physlr project.
physlr_path=$(shell dirname $$(dirname $(realpath $(MAKEFILE_LIST))))

# Path to the Physlr executable.
bin=$(physlr_path)/bin

# Python interpreter.
python=env PYTHONPATH=$(physlr_path) pypy3

# Reference genome.
ref=Error_missing_parameter_ref

# Linked reads.
lr=Error_missing_parameter_lr

# Name of the dataset.
ifneq ($(filter z11%, $(ref)),)
name=fish
else
ifneq ($(filter b73%, $(ref)),)
name=maize
else
name=$(ref)
endif
endif

# Draft genome assembly.
draft=$(name).supernova.scaftigs

# Minimum number of overlapping minimizsers
n=50

# Minimum component size for backbones
min_component_size=50

# Level of verbosity
# 0: silent
# 1: periodic messages
# 2: progress bar
# 3: verbose messages
V=1

SHELL=bash -e -o pipefail
ifeq ($(shell zsh -e -o pipefail -c 'true' 2>/dev/null; echo $$?), 0)
# Set pipefail to ensure that all commands of a pipe succeed.
SHELL=zsh -e -o pipefail
# Report run time and memory usage with zsh.
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J
endif

# Record run time and memory usage in a file using GNU time.
ifneq ($(shell command -v gtime),)
time=command gtime -v -o $@.time
else
time=command time -v -o $@.time
endif

# Compress in parallel.
ifneq ($(shell command -v pigz),)
gzip=pigz -p$t
zcat=pigz -p$t -cd
else
gzip=gzip
zcat=gunzip -c
endif

.DELETE_ON_ERROR:
.SECONDARY:

.PHONY: f1chr4 f1chr2R f1 fishchr25 fish
all: f1chr4 f1chr2R f1 fishchr25 fish

################################################################################
# Fly
# See https://support.10xgenomics.com/de-novo-assembly/datasets/2.1.0/fly
# and https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/performance

# Assemble a physical map of fly chromosome 4.
# ref=fly lr=f1chr4 draft=f1.supernova.scaftigs n=50
f1chr4: \
	f1chr4.k32-w32.n100-1000.physlr.mxperbx.tsv.pdf \
	f1chr4.k32-w32.n100-1000.physlr.depth.tsv.pdf \
	f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.path.$(ref).molecule.bed.pdf \
	f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(ref).n10.paf.gz.pdf \
	f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(draft).n10.sort.best.bed.pdf \
	f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(draft).n10.sort.best.bed.path.$(ref).minidot.pdf \
	f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.label.gv.pdf \
	f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(ref).n10.ann.gv.squish.pdf

# Assemble a physical map of fly chromosome 2R.
# ref=fly lr=f1chr2R draft=f1.supernova.scaftigs n=50
f1chr2R: \
	f1chr2R.k32-w32.n100-1000.physlr.mxperbx.tsv.pdf \
	f1chr2R.k32-w32.n100-1000.physlr.depth.tsv.pdf \
	f1chr2R.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.path.$(ref).molecule.bed.pdf \
	f1chr2R.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(ref).n10.paf.gz.pdf \
	f1chr2R.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(draft).n10.sort.best.bed.pdf \
	f1chr2R.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(draft).n10.sort.best.bed.path.$(ref).minidot.pdf \
	f1chr2R.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.label.gv.pdf \
	f1chr2R.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(ref).n10.ann.gv.squish.pdf

# Assemble a physical map of the fly genome.
# ref=fly lr=f1 draft=f1.supernova.scaftigs n=25
f1: \
	f1.k32-w32.n100-1000.physlr.mxperbx.tsv.pdf \
	f1.k32-w32.n100-1000.physlr.depth.tsv.pdf \
	f1.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.path.$(ref).molecule.bed.pdf \
	f1.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(ref).n10.paf.gz.pdf \
	f1.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(draft).n10.sort.best.bed.pdf \
	f1.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(draft).n10.sort.best.bed.path.$(ref).minidot.pdf \
	f1.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(draft).n10.sort.best.bed.path.quast.tsv \
	f1.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map-mkt.$(draft).n10.sort.best.bed.path.quast.tsv \
	f1.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.label.gv.pdf \
	f1.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(ref).n10.ann.gv.squish.pdf

# Scaffold the ABySS assembly of the fly genome.
# ref=fly lr=f1 draft=f1.abyss n=25
f1-abyss: \
	f1.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(draft).n10.sort.best.bed.pdf \
	f1.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(draft).n10.sort.best.bed.path.$(ref).minidot.pdf \
	f1.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(draft).n10.sort.best.bed.path.quast.tsv \
	f1.k32-w32.n100-1000.c2-x.physlr.overlap.n$n.mol.backbone.map-mkt.$(draft).n10.sort.best.bed.path.quast.tsv

# Aggregate QUAST metrics.
f1.quast.tsv: \
		f1.abyss.quast.tsv \
		f1.abyss.scaftigs.quast.tsv \
		f1.supernova.quast.tsv \
		f1.supernova.scaftigs.quast.tsv \
		f1.k32-w32.n100-1000.c2-x.physlr.overlap.n25.mol.backbone.map.f1.abyss.n10.sort.best.bed.path.quast.tsv \
		f1.k32-w32.n100-1000.c2-x.physlr.overlap.n25.mol.backbone.map-mkt.f1.abyss.n10.sort.best.bed.path.quast.tsv \
		f1.k32-w32.n100-1000.c2-x.physlr.overlap.n25.mol.backbone.map.f1.supernova.scaftigs.n10.sort.best.bed.path.quast.tsv \
		f1.k32-w32.n100-1000.c2-x.physlr.overlap.n25.mol.backbone.map-mkt.f1.supernova.scaftigs.n10.sort.best.bed.path.quast.tsv
	mlr --tsvlite cut -x -f NG75,NGA75,LG75,LGA75 $^ >$@

# Download the fly genome from NCBI.
fly/fly.all.fa:
	mkdir -p $(@D)
	curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz | seqtk seq >$@

# Filter out the unplaced contigs from the fly genome.
fly/fly.fa: fly/fly.all.fa
	paste - - <$< | grep -v NW_ | tr '\t' '\n' >$@

# Download the fly annotation from NCBI.
fly/fly.gff:
	mkdir -p $(@D)
	curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz | $(zcat) >$@

# Download the Supernova assembly of the linked reads from 10x Genomics.
f1.supernova.fa:
	curl http://cf.10xgenomics.com/samples/assembly/2.1.0/fly/fly_pseudohap.fasta.gz | $(zcat) >$@

# Download the fly linked reads from 10x Genomics.
fly/f1.tar:
	mkdir -p $(@D)
	curl -o $@ http://s3-us-west-2.amazonaws.com/10x.files/samples/assembly/2.1.0/fly/fly_fastqs.tar

# Extract the tar file of fly FASTQ reads.
fly/f1.fq.gz: fly/f1.tar
	tar --wildcards -Oxf fly/f1.tar 'fly/H3C7LDMXX/read-RA*.fastq.gz' >$@

# Symlink the fly reads.
f1.fq.gz: fly/f1.bx.trimadap.fq.gz
	ln -s $< $@

# Download the fly linked reads from 10x Genomics.
# Extract the reads that map to chromosome 4.
%.chr4.sortbxn.bam: %.sortbxn.bam
	samtools view -h $< | awk '/^@/ || $$3 == "NC_004353.4"' | samtools view -@$t -o $@

# Symlink the chromosome 4 reads.
f1chr4.fq.gz: fly/fly.f1.chr4.sortbxn.dropse.fq.gz
	ln -sf $< $@

# Extract the reads that map to chromosome 2R.
%.chr2R.sortbxn.bam: %.sortbxn.bam
	samtools view -h $< | awk '/^@/ || $$3 == "NT_033778.4"' | samtools view -@$t -o $@

# Symlink the chromosome 2R reads.
f1chr2R.fq.gz: fly/fly.f1.chr2R.sortbxn.dropse.fq.gz
	ln -sf $< $@

# Symlink the subsampled chromosome 2R reads.
f1chr2R-bx%.fq.gz: fly/fly.f1.chr2R.sortbxn.dropse.bx%.fq.gz
	ln -sf $< $@

################################################################################
# Zebrafish
# See https://support.10xgenomics.com/de-novo-assembly/datasets/2.1.0/fish
# and https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/performance

# Assemble a physical map of fish.
# ref=z11 lr=fish draft=fish.supernova n=40 min_component_size=200
fish: \
	fish.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.k3.mol.backbone.path.z11.molecule.bed.pdf \
	fish.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.k3.mol.backbone.map.z11.n10.paf.gz.pdf \
	fish.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.k3.mol.backbone.map.fish.supernova.n10.sort.best.bed.pdf \
	fish.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.k3.mol.backbone.map.fish.supernova.n10.sort.best.bed.path.z11.minidot.pdf \
	fish.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.k3.mol.backbone.map.fish.supernova.n10.sort.best.bed.path.quast.tsv \
	fish.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.k3.mol.backbone.map-mkt.fish.supernova.n10.sort.best.bed.path.quast.tsv \
	fish.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.k3.mol.backbone.map.z11.n10.ann.gv.squish.pdf

# Aggregate QUAST metrics.
fish.quast.tsv: \
		fish.supernova.scaftigs.quast.tsv \
		fish.supernova.quast.tsv \
		fish.k32-w32.n100-2000.c2-x.physlr.overlap.n40.k3.mol.backbone.map.fish.supernova.n10.sort.best.bed.path.quast.tsv \
		fish.k32-w32.n100-2000.c2-x.physlr.overlap.n40.k3.mol.backbone.map-mkt.fish.supernova.n10.sort.best.bed.path.quast.tsv
	mlr --tsvlite cut -x -f NG75,NGA75,LG75,LGA75 $^ >$@

# Assemble a physical map of fish chromosome 25.
# ref=z11chr25 lr=fishchr25 draft=fish.supernova.chr25 n=50
fishchr25: \
	fishchr25.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.mol.backbone.label.gv.pdf \
	fishchr25.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.mol.backbone.path.z11chr25.molecule.bed.pdf \
	fishchr25.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.mol.backbone.map.z11chr25.n10.paf.gz.pdf \
	fishchr25.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.mol.backbone.map.fish.supernova.chr25.n10.sort.best.bed.pdf \
	fishchr25.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.mol.backbone.map.fish.supernova.chr25.n10.sort.best.bed.path.z11chr25.minidot.pdf \
	fishchr25.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.mol.backbone.map.fish.supernova.chr25.n10.sort.best.bed.path.quast.tsv \
	fishchr25.k32-w32.n100-2000.c2-x.physlr.overlap.n$n.mol.backbone.map.z11chr25.n10.ann.gv.squish.pdf

# Download the zebrafish genome from NCBI.
fish/z11.ncbi.fa:
	mkdir -p $(@D)
	curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.fna.gz | seqtk seq >$@

# Download the zebrafish genome from Ensembl.
fish/z11.all.fa:
	mkdir -p $(@D)
	curl ftp://ftp.ensembl.org/pub/release-95/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz | seqtk seq >$@

# Filter out the unplaced contigs from the reference genome.
fish/z11.fa: fish/z11.all.fa
	paste - - <$< | grep -v scaffold | gsort -V | tr '\t' '\n' >$@

# Extract chromosome 25 from the fish z11 reference genome.
fish/z11chr25.fa: fish/z11.fa
	samtools faidx -n999999999 $< 25 >$@

# Download the Supernova assembly of the linked reads from 10x Genomics.
fish.supernova.fa:
	curl http://cf.10xgenomics.com/samples/assembly/2.1.0/fish/fish_pseudohap.fasta.gz | seqtk seq >$@

# Extract chromosome 25 from the Supernova fish assembly.
fish.supernova.chr25.fa: fish.supernova.z11.chr25.sam.fa
	ln -s $< $@

# Download the fish linked reads from 10x Genomics.
fish/fish.tar:
	mkdir -p $(@D)
	curl -o $@ http://s3-us-west-2.amazonaws.com/10x.files/samples/assembly/2.1.0/fish/fish_fastqs.tar

# Extract the tar file of fish FASTQ reads.
fish/fish.fq.gz: fish/fish.tar
	tar --wildcards -Oxf fish/fish.tar 'fish/H3*/read-RA*.fastq.gz' >$@

# Symlink the fly reads.
fish.fq.gz: fish/fish.bx.trimadap.fq.gz
	ln -s $< $@

# Extract the reads that map to chromosome 25.
%.chr25.sortbxn.bam: %.sortbxn.bam
	samtools view -h $< | awk '/^@/ || $$3 == "25"' | samtools view -@$t -o $@

# Symlink the chromosome 25 reads.
fishchr25.fq.gz: fish/z11.fish.chr25.sortbxn.dropse.fq.gz
	ln -sf $< $@

################################################################################
# Maize (Zea mays)
# See https://support.10xgenomics.com/de-novo-assembly/datasets/2.1.0/maize
# and https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/performance

# Assemble a physical map of maize.
# ref=b73 lr=maize draft=maize.supernova n=50 min_component_size=200
maize: \
	maize.k32-w32.n100-5000.physlr.mxperbx.tsv.pdf \
	maize.k32-w32.n100-5000.physlr.depth.tsv.pdf \
	maize.k32-w32.n100-5000.c2-x.physlr.overlap.n$n.k3.mol.backbone.map.b73.n10.paf.gz.pdf \
	maize.k32-w32.n100-5000.c2-x.physlr.overlap.n$n.k3.mol.backbone.map.maize.supernova.n10.sort.best.bed.pdf \
	maize.k32-w32.n100-5000.c2-x.physlr.overlap.n$n.k3.mol.backbone.map.maize.supernova.n10.sort.best.bed.path.b73.minidot.pdf \
	maize.k32-w32.n100-5000.c2-x.physlr.overlap.n$n.k3.mol.backbone.map.maize.supernova.n10.sort.best.bed.path.quast.tsv \
	maize.k32-w32.n100-5000.c2-x.physlr.overlap.n$n.k3.mol.backbone.map-mkt.maize.supernova.n10.sort.best.bed.path.quast.tsv \
	maize.k32-w32.n100-5000.c2-x.physlr.overlap.n$n.k3.mol.backbone.map.b73.n10.ann.gv.squish.pdf

# Aggregate QUAST metrics.
maize.quast.tsv: \
		maize.supernova.scaftigs.quast.tsv \
		maize.supernova.quast.tsv \
		maize.k32-w32.n200-10000.c2-x.physlr.overlap.n50.k3.mol.backbone.map.maize.supernova.n10.sort.best.bed.path.quast.tsv \
		maize.k32-w32.n200-10000.c2-x.physlr.overlap.n50.k3.mol.backbone.map-mkt.maize.supernova.n10.sort.best.bed.path.quast.tsv
	mlr --tsvlite cut -x -f NG75,NGA75,LG75,LGA75 $^ >$@

# Download the maize genome from NCBI.
maize/b73.ncbi.fa:
	mkdir -p $(@D)
	curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/005/GCF_000005005.2_B73_RefGen_v4/ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/005/GCF_000005005.2_B73_RefGen_v4/GCF_000005005.2_B73_RefGen_v4_genomic.fna.gz | seqtk seq >$@

# Download the maize genome from Ensembl.
maize/b73.all.fa:
	mkdir -p $(@D)
	curl ftp://ftp.ensemblgenomes.org/pub/plants/release-43/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz | seqtk seq >$@

# Filter out the unplaced contigs from the reference genome.
maize/b73.fa: maize/b73.all.fa
	paste - - <$< | grep -v contig | gsort -V | tr '\t' '\n' >$@

# Download the Supernova assembly of the linked reads from 10x Genomics.
maize.supernova.fa:
	curl http://cf.10xgenomics.com/samples/assembly/2.1.0/maize/maize_pseudohap.fasta.gz | seqtk seq >$@

# Download the maize linked reads from 10x Genomics.
maize/maize.tar:
	mkdir -p $(@D)
	curl -o $@ http://s3-us-west-2.amazonaws.com/10x.files/samples/assembly/2.1.0/maize/maize_fastqs.tar

# Extract the tar file of maize FASTQ reads.
maize/maize.fq.gz: maize/maize.tar
	tar --wildcards -Oxf maize/maize.tar 'maize/H3*/read-RA*.fastq.gz' >$@

# Symlink the maize reads.
maize.fq.gz: maize/maize.bx.trimadap.fq.gz
	ln -s $< $@

################################################################################
# Human (HG004)

# Assemble a physical map of hg004.
# ref=grch38 lr=hg004 draft=hg004.supernova.scaftigs k=40 w=32 n=62 min_component_size=200
hg004: \
	hg004.k40-w32.n100-5000.physlr.mxperbx.tsv.pdf \
	hg004.k40-w32.n100-5000.physlr.depth.tsv.pdf \
	hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map.grch38.n10.paf.gz.pdf \
	hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map.hg004.supernova.scaftigs.n10.sort.best.bed.path.quast.tsv \
	hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map.grch38.n10.ann.gv.squish.pdf \
	hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map-mkt.hg004.supernova.scaftigs.n10.sort.best.bed.pdf \
	hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map-mkt.hg004.supernova.scaftigs.n10.sort.best.bed.path.grch38.minidot.pdf \
	hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map-mkt.hg004.supernova.scaftigs.n10.sort.best.bed.path.quast.tsv

# Scaffold the ABySS assembly with the physical map of hg004.
# ref=grch38 lr=hg004 draft=hg004.abyss k=40 w=32 n=62 min_component_size=200
hg004-abyss: \
	hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map.hg004.abyss.n10.sort.best.bed.pdf \
	hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map.hg004.abyss.n10.sort.best.bed.path.grch38.minidot.pdf \
	hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map.hg004.abyss.n10.sort.best.bed.path.quast.tsv \
	hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map-mkt.hg004.abyss.n10.sort.best.bed.path.quast.tsv

# Download the human genome from NCBI.
grch38/grch38.all.fa:
	mkdir -p $(@D)
	curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
		| seqtk seq >$@

# Filter out the unplaced contigs from the reference genome.
grch38/grch38.fa: grch38/grch38.all.fa
	paste - - <$< | egrep -v 'chrEBV|chrUn|random' | tr '\t' '\n' >$@

# Download coordinates of modeled centromeres and heterochromatin regions.
grch38/grch38.gap.tsv:
	mkdir -p $(@D)
	curl https://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/38/Modeled_regions_for_GRCh38.tsv \
		| sed 's/^#//' >$@

# Convert coordinates of centromeres and gaps to BED format.
grch38/grch38.gap.bed: grch38/grch38.gap.tsv grch38/grch38.fa.gap.bed
	(awk 'NR>1 { print "chr" $$2 "\t" $$3 "\t" $$4 "\t" $$1 }' $<; \
		awk '{ print $$1 "\t" $$2 "\t" $$3 "\tgap" NR }' grch38/grch38.fa.gap.bed) >$@

# Symlink the reads.
hg004.fq.gz: hg004.lrbasic.trimadap.fq.gz
	ln -s $< $@

# Aggregate QUAST metrics.
hg004.quast.tsv: \
		hg004.abyss.quast.tsv \
		hg004.supernova.quast.tsv \
		hg004.supernova.scaftigs.quast.tsv \
		hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map.hg004.abyss.n10.sort.best.bed.path.quast.tsv \
		hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map-mkt.hg004.abyss.n10.sort.best.bed.path.quast.tsv \
		hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map.hg004.supernova.scaftigs.n10.sort.best.bed.path.quast.tsv \
		hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map-mkt.hg004.supernova.scaftigs.n10.sort.best.bed.path.quast.tsv \
		hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n62.k3.mol.backbone.map.hg004.supernova.n10.sort.best.bed.path.quast.tsv
	mlr --tsvlite cut -x -f NG75,NGA75,LG75,LGA75 $^ >$@

################################################################################
# Western redcedar (Thuja plicata)
# See https://github.com/sjackman/redcedar
# and https://github.com/sjackman/thesis/blob/master/redcedar.md

# Assemble a physical map of western redcedar.
# lr=redcedar draft=redcedar.abyss n=10 min_component_size=200
redcedar: \
	redcedar.k32-w32.n100-5000.physlr.mxperbx.tsv.pdf \
	redcedar.k32-w32.n100-5000.physlr.depth.tsv.pdf \
	redcedar.k32-w32.n100-5000.c2-x.physlr.overlap.n$n.mol.backbone.path \
	redcedar.k32-w32.n100-5000.c2-x.physlr.overlap.n$n.mol.backbone.map.redcedar.abyss.n10.sort.best.bed.path.fa \
	redcedar.k32-w32.n100-5000.c2-x.physlr.overlap.n$n.mol.backbone.map.redcedar.abyss.n10.sort.best.bed.pdf

# Symlink the linked reads.
redcedar.fq.gz: redcedar.lrbasic.trimadap.fq.gz
	ln -s $< $@

################################################################################
# White spruce (Picea glauca genotype WS77111)

# Assemble a physical map of white spruce.
# lr=ws77111 draft=ws77111.abyss n=50 min_component_size=200
ws77111: \
	ws77111.k32-w32.n100-5000.physlr.mxperbx.tsv.pdf \
	ws77111.k32-w32.n100-5000.physlr.depth.tsv.pdf \
	ws77111.k32-w32.n100-5000.c2-x.physlr.overlap.n$n.mol.backbone.path \
	ws77111.k32-w32.n100-5000.c2-x.physlr.overlap.n$n.mol.backbone.map.ws77111.abyss.n10.sort.best.bed.path.fa \
	ws77111.k32-w32.n100-5000.c2-x.physlr.overlap.n$n.mol.backbone.map.ws77111.abyss.n10.sort.best.bed.pdf

# Trim adapter sequences using trimadap.
ws77111.lrbasic.trimadap.fq.gz: path=/projects/spruceup/pglauca/WS77111/data/reads/DNA/fastq/chromium
ws77111.lrbasic.trimadap.fq.gz:
	( \
		seqtk mergepe $(path)/A64454_first.fastq.gz $(path)/A64454_second.fastq.gz; \
		seqtk mergepe $(path)/A64454-2_first.fastq.gz $(path)/A64454-2_second.fastq.gz; \
		seqtk mergepe $(path)/A64455_first.fastq.gz $(path)/A64455_second.fastq.gz; \
		seqtk mergepe $(path)/A64455-2_first.fastq.gz $(path)/A64455-2_second.fastq.gz; \
		seqtk mergepe $(path)/A64456_first.fastq.gz $(path)/A64456_second.fastq.gz; \
		seqtk mergepe $(path)/A64456-2_first.fastq.gz $(path)/A64456-2_second.fastq.gz; \
		seqtk mergepe $(path)/A64457_first.fastq.gz $(path)/A64457_second.fastq.gz; \
		seqtk mergepe $(path)/A64457-2_first.fastq.gz $(path)/A64457-2_second.fastq.gz \
	) | trimadap-mt -p$t -t1 - | sed 's/^X$$/N/' | $(gzip) >$@

# Symlink the linked reads.
ws77111.fq.gz: ws77111.lrbasic.trimadap.fq.gz
	ln -s $< $@

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

# Align molecule-split linked reads to a target genome.
%.$(lr).n100-2000.overlap.n$n.split.sortn.bam: %.fa.bwt $(lr).n100-2000.physlr.overlap.n$n.split.fq.gz
	bwa mem -t$t -pC $*.fa $(lr).n100-2000.physlr.overlap.n$n.split.fq.gz | samtools view -@$t -F4 -o $@

################################################################################
# EMA

# Map linked reads to the draft genome using EMA.
# Filter out reads without barcodes.
%.$(lr).ema.sortn.bam: $(lr).fq.gz %.fa.bwt
	$(zcat) $< | paste - - - - - - - - | grep "BX:Z:" | tr '\t' '\n' \
	| $(time) ema align -t$t -r $*.fa -1 /dev/stdin | samtools view -@$t -h -F4 -o $@

################################################################################
# minimap2

# Align linked reads to a target genome.
%.$(lr).minimap2.sortn.bam: %.fa.bwt $(lr).fq.gz
	minimap2 -t$t -a -xsr -y $*.fa $(lr).fq.gz | samtools view -@$t -F4 -o $@

# Align a FASTA file to the reference genome and produce a PAF file.
%.$(ref).paf.gz: $(name)/$(ref).fa %.fa
	$(time) minimap2 -t$t -xasm10 $^ | $(gzip) >$@

# Align a FASTA file to the reference genome and produce a SAM file.
%.$(ref).sam.gz: $(name)/$(ref).fa %.fa
	$(time) minimap2 -t$t -xasm10 -a $^ | $(gzip) >$@

################################################################################
# Filter alignments

# Extract the sequences that map to chromosome 25.
%.chr25.paf.gz: %.paf.gz
	$(zcat) $< | awk '$$6 == "25"' | $(gzip) >$@

# Extract the sequences that map to chromosome 25.
%.chr25.sam.gz: %.sort.bam %.sort.bam.bai
	samtools view -h $< 25 | $(gzip) >$@

################################################################################
# miniasm

# Draw a dot plot of a PAF file.
# Skip alignments to non-chromosomal sequences.
%.minidot.pdf: %.paf.gz
	$(zcat) $< | grep -v NW_ | minidot /dev/stdin | ps2pdf - $@

################################################################################
# samtools

# Index a FASTA file.
%.fa.fai: %.fa
	samtools faidx $<

# Sort a SAM file by position.
%.sort.bam: %.sam.gz
	samtools sort -@$t -T$$(mktemp -u -t $(@F).XXXXXX) -o $@ $<

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

# Convert a SAM file to FASTA.
%.sam.fa: %.sam.gz
	samtools fasta -@$t -F0x904 $< >$@

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
	$(zcat) $< | ema preproc -t$t -b -n1 -w 4M-with-alts-february-2016.txt -o $*.ema $*.ema-ncnt
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
dist=50000
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
	ntcard -t$t -c1000 -k 32,48,64,80,96,112,128 -p $*.ntcard $<

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

# Extract the coordinates in BED format of gaps in a FASTA file.
%.fa.gap.bed: %.fa
	$(python) $(bin)/physlr fasta-gaps $< >$@

# Index a FASTA file.
%.k$k-w$w.physlr.tsv: %.fa
	$(time) $(physlr_path)/src/physlr-indexlr -t$t -k$k -w$w --pos -o $@ $<

# Index a compressed FASTQ file of linked reads.
%.k$k-w$w.physlr.tsv.gz: %.fq.gz
	$(zcat) $< | $(time) $(physlr_path)/src/physlr-indexlr -t$t -k$k -w$w - | $(gzip) >$@

# Identify the overlapping minimizers of each pair of barcodes.
%.physlr.intersect.tsv: %.physlr.tsv
	$(python) $(bin)/physlr intersect $< >$@

# Filter barcodes by number of minimizers.
%.n100-1000.physlr.tsv.gz: %.physlr.tsv.gz
	$(zcat) $< | $(time) $(physlr_path)/src/physlr-filter-barcodes -n100 -N1000 - | $(gzip) >$@

# Filter barcodes by number of minimizers.
%.n100-2000.physlr.tsv.gz: %.physlr.tsv.gz
	$(zcat) $< | $(time) $(physlr_path)/src/physlr-filter-barcodes -n100 -N2000 - | $(gzip) >$@

# Filter barcodes by number of minimizers.
%.n100-5000.physlr.tsv.gz: %.physlr.tsv.gz
	$(zcat) $< | $(time) $(physlr_path)/src/physlr-filter-barcodes -n100 -N5000 - | $(gzip) >$@

# Filter barcodes by number of minimizers.
%.n100-10000.physlr.tsv.gz: %.physlr.tsv.gz
	$(zcat) $< | $(time) $(physlr_path)/src/physlr-filter-barcodes -n100 -N10000 - | $(gzip) >$@

# Count the number of minimizers per barcode.
%.physlr.mxperbx.tsv: %.physlr.tsv.gz
	$(zcat) $< | awk 'BEGIN { print "Minimizers" } { print NF - 1 }' >$@

# Count the frequency of the minimizers.
%.physlr.depth.tsv: %.physlr.tsv.gz
	$(zcat) $< | $(time) $(python) $(bin)/physlr count-minimizers -V$V - >$@

# Filter minimizers by their frequency.
%.c2-x.physlr.tsv: %.physlr.tsv.gz
	$(zcat) $< | $(time) $(physlr_path)/src/physlr-filter-bxmx -o $@ -

# Determine overlaps and output the graph in TSV.
%.physlr.overlap.tsv: %.physlr.tsv
	$(time) $(physlr_path)/src/physlr-overlap -t1 -n10 $< >$@

# Determine the degree of each vertex.
%.deg.tsv: %.tsv
	$(python) $(bin)/physlr degree $< >$@

# Determine the maximum spanning tree.
%.mst.tsv: %.tsv
	$(python) $(bin)/physlr mst -V$V --prune-branches=10 --prune-bridges=10 $< >$@

# Separate a graph into its biconnected components by removing its cut vertices.
%.bic.tsv: %.tsv
	$(python) $(bin)/physlr biconnected-components $< >$@

# Determine the backbone graph from the overlap TSV.
%.backbone.tsv: %.tsv
	$(time) $(python) $(bin)/physlr backbone-graph --prune-branches=10 --prune-bridges=10 -s0 $< >$@

# Determine the corrected backbone graph from the overlap TSV.
min_path_size=200
%.s3.backbone.tsv: %.tsv
	$(time) $(python) $(bin)/physlr backbone-graph --prune-branches=10 --prune-bridges=10 -s3 --min-path-size=$(min_path_size) $< >$@

# Determine the backbone path of the backbone graph.
%.path: %.tsv
	$(time) $(python) $(bin)/physlr backbone --prune-branches=0 -s0 $< >$@

# Cut chimeric paths to correct misassemblies.
%.backbone.depth.tsv: %.tsv %.backbone.path
	$(time) $(python) $(bin)/physlr cut-chimera -V$V -s3 --min-path-size=0 $^ -o $*.backbone.depth.tsv >$*.backbone.s3.path

# Cut chimeric paths to correct misassemblies.
%.s3.backbone.depth.tsv: %.tsv %.s3.backbone.path
	$(time) $(python) $(bin)/physlr cut-chimera -V$V -s3 --min-path-size=0 $^ -o $*.s3.backbone.depth.tsv >$*.s3.backbone.s3.path

# Flesh out the backbone path
%.backbone.fleshed.path: %.tsv %.backbone.path
	$(python) $(bin)/physlr flesh-backbone --min-component-size=$(min_component_size) $< $*.backbone.path >$@

# Split the minimizers to molecules
%.overlap.n$n.mol.split.tsv: %.overlap.n$n.mol.tsv %.tsv
	$(python) $(bin)/physlr split-minimizers -t$t $< $*.tsv >$@

# Split the reads into molecules
%.split.fq.gz: $(lr).fq.gz %.split.tsv $(lr).k$k-w$w.physlr.tsv
	$(zcat) $< | $(python) $(bin)/physlr split-reads-molecules $*.split.tsv \
		$(lr).k$k-w$w.physlr.tsv - | $(gzip) >$@

# Determine the minimum tiling graph of the backbone graph.
%.backbone.tiling.tsv: %.backbone.tsv
	$(python) $(bin)/physlr tiling-graph $< >$@

# Map the barcodes to the backbone graph and output BED.
%.backbone.map.$(lr).n100-2000.n10.bed: %.backbone.path $(lr).n100-2000.physlr.tsv
	$(time) $(python) $(bin)/physlr map -V$V -n10 $^ $(lr).n100-2000.physlr.tsv >$@

# Map the draft assembly to the backbone graph and output BED.
%.overlap.n$n.mol.backbone.map.$(draft).n10.bed: %.overlap.n$n.mol.backbone.path %.tsv $(draft).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map -V$V -n10 $^ >$@

# Map the draft assembly to the backbone graph and output BED (K3-cliques).
%.overlap.n$n.k3.mol.backbone.map.$(draft).n10.bed: %.overlap.n$n.k3.mol.backbone.path %.tsv $(draft).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map -V$V -n10 $^ >$@

# Map the draft assembly to the Louvain backbone graph and output BED.
%.overlap.n$n.louvain.mol.backbone.map.$(draft).n10.bed: %.overlap.n$n.louvain.mol.backbone.path %.tsv $(draft).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map -V$V -n10 $^ >$@

# Map the reference assembly to the backbone graph and output BED.
%.overlap.n$n.mol.backbone.map.$(ref).n10.bed: %.overlap.n$n.mol.backbone.path %.tsv $(name)/$(ref).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map -V$V -n10 $^ >$@

# Map the reference assembly to the backbone graph and output BED (K3-cliques).
%.overlap.n$n.k3.mol.backbone.map.$(ref).n10.bed: %.overlap.n$n.k3.mol.backbone.path %.tsv $(name)/$(ref).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map -V$V -n10 $^ >$@

# Map the reference assembly to the Louvain backbone graph and output BED.
%.overlap.n$n.louvain.mol.backbone.map.$(ref).n10.bed: %.overlap.n$n.louvain.mol.backbone.path %.tsv $(name)/$(ref).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map -V$V -n10 $^ >$@

# Map the draft assembly to the backbone graph using physlr map-mkt and output BED.
%.overlap.n$n.mol.backbone.map-mkt.$(draft).n10.bed: %.overlap.n$n.mol.backbone.path %.tsv $(draft).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map-mkt -V$V -n10 $^ >$@

# Map the draft assembly to the k3-cliques backbone graph using physlr map-mkt and output BED.
%.overlap.n$n.k3.mol.backbone.map-mkt.$(draft).n10.bed: %.overlap.n$n.k3.mol.backbone.path %.tsv $(draft).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map-mkt -V$V -n10 $^ >$@

# Map the draft assembly to the Louvain backbone graph using physlr map-mkt and output BED.
%.overlap.n$n.louvain.mol.backbone.map-mkt.$(draft).n10.bed: %.overlap.n$n.louvain.mol.backbone.path %.tsv $(draft).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map-mkt -V$V -n10 $^ >$@

# Map the reference to the backbone graph and output PAF.
%.overlap.n$n.mol.backbone.map.$(ref).n10.paf.gz: %.overlap.n$n.mol.backbone.path %.tsv $(name)/$(ref).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map-paf -V$V -n10 $^ | $(gzip) >$@

# Map the reference to the k3-cliques backbone graph and output PAF.
%.overlap.n$n.k3.mol.backbone.map.$(ref).n10.paf.gz: %.overlap.n$n.k3.mol.backbone.path %.tsv $(name)/$(ref).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map-paf -V$V -n10 $^ | $(gzip) >$@

# Map the reference to the Louvain backbone graph and output PAF.
%.overlap.n$n.louvain.mol.backbone.map.$(ref).n10.paf.gz: %.overlap.n$n.louvain.mol.backbone.path %.tsv $(name)/$(ref).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map-paf -V$V -n10 $^ | $(gzip) >$@

# Map the reference to the corrected backbone graph and output PAF.
%.overlap.n$n.mol.s3.backbone.map.$(ref).n10.paf.gz: %.overlap.n$n.mol.s3.backbone.path %.tsv $(name)/$(ref).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map-paf -V$V -n10 $^ | $(gzip) >$@

# Lift over query coordinates of a PAF file from minimzer index to nucleotide coordinate.
%.qpos.paf.gz: $(name)/$(ref).k$k-w$w.physlr.tsv %.paf.gz
	$(zcat) $*.paf.gz | $(time) $(python) $(bin)/physlr liftover-paf -V$V $< - | $(gzip) >$@

# Chain alignments in PAF format.
%.chain.paf.gz: %.paf.gz
	bin/chainpaf $< $@

# Reference genome size.
g=$(shell awk '{x += $$2} END{print x}' $(name)/$(ref).fa.fai)

# Report assembly metrics of a PAF file.
%.qpos.chain.metrics.tsv: %.qpos.chain.paf.gz $(name)/$(ref).fa.fai
	$(zcat) $< | $(time) $(python) $(bin)/physlr paf-metrics -V$V -g$g - \
		| sed 's/\/dev\/stdin/$</' >$@

# Annotate a backbone graph with a BED file of mappings.
# Filter out small components.
%.map.$(ref).n10.sort.best.ann.gv: %.tsv %.path %.map.$(ref).n10.sort.best.bed
	$(python) $(bin)/physlr annotate-graph --min-component-size=$(min_component_size) -Ogv $^ >$@

# Annotate a backbone graph with a PAF file of mappings.
# Filter out small components.
%.backbone.map.$(ref).n10.ann.gv: %.backbone.map.$(ref).n10.paf.gz %.backbone.tsv %.backbone.path
	gunzip -c $< | $(python) $(bin)/physlr annotate-graph -V$V --min-component-size=$(min_component_size) -Ogv $*.backbone.tsv $*.backbone.path - >$@

# Annotate a MST graph with a PAF file of mappings.
# Filter out small components.
%.mst.map.$(ref).n10.ann.gv: %.backbone.map.$(ref).n10.paf.gz %.mst.tsv %.backbone.path
	gunzip -c $< | $(python) $(bin)/physlr annotate-graph -V$V --min-component-size=$(min_component_size) -Ogv $*.mst.tsv $*.backbone.path - >$@

# Filter a BED file by score.
%.n100.bed: %.n10.bed
	awk '$$5 >= 100' $< >$@

# Estimate the number of molecules per barcode.
%.physlr.overlap.n20.countmol.tsv: %.physlr.overlap.tsv
	$(python) $(bin)/physlr count-molecules -n20 $< >$@

# Remove barcodes with more than one molecule.
%.physlr.overlap.molecules.M2.tsv: %.physlr.overlap.molecules.tsv
	$(python) $(bin)/physlr filter -M2 $< >$@

# Count the number of common neighbours for each edge.
%.cn.tsv: %.tsv
	$(python) $(bin)/physlr common-neighbours $< >$@

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

# Filter edges by n using Miller.
%.n$n.tsv: %.tsv
	mlr --tsvlite filter '$$n >= $n' $< >$@

# Keep the best edges of each vertex.
%.best$(bestn).tsv: %.tsv
	$(python) $(bin)/physlr best-edges --bestn=$(bestn) $< >$@

# Separate barcodes into molecules.
mol_strategy=bc
%.mol.tsv: %.tsv
	$(time) $(python) $(bin)/physlr molecules -V$V -t$t --separation-strategy=$(mol_strategy) $< >$@

# Separate barcodes into molecules using K3-clique community detection.
%.k3.mol.tsv: %.tsv
	$(time) $(python) $(bin)/physlr molecules -V$V -t$t --separation-strategy=bc+k3 $< >$@

# Separate barcodes into molecules using Louvain community detection.
%.louvain.mol.tsv: %.tsv
	$(time) $(python) $(bin)/physlr molecules -V$V -t$t --separation-strategy=bc+louvain $< >$@

# Make a vertex-induced subgraph
%.subgraph.tsv: %.tsv
	$(python) $(bin)/physlr subgraph -v$v -d$d $< >$@

# Make multiple vertex-induced subgraphs
%.subgraphs: %.tsv
	$(python) $(bin)/physlr subgraphs -v$v -d$d -o $*.subgraphs $<

# Extract statistics of the subgraphs
%.subgraphs-stats.tsv: %.tsv
	$(python) $(bin)/physlr subgraphs-stats -t$t $< >$@

# Convert a graph from TSV to GraphViz.
# Filter out small components.
%.gv: %.tsv
	$(python) $(bin)/physlr filter --min-component-size=$(min_component_size) -Ogv $< >$@

# Extract a BED file of the backbone barcodes.
# Filter out small components.
%.path.$(ref).molecule.bed: $(name)/$(ref).$(lr).a0.65.d$(dist).n5.q1.s2000.molecule.bed %.path
	$(python) $(bin)/physlr filter-bed --min-component-size=$(min_component_size) $^ >$@

# Extract a BED file of the fleshed-out backbone barcodes.
# Filter out small components.
%.backbone.fleshed.path.$(ref).molecule.bed: $(name)/$(ref).$(lr).a0.65.d$(dist).n5.q1.s2000.molecule.bed %.backbone.fleshed.path
	$(python) $(bin)/physlr filter-bed --min-component-size=$(min_component_size) $^ >$@

# Extract a BED file of fleshed-out backbone barcodes from molecule-split reads
# Filter out small components
%.backbone.fleshed.path.$(ref).molecule.split.bed:  $(name)/$(ref).$(lr).n100-2000.overlap.n$n.split.a0.65.d$(dist).n5.q1.s2000.molecule.bed %.backbone.fleshed.path
	$(python) $(bin)/physlr filter-bed --min-component-size=$(min_component_size) --molecule-bed=1 $^ >$@

# Sort a BED file.
%.sort.bed: %.bed
	sort -k1,1n -k1,1 -k2,2n -k3,3n -k5,5nr -k4,4 $< >$@

# Keep the best record at each position.
%.sort.best.bed: %.sort.bed
	awk '{ keep = $$1 " " $$2 " " $$3 != x; x = $$1 " " $$2 " " $$3 } keep' $< >$@

# Extract scaffolds paths from a BED file.
%.bed.path: %.bed
	$(time) $(python) $(bin)/physlr bed-to-path $^ >$@

# Sort a PAF file by target position.
%.sort.paf.gz: %.paf.gz
	gunzip -c $< | sort -k6,6n -k8,8n -k9,9n -k10,10nr | $(gzip) >$@

# Keep the best alignment for each PAF (Tname,Tstart,Tend) tuple.
%.sort.best.paf.gz: %.sort.paf.gz
	gunzip -c $< | awk '{ keep = $$6 " " $$8 " " $$9 != x; x = $$6 " " $$8 " " $$9 } keep' | $(gzip) >$@

# Produce sequences in FASTA format from paths.
%.$(draft).n10.sort.best.bed.path.fa: $(draft).fa %.$(draft).n10.sort.best.bed.path
	$(time) $(python) $(bin)/physlr path-to-fasta --min-length=3000 $^ >$@

# Convert a PAF file to a BED file.
%.paf.bed: %.paf.gz
	gunzip -c $< | awk -F'\t' -vOFS='\t' '{ print $$1, $$3, $$4, $$6 "_" $$8, $$10, $$5}' >$@

################################################################################
# Bedtools

# Compute genome coverage.
%.bed.$(ref).cov.tsv: %.bed $(name)/$(ref).fa.fai
	grep -v NA $< | sort -k1,1 -k2,2n -k3,3n | bedtools genomecov -max 1 -g $(name)/$(ref).fa.fai -i - | awk '$$2 != 0 || $$5 != 1' >$@

################################################################################
# QUAST

# Calculate assembly contiguity and correctness metrics using QUAST.
%.quast.tsv: %.fa $(name)/$(ref).fa
	quast-lg -t$t -es --fast --large --scaffold-gap-max-size 100000 --min-identity 95 -R $(name)/$(ref).fa -o $*.quast $<
	cp $*.quast/transposed_report.tsv $@

# Aggregate QUAST metrics.
%.quast.tsv: \
		%.supernova.scaftigs.quast.tsv \
		%.supernova.quast.tsv \
		%.n100-2000.c2-x.physlr.overlap.n$n.mol.backbone.map.$(name).supernova.scaftigs.n10.sort.best.bed.path.quast.tsv
	mlr --tsvlite cut -x -f NG75,NGA75,LG75,LGA75 $^ >$@

################################################################################
# GraphViz

# Extract the largest component.
%.comp1.gv: %.gv
	ccomps -z -X'#0' $< >$@ || true

# Extract components of at least 100 vertices.
%.X100.gv: %.gv
	(ccomps -z -X'%100-' $< || true) \
	| awk 'BEGIN { print "strict graph {" } END { print "}" } { sub("^strict graph ", "subgraph ") } 1' >$@

# Label the edges with edge weight.
%.label.gv: %.gv
	gvpr -c 'E { label = n }' $< >$@

# Remove text labels.
%.nolabel.gv: %.gv
	gvpr -c 'BEG_G { setDflt ($$, "N", "shape", "point") } E { label = "" }' $< >$@

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
# Do not prevent overlapping vertices.
%.gv.squish.pdf: %.gv
	neato -Gsize=100 -Tpdf -o $@ $<

# Layout and render an undirected graph to PDF.
# Prevent overlapping vertices by scaling the layout.
%.gv.pdf: %.gv
	neato -Gsize=100 -Goverlap=scale -Tpdf -o $@ $<

# Layout and render an undirected graph to PDF using sfdp.
%.gv.sfdp.pdf: %.gv
	sfdp -Gsize=100 -Goverlap_scaling=200 -Tpdf -o $@ $<

################################################################################
# RMarkdown reports

# Plot number of markers per barcode.
%.mxperbx.tsv.pdf: %.mxperbx.tsv
	Rscript -e 'rmarkdown::render("$(physlr_path)/data/plot-mxperbx.rmd", "html_document", "$*.mxperbx.tsv.html", knit_root_dir="$(PWD)", output_dir="$(PWD)", params=list(input_tsv="$<"))'

# Plot depth of coverage.
%.depth.tsv.pdf: %.depth.tsv
	Rscript -e 'rmarkdown::render("$(physlr_path)/data/plot-depth.rmd", "html_document", "$*.depth.tsv.html", knit_root_dir="$(PWD)", output_dir="$(PWD)", params=list(input_tsv="$<"))'

# Plot edge property n.
%.tsv.n.pdf: %.tsv
	Rscript -e 'rmarkdown::render("$(physlr_path)/data/plot-edge-n.rmd", "html_document", "$*.tsv.n.html", knit_root_dir="$(PWD)", output_dir="$(PWD)", params=list(input_tsv="$<"))'

# Plot vertex degree.
%.deg.tsv.pdf: %.deg.tsv
	Rscript -e 'rmarkdown::render("$(physlr_path)/data/plot-degree.rmd", "html_document", "$*.deg.tsv.html", knit_root_dir="$(PWD)", output_dir="$(PWD)", params=list(input_tsv="$<"))'

# Plot a BED file.
%.bed.pdf: %.bed
	Rscript -e 'rmarkdown::render("$(physlr_path)/data/plotbed.rmd", "html_document", "$*.bed.html", knit_root_dir="$(PWD)", output_dir="$(PWD)", params=list(input_bed="$<"))'

# Plot a PAF file.
%.paf.gz.pdf: %.paf.gz $(name)/$(ref).fa.gap.bed
	Rscript -e 'rmarkdown::render("$(physlr_path)/data/plotpaf.rmd", "html_document", "$*.paf.gz.html", knit_root_dir="$(PWD)", output_dir="$(PWD)", params=list(input_paf="$<", gap_bed="$(name)/$(ref).fa.gap.bed"))'

# Compare assembly metrics.
%.quast.html: %.quast.tsv
	Rscript -e 'rmarkdown::render("$(physlr_path)/data/quast.rmd", "html_document", "$*.quast.html", knit_root_dir="$(PWD)", output_dir="$(PWD)", params=list(input_tsv="$<", output_tsv="$*.quast.table.tsv"))'
