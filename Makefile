.DELETE_ON_ERROR:
.SECONDARY:

all: lint

# Check the source code for errors with Pylint.
lint:
	pylint physlr

# Assemble the test data.
check:
	make test m=1 prune_branches=1 prune_bridges=1 -C data
	diff -q data/test.k1-w1.n1-2.c2-x.physlr.overlap.m1.mol.backbone.path data/test.k1-w1.n1-2.c2-x.physlr.overlap.m1.mol.backbone.path.good

# Render the diagram of the pipeline.
pipeline.pdf: pipeline.gv
	dot -Tpdf -o $@ $<
