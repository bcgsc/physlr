.DELETE_ON_ERROR:
.SECONDARY:

all: lint

# Check the source code for errors with Pylint.
lint:
	pylint physlr

# Assemble the test data.
check:
	make test m=1 prune_branches=1 prune_bridges=1 -C data
	diff -q data/test.k1-w1.n1-2.c2-x.physlr.overlap.m0.mol.backbone0.path data/test.k1-w1.n1-2.c2-x.physlr.overlap.m0.mol.backbone0.path.good

# Render the diagram of the pipeline.
pipeline.pdf: pipeline.gv
	dot -Tpdf -o $@ $<
