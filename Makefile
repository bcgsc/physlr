.DELETE_ON_ERROR:
.SECONDARY:

all: lint

# Check the source code for errors with Pylint.
lint:
	pylint physlr

# Assemble the test data.
check:
	make -C data

# Render the diagram of the pipeline.
pipeline.pdf: pipeline.gv
	dot -Tpdf -o $@ $<
