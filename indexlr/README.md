    Compile
        Release: make

    Usage
	    ./physlr-indexlr -k K -w W [-v] file...

    Run
        ./physlr-indexlr -k100 -w5 data/tiny.fq >out

    Output
	    Each line of output is a barcode followed by a tab, and then a list of space-separated minimizers.
        <barcode>\t<minimizer1> <minimizer2> ... <minimizerN>
