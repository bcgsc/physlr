    Compile
        Release: make
        Debug: make DEBUG=1

    Usage
	    minimizereads -k K -w W [-v] file...

    Run
        ./minimizereads -k 100 -w 5 data/tiny.fq >out

    Output
	    Each line of output is a barcode followed by a list of minimizers.
            <barcode>\t<minimizer1> <minimizer2> ... <minimizerN>
