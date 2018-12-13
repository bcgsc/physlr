    Compile
        Release: make
        Debug: make DEBUG=1

    Run
        ./minimizereads -k 100 -w 5 -i data/tiny.fq >outfile

    Usage:  ./minimizereads -w <window size> -k <kmer length> [-i file] [-v]
    
    Output: Each line of output is a barcode followed by a list of minimizers.
            <barcode>\t<minimizer1> <minimizer2> ... <minimizerN>
