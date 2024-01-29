####################################################################################################
# Reads a chunks of fasta sequences or fastq reads (-b reads at a time for fastq) and performs RLE (without keeping track of counts), in other words homopolymer compression.
# example command:
# pypy3 hpc.py -i input.fq -o input.rle.fq -b 100000
# (automnatcially detects if input is fasta or fastq)
# Author: Amirhossein Afshinfard aafshinfard@gmail.com

import argparse
import sys
import datetime
import re

def rle_fasta(seq):
    """
    HPC encoding for a sequence
    """
    compressed_seq = []
    prev_base = None
    for base in seq:
        if base != prev_base:
            compressed_seq.append(base)
            prev_base = base
    return ''.join(compressed_seq)

### Faster version of RLE
def rle(seq, qual):
    """
    HPC encoding for reads (FAKE quality score)
    """
    compressed_seq = []
    prev_base = None
    for base in seq:
        if base != prev_base:
            compressed_seq.append(base)
            prev_base = base
    # fake quality scores: to perform a fast conversion.
    compressed_qual = qual[:len(compressed_seq)]
    return ''.join(compressed_seq), compressed_qual
### Another version of RLE
def rle_2(seq, qual):
    """
    HPC encoding for a sequence  (FAKE quality score)
    """
    compressed_seq = ""
    # compressed_qual = ""
    prev_base = None
    for base in seq:
        if base != prev_base:
            compressed_seq += base
            # compressed_qual += q
            prev_base = base
    # fake quality scores: to perform a fast conversion.
    # Checkout rle-q.py for code that preservers quality scores.
    compressed_qual = qual[:len(compressed_seq)]
    return compressed_seq, compressed_qual

def process_fasta(input_file, output_file, batch_size):
    """
    Read in a chunk of fasta sequences and perform HPC encoding
    """
    with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
        counter = 0
        while True:
            batch = [in_file.readline() for i in range(batch_size * 2)]
            batch = [line.strip() for line in batch if line.strip()]
            if not batch:
                break
            for i in range(0, len(batch), 2):
                header = batch[i].strip()
                seq = batch[i + 1].strip()
                compressed_seq = rle(seq)
                out_file.write("{}\n{}\n".format(header, compressed_seq))
            print("processed {} sequences".format(counter * batch_size + int(len(batch) / 2)), file=sys.stderr)
            if len(batch) != batch_size * 2:
                break
            counter += 1

def process_reads(input_file, output_file, batch_size):
    """
    Read in a chunk of reads and do HPC encoding
    """
    # make a function map between open ang gzip.open
    import gzip
    open_map = {'gz': gzip.open, 'fq': open, 'fastq': open, 'fa': open, 'fasta': open}
    opener = open_map[input_file.split('.')[-1]]
    
    with opener(input_file, 'rt') as in_file, open(output_file, 'w') as out_file:
        counter = 0
        while True:
            print("Working on batch {} - time stamp {}".format(counter, datetime.datetime.now()), file=sys.stderr)
            batch = [in_file.readline().strip() for i in range(batch_size * 4)]
            if len(batch[-1]) == 0:
                print("-- refining last batch - time stamp {}".format(counter, datetime.datetime.now()), file=sys.stderr)
                batch = [line for line in batch if len(line) > 0]
            print("-- read batch {} - time stamp {}".format(counter, datetime.datetime.now()), file=sys.stderr)
            batch_len = len(batch)
            print("-- batch length: {} lines - time stamp {}".format(batch_len, datetime.datetime.now()), file=sys.stderr)
            # batch = [line.strip() for line in batch if line.strip()]
            if not batch:
                break
            print("-- HomoPolymer Compression - time stamp {}".format(datetime.datetime.now()), file=sys.stderr)
            for i in range(0, batch_len, 4):
                header = batch[i]
                seq = batch[i + 1]
                plus = batch[i + 2]
                qual = batch[i + 3]
                # compressed_seq, compressed_qual = rle(seq, qual)
                compressed_seq = re.sub(r"(.)\1+", r"\1", seq)
                compressed_qual = qual[:len(compressed_seq)]
                out_file.write("{}\n{}\n{}\n{}\n".format(header, compressed_seq, plus, compressed_qual))
            print("-- Compressed {} reads - time stamp {}".format(counter * batch_size + int(batch_len / 4), datetime.datetime.now()), file=sys.stderr)
            if batch_len != batch_size * 4:
                break
            counter += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform HPC encoding on FASTQ sequencing reads')
    parser.add_argument('-i', '--input', type=str, required=True, help='input FASTQ file path')
    parser.add_argument('-o', '--output', type=str, required=True, help='output FASTQ file path')
    parser.add_argument('-b', '--batch', type=int, default=1000, help='batch size (default: 1000)')
    
    args = parser.parse_args()
    
    # if args.input.endswith('.fasta'): run with fasta otherwise run with fastq
    if args.input.endswith('.fasta'):
        args.batch = 1
        # print that we are working with fasta and batch size is 1
        print("Working with fasta, batch size is set to 1")
        process_fasta(args.input, args.output, args.batch)
    else:
        print("Working with fastq, batch size is set to {}".format(args.batch))
        process_reads(args.input, args.output, args.batch)