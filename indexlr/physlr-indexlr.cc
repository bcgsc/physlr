/* Convert linked-reads to minimizers using ntHash-2.0.0
   Usage:  physlr-indexlr -k K -w W [-v] file...
   Output: Each line of output is a barcode followed by a list of minimzers.
*/

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// ntHash 2.0.0
#include "nthash.h"
#include "ntHashIterator.h"

using namespace std;

// Hash the k-mers of a read
static vector<uint64_t> hashKmers(const string &readstr, const size_t k) {
	vector<uint64_t> hashes;
	hashes.reserve(readstr.size() - k + 1);
	for (ntHashIterator iter(readstr, 1, k);
	     iter != iter.end(); ++iter) {
		hashes.push_back((*iter)[0]);
	}
	return hashes;
}

// Find all the minimizers from a list of hash values. A minimizer is the
// minimum hash value in a window of size w sliding along that list.
static vector<uint64_t> getMinimizers(const vector<uint64_t> &hashes,
				      const size_t w) {
	vector<uint64_t> minimizerSet;
	minimizerSet.reserve(w);
	int prevIndex = -1;
	vector<uint64_t>::const_iterator windowIter = hashes.begin();
	for (unsigned i = 1; i < hashes.size() - w + 1; ++i, ++windowIter) {
		vector<uint64_t>::const_iterator min =
			std::min_element(windowIter, windowIter + w);
		int minIndex = std::distance(windowIter, min);
		minIndex += i;
		if (minIndex > prevIndex) {
			prevIndex = minIndex;
			minimizerSet.push_back(*min);
		}
	}
	return minimizerSet;
}

// Print the barcode and minimzers of a read
static void printMinimizedRead(const string &barcode,
			       vector<uint64_t> &minimizerSet) {
	cout << barcode << "\t";
	for (auto &m : minimizerSet)
		cout << m << " ";
	cout << "\n";
}

// Read a FASTQ file and reduce each read to a set of minimizers
static void minimizeReads(std::istream *is, const size_t k, const size_t w, bool verbose) {
	size_t nread = 0, nline = 0;
	string sequence, header, comment, barcode;
	while (getline(*is, header)
	       && getline(*is, sequence)
	       && (*is).ignore(std::numeric_limits<std::streamsize>::max(), '\n')
	       && (*is).ignore(std::numeric_limits<std::streamsize>::max(), '\n')) {
		// Ignore the next two lines; '+' and quality string.
		nline += 4;
		nread += 1;
		// Extract the barcode.
		istringstream iss(header);
		iss >> comment >> barcode;
		string prefix = barcode.substr(0, 5);
		if (prefix != "BX:Z:") {
			std::cerr << "physlr-indexlr: error: expected BX:Z: and saw "
				  << prefix << " at line " << nline - 3 << "\n";
			exit(EXIT_FAILURE);
		}
		barcode = barcode.substr(5);
		// Validate parameters.
		if (k > sequence.size()) {
			if (verbose) {
				std::cerr << "warning: skip read " << nread
					  << " on line " << nline - 2
					  << "; k > read length " << "(k = " << k
					  << ", read length = " << sequence.size() << ")\n";
			}
			continue;
		}
		// Hash the kmers.
		vector<uint64_t> hashes = hashKmers(sequence, k);
		// NOTE: The predicate P(#kmers != #hashes) will be true when
		// reads contains Ns, so check with the number of hashes ntHash
		// gives after hashing a read. ntHash takes care to skip Ns.
		if (w > hashes.size()) {
			if (verbose) {
				std::cerr << "warning: skip read " << nread
					  << " on line " << nline - 2
					  << "; window size > #hashes (w = " << w
					  << ", #hashes = " << hashes.size() << ")\n";
			}
			continue;
		}
		// Minimerize.
		vector<uint64_t> minimizerSet = getMinimizers(hashes, w);
		printMinimizedRead(barcode, minimizerSet);
		hashes.clear();
		minimizerSet.clear();
	}
}

static void printUsage(const string &progname) {
	cout << "Usage:  " << progname <<
		"  -k K -w W [-v] file...\n\n"
		"  -k K     use K as k-mer size\n"
		"  -w W     use W as sliding-window size\n"
		"  -v       enable verbose output\n"
		"  --help   display this help and exit\n"
		"  file     space separated list of FASTQ files\n";
}

static void printErrorMsg(const string &progname, const string &msg) {
	cerr << progname << ": " << msg << "\n"
		"Try 'physlr-indexlr --help' for more information.\n";
}

int main(int argc, char *argv[]) {
	string progname = "physlr-indexlr";
	int c;
	int optindex = 0;
	int help     = 0;
	unsigned k   = 0;
	unsigned w   = 0;
	bool verbose = false;
	bool failed  = false;
	bool w_set   = false;
	bool k_set   = false;
	static const struct option longopts[]= {
				   {"help", no_argument, &help, 1},
				   {0, 0, 0, 0}
	};
	while((c = getopt_long(argc, argv, "k:w:v",
			       longopts, &optindex))!= -1) {
		switch (c) {
		case 0:
			break;
		case 'k':
			k_set = true;
			k = atoi(optarg);
			break;
		case 'w':
			w_set = true;
			w = atoi(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		default:
			exit(EXIT_FAILURE);
		}
	}
	std::vector<string> infile(&argv[optind], &argv[argc]);
	if (argc < 2) {
		printUsage(progname);
		exit(EXIT_FAILURE);
	}
	if (help) {
		printUsage(progname);
		exit(EXIT_SUCCESS);
	}
	else if (!k_set) {
		printErrorMsg(progname, "missing option -- 'k'");
		failed = true;
	}
	else if (!w_set) {
		printErrorMsg(progname, "missing option -- 'w'");
		failed = true;
	}
	else if (k == 0) {
		printErrorMsg(progname,
			      "option has incorrect argument -- 'k'");
		failed = true;
	}
	else if (w == 0) {
		printErrorMsg(progname,
			      "option has incorrect argument -- 'w'");
		failed = true;
	}
	else if (infile.empty()) {
		printErrorMsg(progname, "missing file operand");
		failed = true;
	}
	if (failed)
		exit(EXIT_FAILURE);
	for (auto &f : infile) {
		ifstream ifs(f);
		if (!ifs) {
			std::cerr << "phslyr-indexlr: error: failed to open: " << f << "\n";
			exit(EXIT_FAILURE);
		}
		minimizeReads(&ifs, k, w, verbose);
	}
	return 0;
}
