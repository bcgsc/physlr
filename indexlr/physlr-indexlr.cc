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

// Find the minimum hash values in a sliding window of size w.
//
// A window of hash values has one new hash value with respect to the previous
// window (and one less). If the new hash value is less than or equal to the
// minimum of the previous window, it is also the minimum of the current
// window. The if statement within the outer for loop checks for this condition
// and skips the inner for loop which finds the minimum in the current window.
static vector<uint64_t> getMinimizers(const vector<uint64_t> &hashes,
				      const size_t w) {
	vector<uint64_t> minSet;
	minSet.reserve(w);
	uint64_t oldMin = 0;
	for (size_t i = 0; i < hashes.size() - w + 1; ++i) {
		uint64_t newMin = hashes[i + w - 1];
		if (newMin > oldMin) {
			for (size_t j = i; j < i + w - 1; ++j) {
				if (hashes[j] < newMin)
					newMin = hashes[j];
			}
		}
		vector<uint64_t>::iterator iter =
			std::lower_bound(minSet.begin(), minSet.end(), newMin);
		if (iter == minSet.end() || newMin < *iter)
			minSet.insert(iter, newMin);
		oldMin = newMin;
	}
	return minSet;
}

// Print the barcode and minimzers of a read
static void printMinimizedRead(const string &barcode,
			       vector<uint64_t> &minSet) {
	cout << barcode << "\t";
	for (vector<uint64_t>::iterator iter = minSet.begin();
	     iter != minSet.end(); ++iter)
		cout << *iter << " ";
	cout << "\n";
}

// Read a FASTQ file and reduce each read to a set of minimizers
static void minimizeReads(std::istream *is, const size_t k, const size_t w,
			  bool verbose) {
	size_t nread = 0, nline = 0, readlen;
	bool io_error = false;
	while (!is->eof() && !io_error) {
		bool barcode_error = false;
		string line, comment, barcode;
		for (int i = 1; i <= 4; ++i) {
			if (!getline(*is, line)) {
				io_error = true;
				break;
			}
			++nline;
			if (i == 1) {
				istringstream iss(line);
				iss >> comment >> barcode;
				string prefix = barcode.substr(0, 5);
				if (prefix != "BX:Z:") {
					if (verbose) {
						cerr << "WARN: Unknown "
						     << "barcode prefix '"
						     << prefix
						     << "' on line " << nline
						     << "\n";
					}
					barcode_error = true;
					continue;
				}
				barcode = barcode.substr(5);
			}
			else if (i == 2) {
				++nread;
				if (barcode_error) {
					if (verbose) {
						cerr << "WARN: Skip read "
						     << nread
						     << " on line " << nline
						     << "; malformed barcode "
						     << "on line " << nline - 1
						     << "\n";
					}
					continue;
				}
				readlen = line.size();
				if (k > readlen) {
					if (verbose) {
						cerr << "WARN: Skip read "
						     << nread
						     << " on line " << nline
						     << "; k > read length "
						     << "(k = " << k
						     << ", read length = "
						     << readlen << ")\n";
					}
					continue;
				}
				vector<uint64_t> hashes = hashKmers(line, k);
				// NOTE: The predicate P(#kmers != #hashes) will
				// be true when reads contains Ns, so check with
				// the number of hashes ntHash gives after
				// hashing a read. ntHash takes care to skip Ns.
				if (w > hashes.size()) {
					if (verbose) {
						cerr << "WARN: Skip read "
						     << nread
						     << " on line " << nline
						     << "; window size > #hashes"
						     << "(w = " << w
						     << ", #hashes = "
						     << hashes.size() << ")\n";
					}
					continue;
				}
				vector<uint64_t> minSet =
					getMinimizers(hashes, w);
				printMinimizedRead(barcode, minSet);
				hashes.clear();
				minSet.clear();
			}
		}
	}
}

void printUsage(const string &progname) {
	cout << "Usage:  " << progname <<
		"  -k K -w W [-v] file...\n\n"
		"  -k K     use K as k-mer size\n"
		"  -w W     use W as sliding-window size\n"
		"  -v       enable verbose output\n"
		"  --help   display this help and exit\n"
		"  file     space separated list of FASTQ files\n";
}

void printErrorMsg(const string &progname, const string &msg) {
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
	struct option longopts[]= {
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
	vector<string> infile;
	infile.reserve(argc - optind);
	for (int i = optind; i < argc; ++i)
		infile.push_back(argv[i]);
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
	for (unsigned i = 0; i < infile.size(); ++i) {
		ifstream ifs(infile[i]);
		if (!ifs) {
			cerr << "ERROR: Failed to open: " << infile[i] << "\n";
			continue;
		}
		minimizeReads(&ifs, k, w, verbose);
	}
	return 0;
}
