/* Convert linked-reads to minimizers using ntHash-2.0.0
   usage:  physlr-indexlr -k K -w W [-v] file...
   Output: Each line of output is a barcode followed by a list of minimzers.
*/

#include <cassert>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <limits>

// ntHash 2.0.0
#include "nthash.h"
#include "ntHashIterator.h"

using namespace std;

// -------------------------- Prototypes ------------------------------
vector<uint64_t> hashKmers(string &readstr, size_t k);
set<uint64_t> getMinimizers1(vector<uint64_t> hashes, size_t w);
set<uint64_t> getMinimizers(vector<uint64_t> hashes, size_t w);
void printRead(string &barcode, string &readstr);
void printMinimizedRead(string barcode, set<uint64_t> &minSet);
void minimizeReads(std::istream*, size_t k, size_t w, bool verbose);
void printUsage(string);
void printErrorMsg(string, string);
// --------------------------------------------------------------------

int main(int argc, char *argv[]) {
	string progname = "physlr-indexlr";
	int c, help, optindex;
	unsigned k, w;
	string path;
	vector<string> infile;
	bool verbose, abort, w_set, k_set;
	k = w = 0;
	verbose = abort = w_set = k_set = false;
	optindex = help = 0;
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
			return EXIT_FAILURE;
		}
	}
	for (int i = optind; i < argc; ++i)
		infile.push_back(argv[i]);
	if (argc == 1 || help) {
		printUsage(progname);
		abort = true;
	}
	else if (!k_set) {
		printErrorMsg(progname, "missing option -- 'k'");
		abort = true;
	}
	else if (!w_set) {
		printErrorMsg(progname, "missing option -- 'w'");
		abort = true;
	}
	else if (k == 0) {
		printErrorMsg(progname,
				"option has incorrect argument -- 'k'");
		abort = true;
	}
	else if (w == 0) {
		printErrorMsg(progname,
				"option has incorrect argument -- 'w'");
		abort = true;
	}
	else if (infile.empty()) {
		printErrorMsg(progname, "missing file operand");
		abort = true;
	}
	if (abort)
		return EXIT_FAILURE;
	for (unsigned i = 0; i < infile.size(); ++i) {
		ifstream ifs;
		ifs.open(infile[i]);
		if (!ifs.is_open()) {
			cerr << "ERROR: Failed to open: " << infile[i] << "\n";
			continue;
		}
		minimizeReads(&ifs, k, w, verbose);
		ifs.close();
	}
	return EXIT_SUCCESS;
}

// ------------------------ Function declarations ---------------------

// Read a FASTQ file and reduce each read to a set of minimizers
void minimizeReads(std::istream *is, size_t k, size_t w, bool verbose) {
	size_t nread = 0, nline = 0, readlen;
	string line, comment, barcode;
	while (getline(*is, line)) {
		nline++;
		if (line[0] == '@') {
			istringstream iss(line);
			iss >> comment >> barcode;
			barcode = barcode.substr(5);
		}
		else if(line[0] == 'A' || line[0] == 'C' ||
			line[0] == 'G' || line[0] == 'T') {
			nread++;
			readlen = line.size();
			if (k > readlen) {
				if (verbose) {
					cerr << "WARN: Skip read " << nread;
					cerr << " on line " << nline << "; ";
					cerr << "k > read length ";
					cerr << "(k = " << k;
					cerr << ", read length = " << readlen << ")\n";
				}
				continue;
			}
			vector<uint64_t> hashes = hashKmers(line, k);
			// NOTE: The predicate P(#kmers != #hashes) will be true
			// when reads contains Ns, so check with the number of
			// hashes ntHash returns after hashing a read. ntHash
			// takes care to skip Ns.
			if (w > hashes.size()) {
				if (verbose) {
					cerr << "WARN: Skip read " << nread;
					cerr << " on line " << nline << "; ";
					cerr << "window size > #hashes ";
					cerr << "(w = " << w;
					cerr << ", #hashes = " << hashes.size() << ")\n";
				}
				continue;
			}			
			set<uint64_t> minSet = getMinimizers(hashes, w);
			printMinimizedRead(barcode, minSet);
			hashes.clear();
			minSet.clear();
		}
	}
}

// Hash the k-mers of a read
vector<uint64_t> hashKmers(string &readstr, size_t k) {
	vector<uint64_t> hashes;
	for (ntHashIterator iter(readstr, 1, k);
	     iter != iter.end(); ++iter) {
		hashes.push_back((*iter)[0]);
	}
	return hashes;
}

// Find the minimum hash values in a sliding window of size w. This is
// straight-forward: for each window, find the minimum.
set<uint64_t> getMinimizers(vector<uint64_t> hashes, size_t w) {
	set<uint64_t> minSet;
	for (size_t i = 0; i < hashes.size() - w + 1; ++i) {
		uint64_t min = hashes[i];
		for (size_t j = i + 1; j < i + w; ++j) {
			if (hashes[j] < min)
				min = hashes[j];
		}
		minSet.insert(min);
	}
	return minSet;
}

/* Find the minimum hash values in a sliding window of size w.

   This is slightly different from the straight-forward implementation in
   getMinimizers() but may not be faster, and, both have been tested to produce
   identical output.

   A window of hash values has one new hash value with respect to the previous
   window (and one less). If the new hash value is less than or equal to the
   minimum of the previous window, it is also the minimum of the current
   window. The if statement within the outer for loop checks for this condition
   and skips the inner for loop which finds the minimum in the current window.
 */
set<uint64_t> getMinimizers1(vector<uint64_t> hashes, size_t w) {
	set<uint64_t> minSet;
	uint64_t oldMin = 0;
	for (size_t i = 0; i < hashes.size() - w + 1; ++i) {
		uint64_t newMin = hashes[i + w - 1];
		if (newMin > oldMin) {
			for (size_t j = i; j < i + w - 1; ++j) {
				if (hashes[j] < newMin)
					newMin = hashes[j];
			}
		}
		minSet.insert(newMin);
		oldMin = newMin;
	}
	return minSet;
}

// Print the barcode, read and read length
void printRead(string &barcode, string &readstr) {
	cout << barcode << ' ' << readstr << ' ' << readstr.size() << '\n';
}

// Print the barcode and minimzers of a read
void printMinimizedRead(string barcode, set<uint64_t> &minSet) {
	cout << barcode << '\t';
	for (set<uint64_t>::iterator iter = minSet.begin(); iter != minSet.end();
	     ++iter)
		cout << *iter << ' ';
	cout << '\n';
}

void printUsage(string progname) {
	cout << "usage:  " << progname;
	cout << " -k K -w W [-v] file...\n\n";
	cout << " -k K     use K as k-mer size\n";
	cout << " -w W     use W as sliding-window size\n";
	cout << " -v       enable verbose output\n";
	cout << " file     space separated list of FASTQ files\n";
	cout << " --help   display this help and exit\n";
}

void printErrorMsg(string progname, string msg) {
	cerr << progname << ": " << msg << "\n";
	cerr << "Try 'physlr-indexlr --help' for more information.\n";
}

