/* Convert linked-reads to minimizers using ntHash.
   Usage: ./minimizereads -k 100 -w 5 [-i <FASTQ file>] [-v] >outfile
   Output: Each line of output is a barcode followed by the list of minimzers.
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
void printUsage(char *argv_0);
// --------------------------------------------------------------------

int main(int argc, char *argv[]) {
	ifstream ifs;
	istream *is = &std::cin;
	int c;
	size_t k = 0, w = 0;
	string path;
	bool verbose = false;
	while((c = getopt(argc, argv, "k:w:i:v")) != -1) {
		switch (c) {
		case 'k':
			k = atoi(optarg);
			break;
		case 'w':
			w = atoi(optarg);
			break;
		case 'i':
			path = string(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		default:
			printUsage(argv[0]);
			return EXIT_FAILURE;
		}
	}
	// Switch input stream to a file
	if (!path.empty()) {
		ifs.open(path);
		if (!ifs.is_open()) {
			cerr << "ERROR: Failed to open: " << path << "\n";
			return EXIT_FAILURE;
		}
		is = &ifs;
	}
	if (k == 0 || w == 0) {
		cerr << "ERROR: Missing option\n";
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}
	minimizeReads(is, k, w, verbose);
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

void printUsage(char *argv_0) {
	cout << "Usage:  " << argv_0;
	cout << " -w <window> -k <kmer> [-i <file>] [-v]\n";
	cout << "\t <window> = Size of sliding window\n";
	cout << "\t <kmer>   = k-mer length\n";
	cout << "\t <file>   = Input FASTQ file\n";
}
/* -------------------------------------------------------------------- */
