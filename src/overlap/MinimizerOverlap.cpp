/*
 * MinimizerOverlap.cpp
 *
 *  Created on: Feb 11, 2019
 *      Author: cjustin
 */

#include "config.h"
#include "iostream"
#include "stdlib.h"
#include "cassert"
#include "sstream"
#include <getopt.h>
#include <omp.h>
#include <limits>
#include <math.h>
#include <string>
//#include "Options.h"
#include "Options.cpp"
#include <google/dense_hash_map>
//#include <google/dense_hash_set>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <fstream>
using namespace std;

#define PROGRAM "physlr-overlap"

void printVersion() {
	const char VERSION_MESSAGE[] =
	GIT_REVISION " (" PACKAGE_NAME ") " GIT_REVISION
	"\n"
	"Written by Justin Chu.\n"
	"\n"
	"Copyright 2019 Canada's Michael Smith Genome Science Centre\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog() {
	static const char dialog[] = "Usage: physlr-overlap [OPTION]... [MINIZERS.tsv]\n"
			"Read a sketch of linked reads and find overlapping barcodes.\n"
			"  -n, --min-n=INT   remove edges with fewer than n shared markers [0].\n"
			"  -t, --threads=INT Number of threads [1].\n"
			"Report bugs to <cjustin@bcgsc.ca>.";
	cerr << dialog << endl;
	exit(0);
}

// Only for pairs of std::hash-able types for simplicity.
// You can of course template this struct to allow other hash functions
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;
    }
};

int main(int argc, char *argv[]) {
	bool die = false;

	//switch statement variable
	int c;

	//long form arguments
	static struct option long_options[] = {
			{ "min-n", required_argument, NULL,'n' },
			{ "threads", required_argument, NULL,'t' },
			{ NULL, 0, NULL, 0 }};

	int i = 0;
	while ((c = getopt_long(argc, argv, "n:t:", long_options, &i))
			!= -1) {
		switch (c) {
		case 't': {
			stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				cerr << "Error - Invalid parameters! t: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'n': {
			stringstream convert(optarg);
			if (!(convert >> opt::minN)) {
				cerr << "Error - Invalid parameters! n: " << optarg << endl;
				return 0;
			}
			break;
		}
		default: {
			die = true;
			break;
		}
		}
	}
	omp_set_num_threads(opt::threads);

	//Stores fasta input file names
	vector<string> inputFiles;

	while (optind < argc) {
		inputFiles.push_back(argv[optind]);
		optind++;
	}

	if (inputFiles.size() == 0) {
		cerr << "Missing Input Files" << endl;
		die = true;
	}

	if (die) {
		printHelpDialog();
		exit(EXIT_FAILURE);
	}


	//constuct minimizers to barcodes
	//barcode to ID table (index in vector)
	google::dense_hash_map<string,BarcodeID> barcodes;
	barcodes.set_empty_key("");

	//vector of barcodes
	//TODO preallocate smaller number?
//	unsigned barcodeMaxCount = 4000000;
//	vector<unordered_set<Minimizer>> barcodeToMinimizer(barcodeMaxCount, unordered_set<Minimizer>());

	unordered_map<Minimizer, unordered_set<BarcodeID>> minimizerToBarcode;

	BarcodeID numBarcodes = 0;

	//read in minimizer file
	//format: GAGGTCCGTGGAGAGG-1	472493953667297251 1168973555595507959 342455687043295195 283275954102976652
	for(vector<string>::const_iterator itr = inputFiles.begin(); itr != inputFiles.end(); itr++){
		ifstream fh;
		fh.open(itr->c_str());
		string line;
		while (getline(fh, line)) {
			stringstream ss(line);
			string barcode;
			ss >> barcode;
			Minimizer minimizer;
			vector<Minimizer> minimizerList;
			if (barcodes.find(barcode) != barcodes.end()) {
				barcodes[barcode] = numBarcodes++;
			}
			while(ss >> minimizer){
				minimizerToBarcode[minimizer].insert(barcodes[barcode]);
			}
		}
	}

	//store into 2d matrix / hash table
	//todo revisit Counts? -> can be smaller
	typedef google::dense_hash_map< pair<BarcodeID,BarcodeID>, Count, pair_hash> SimMat;
	SimMat barcodeSimMat;
	barcodeSimMat.set_empty_key(pair<BarcodeID, BarcodeID>(numBarcodes, numBarcodes));

	//counts of vector
	vector<Count> barcodeCount(barcodes.size(), 0);

	//todo some how parallelize me?
	for (unordered_map<Minimizer, unordered_set<BarcodeID>>::const_iterator itr =
			minimizerToBarcode.begin(); itr != minimizerToBarcode.end();
			itr++) {
		for (unordered_set<BarcodeID>::const_iterator i = itr->second.begin();
				i != itr->second.end(); i++) {
			BarcodeID barcode1 = *i;
			barcodeCount[barcode1]++;
			for (unordered_set<BarcodeID>::const_iterator j =
					itr->second.begin(); j != itr->second.end(); j++) {
				BarcodeID barcode2 = *j;
				if (barcode1 < barcode2) {
					pair<BarcodeID, BarcodeID> tempPair(barcode1, barcode2);
					SimMat::iterator result = barcodeSimMat.find(tempPair);
					if (result == barcodeSimMat.end()) {
						barcodeSimMat[tempPair] = 1;
					} else {
						result->second++;
					}
				}
			}
		}
	}

	cout << "U\tn\n";
	string bufferString = "";
	//print out vertexes + counts
#pragma omp parallel
	for (google::dense_hash_map<string,BarcodeID>::const_iterator itr =
			barcodes.begin(); itr != barcodes.end();
			itr++) {
		bufferString.clear();
		bufferString += itr->first;
		bufferString += "\t";
		bufferString += barcodeCount[itr->second].size();
		bufferString += "\n";
		cout << bufferString;
	}

	cout << "U\tV\tn\n";
#pragma omp parallel
	for (SimMat::const_iterator itr =
			barcodes.begin(); itr != barcodes.end();
			itr++) {
		bufferString.clear();
		bufferString += itr->first.first;
		bufferString += "\t";
		bufferString += itr->first.second;
		bufferString += "\t";
		bufferString += itr->second;
		bufferString += "\n";
		cout << bufferString;
	}
	return 0;
}

