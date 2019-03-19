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
#include <google/sparse_hash_map>
#include "tsl/robin_map.h"
//#include <google/dense_hash_set>
//#include "SimMat.h"
#include "city.h"
#include "city.cc"
#include <set>
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
//    template <class T1, class T2>
    std::size_t operator () (const std::pair<string,string> &p) const {
        auto h1 = CityHash64(p.first.c_str(), p.first.length());
        auto h2 = CityHash64(p.second.c_str(), p.first.length());

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;
    }
};

struct fastHash {
    std::size_t operator () (uint64_t key) const {
		key = (~key + (key << 21)); // key = (key << 21) - key - 1;
		key = key ^ key >> 24;
		key = ((key + (key << 3)) + (key << 8)); // key * 265
		key = key ^ key >> 14;
		key = ((key + (key << 2)) + (key << 4)); // key * 21
		key = key ^ key >> 28;
		key = (key + (key << 31));
		return key;
    }
};

//returns memory of program in kb
int memory_usage() {
	int mem = 0;
	ifstream proc("/proc/self/status");
	string s;
	while (getline(proc, s), !proc.fail()) {
		if (s.substr(0, 6) == "VmSize") {
			stringstream convert(
					s.substr(s.find_last_of('\t'), s.find_last_of('k') - 1));
			if (!(convert >> mem)) {
				return 0;
			}
			return mem;
		}
	}
	return mem;
}

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

//	unordered_map<Minimizer, unordered_set<BarcodeID>> minimizerToBarcode;
	//TODO optimize pre-set size
	tsl::robin_map<Minimizer, set<BarcodeID>> minimizerToBarcode;
	vector<string> barcodeToStr;

	double sTime = omp_get_wtime();
	string barcodeBuffer;
	Minimizer minimizerBuffer;

	//read in minimizer file
	//format: GAGGTCCGTGGAGAGG-1	472493953667297251 1168973555595507959 342455687043295195 283275954102976652
	for(vector<string>::const_iterator itr = inputFiles.begin(); itr != inputFiles.end(); itr++){
		ifstream fh;
		fh.open(itr->c_str());
		string line;
		cerr << "Loading file" << *itr << endl;
		while (getline(fh, line)) {
			stringstream ss(line);
			ss >> barcodeBuffer;
			const auto &barcode = barcodes.find(barcodeBuffer);
			if (barcode == barcodes.end()) {
				barcodeToStr.push_back(barcodeBuffer);
				barcodes[barcodeBuffer] = barcodeToStr.size() - 1;
				while(ss >> minimizerBuffer){
					minimizerToBarcode[minimizerBuffer].insert(barcodeToStr.size() - 1);
				}
			}
			else{
				while(ss >> minimizerBuffer){
					minimizerToBarcode[minimizerBuffer].insert(barcode->second);
				}
			}
		}
	}

	cerr << "Finished constructing minimizerToBarcodes in sec: "
			<< omp_get_wtime() - sTime << endl;
	cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB"
			<< endl;
	sTime = omp_get_wtime();

	//store into 2d matrix / hash table
	//todo revisit Counts? -> can be smaller
	typedef tsl::robin_map<uint64_t, Count, fastHash> SimMat;
//	typedef unordered_map< uint64_t, Count, std::hash<uint64_t>> SimMat;
	SimMat barcodeSimMat;
//	barcodeSimMat.set_empty_key((size_t) barcodeToStr.size() << 32 | barcodeToStr.size());
//	barcodeSimMat.set_empty_key(pair<string,string>("",""));

	//counts of vector
	vector<Count> barcodeCount(barcodes.size(), 0);

	cerr << "Populating Overlaps" << endl;
	cerr << "Total Minimizers: " << minimizerToBarcode.size() << endl;
	cerr << "Total Barcodes: " << barcodes.size() << endl;

//	11706234 edges, 68904 barcodes, 2615306 minimizers
	//todo some how parallelize this
	for (const auto &itr : minimizerToBarcode) {
		for (auto barcode_i = itr.second.begin();
				barcode_i != itr.second.end(); barcode_i++) {
			barcodeCount[*barcode_i]++;
			auto barcode_j = std::set<BarcodeID>::const_iterator(barcode_i);
			for (barcode_j++; barcode_j != itr.second.end(); barcode_j++) {
				barcodeSimMat[(size_t) *barcode_i << 32 | *barcode_j]++;
//				if(barcodeSimMat.find((size_t) *barcode_i << 32 | *barcode_j) == barcodeSimMat.end())
//				{
//					barcodeSimMat[(size_t) *barcode_i << 32 | *barcode_j] = 0;
//				}
//				else{
//					barcodeSimMat[(size_t) *barcode_i << 32 | *barcode_j]++;
//				}
		    }
		}
	}

//	SimMat simMat(barcodes.size());
//	for (const auto &itr : minimizerToBarcode) {
//		for (std::set<BarcodeID>::const_iterator barcode_i = itr.second.begin();
//				barcode_i != itr.second.end(); barcode_i++) {
//			auto barcode_j = std::set<BarcodeID>::const_iterator(barcode_i);
//			for (barcode_j++; barcode_j != itr.second.end(); barcode_j++) {
//				simMat.updateCount(*barcode_i , *barcode_j);
//		    }
//		}
//	}

	cerr << "Finished computing overlaps in sec: " << omp_get_wtime() - sTime << endl;
	cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB"
			<< endl;

	cerr << "Total number of unfiltered edges: " << barcodeSimMat.size() << endl;

	cout << "U\tn\n";
	string bufferString = "";
	//print out vertexes + counts
#pragma omp parallel
	for (google::dense_hash_map<string, BarcodeID>::const_iterator itr =
			barcodes.begin(); itr != barcodes.end(); itr++) {

		bufferString.clear();
		bufferString += itr->first;
		bufferString += "\t";
		bufferString += to_string(barcodeCount[itr->second]);
		bufferString += "\n";
		cout << bufferString;
	}

	size_t edgeCount = 0;
	cout << "\nU\tV\tn\n";
#pragma omp parallel
	for (SimMat::const_iterator itr =
			barcodeSimMat.begin(); itr != barcodeSimMat.end();
			itr++) {
		//filter by n
		if (opt::minN <= itr->second) {
			bufferString.clear();
			bufferString += barcodeToStr[itr->first >> 32];
			bufferString += "\t";
			bufferString += barcodeToStr[(uint32_t) itr->first];
			bufferString += "\t";
			bufferString += to_string(itr->second);
			bufferString += "\n";
			cout << bufferString;
			++edgeCount;
		}
	}
	cerr << "Total number of filtered edges: " << edgeCount << endl;

	return 0;
}

