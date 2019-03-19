/*
 * MinimizerOverlap.cpp
 *
 *  Created on: Feb 11, 2019
 *      Author: cjustin
 */

#include "config.h"
#include "Options.cpp"
#include "tsl/robin_map.h"
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <omp.h>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#define PROGRAM "physlr-overlap"

void printVersion() {
	const char VERSION_MESSAGE[] =
	GIT_REVISION " (" PACKAGE_NAME ") " GIT_REVISION
	"\n"
	"Written by Justin Chu.\n"
	"\n"
	"Copyright 2019 Canada's Michael Smith Genome Science Centre\n";
	std::cerr << VERSION_MESSAGE << std::endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog() {
	static const char dialog[] =
			"Usage: physlr-overlap [OPTION]... [MINIZERS.tsv]\n"
					"Read a sketch of linked reads and find overlapping barcodes.\n"
					"  -n, --min-n=INT   remove edges with fewer than n shared markers [0].\n"
					"  -t, --threads=INT Number of threads [1].\n"
					"Report bugs to <cjustin@bcgsc.ca>.";
	std::cerr << dialog << std::endl;
	exit(0);
}

//invertible hash function
struct fastHash {
	std::size_t operator ()(uint64_t key) const {
		key = (~key + (key << 21u)); // key = (key << 21) - key - 1;
		key = key ^ key >> 24u;
		key = ((key + (key << 3u)) + (key << 8u)); // key * 265
		key = key ^ key >> 14u;
		key = ((key + (key << 2u)) + (key << 4u)); // key * 21
		key = key ^ key >> 28u;
		key = (key + (key << 31u));
		return key;
	}
};

//returns memory of program in kb
int memory_usage() {
	int mem = 0;
	std::ifstream proc("/proc/self/status");
	std::string s;
	while (getline(proc, s), !proc.fail()) {
		if (s.substr(0, 6) == "VmSize") {
			std::stringstream convert(
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
	static struct option long_options[] = { { "min-n", required_argument,
			nullptr, 'n' }, { "threads", required_argument, nullptr, 't' }, {
			nullptr, 0, nullptr, 0 } };

	int i = 0;
	while ((c = getopt_long(argc, argv, "n:t:", long_options, &i)) != -1) {
		switch (c) {
		case 't': {
			std::stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				std::cerr << "Error - Invalid parameters! t: " << optarg
						<< std::endl;
				return 0;
			}
			break;
		}
		case 'n': {
			std::stringstream convert(optarg);
			if (!(convert >> opt::minN)) {
				std::cerr << "Error - Invalid parameters! n: " << optarg
						<< std::endl;
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
	std::vector<std::string> inputFiles;

	while (optind < argc) {
		inputFiles.emplace_back(argv[optind]);
		optind++;
	}

	if (inputFiles.empty()) {
		std::cerr << "Missing Input Files" << std::endl;
		die = true;
	}

	if (die) {
		printHelpDialog();
		exit(EXIT_FAILURE);
	}

	//constuct minimizers to barcodes
	//barcode to ID table (index in vector)
	tsl::robin_map<std::string, BarcodeID> barcodes;

	//vector of barcodes
	tsl::robin_map<Minimizer, std::set<BarcodeID>> minimizerToBarcode;
	std::vector<std::string> barcodeToStr;

	double sTime = omp_get_wtime();
	std::string barcodeBuffer;
	Minimizer minimizerBuffer;

	//read in minimizer file
	//format: GAGGTCCGTGGAGAGG-1	472493953667297251 1168973555595507959 342455687043295195 283275954102976652
	for (std::vector<std::string>::const_iterator itr = inputFiles.begin();
			itr != inputFiles.end(); itr++) {
		std::ifstream fh;
		fh.open(itr->c_str());
		std::string line;
		std::cerr << "Loading file" << *itr << std::endl;
		while (getline(fh, line)) {
			std::stringstream ss(line);
			ss >> barcodeBuffer;
			const auto &barcode = barcodes.find(barcodeBuffer);
			if (barcode == barcodes.end()) {
				barcodeToStr.emplace_back(barcodeBuffer);
				barcodes[barcodeBuffer] = barcodeToStr.size() - 1;
				while (ss >> minimizerBuffer) {
					minimizerToBarcode[minimizerBuffer].insert(
							barcodeToStr.size() - 1);
				}
			} else {
				while (ss >> minimizerBuffer) {
					minimizerToBarcode[minimizerBuffer].insert(barcode->second);
				}
			}
		}
	}

	std::cerr << "Finished constructing minimizerToBarcodes in sec: "
			<< omp_get_wtime() - sTime << std::endl;
	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576)
			<< "GB" << std::endl;
	sTime = omp_get_wtime();

	//store into 2d matrix / hash table
	//todo revisit Counts? -> can be smaller
	using SimMat = tsl::robin_map<uint64_t, Count, fastHash>;
	SimMat barcodeSimMat;

	//counts of vector
	std::vector<Count> barcodeCount(barcodes.size(), 0);

	std::cerr << "Populating Overlaps" << std::endl;
	std::cerr << "Total Minimizers: " << minimizerToBarcode.size() << std::endl;
	std::cerr << "Total Barcodes: " << barcodes.size() << std::endl;

//	11706234 edges, 68904 barcodes, 2615306 minimizers
	//todo some how parallelize this
	for (const auto &itr : minimizerToBarcode) {
		for (auto barcode_i = itr.second.begin(); barcode_i != itr.second.end();
				barcode_i++) {
			barcodeCount[*barcode_i]++;
			auto barcode_j = std::set<BarcodeID>::const_iterator(barcode_i);
			for (barcode_j++; barcode_j != itr.second.end(); barcode_j++) {
				barcodeSimMat[static_cast<size_t>(*barcode_i) << 32u | *barcode_j]++;
			}
		}
	}

	std::cerr << "Finished computing overlaps in sec: "
			<< omp_get_wtime() - sTime << std::endl;
	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576)
			<< "GB" << std::endl;
	std::cerr << "Total number of unfiltered edges: " << barcodeSimMat.size()
			<< std::endl;

	std::cout << "U\tn\n";
	std::string bufferString;
	//print out vertexes + counts
#pragma omp parallel
	for (const auto &itr : barcodes) {
		bufferString.clear();
		bufferString += itr.first;
		bufferString += "\t";
		bufferString += std::to_string(barcodeCount[itr.second]);
		bufferString += "\n";
		std::cout << bufferString;
	}

	size_t edgeCount = 0;
	std::cout << "\nU\tV\tn\n";
#pragma omp parallel
	for (const auto &itr : barcodeSimMat) {
		//filter by n
		if (opt::minN <= itr.second) {
			bufferString.clear();
			bufferString += barcodeToStr[itr.first >> 32u];
			bufferString += "\t";
			bufferString += barcodeToStr[static_cast<uint32_t>(itr.first)];
			bufferString += "\t";
			bufferString += std::to_string(itr.second);
			bufferString += "\n";
			std::cout << bufferString;
			++edgeCount;
		}
	}
	std::cerr << "Total number of filtered edges: " << edgeCount << std::endl;

	return 0;
}

