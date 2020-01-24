#include "tsl/robin_map.h"
#include "tsl/robin_set.h"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/subgraph.hpp>

#if _OPENMP
#include <omp.h>
#endif

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt {
static unsigned threads = 1;
} // namespace opt

#define PROGRAM "physlr-split-minimizers"
#define PACKAGE_NAME "physlr"
#define GIT_REVISION "pre-autotools"

static uint64_t
memory_usage()
{
	int mem = 0;
	std::ifstream proc("/proc/self/status");
	for (std::string s; std::getline(proc, s);) {
		if (s.substr(0, 6) == "VmSize") {
			std::stringstream convert(s.substr(s.find_last_of('\t'), s.find_last_of('k') - 1));
			if (!(convert >> mem)) {
				return 0;
			}
			return mem;
		}
	}
	return mem;
}

struct vertexProperties
{
	std::string name = "";
	int weight = 0;
	size_t indexOriginal = 0;
};

struct edgeProperties
{
	int weight = 0;
};

struct edgeComponent_t
{
	enum
	{
		num = INT_MAX
	};
	using kind = boost::edge_property_tag;
} edgeComponent;

using graph_t = boost::subgraph<boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    vertexProperties,
    boost::property<
        boost::edge_index_t,
        int,
        boost::property<edgeComponent_t, std::size_t, edgeProperties>>>>;
using vertex_t = graph_t::vertex_descriptor;
using edge_t = graph_t::edge_descriptor;
using barcodeToIndex_t = std::unordered_map<std::string, vertex_t>;
using indexToBarcode_t = std::unordered_map<vertex_t, std::string>;
using bxToMolIdx_t = std::unordered_map<std::string, std::vector<vertex_t>>;
using BarcodeID = uint32_t;
using Minimizer = uint64_t;

static void
printVersion()
{
	const char VERSION_MESSAGE[] =
	    PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	            "Written by Johnathan Wong.\n"
	            "\n"
	            "Copyright 2020 Canada's Michael Smith Genome Science Centre\n";
	std::cerr << VERSION_MESSAGE << std::endl;
	exit(EXIT_SUCCESS);
}

static void
printHelpDialog()
{
	static const char dialog[] =
	    "Usage: physlr-split-minimizers [OPTION]... [GRAPH.tsv] [MINIMIZERS.tsv]\n"
	    "Split minimizers based on the molecule overlap graph.\n"
	    "  -t, --threads     threads [1]\n"
	    "  -v         enable verbose output\n"
	    "  --version     Print version\n"
	    "  --help     display this help and exit\n"
	    "Report bugs to <jowong@bcgsc.ca>.";
	std::cerr << dialog << std::endl;
	exit(0);
}

static void
printErrorMsg(const std::string& progname, const std::string& msg)
{
	std::cerr << progname << ": " << msg << "\nTry '" << progname
	          << " --help' for more information.\n";
}

void
readTSV(graph_t& g, const std::string& infile, bool verbose)
{
	std::cerr << "Loading graph" << std::endl;
#if _OPENMP
	double sTime = omp_get_wtime();
#endif
	barcodeToIndex_t barcodeToIndex;
	indexToBarcode_t indexToBarcode;

	infile == "-" ? "/dev/stdin" : infile;
	std::ifstream infileStream(infile);
	for (std::string line; std::getline(infileStream, line);) {
		if (line == "U\tm") {
			continue;
		}
		if (line.empty()) {
			break;
		}
		std::string node1;
		int weight;
		std::istringstream ss(line);
		if (ss >> node1 >> weight) {
			auto u = boost::add_vertex(g);
			g[u].name = node1;
			g[u].weight = weight;
			g[u].indexOriginal = u;
			barcodeToIndex[node1] = u;
			indexToBarcode[u] = node1;
		} else {
			printErrorMsg(PROGRAM, "unknown graph format");
			exit(EXIT_FAILURE);
		}
	}

	if (verbose) {
		std::cerr << "Loaded vertices to graph ";
#if _OPENMP
		std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
		sTime = omp_get_wtime();
#endif
	}
	for (std::string line; std::getline(infileStream, line);) {
		if (line == "U\tV\tm") {
			continue;
		}
		if (line.empty()) {
			printErrorMsg(PROGRAM, "unknown graph format");
			exit(EXIT_FAILURE);
		}
		std::string node1, node2;
		int weight;
		std::istringstream ss(line);
		if (ss >> node1 >> node2 >> weight) {
			auto E = boost::add_edge(barcodeToIndex[node1], barcodeToIndex[node2], g).first;
			g[E].weight = weight;
		} else {
			printErrorMsg(PROGRAM, "unknown graph format");
			exit(EXIT_FAILURE);
		}
	}

	if (verbose) {
		std::cerr << "Loaded edges to graph ";
	} else {
		std::cerr << "Loaded graph ";
	}
#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif
	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB" << std::endl;
}

void
findMoleculesPerBarcode(bxToMolIdx_t& bxToMolIdx, const graph_t& g)
{
	auto vertexItRange = boost::vertices(g);
	for (auto vertexIt = vertexItRange.first; vertexIt != vertexItRange.second; ++vertexIt) {
		std::string pattern = R"((\S+)_\d+_\d+$)";
		std::regex rgx(pattern);
		std::smatch matches;

		if (std::regex_search(g[*vertexIt].name, matches, rgx)) {
			auto bx = matches[1].str();
			bxToMolIdx[bx].emplace_back(*vertexIt);
		} else {
			std::cerr << "Unknown vertex Format" << std::endl;
			exit(1);
		}
	}
	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB" << std::endl;
}

void
splitMinimizers(
    bxToMolIdx_t& bxToMolIdx,
    std::vector<std::vector<Minimizer>>& barcodeToMinimizer,
    graph_t& g,
    tsl::robin_map<std::string, BarcodeID>& barcodes)
{
	size_t numBx = bxToMolIdx.size();
	// Canonical for loop for openMP
#if _OPENMP
#pragma omp parallel for
#endif
	for (size_t i = 0; i < numBx; ++i) {
		auto it = bxToMolIdx.begin();
		for (size_t j = 0; j < i; ++j) {
			++it;
		}

		std::stringstream ssOut;
		std::stringstream ssErr;

		auto& bx = (*it).first;
		auto& bxIdx = barcodes[bx];
		auto& minimizerSet = barcodeToMinimizer[bxIdx];

		for (auto& mol : (*it).second) {
			// Get Union of minimizers of neighbours
			auto neighbours = boost::adjacent_vertices(mol, g);
			std::unordered_set<Minimizer> neighbourMxsUnion;
			for (auto neighbourItr = neighbours.first; neighbourItr != neighbours.second;
			     ++neighbourItr) {
				std::string pattern = R"((\S+)_\d+_\d+$)";
				std::regex rgx(pattern);
				std::smatch matches;
				if (std::regex_search(g[*neighbourItr].name, matches, rgx)) {
					auto neighbourBx = matches[1].str();
					auto& neighbourBxIdx = barcodes[neighbourBx];
					auto& neighbourMxs = barcodeToMinimizer[neighbourBxIdx];
					for (auto& mx : neighbourMxs) {
						neighbourMxsUnion.insert(mx);
					}
				}
			}
			// Intersect minimizers of barcode with union
			std::vector<Minimizer> splitMinimizers;
			for (auto& mx : minimizerSet) {
				if (neighbourMxsUnion.find(mx) != neighbourMxsUnion.end()) {
					splitMinimizers.emplace_back(mx);
				}
			}

			if (splitMinimizers.empty()) {
				ssErr << "Warning: " << g[mol].name << " has no associated minimizers\n";
				ssOut << g[mol].name << "\t\n";
			} else {
				ssOut << g[mol].name << "\t";
				for (unsigned i = 0; i < splitMinimizers.size(); ++i) {
					if (i == splitMinimizers.size() - 1) {
						ssOut << splitMinimizers[i] << "\n";
					} else {
						ssOut << splitMinimizers[i] << " ";
					}
				}
			}
		}
#if _OPENMP
#pragma omp critical
#endif
		{
			std::cout << ssOut.str();
			std::cerr << ssErr.str();
		}
	}
	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB" << std::endl;
}

void
getBarcodeToMinimizer(
    tsl::robin_map<std::string, BarcodeID>& barcodes,
    std::vector<std::vector<Minimizer>>& barcodeToMinimizer,
    std::string& inputFile)
{
	// vector of barcodes
	std::vector<std::string> barcodeToStr;
	std::string barcodeBuffer;
	Minimizer minimizerBuffer;

	// read in minimizer file
	// format: GAGGTCCGTGGAGAGG-1	472493953667297251 1168973555595507959 342455687043295195
	// 283275954102976652

	std::ifstream fh;
	fh.open(inputFile);
	std::string line;
	std::cerr << "Loading file " << inputFile << std::endl;
	while (getline(fh, line)) {
		std::stringstream ss(line);
		ss >> barcodeBuffer;
		const auto& barcode = barcodes.find(barcodeBuffer);
		if (barcode == barcodes.end()) {
			barcodeToStr.emplace_back(barcodeBuffer);
			barcodes[barcodeBuffer] = barcodeToStr.size() - 1;
			barcodeToMinimizer.emplace_back(std::vector<Minimizer>());
			while (ss >> minimizerBuffer) {
				barcodeToMinimizer[barcodeToStr.size() - 1].emplace_back(minimizerBuffer);
			}
		} else {
			while (ss >> minimizerBuffer) {
				barcodeToMinimizer[barcode->second].emplace_back(minimizerBuffer);
			}
		}
	}
}

int
main(int argc, char* argv[])
{

	bool die = false;
	bool verbose = false;
	static int help = 0;
	static int version = 0;
	int optindex = 0;

	// long form arguments
	static struct option longopts[] = { { "help", no_argument, &help, 1 },
		                                { "threads", required_argument, nullptr, 't' },
		                                { "version", no_argument, &version, 1 },
		                                { nullptr, 0, nullptr, 0 } };

	for (int c; (c = getopt_long(argc, argv, "t:", longopts, &optindex)) != -1;) {
		switch (c) {
		case 't': {
			std::stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				std::cerr << "Error - Invalid parameters! t: " << optarg << std::endl;
				return 0;
			}
			break;
		}
		case 'v': {
			verbose = true;
			break;
		}
		default: {
			die = true;
			break;
		}
		}
	}
	if (help != 0) {
		printVersion();
		exit(EXIT_SUCCESS);
	} else if (version != 0) {
		printHelpDialog();
		exit(EXIT_SUCCESS);
	}

#if _OPENMP
	omp_set_num_threads(opt::threads);
#endif

	// Stores input file names
	std::vector<std::string> inputFiles;

	while (optind < argc) {
		inputFiles.emplace_back(argv[optind]);
		optind++;
	}

	if (inputFiles.empty()) {
		printErrorMsg(PROGRAM, "missing file operand");
		die = true;
	}

	if (die) {
		printHelpDialog();
		exit(EXIT_FAILURE);
	}

	// barcode to ID table (index in vector)
	tsl::robin_map<std::string, BarcodeID> barcodes;

	// barcodeID (index) to minimizer vector of vector
	// Note: Because a vector isn't a set, the input cannot have duplicates.
	std::vector<std::vector<Minimizer>> barcodeToMinimizer;

#if _OPENMP
	double sTime = omp_get_wtime();
#endif
	getBarcodeToMinimizer(barcodes, barcodeToMinimizer, inputFiles[1]);

#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif
	graph_t g;
	readTSV(g, inputFiles[0], verbose);
#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif

	bxToMolIdx_t bxToMolIdx;
	findMoleculesPerBarcode(bxToMolIdx, g);
	splitMinimizers(bxToMolIdx, barcodeToMinimizer, g, barcodes);
}
