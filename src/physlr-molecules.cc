#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "tsl/robin_set.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/subgraph.hpp>

#if _OPENMP
#include <omp.h>
#endif

#define PROGRAM "physlr-molecules"
#define PACKAGE_NAME "physlr"
#define GIT_REVISION "pre-autotools"

using namespace boost;

struct vertexProperties
{
	std::string name;
	int weight;
	int indexOriginal;
};

struct edgeProperties
{
	int weight;
};

namespace boost {
struct edge_component_t
{
	enum
	{
		num = 555
	};
	typedef edge_property_tag kind;
} edge_component;
}

using graph_t = subgraph<adjacency_list<
    vecS,
    vecS,
    undirectedS,
    vertexProperties,
    property<edge_index_t, int, property<edge_component_t, std::size_t, edgeProperties>>>>;
using vertex_t = graph_t::vertex_descriptor;
using edge_t = graph_t::edge_descriptor;
using barcodeToIndex_t = std::unordered_map<std::string, vertex_t>;
using indexToBarcode_t = std::unordered_map<vertex_t, std::string>;

static void
printVersion()
{
	const char VERSION_MESSAGE[] =
	    PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	            "Written by Johnathan Wong.\n"
	            "\n"
	            "Copyright 2019 Canada's Michael Smith Genome Science Centre\n";
	std::cerr << VERSION_MESSAGE << std::endl;
	exit(EXIT_SUCCESS);
}

static void
printErrorMsg(const std::string& progname, const std::string& msg)
{
	std::cerr << progname << ": " << msg
	          << "\nTry 'physlr-molecules --help' for more information.\n";
}

static void
printUsage(const std::string& progname)
{
	std::cout << "Usage:  " << progname
	          << "  [-s SEPARATION-STRATEGY] [-v -t T] FILE...\n\n"
	             "  -v         enable verbose output\n"
	             "  -s --separation-strategy   \n"
	             "  -t N       use N number of threads [1]\n"
	             "  --help     display this help and exit\n"
	             "  SEPARATION-STRATEGY      space separated list of ntHits tsv output files\n";
}

void
printGraph(graph_t& g, indexToBarcode_t& indexToBarcode)
{
	std::cout << "U\tn" << std::endl;
	std::string node1, node2;
	int weight;
	auto vertexItRange = vertices(g);
	for (auto vertexIt = vertexItRange.first; vertexIt != vertexItRange.second; ++vertexIt) {
		node1 = indexToBarcode[*vertexIt];
		weight = g[*vertexIt].weight;
		std::cout << node1 << "\t" << weight << "\n";
	}
	std::cout << "\nU\tV\tn" << std::endl;
	auto edgeItRange = edges(g);
	for (auto edgeIt = edgeItRange.first; edgeIt != edgeItRange.second; ++edgeIt) {
		weight = g[*edgeIt].weight;
		node1 = indexToBarcode[boost::source(*edgeIt, g)];
		node2 = indexToBarcode[boost::target(*edgeIt, g)];
		std::cout << node1 << "\t" << node2 << "\t" << weight << "\n";
	}
}

int
main(int argc, char* argv[])
{

	auto progname = "physlr-molecules";
	int c;
	int optindex = 0;
	char* end = nullptr;
	static int help = 0;
	std::string separationStrategy;
	// bool verbose = false;
	unsigned t = 1;
	// bool failed = false;
	static const struct option longopts[] = {
		{ "help", no_argument, &help, 1 },
		{ "separation-strategy", required_argument, nullptr, 's' },
		{ nullptr, 0, nullptr, 0 }
	};
	while ((c = getopt_long(argc, argv, "s:vt:", longopts, &optindex)) != -1) {
		switch (c) {
		case 0:
			break;
		case 's':
			separationStrategy.assign(optarg);
			break;
		case 'v':
			// verbose = true;
			break;
		case 't':
			t = strtoul(optarg, &end, 10);
			break;
		default:
			exit(EXIT_FAILURE);
		}
	}
	std::vector<std::string> infiles(&argv[optind], &argv[argc]);
	if (argc < 1) {
		printUsage(progname);
		exit(EXIT_FAILURE);
	}
	if (help != 0) {
		printVersion();
		printUsage(progname);
		exit(EXIT_SUCCESS);
	} else if (infiles.empty()) {
		printErrorMsg(progname, "missing file operand");
		// failed = true;
	}
	t = t + 1;
	// TODO add checks for separation strategy
	/*if (failed) {
	    printUsage(progname);
	    exit(EXIT_FAILURE);
	}*/

	/*enum
	{
	    A,
	    B,
	    C,
	    D,
	    E,
	    F
	};
*/

	graph_t g;

	barcodeToIndex_t barcodeToIndex;
	indexToBarcode_t indexToBarcode;
	for (auto& infile : infiles) {
		infile == "-" ? "/dev/stdin" : infile;
		vertex_t U;
		edge_t E;
		std::ifstream infileStream(infile);
		std::string line;
		bool atEdges = false;
		while (std::getline(infileStream, line)) {
			if (line == "" or line == "U\tn") {
				continue;
			}
			if (line == "" or line == "U\tV\tn") {
				atEdges = true;
				continue;
			}

			std::string node1, node2;
			int weight;
			std::istringstream ss(line);
			ss >> node1;
			if (!atEdges) {
				ss >> weight;
				U = add_vertex(g);
				g[U].name = node1;
				g[U].weight = weight;
				g[U].indexOriginal = U;
				barcodeToIndex[node1] = U;
				indexToBarcode[U] = node1;
			} else {
				ss >> node2 >> weight;
				E = add_edge(barcodeToIndex[node1], barcodeToIndex[node2], g).first;
				g[E].weight = weight;
			}
		}
	}
	// barcodeToIndex.clear();
	// boost::print_graph(g);
	auto vertexItRange = vertices(g);
	for (auto vertexIt = vertexItRange.first; vertexIt != vertexItRange.second; ++vertexIt) {
		std::cout << "Index: " << *vertexIt << std::endl;
		auto neighbours = boost::adjacent_vertices(*vertexIt, g);
		/*graph_t& g1 = g.create_subgraph();
		for (auto vd : make_iterator_range(neighbours))
		    boost::add_vertex(vd, g1);
		boost::print_graph(g1);
		std::cout << "Index: " << *vertexIt << std::endl;*/
		graph_t& g2 = g.create_subgraph(neighbours.first, neighbours.second);
		// boost::print_graph(g2);
		property_map<graph_t, edge_component_t>::type component = get(edge_component, g2);
		std::size_t num_comps = biconnected_components(g2, component);
		std::cerr << "Found " << num_comps << " biconnected components.\n";

		std::vector<vertex_t> art_points;
		articulation_points(g2, std::back_inserter(art_points));
		std::cerr << "Found " << art_points.size() << " articulation points.\n";
		graph_traits<graph_t>::edge_iterator ei, ei_end;
		for (boost::tie(ei, ei_end) = edges(g2); ei != ei_end; ++ei)
			std::cout << g2[source(*ei, g2)].indexOriginal << " -- "
			          << g2[target(*ei, g2)].indexOriginal << "[label=\"" << component[*ei]
			          << "\"]\n";
		std::cout << "}\n";
		// std::cout << "len of component: " << component.size() << std::endl;
	}

	printGraph(g, indexToBarcode);
}
