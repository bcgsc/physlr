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

static int
memory_usage()
{
	int mem = 0;
	std::ifstream proc("/proc/self/status");
	std::string s;
	while (getline(proc, s), !proc.fail()) {
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

using namespace boost;

struct vertexProperties
{
	std::string name;
	int weight;
	size_t indexOriginal;
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
printGraph(graph_t& g)
{
	std::cout << "U\tn" << std::endl;
	std::string node1, node2;
	int weight;
	auto vertexItRange = vertices(g);
	for (auto vertexIt = vertexItRange.first; vertexIt != vertexItRange.second; ++vertexIt) {
		node1 = g[*vertexIt].name;
		weight = g[*vertexIt].weight;
		std::cout << node1 << "\t" << weight << "\n";
	}
	std::cout << "\nU\tV\tn" << std::endl;
	auto edgeItRange = edges(g);
	for (auto edgeIt = edgeItRange.first; edgeIt != edgeItRange.second; ++edgeIt) {
		weight = g[*edgeIt].weight;
		node1 = g[boost::source(*edgeIt, g)].name;
		node2 = g[boost::target(*edgeIt, g)].name;
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
	static int benchmark = 0;
	static int help = 0;
	std::string separationStrategy;
	// bool verbose = false;
	unsigned t = 1;
	// bool failed = false;
	static const struct option longopts[] = {
		{ "help", no_argument, &help, 1 },
		{ "benchmark", no_argument, &benchmark, 1 },
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
#if _OPENMP
	double sTime = omp_get_wtime();
#endif
	std::cerr << "Using bgl" << std::endl;
	std::cerr << "Loading Graph" << std::endl;

	barcodeToIndex_t barcodeToIndex;
	indexToBarcode_t indexToBarcode;
	for (auto& infile : infiles) {
		infile == "-" ? "/dev/stdin" : infile;
		vertex_t u;
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
				std::cerr << "Added vertices to graph" << std::endl;
#if _OPENMP
				std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
				sTime = omp_get_wtime();
#endif
				continue;
			}

			std::string node1, node2;
			int weight;
			std::istringstream ss(line);
			ss >> node1;
			if (!atEdges) {
				ss >> weight;
				u = add_vertex(g);
				g[u].name = node1;
				g[u].weight = weight;
				g[u].indexOriginal = u;
				barcodeToIndex[node1] = u;
				indexToBarcode[u] = node1;
			} else {
				ss >> node2 >> weight;
				E = add_edge(barcodeToIndex[node1], barcodeToIndex[node2], g).first;
				g[E].weight = weight;
			}
		}
	}
	std::cerr << "Added edges to graph" << std::endl;
#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif
	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB" << std::endl;
	// barcodeToIndex.clear();
	// boost::print_graph(g);
	using vertexSet_t = std::unordered_set<vertex_t>;
	using componentVertexSet_t = std::vector<vertexSet_t>;
	uint64_t vertexNum = indexToBarcode.size();

	std::vector<std::unordered_map<vertex_t, size_t>> componentsOfVertices;
	componentsOfVertices.resize(indexToBarcode.size());
	// double subgraphTime = 0, biconnectedTime = 0, cleanupTime = 0, componentTime = 0,
	//       insertMapTime = 0;
	if (benchmark != 0) {
		if (vertexNum > 100000) {
			vertexNum = 100000;
		}
		sTime = omp_get_wtime();
	}
#if _OPENMP
	sTime = omp_get_wtime();
#endif
	double totalBiconnectedTime = 0, totalSubgraphTime = 0;
	// auto vertexItRange = vertices(g);
	for (uint64_t vertexId = 0; vertexId < vertexNum; ++vertexId) {
		// for (auto vertexIt = vertexItRange.first; vertexIt != vertexItRange.second; ++vertexIt) {
		// std::cout << "Index: " << *vertexIt << std::endl;
		// auto neighbours = boost::adjacent_vertices(*vertexIt, g);
#if _OPENMP
		double subgraphTime = omp_get_wtime();
#endif
		auto neighbours = boost::adjacent_vertices(vertexId, g);
		/*graph_t& g1 = g.create_subgraph();
		for (auto vd : make_iterator_range(neighbours))
		    boost::add_vertex(vd, g1);
		boost::print_graph(g1);
		std::cout << "Index: " << *vertexIt << std::endl;*/
		graph_t& g1 = g.create_subgraph(neighbours.first, neighbours.second);
#if _OPENMP
		totalSubgraphTime += (omp_get_wtime() - subgraphTime);
#endif
		// boost::print_graph(g1);
#if _OPENMP
		double biconnectedTime = omp_get_wtime();
#endif
		property_map<graph_t, edge_component_t>::type component = get(edge_component, g1);
		biconnected_components(g1, component);
		// std::cerr << "Found " << num_comps << " biconnected components.\n";

		std::vector<vertex_t> art_points_vec;
		articulation_points(g1, std::back_inserter(art_points_vec));
		std::unordered_set<vertex_t> art_points(art_points_vec.begin(), art_points_vec.end());
		// std::cerr << "Found " << art_points.size() << " articulation points.\n";
		// for (auto&& x : art_points)
		//	std::cerr << g1[x].indexOriginal << std::endl;
#if _OPENMP
		totalBiconnectedTime += (omp_get_wtime() - biconnectedTime);
#endif
		graph_traits<graph_t>::edge_iterator ei, ei_end;
		componentVertexSet_t componentVertices;
		for (boost::tie(ei, ei_end) = edges(g1); ei != ei_end; ++ei) {
			size_t componentNum = component[*ei];
			if (componentNum + 1 > componentVertices.size()) {
				componentVertices.resize(componentNum + 1);
			}
			vertex_t node1 = source(*ei, g1);
			vertex_t node2 = target(*ei, g1);
			if (art_points.find(node1) == art_points.end()) {
				// std::cout << g1[node1].indexOriginal << " " << componentNum << std::endl;
				componentVertices[componentNum].insert(g1[node1].indexOriginal);
			}
			if (art_points.find(node2) == art_points.end()) {
				// std::cout << g1[node2].indexOriginal << " " << componentNum << std::endl;
				componentVertices[componentNum].insert(g1[node2].indexOriginal);
			}
		}
		size_t moleculeNum = 0;
		std::unordered_map<vertex_t, size_t> vertexToComponent;
		for (auto&& vertexSet : componentVertices) {
			// std::cout << "size: " << vertexSet.size() << std::endl;
			if (vertexSet.size() <= 1) {
				continue;
			}
			for (auto&& vertex : vertexSet) {
				vertexToComponent[vertex] = moleculeNum;
				// std::cout << vertex << std::endl;
			}
			moleculeNum++;
		}
		// std::cerr << g.m_children.size() << std::endl;
		for (auto i = g.m_children.begin(); i != g.m_children.end(); ++i) {
			delete *i;
		}
		g.m_children.clear();
		// auto edgeItRange = edges(g1);
		// for (auto edgeIt = edgeItRange.first; edgeIt != edgeItRange.second; ++edgeIt) {
		//		remove_edge(*edgeIt, g1);
		//	}

		componentsOfVertices[vertexId] = vertexToComponent;
		// std::cout << "len of component: " << component.size() << std::endl;
	}
	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB" << std::endl;
	std::cerr << "Finish mol sep " << std::endl;
#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif

	if (benchmark != 0) {
		std::cerr << "biconnected time: " << totalBiconnectedTime << std::endl;
		std::cerr << "subgraph time: " << totalSubgraphTime << std::endl;
		exit(0);
	}
	// printGraph(g);
	graph_t outG;
	barcodeToIndex_t outGBarcodeToIndex;
	for (size_t i = 0; i < componentsOfVertices.size(); i++) {
		/*std::pair<vertex_t, size_t> result;
		result = *std::max_element(
		    componentsOfVertices[i].begin(),
		    componentsOfVertices[i].end(),
		    [](const std::pair<vertex_t, size_t>& p1, const std::pair<vertex_t, size_t>& p2) {
		        return p1.second < p2.second;
		    });*/

		size_t maxVal = 0;
		for (auto&& val : componentsOfVertices[i]) {
			if (val.second > maxVal) {
				maxVal = val.second;
			}
		}
		for (size_t j = 0; j < maxVal + 1; j++) {

			vertex_t u = add_vertex(outG);
			outG[u].name = g[i].name + "_" + std::to_string(j);
			outG[u].weight = g[i].weight;
			outG[u].indexOriginal = u;
			outGBarcodeToIndex[outG[u].name] = u;
		}
	}
	auto edgeItRange = edges(g);
	for (auto edgeIt = edgeItRange.first; edgeIt != edgeItRange.second; ++edgeIt) {
		vertex_t u, v;
		u = g[boost::source(*edgeIt, g)].indexOriginal;
		v = g[boost::target(*edgeIt, g)].indexOriginal;
		if (componentsOfVertices[u].find(v) == componentsOfVertices[u].end() ||
		    componentsOfVertices[v].find(u) == componentsOfVertices[v].end()) {
			continue;
		}
		size_t uMolecule = componentsOfVertices[u][v];
		size_t vMolecule = componentsOfVertices[v][u];
		std::string uName = g[u].name + "_" + std::to_string(uMolecule);
		std::string vName = g[v].name + "_" + std::to_string(vMolecule);
		edge_t e = add_edge(outGBarcodeToIndex[uName], outGBarcodeToIndex[vName], outG).first;
		outG[e].weight = g[*edgeIt].weight;
	}
	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB" << std::endl;
	printGraph(outG);
}
