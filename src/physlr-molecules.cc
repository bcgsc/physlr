#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
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

struct edge_component_t
{
	enum
	{
		num = 555
	};
	using kind = boost::edge_property_tag;
} edge_component;

using graph_t = boost::subgraph<boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    vertexProperties,
    boost::property<
        boost::edge_index_t,
        int,
        boost::property<edge_component_t, std::size_t, edgeProperties>>>>;
using vertex_t = graph_t::vertex_descriptor;
using edge_t = graph_t::edge_descriptor;
using barcodeToIndex_t = std::unordered_map<std::string, vertex_t>;
using indexToBarcode_t = std::unordered_map<vertex_t, std::string>;
using vertexSet_t = std::unordered_set<vertex_t>;
using componentToVertexSet_t = std::vector<vertexSet_t>;
using vertexToComponent_t = std::unordered_map<vertex_t, size_t>;
using vecVertexToComponent_t = std::vector<vertexToComponent_t>;

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
printGraph(const graph_t& g)
{
	std::cout << "U\tn" << std::endl;
	std::string node1, node2;
	int weight;
	auto vertexItRange = boost::vertices(g);
	for (auto vertexIt = vertexItRange.first; vertexIt != vertexItRange.second; ++vertexIt) {
		node1 = g[*vertexIt].name;
		weight = g[*vertexIt].weight;
		std::cout << node1 << "\t" << weight << "\n";
	}
	std::cout << "\nU\tV\tn" << std::endl;
	auto edgeItRange = boost::edges(g);
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
	std::string separationStrategy = "bc";
	bool verbose = false;
	unsigned t = 1;
	bool failed = false;
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
			verbose = true;
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
		failed = true;
	}
	if (t > 1) {
		std::cerr << "physlr-molecules does not support multithreading currently." << std::endl;
		t = 1;
	}
	if (separationStrategy.compare("bc")) {
		std::cerr << "physlr-molecules only supports biconnected components currently."
		          << std::endl;
		separationStrategy = "bc";
	}
	// TODO add checks for separation strategy
	if (failed) {
		printUsage(progname);
		exit(EXIT_FAILURE);
	}

	graph_t g;
	std::cerr << "Loading graph" << std::endl;

#if _OPENMP
	double sTime = omp_get_wtime();
#endif

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
			if (line.empty() or line == "U\tn") {
				continue;
			}
			if (line.empty() or line == "U\tV\tn") {
				atEdges = true;
				if (verbose) {
					std::cerr << "Loaded vertices to graph ";
#if _OPENMP
					std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
					sTime = omp_get_wtime();
#endif
				}
				continue;
			}

			std::string node1, node2;
			int weight;
			std::istringstream ss(line);
			ss >> node1;
			if (!atEdges) {
				ss >> weight;
				u = boost::add_vertex(g);
				g[u].name = node1;
				g[u].weight = weight;
				g[u].indexOriginal = u;
				barcodeToIndex[node1] = u;
				indexToBarcode[u] = node1;
			} else {
				ss >> node2 >> weight;
				E = boost::add_edge(barcodeToIndex[node1], barcodeToIndex[node2], g).first;
				g[E].weight = weight;
			}
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

	uint64_t vertexNum = indexToBarcode.size();

	vecVertexToComponent_t vecVertexToComponent;
	vecVertexToComponent.resize(vertexNum);

#if _OPENMP
	sTime = omp_get_wtime();
#endif

	for (uint64_t vertexId = 0; vertexId < vertexNum; ++vertexId) {
		// Find neighbour of vertex and generate neighbour induced subgraph
		auto neighbours = boost::adjacent_vertices(vertexId, g);
		graph_t& subgraph = g.create_subgraph(neighbours.first, neighbours.second);

		// Find biconnected components
		boost::property_map<graph_t, edge_component_t>::type component =
		    boost::get(edge_component, subgraph);
		boost::biconnected_components(subgraph, component);

		std::vector<vertex_t> art_points_vec;
		articulation_points(subgraph, std::back_inserter(art_points_vec));
		std::unordered_set<vertex_t> art_points(art_points_vec.begin(), art_points_vec.end());

		// Remove articulation points from biconnected components
		boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
		componentToVertexSet_t componentToVertexSet;

		for (boost::tie(ei, ei_end) = boost::edges(subgraph); ei != ei_end; ++ei) {
			size_t componentNum = component[*ei];
			if (componentNum + 1 > componentToVertexSet.size()) {
				componentToVertexSet.resize(componentNum + 1);
			}

			vertex_t node1 = source(*ei, subgraph);
			vertex_t node2 = target(*ei, subgraph);

			if (art_points.find(node1) == art_points.end()) {
				componentToVertexSet[componentNum].insert(subgraph[node1].indexOriginal);
			}
			if (art_points.find(node2) == art_points.end()) {
				componentToVertexSet[componentNum].insert(subgraph[node2].indexOriginal);
			}
		}

		size_t moleculeNum = 0;
		vertexToComponent_t vertexToComponent;

		// Remove components with size less than 1
		for (auto&& vertexSet : componentToVertexSet) {
			if (vertexSet.size() <= 1) {
				continue;
			}
			for (auto&& vertex : vertexSet) {
				vertexToComponent[vertex] = moleculeNum;
			}
			moleculeNum++;
		}

		// Delete subgraph to keep memory in control
		for (auto i = g.m_children.begin(); i != g.m_children.end(); ++i) {
			delete *i;
		}
		g.m_children.clear();

		vecVertexToComponent[vertexId] = vertexToComponent;
	}

	std::cerr << "Finished molecule separation ";
#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif
	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB" << std::endl;

	std::cerr << "Generating new graph" << std::endl;

	graph_t outG;
	barcodeToIndex_t outGBarcodeToIndex;
	for (size_t i = 0; i < vecVertexToComponent.size(); i++) {

		size_t maxVal = 0;
		for (auto&& val : vecVertexToComponent[i]) {
			if (val.second > maxVal) {
				maxVal = val.second;
			}
		}

		for (size_t j = 0; j < maxVal + 1; j++) {
			vertex_t u = boost::add_vertex(outG);
			outG[u].name = g[i].name + "_" + std::to_string(j);
			outG[u].weight = g[i].weight;
			outG[u].indexOriginal = u;
			outGBarcodeToIndex[outG[u].name] = u;
		}
	}

	auto edgeItRange = boost::edges(g);
	for (auto edgeIt = edgeItRange.first; edgeIt != edgeItRange.second; ++edgeIt) {
		vertex_t u, v;
		u = g[boost::source(*edgeIt, g)].indexOriginal;
		v = g[boost::target(*edgeIt, g)].indexOriginal;

		if (vecVertexToComponent[u].find(v) == vecVertexToComponent[u].end() ||
		    vecVertexToComponent[v].find(u) == vecVertexToComponent[v].end()) {
			continue;
		}

		size_t uMolecule = vecVertexToComponent[u][v];
		size_t vMolecule = vecVertexToComponent[v][u];
		std::string uName = g[u].name + "_" + std::to_string(uMolecule);
		std::string vName = g[v].name + "_" + std::to_string(vMolecule);
		edge_t e =
		    boost::add_edge(outGBarcodeToIndex[uName], outGBarcodeToIndex[vName], outG).first;
		outG[e].weight = g[*edgeIt].weight;
	}

	std::cerr << "Generated new graph ";
#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif

	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB" << std::endl;

	printGraph(outG);
	if (verbose) {
		std::cerr << "Printed graph" << std::endl;
#if _OPENMP
		std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
#endif
	}
}
