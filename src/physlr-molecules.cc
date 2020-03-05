#include "tsl/robin_map.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
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
	uint64_t indexOriginal = 0;
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

// this definition assumes there is no redundant edge in the undirected graph
using graph_t = boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    vertexProperties,
    boost::property<
        boost::edge_index_t,
        int,
        boost::property<edgeComponent_t, std::uint64_t, edgeProperties>>>;
using vertex_t = graph_t::vertex_descriptor;
using edge_t = graph_t::edge_descriptor;
using barcodeToIndex_t = std::unordered_map<std::string, vertex_t>;
using indexToBarcode_t = std::unordered_map<vertex_t, std::string>;
using vertexSet_t = std::unordered_set<vertex_t>;
using componentToVertexSet_t = std::vector<vertexSet_t>;
using vertexToComponent_t = std::unordered_map<vertex_t, uint64_t>;
using vecVertexToComponent_t = std::vector<vertexToComponent_t>;
using vertexToIndex_t =
    std::unordered_map<vertex_t, uint64_t>; // wanna improve this? checkout boost::bimap
using indexToVertex_t =
    std::unordered_map<uint64_t, vertex_t>; // wanna improve this? checkout boost::bimap
using adjacencyMatrix_t = std::vector<std::vector<uint_fast32_t>>;
using adjacencyVector_t = std::vector<uint_fast32_t>;
using Clique_type = std::unordered_map<vertex_t, uint64_t>;

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
	          << "  [-s SEPARATION-STRATEGY] [-v] FILE...\n\n"
	             "  -v         enable verbose output\n"
	             "  -s --separation-strategy   \n"
	             "  SEPARATION-STRATEGY      `+` separated list of molecule separation strategies "
	             "[bc]\n"
	             "  --help     display this help and exit\n";
}

void
printGraph(const graph_t& g)
{
	std::cout << "U\tm" << std::endl;
	auto vertexItRange = boost::vertices(g);
	for (auto vertexIt = vertexItRange.first; vertexIt != vertexItRange.second; ++vertexIt) {
		auto& node1 = g[*vertexIt].name;
		auto& weight = g[*vertexIt].weight;
		std::cout << node1 << "\t" << weight << "\n";
	}
	std::cout << "\nU\tV\tm" << std::endl;
	auto edgeItRange = boost::edges(g);
	for (auto edgeIt = edgeItRange.first; edgeIt != edgeItRange.second; ++edgeIt) {
		auto& weight = g[*edgeIt].weight;
		auto& node1 = g[boost::source(*edgeIt, g)].name;
		auto& node2 = g[boost::target(*edgeIt, g)].name;
		std::cout << node1 << "\t" << node2 << "\t" << weight << "\n";
	}
}

void
readTSV(graph_t& g, const std::vector<std::string>& infiles, bool verbose)
{
	auto progname = "physlr-molecules";
	std::cerr << "Loading graph" << std::endl;
#if _OPENMP
	double sTime = omp_get_wtime();
#endif
	barcodeToIndex_t barcodeToIndex;
	indexToBarcode_t indexToBarcode;
	for (auto& infile : infiles) {
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
				printErrorMsg(progname, "unknown graph format");
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
				printErrorMsg(progname, "unknown graph format");
				exit(EXIT_FAILURE);
			}
			std::string node1, node2;
			int weight;
			std::istringstream ss(line);
			if (ss >> node1 >> node2 >> weight) {
				auto E = boost::add_edge(barcodeToIndex[node1], barcodeToIndex[node2], g).first;
				g[E].weight = weight;
			} else {
				printErrorMsg(progname, "unknown graph format");
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
		std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB"
		          << std::endl;
	}
}

/* Generate a molecule separated graph (molSepG) using component/community information from
molecule separation (vecVertexToComponent). The input graph (inG) is the barcode overlap graph
or a molecule separated graph from the previous round of molecule separation.*/
void
componentsToNewGraph(
    const graph_t& inG,
    graph_t& molSepG,
    vecVertexToComponent_t& vecVertexToComponent)
{
	barcodeToIndex_t molSepGBarcodeToIndex;
#if _OPENMP
	double sTime = omp_get_wtime();
#endif
	for (uint64_t i = 0; i < vecVertexToComponent.size(); ++i) {

		uint64_t maxVal = 0;
		if (!vecVertexToComponent[i].empty()) {
			maxVal =
			    std::max_element(
			        vecVertexToComponent[i].begin(),
			        vecVertexToComponent[i].end(),
			        [](const vertexToComponent_t::value_type& p1,
			           const vertexToComponent_t::value_type& p2) { return p1.second < p2.second; })
			        ->second;
		}

		for (uint64_t j = 0; j < maxVal + 1; ++j) {
			auto u = boost::add_vertex(molSepG);
			molSepG[u].name = inG[i].name + "_" + std::to_string(j);
			molSepG[u].weight = inG[i].weight;
			molSepG[u].indexOriginal = u;
			molSepGBarcodeToIndex[molSepG[u].name] = u;
		}
	}

	auto edgeItRange = boost::edges(inG);
	for (auto edgeIt = edgeItRange.first; edgeIt != edgeItRange.second; ++edgeIt) {
		auto& u = inG[boost::source(*edgeIt, inG)].indexOriginal;
		auto& v = inG[boost::target(*edgeIt, inG)].indexOriginal;

		if (vecVertexToComponent[u].find(v) == vecVertexToComponent[u].end() ||
		    vecVertexToComponent[v].find(u) == vecVertexToComponent[v].end()) {
			continue;
		}

		auto& uMolecule = vecVertexToComponent[u][v];
		auto& vMolecule = vecVertexToComponent[v][u];
		auto uName = inG[u].name + "_" + std::to_string(uMolecule);
		auto vName = inG[v].name + "_" + std::to_string(vMolecule);
		auto e =
		    boost::add_edge(molSepGBarcodeToIndex[uName], molSepGBarcodeToIndex[vName], molSepG)
		        .first;
		molSepG[e].weight = inG[*edgeIt].weight;
	}

	std::cerr << "Generated new graph ";
#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif

	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB" << std::endl;
}

uint64_t
biconnectedComponents(
    graph_t& subgraph,
    vertexToComponent_t& vertexToComponent,
    uint64_t initial_community_id = 0)
{
	// Find biconnected components
	boost::property_map<graph_t, edgeComponent_t>::type component =
	    boost::get(edgeComponent, subgraph);

	std::vector<vertex_t> artPointsVec;
	boost::biconnected_components(subgraph, component, std::back_inserter(artPointsVec));

	vertexSet_t artPoints(artPointsVec.begin(), artPointsVec.end());

	// Remove articulation points from biconnected components
	boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
	componentToVertexSet_t componentToVertexSet;

	for (boost::tie(ei, ei_end) = boost::edges(subgraph); ei != ei_end; ++ei) {
		uint64_t componentNum = component[*ei];
		if (componentNum + 1 > componentToVertexSet.size()) {
			componentToVertexSet.resize(componentNum + 1);
		}

		auto node1 = source(*ei, subgraph);
		auto node2 = target(*ei, subgraph);

		if (artPoints.find(node1) == artPoints.end()) {
			componentToVertexSet[componentNum].insert(subgraph[node1].indexOriginal);
		}
		if (artPoints.find(node2) == artPoints.end()) {
			componentToVertexSet[componentNum].insert(subgraph[node2].indexOriginal);
		}
	}

	uint64_t moleculeNum = initial_community_id;

	// Remove components with size less than 1
	for (auto&& vertexSet : componentToVertexSet) {
		if (vertexSet.size() <= 1) {
			continue;
		}
		for (auto&& vertex : vertexSet) {
			vertexToComponent[vertex] = moleculeNum;
		}
		++moleculeNum;
	}
	return moleculeNum;
}

void
bin_components(
    componentToVertexSet_t& source,
    componentToVertexSet_t& binned_neighbours,
    uint64_t bin_size = 50)
{
	// //   Iterate over each component and if its bigger than bin_size:
	// //   randomly split the component (set of vertices) into smaller even bins

	std::vector<uint64_t> components_size;
	uint64_t neighborhood_size;
	uint64_t components_count;
	for (uint64_t i = 0; i < source.size(); i++) { // NOLINT
		neighborhood_size = source[i].size();
		components_count = ((neighborhood_size - 1) / bin_size) + 1;
		components_size.push_back(components_count);
	}
	uint64_t new_size =
	    std::accumulate(components_size.begin(), components_size.end(), uint64_t(0));
	binned_neighbours.resize(new_size);
	uint64_t counter_new = 0;
	uint64_t base_com_size;
	uint64_t leftover;

	for (uint64_t i = 0; i < source.size(); i++) { // NOLINT
		// Using unordered_set, we make use of its random nature and we do not shuffle randomly
		base_com_size = source[i].size() / components_size[i];
		leftover = source[i].size() % components_size[i];
		uint64_t yet_leftover = (leftover ? 1 : 0);

		auto elementIt = source[i].begin();
		while (elementIt != source[i].end()) {
			uint64_t length = base_com_size + yet_leftover;
			if (--leftover == 0) {
				yet_leftover = 0;
			}

			for (uint64_t j = 0; j < length; j++) {
				binned_neighbours[counter_new].insert(*elementIt);
				++elementIt;
			}
			counter_new++;
		}
	}
}

template<class Neighbours_Type>
void
bin_neighbours(
    Neighbours_Type neighbours,
    componentToVertexSet_t& binned_neighbours,
    uint64_t bin_size = 50)
{
	// //   Randomly split the set of vertices (neighbours) into bins

	componentToVertexSet_t compToVertset(1, vertexSet_t(neighbours.first, neighbours.second));
	if (compToVertset[0].size() > bin_size) {
		bin_components(compToVertset, binned_neighbours, bin_size);
	} else {
		binned_neighbours = compToVertset;
	}
}

template<class Graph, class vertexIter, class edgeSet>
void
make_subgraph(Graph& g, Graph& subgraph, edgeSet& edge_set, vertexIter vBegin, vertexIter vEnd)
{
	// //   Make a vertex-induced subgraph of graph g, based on vertices from vBegin to vEnd
	// //   track the source node by indexOriginal

	// Add vertices into the subgraph, but set `indexOriginal` the index of it in the source graph.
	for (auto& vIter = vBegin; vIter != vEnd; ++vIter) {
		auto u = boost::add_vertex(subgraph);

		subgraph[u].name = g[*vIter].name;
		subgraph[u].weight = g[*vIter].weight;
		subgraph[u].indexOriginal = g[*vIter].indexOriginal;
	}

	// Iterate over all pairs of vertices in the subgraph:
	// check whether there exist an edge between their corresponding vertices in the source graph.
	graph_t::vertex_iterator vIter1, vIter2, vend1, vend2;
	for (boost::tie(vIter1, vend1) = vertices(subgraph); vIter1 != vend1; ++vIter1) {
		for (boost::tie(vIter2, vend2) = vertices(subgraph); vIter2 != vend2; ++vIter2) {
			if (vIter1 != vIter2) {
				auto got = edge_set.find(std::make_pair(
				    subgraph[*vIter1].indexOriginal, subgraph[*vIter2].indexOriginal));

				if (got != edge_set.end()) {
					auto new_edge = boost::add_edge(*vIter1, *vIter2, subgraph).first;
					subgraph[new_edge].weight = got->second;
				}
			}
		}
	}
}

adjacencyMatrix_t
convert_adj_list_adj_mat(graph_t& subgraph, vertexToIndex_t& vertexToIndex)
{
	// Inputs:
	// - subgraph: adjacency list to convert to adjacency list
	// - vertexToIndex: (empty, to be filled in)
	//      Dictionary of (vertex name) -> (index in temporary adjacency matrix)
	// Ouput(s):
	// - adj_mat: the adjacency matrix for subgraph
	// - vertexToIndex (referenced input)

	int N = boost::num_vertices(subgraph);
	adjacencyVector_t tempVector(N, 0);
	adjacencyMatrix_t adj_mat(N, tempVector);

	typedef boost::graph_traits<graph_t>::edge_iterator edge_iterator;

	pair<edge_iterator, edge_iterator> ei = edges(subgraph);

	vertexToIndex_t::iterator got_a;
	vertexToIndex_t::iterator got_b;
	uint64_t adj_mat_index = 0;
	for (edge_iterator edge_iter = ei.first; edge_iter != ei.second; ++edge_iter) {
		vertex_t a = source(*edge_iter, subgraph);
		vertex_t b = target(*edge_iter, subgraph);
		// if not visited a or b
		//      add to dictionary
		// Could be more efficient by adding a "visited" property to vertices of the graph
		// Now we implement by hash table lookup:

		got_a = vertexToIndex.find(a);
		uint64_t index_a;
		if (got_a == vertexToIndex.end()) {
			vertexToIndex.insert(std::pair<vertex_t, uint64_t>(a, adj_mat_index));
			index_a = adj_mat_index++;
		} else {
			index_a = got_a->second;
		}

		got_b = vertexToIndex.find(b);
		uint64_t index_b;
		if (got_b == vertexToIndex.end()) {
			vertexToIndex.insert(std::pair<vertex_t, uint64_t>(b, adj_mat_index));
			index_b = adj_mat_index++;
		} else {
			index_b = got_b->second;
		}

		adj_mat[index_a][index_b] = (int)subgraph[*edge_iter].weight;
		adj_mat[index_b][index_a] = adj_mat[index_a][index_b];
	}
	return adj_mat;
}




uint64_t
community_detection_cosine_similarity(
    graph_t& subgraph,
    vertexToComponent_t& vertexToComponent,
    uint64_t initial_community_id = 0,
    bool squaring = true,
    double threshold=0.3)
{
    // Detect communities using cosine similarity of vertices

    // 0- Map indices and vertex names

    vertexToIndex_t vertexToIndex;
    uint64_t subgraph_size = boost::num_vertices(subgraph);
    vertexToIndex.reserve(subgraph_size);

    if (subgraph_size < 10)
        // Do nothing on subgraphs smaller than a certain size
        threshold = 0;

    adjacencyMatrix_t adj_mat(convert_adj_list_adj_mat(subgraph, vertexToIndex));
    indexToVertex_t indexToVertex = inverse_map(vertexToIndex);

    // 1- Calculate the cosine similarity:

    int size_adj_mat = adj_mat.size();
    vector<double> tempVector(size_adj_mat, 0);
    vector<vector<double>> cosSimilarity2d(size_adj_mat, tempVector);

    if (squaring){
        calculate_cosine_similarity_2d(square_matrix_ikj(adj_mat, true),
                                        // may need some change
                                        //square_matrix_ijk(adj_mat, true),
                                        //square_matrix_boost(adj_mat),
                                        cosSimilarity2d);
    }else{
        calculate_cosine_similarity_2d(adj_mat, cosSimilarity2d);
        }

    // 2- Determine the threshold:
    // not implemented yet; so use a predefined universal threshold.

    threshold = threshold;

    // 3- Filter out edges:

    auto start3 = timeNow();
    for (int i = 0; i < adj_mat.size() ; i++)
    {
        for (int j = i+1; j < adj_mat.size() ; j++)
            {
                if (cosSimilarity2d[i][j] < threshold)
                {
                    adj_mat[i][j] = 0;
                    adj_mat[j][i] = 0;
                }
            }
    }

    // 4- Detect Communities (find connected components - DFS)
    //      Alternative implementation: convert to adjacency list and use boost to find cc

    uint64_t community_id = initial_community_id;
    stack<uint64_t> toCheck;
    stack<uint64_t> toAdd;
    vector<int> zeros(adj_mat.size(),0);
    vector<int> isDetected(adj_mat.size(),0);
    bool isSingleton = false;
    for (uint64_t i = 0 ; i < adj_mat.size(); i++)
    {
        // DFS traversal
        if (isDetected[i])
            continue; // this node is included in a community already.
        toCheck.push(i);
        isDetected[i] = 1;
        isSingleton = true;
        uint64_t ii;
        uint64_t node_to_add;

        while(!toCheck.empty()){

            ii = toCheck.top();
            toCheck.pop();
            toAdd.push(ii);
            for (uint64_t j = 0 ; j < adj_mat.size(); j++)
            {
                if (isDetected[j])
                    continue; // this node is included in a community already.
                if (adj_mat[ii][j] > 0){
                    toCheck.push(j);
                    isDetected[j]=1;
                    isSingleton = false;
                }
            }
        }
        if (toAdd.size() < 2)
        {
            while(!toAdd.empty())
                toAdd.pop();
        }
        else
        {
            while(!toAdd.empty())
            {
                node_to_add = toAdd.top();
                toAdd.pop();
                auto vt = indexToVertex.find(node_to_add);
                if (vt != indexToVertex.end() ){
                    vertexToComponent[subgraph[vt->second].indexOriginal] = community_id;
                }
                else{
                    std::cerr<<"BUG: not found in the hash table!"<<std::endl;
                    }
            }
            community_id++;
        }
    }

    return community_id;
}

int
main(int argc, char* argv[])
{
	auto progname = "physlr-molecules";
	int optindex = 0;
	static int help = 0;
	std::string separationStrategy = "bc";
	uint64_t threads = 1;
	bool verbose = false;
	bool failed = false;
	static const struct option longopts[] = {
		{ "help", no_argument, &help, 1 },
		{ "separation-strategy", required_argument, nullptr, 's' },
		{ nullptr, 0, nullptr, 0 }
	};

	for (int c; (c = getopt_long(argc, argv, "s:vt:", longopts, &optindex)) != -1;) {
		switch (c) {
		case 0:
			break;
		case 's':
			separationStrategy.assign(optarg);
			break;
		case 't':
			threads = std::stoi(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		default:
			exit(EXIT_FAILURE);
		}
	}

	std::cerr << " using " << threads << " thread(s)." << std::endl;
#if _OPENMP
	omp_set_num_threads(threads);
#endif

	std::vector<std::string> infiles(&argv[optind], &argv[argc]);
	if (argc < 1) {
		failed = true;
	}
	if (help != 0) {
		printVersion();
		printUsage(progname);
		exit(EXIT_SUCCESS);
	} else if (infiles.empty()) {
		printErrorMsg(progname, "missing file operand");
		failed = true;
	}
	if (separationStrategy != "bc") {
		printErrorMsg(progname, "unsupported molecule separation strategy");
		failed = true;
	}
	if (failed) {
		printUsage(progname);
		exit(EXIT_FAILURE);
	}

	graph_t g;
	readTSV(g, infiles, verbose);

	barcodeToIndex_t barcodeToIndex;
	indexToBarcode_t indexToBarcode;

	vecVertexToComponent_t vecVertexToComponent;
	vecVertexToComponent.resize(boost::num_vertices(g));

#if _OPENMP
	double sTime = omp_get_wtime();
#endif

	uint64_t initial_community_id = 0;

	// // auxiliary data-set: set of edges for faster lookup
	tsl::robin_map<
	    std::pair<std::uint64_t, uint64_t>,
	    int,
	    boost::hash<std::pair<uint64_t, uint64_t>>>
	    edge_set;
	edge_set.reserve(num_edges(g));

	auto edgeItRange = boost::edges(g);
	for (auto edgeIt = edgeItRange.first; edgeIt != edgeItRange.second; ++edgeIt) {
		auto& weight = g[*edgeIt].weight;
		auto& node1 = g[boost::source(*edgeIt, g)].indexOriginal;
		auto& node2 = g[boost::target(*edgeIt, g)].indexOriginal;
		edge_set[std::pair<uint64_t, uint64_t>(node1, node2)] = weight;
	}

	auto vertexItRange = vertices(g);

#if !_OPENMP
	if (threads > 1) {
		threads = 1;
		std::cerr << "Setting threads to 1." << std::endl;
	}
#endif

	if (threads > 1) {
		const uint64_t array_size = boost::num_vertices(g);
		std::vector<boost::graph_traits<graph_t>::vertex_iterator> iterators_array;
		iterators_array.resize(array_size);
		boost::graph_traits<graph_t>::vertex_iterator allocate_it = vertexItRange.first;
		for (uint64_t j = 0; j < array_size; ++j) {
			iterators_array[j] = allocate_it++;
		}

#pragma omp parallel for
		for (uint64_t j = 0; j < array_size; ++j) {
			initial_community_id = 0;
			componentToVertexSet_t componentsVec;
			vertexToComponent_t vertexToComponent;

			// Find neighbour of vertex and generate neighbour induced subgraph
			auto neighbours = boost::adjacent_vertices(*(iterators_array[j]), g);

			// binning
			bin_neighbours(neighbours, componentsVec, 50);

			// for (uint64_t comp_i = 0; comp_i < componentsVec.size(); comp_i++) {
			for (auto& comp_i : componentsVec) {
				graph_t subgraph;
				make_subgraph(g, subgraph, edge_set, comp_i.begin(), comp_i.end());

				initial_community_id =
				    biconnectedComponents(subgraph, vertexToComponent, initial_community_id);
			}
			vecVertexToComponent[*(iterators_array[j])] = vertexToComponent;
		}
	} else {
		for (auto vertexIt = vertexItRange.first; vertexIt != vertexItRange.second; ++vertexIt) {
			initial_community_id = 0;
			componentToVertexSet_t componentsVec;
			vertexToComponent_t vertexToComponent;

			// Find neighbour of vertex and generate neighbour induced subgraph
			auto neighbours = boost::adjacent_vertices(*vertexIt, g);

			// binning
			bin_neighbours(neighbours, componentsVec, 50);

			// for (uint64_t comp_i = 0; comp_i < componentsVec.size(); comp_i++) {
			for (auto& comp_i : componentsVec) {
				graph_t subgraph;
				make_subgraph(g, subgraph, edge_set, comp_i.begin(), comp_i.end());

				initial_community_id =
				    biconnectedComponents(subgraph, vertexToComponent, initial_community_id);
			}
			vecVertexToComponent[*vertexIt] = vertexToComponent;
		}
	}

	std::cerr << "Finished molecule separation ";
#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif
	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB" << std::endl;

	std::cerr << "Generating molecule overlap graph" << std::endl;

	graph_t molSepG;
	componentsToNewGraph(g, molSepG, vecVertexToComponent);
	printGraph(molSepG);
	if (verbose) {
		std::cerr << "Printed graph" << std::endl;
#if _OPENMP
		std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
#endif
	}
}
