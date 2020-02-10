#include "tsl/robin_map.h"
#include "tsl/robin_set.h"

#include <cfenv>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/functional/hash.hpp>

#if _OPENMP
#include <omp.h>
#endif

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt {
static unsigned threads = 1;
static unsigned scoreThreshold = 10;
static unsigned mapPositions = 10;
} // namespace opt

uint64_t max = std::numeric_limits<std::uint64_t>::max();

#define PROGRAM "physlr-map"
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

using BarcodeID = uint32_t;
using BarcodeID = uint32_t;
using Minimizer = uint64_t;
using pair = std::pair<uint64_t, uint64_t>;

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
	    "Usage: physlr-map [OPTION]... [TargetPATHS.path] [TargetMINIMIZERS.tsv] "
	    "[QueryMINIMIZERS.tsv]\n"
	    "Map sequences to a physical map.\n"
	    "  -t, --threads     threads [1]\n"
	    "  -n, --scoreThreshold     minimum score to map [10]\n"
	    "  -m, --mapPositions     number of positions used to orient [10]\n"
	    "  -v         enable verbose output\n"
	    "  --version     Print version\n"
	    "  --help     display this help and exit\n"
	    "Report bugs to <jowong@bcgsc.ca>.";
	std::cerr << dialog << std::endl;
	exit(EXIT_SUCCESS);
}

static void
printErrorMsg(const std::string& progname, const std::string& msg)
{
	std::cerr << progname << ": " << msg << "\nTry '" << progname
	          << " --help' for more information.\n";
}

uint64_t
lowerMedian(std::vector<uint64_t> scores)
{
	uint64_t size = scores.size();

	if (size == 0) {
		return 0;
	}
	sort(scores.begin(), scores.end());
	if (size % 2 == 0) {
		return (scores[size / 2 - 1]);
	}
	return scores[size / 2];
}

uint64_t
upperMedian(std::vector<uint64_t> scores)
{
	uint64_t size = scores.size();

	if (size == 0) {
		return 0;
	}
	sort(scores.begin(), scores.end());
	return scores[size / 2];
}

std::string
determineOrientation(uint64_t prev, uint64_t curr, uint64_t next)
{

	if (prev != max && next != max) {
		if (prev == curr && curr == next) {
			return ".";
		}
		if (prev <= curr && curr <= next) {
			return "+";
		}
		if (prev >= curr && curr >= next) {
			return "-";
		}
		return ".";
	}

	if (prev != max) {
		if (prev < curr) {
			return "+";
		}
		if (prev > curr) {
			return "-";
		}
		return ".";
	}

	if (next != max) {
		if (curr < next) {
			return "+";
		}
		if (curr > next) {
			return "-";
		}
		return ".";
	}

	return ".";
}

void
readPaths(std::vector<std::vector<std::string>>& paths, const std::string& pathFile)
{
#if _OPENMP
	double sTime = omp_get_wtime();
#endif
	std::ifstream fh;
	std::cerr << "Loading file " << pathFile << std::endl;
	fh.open(pathFile);
	if (!fh) {
		std::cerr << "Invalid file: " << pathFile << std::endl;
		exit(EXIT_FAILURE);
	}
	std::string line;
	while (getline(fh, line)) {
		std::stringstream ss(line);
		std::string molecule;
		paths.emplace_back(std::vector<std::string>());
		while (ss >> molecule) {
			paths[paths.size() - 1].emplace_back(molecule);
		}
	}
	std::cerr << "Loaded" << std::endl;
#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif
}

void
getMoleculeToMinimizer(
    tsl::robin_map<std::string, std::vector<Minimizer>>& moleculeToMinimizer,
    std::string& inputFile)
{
#if _OPENMP
	double sTime = omp_get_wtime();
#endif
	std::string molecule;
	std::string minimizerAndLoc;
	Minimizer minimizer;

	// read in minimizer file
	// format: GAGGTCCGTGGAGAGG-1	472493953667297251 1168973555595507959 342455687043295195
	// 283275954102976652

	std::ifstream fh;
	std::cerr << "Loading file " << inputFile << std::endl;
	fh.open(inputFile);
	if (!fh) {
		std::cerr << "Invalid file: " << pathFile << std::endl;
		exit(EXIT_FAILURE);
	}
	std::string line;
	while (getline(fh, line)) {
		std::stringstream ss(line);
		ss >> molecule;

		moleculeToMinimizer[molecule] = std::vector<Minimizer>();
		while (ss >> minimizerAndLoc) {
			uint64_t colonLoc = minimizerAndLoc.find(':');
			if (colonLoc == std::string::npos) {
				minimizer = std::stoull(minimizerAndLoc);
				moleculeToMinimizer[molecule].emplace_back(minimizer);
			} else {
				minimizer = std::stoull(minimizerAndLoc.substr(0, colonLoc));
				moleculeToMinimizer[molecule].emplace_back(minimizer);
			}
		}
	}
	std::cerr << "Loaded" << std::endl;
#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif
}

void
getMinimizerToPos(
    const std::vector<std::vector<std::string>>& paths,
    const tsl::robin_map<std::string, std::vector<Minimizer>>& moleculeToMinimizer,
    tsl::robin_map<Minimizer, tsl::robin_set<pair, boost::hash<pair>>>& minimizerToPos)
{
	std::cerr << "Mapping Minimizers to positions" << std::endl;
#if _OPENMP
	double sTime = omp_get_wtime();
#endif
	for (uint64_t targetId = 0; targetId < paths.size(); ++targetId) {
		auto& path = paths[targetId];
		for (uint64_t pos = 0; pos < path.size(); ++pos) {
			const auto& molecule = path[pos];
			const auto& minimizers = moleculeToMinimizer.at(molecule);
			for (const auto& minimizer : minimizers) {
				if (minimizerToPos.find(minimizer) == minimizerToPos.end()) {
					minimizerToPos[minimizer] = tsl::robin_set<pair, boost::hash<pair>>();
					minimizerToPos[minimizer].insert(std::make_pair(targetId, pos));
				} else {
					minimizerToPos[minimizer].insert(std::make_pair(targetId, pos));
				}
			}
		}
	}
	std::cerr << "Mapped" << std::endl;
#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif
}

void
mapQueryToTarget(
    const tsl::robin_map<std::string, std::vector<Minimizer>>& queryToMinimizer,
    const tsl::robin_map<Minimizer, tsl::robin_set<pair, boost::hash<pair>>>& minimizerToPos)
{

	std::vector<std::string> queryToMinimizerkeys;
	queryToMinimizerkeys.reserve(queryToMinimizer.size());

	for (const auto& queryToMinimizerKeyVal : queryToMinimizer) {
		queryToMinimizerkeys.push_back(queryToMinimizerKeyVal.first);
	}

	unsigned num_mapped = 0;
#if _OPENMP
#pragma omp parallel for
#endif
	for (uint64_t i = 0; i < queryToMinimizerkeys.size(); ++i) { // NOLINT(modernize-loop-convert)
		const auto& queryId = queryToMinimizerkeys.at(i);
		const auto& minimizers = queryToMinimizer.at(queryId);
		tsl::robin_map<pair, std::vector<uint64_t>, boost::hash<pair>> targetIdPosToQuerypos;

		for (uint64_t queryPos = 0; queryPos < minimizers.size(); ++queryPos) {
			auto& minimizer = minimizers[queryPos];
			if (minimizerToPos.find(minimizer) != minimizerToPos.end()) {
				auto& vectorTargetIdPos = minimizerToPos.at(minimizer);
				for (const auto& targetIdPos : vectorTargetIdPos) {
					if (targetIdPosToQuerypos.find(targetIdPos) == targetIdPosToQuerypos.end()) {
						targetIdPosToQuerypos[targetIdPos] = std::vector<uint64_t>();
						targetIdPosToQuerypos[targetIdPos].emplace_back(queryPos);
					} else {
						targetIdPosToQuerypos[targetIdPos].emplace_back(queryPos);
					}
				}
			}
		}

		for (const auto& targetIdPosToQueryposKeyVal : targetIdPosToQuerypos) {

			auto& targetIdPos = targetIdPosToQueryposKeyVal.first;
			auto& queryPos = targetIdPosToQueryposKeyVal.second;

			targetIdPosToQuerypos[targetIdPos] = { lowerMedian(queryPos) };
		}

		tsl::robin_map<pair, uint64_t, boost::hash<pair>> targetIdPosToCount;
		for (const auto& minimizer : minimizers) {
			if (minimizerToPos.find(minimizer) != minimizerToPos.end()) {
				auto& vectorTargetIdPos = minimizerToPos.at(minimizer);
				for (const auto& targetIdPos : vectorTargetIdPos) {
					if (targetIdPosToCount.find(targetIdPos) == targetIdPosToCount.end()) {
						targetIdPosToCount[targetIdPos] = 1;
					} else {
						++targetIdPosToCount[targetIdPos];
					}
				}
			}
		}

		bool mapped = false;

		for (const auto& targetIdPosToCountKeyVal : targetIdPosToCount) {

			auto& targetIdPos = targetIdPosToCountKeyVal.first;
			auto& score = targetIdPosToCountKeyVal.second;

			if (score >= opt::scoreThreshold) {

				mapped = true;
				auto& targetId = std::get<0>(targetIdPos);
				auto& targetPos = std::get<1>(targetIdPos);
				std::string orientation;
				uint64_t prev;
				uint64_t next;
				uint64_t curr = targetIdPosToQuerypos[targetIdPos][0];

				if (opt::mapPositions == 1) {

					auto prevPair = std::make_pair(targetId, targetPos - 1);
					auto nextPair = std::make_pair(targetId, targetPos + 1);

					if (targetIdPosToQuerypos.find(prevPair) == targetIdPosToQuerypos.end()) {
						prev = max;
					} else {
						prev = targetIdPosToQuerypos[prevPair][0];
					}

					if (targetIdPosToQuerypos.find(nextPair) == targetIdPosToQuerypos.end()) {
						next = max;
					} else {
						next = targetIdPosToQuerypos[nextPair][0];
					}
				} else {

					std::vector<uint64_t> prevVec;

					for (unsigned i = 1; i <= opt::mapPositions; ++i) {
						auto prevPair = std::make_pair(targetId, targetPos - i);
						if (targetIdPosToQuerypos.find(prevPair) != targetIdPosToQuerypos.end()) {
							prevVec.emplace_back(targetIdPosToQuerypos[prevPair][0]);
						}
					}

					if (prevVec.empty()) {
						prevVec.emplace_back(max);
					}
					prev = lowerMedian(prevVec);

					std::vector<uint64_t> nextVec;

					for (unsigned i = 1; i <= opt::mapPositions; ++i) {
						auto nextPair = std::make_pair(targetId, targetPos + i);
						if (targetIdPosToQuerypos.find(nextPair) != targetIdPosToQuerypos.end()) {
							nextVec.emplace_back(targetIdPosToQuerypos[nextPair][0]);
						}
					}

					if (nextVec.empty()) {
						nextVec.emplace_back(max);
					}

					next = upperMedian(nextVec);
				}

				orientation = determineOrientation(prev, curr, next);
#if _OPENMP
#pragma omp critical
#endif
				std::cout << targetId << "\t" << targetPos << "\t" << targetPos + 1 << "\t"
				          << queryId << "\t" << score << "\t" << orientation << std::endl;
			}
		}

		if (mapped) {
#if _OPENMP
#pragma omp atomic
#endif
			++num_mapped;
		}
	}
	std::cerr << "Mapped " << num_mapped << " sequences of " << queryToMinimizer.size() << " ("
	          << std::fixed << std::setprecision(2) << 100 * num_mapped / queryToMinimizer.size()
	          << "%))" << std::endl;
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
		                                { "mapPositions", required_argument, nullptr, 'm' },
		                                { "scoreThreshold", required_argument, nullptr, 'n' },
		                                { "threads", required_argument, nullptr, 't' },
		                                { "version", no_argument, &version, 1 },
		                                { nullptr, 0, nullptr, 0 } };

	for (int c; (c = getopt_long(argc, argv, "t:vn:m:", longopts, &optindex)) != -1;) {
		switch (c) {
		case 't': {
			std::stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				printErrorMsg(PROGRAM, "Invalid parameters! t: ");
				die = true;
			}
			if (opt::threads == 0) {
				printErrorMsg(PROGRAM, "Invalid parameters! t: ");
				die = true;
			}
			break;
		}
		case 'm': {
			std::stringstream convert(optarg);
			if (!(convert >> opt::mapPositions)) {
				printErrorMsg(PROGRAM, "Invalid parameters! m: ");
				die = true;
			}
			if (opt::mapPositions == 0) {
				printErrorMsg(PROGRAM, "Invalid parameters! m: ");
				die = true;
			}
			break;
		}
		case 'n': {
			std::stringstream convert(optarg);
			if (!(convert >> opt::scoreThreshold)) {
				printErrorMsg(PROGRAM, "Invalid parameters! n: ");
				die = true;
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
		printHelpDialog();
		exit(EXIT_SUCCESS);
	} else if (version != 0) {
		printVersion();
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

	if (inputFiles.empty() || inputFiles.size() != 3) {
		printErrorMsg(PROGRAM, "missing file operand");
		die = true;
	}

	if (die) {
		printHelpDialog();
		exit(EXIT_FAILURE);
	}

	std::vector<std::vector<std::string>> paths;
	readPaths(paths, inputFiles[0]);

	tsl::robin_map<std::string, std::vector<Minimizer>> moleculeToMinimizer;
	getMoleculeToMinimizer(moleculeToMinimizer, inputFiles[1]);

	tsl::robin_map<std::string, std::vector<Minimizer>> queryToMinimizer;
	getMoleculeToMinimizer(queryToMinimizer, inputFiles[2]);

	tsl::robin_map<Minimizer, tsl::robin_set<pair, boost::hash<pair>>> minimizerToPos;
	getMinimizerToPos(paths, moleculeToMinimizer, minimizerToPos);

	if (verbose) {
		std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB"
		          << std::endl;
	}

	mapQueryToTarget(queryToMinimizer, minimizerToPos);
}
