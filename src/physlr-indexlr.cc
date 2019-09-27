// Convert linked-reads to minimizers using ntHash-2.0.0.
// Usage:  physlr-indexlr -k K -w W [-v] [-o file] file...
// Output: Each line of output is a barcode followed by a list of minimzers.
// Originally written for Physlr: (https://github.com/bcgsc/physlri
// Written by Vladimir Nikolic (schutzekatze) and Shaun Jackman (@sjackman)

#include "btl_bloomfilter/BloomFilter.hpp"
#include "indexlr-workers.h"

#include <cassert>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace std::chrono;

class Timer
{
  public:
	Timer() { m_startTimepoint = high_resolution_clock::now(); }
	~Timer() { stop(); }
	void stop()
	{
		auto endTimepoint = high_resolution_clock::now();
		auto start = time_point_cast<microseconds>(m_startTimepoint).time_since_epoch().count();
		auto end = time_point_cast<microseconds>(endTimepoint).time_since_epoch().count();
		std::cerr << (end - start) * 0.000001 << std::endl;
	}

  private:
	time_point<high_resolution_clock> m_startTimepoint;
};

// Read a FASTQ file and reduce each read to a set of minimizers
static void
minimizeReads(
    const std::string& ipath,
    const std::string& opath,
    const size_t k,
    const size_t w,
    const size_t t,
    const bool withRBloomFilter,
    const bool withSBloomFilter,
    const bool withPositions,
    const bool withStrands,
    const bool verbose,
    BloomFilter& rBloomFilter,
    BloomFilter& sBloomFilter)
{
	InputWorker inputWorker(ipath);
	OutputWorker outputWorker(opath, inputWorker);

	inputWorker.start();
	outputWorker.start();

	auto minimizeWorkers = std::vector<MinimizeWorker>(
	    t,
	    MinimizeWorker(
	        k,
	        w,
	        withRBloomFilter,
	        withSBloomFilter,
	        withPositions,
	        withStrands,
	        verbose,
	        rBloomFilter,
	        sBloomFilter,
	        inputWorker,
	        outputWorker));
	for (auto& worker : minimizeWorkers) {
		worker.start();
	}
	for (auto& worker : minimizeWorkers) {
		worker.join();
	}

	inputWorker.join();
	outputWorker.join();
}

static void
printErrorMsg(const std::string& progname, const std::string& msg)
{
	std::cerr << progname << ": " << msg << "\nTry 'physlr-indexlr --help' for more information.\n";
}

static void
printUsage(const std::string& progname)
{
	std::cout << "Usage:  " << progname
	          << "  -k K -w W [-b bf_path] [-v] [-o FILE] FILE...\n\n"
	             "  -k K        use K as k-mer size\n"
	             "  -w W        use W as sliding-window size\n"
	             "  -r repetitive_bf_path  use a bloomfilter to filter out bad minimizers\n"
	             "  -s solid_bf_path  use a bloomfilter to filter in solid minimizers\n"
	             "  --pos       include minimizer positions in the output\n"
	             "  --strand    include minimizer strand in the output\n"
	             "  -v          enable verbose output\n"
	             "  -o FILE     write output to FILE, default is stdout\n"
	             "  -t N        use N number of threads (default 1, max 5)\n"
	             "  --help      display this help and exit\n"
	             "  FILE        space separated list of FASTQ files\n";
}

int
main(int argc, char* argv[])
{
	auto progname = "physlr-indexlr";
	int c;
	int optindex = 0;
	static int help = 0;
	unsigned k = 0;
	unsigned w = 0;
	bool verbose = false;
	bool withRBloomFilter = false;
	bool withSBloomFilter = false;
	BloomFilter repetitiveBloomFilter;
	BloomFilter solidBloomFilter;
	unsigned t = 1;
	bool failed = false;
	bool w_set = false;
	bool k_set = false;
	static int withPositions = 0;
	static int withStrands = 0;
	char* end = nullptr;
	std::string outfile("/dev/stdout");
	static const struct option longopts[] = { { "pos", no_argument, &withPositions, 1 },
		                                      { "strand", no_argument, &withStrands, 1 },
		                                      { "help", no_argument, &help, 1 },
		                                      { nullptr, 0, nullptr, 0 } };
	while ((c = getopt_long(argc, argv, "k:w:o:vt:r:s:", longopts, &optindex)) != -1) {
		switch (c) {
		case 0:
			break;
		case 'k':
			k_set = true;
			k = strtoul(optarg, &end, 10);
			break;
		case 'w':
			w_set = true;
			w = strtoul(optarg, &end, 10);
			break;
		case 'o':
			outfile.assign(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		case 't':
			t = strtoul(optarg, &end, 10);
			if (t > 5) {
				t = 5;
				std::cerr << progname
				          << ": Using more than 5 threads does not scale, reverting to 5.\n";
			}
			break;
		case 'r': {
			withRBloomFilter = true;
			std::cerr << "Loading repetitive Bloom filter from " << optarg << std::endl;
			Timer timer;
			try {
				repetitiveBloomFilter.loadFilter(optarg);
			} catch (const std::exception& e) {
				std::cerr << e.what() << '\n';
			}
			std::cerr << "Finished loading repetitive Bloom filter in sec: ";
			timer.~Timer();
			break;
		}
		case 's': {
			withSBloomFilter = true;
			std::cerr << "Loading solid Bloom filter from " << optarg << std::endl;
			Timer timer;
			try {
				solidBloomFilter.loadFilter(optarg);
			} catch (const std::exception& e) {
				std::cerr << e.what() << '\n';
			}
			std::cerr << "Finished loading solid Bloom filter in sec: ";
			timer.~Timer();
			break;
		}
		default:
			exit(EXIT_FAILURE);
		}
	}
	std::vector<std::string> infiles(&argv[optind], &argv[argc]);
	if (argc < 2) {
		printUsage(progname);
		exit(EXIT_FAILURE);
	}
	if (help != 0) {
		printUsage(progname);
		exit(EXIT_SUCCESS);
	} else if (!k_set) {
		printErrorMsg(progname, "missing option -- 'k'");
		failed = true;
	} else if (!w_set) {
		printErrorMsg(progname, "missing option -- 'w'");
		failed = true;
	} else if (k == 0) {
		printErrorMsg(progname, "option has incorrect argument -- 'k'");
		failed = true;
	} else if (w == 0) {
		printErrorMsg(progname, "option has incorrect argument -- 'w'");
		failed = true;
	} else if (infiles.empty()) {
		printErrorMsg(progname, "missing file operand");
		failed = true;
	} else if (withRBloomFilter && withSBloomFilter) {
		printErrorMsg(progname, "missing file operand");
		failed = true;
	}
	if (failed) {
		exit(EXIT_FAILURE);
	}

	for (auto& infile : infiles) {
		minimizeReads(
		    infile == "-" ? "/dev/stdin" : infile,
		    outfile,
		    k,
		    w,
		    t,
		    withRBloomFilter,
		    withSBloomFilter,
		    withPositions,
		    withStrands,
		    verbose,
		    repetitiveBloomFilter,
		    solidBloomFilter);
	}

	return 0;
}
