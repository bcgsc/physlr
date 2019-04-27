#include "tsl/robin_map.h"
#include "tsl/robin_set.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

static std::chrono::time_point<std::chrono::steady_clock> t0; // NOLINT(cert-err58-cpp)

static inline void
assert_good(const std::ios& stream, const std::string& path)
{
	if (!stream.good()) {
		std::cerr << "error: " << strerror(errno) << ": " << path << '\n';
		exit(EXIT_FAILURE);
	}
}

static void
printErrorMsg(const std::string& progname, const std::string& msg)
{
	std::cerr << progname << ": " << msg
	          << "\nTry 'physlr-filter-barcodes --help' for more information.\n";
}

static void
printUsage(const std::string& progname)
{
	std::cout << "Usage:  " << progname
	          << "  -n n -N N [-s] [-o file] file...\n\n"
	             "  -s         silent; disable verbose output\n"
	             "  -o file    write output to file, default is stdout\n"
	             "  -n         minimum number of minimizers per barcode\n"
	             "  -N         maximum number of minimizers per barcode\n"
	             "  --help     display this help and exit\n"
	             "  file       space separated list of FASTQ files\n";
}

using Mx = uint64_t;
using Mxs = tsl::robin_set<Mx>;
using Bx = std::string;
using BxtoMxs = tsl::robin_map<Bx, Mxs>;
using MxtoCounts = tsl::robin_map<Mx, unsigned>;

static BxtoMxs
readMxs(std::istream& is, const std::string& ipath, bool silent)
{
	if (is.peek() == std::ifstream::traits_type::eof()) {
		std::cerr << "physlr-filterbarcodes: error: Empty input file: " << ipath << '\n';
		exit(EXIT_FAILURE);
	}
	BxtoMxs bxtomxs;
	Bx bx;
	std::string mx;
	std::string mx_line;
	while ((is >> bx) && (getline(is, mx_line))) {
		std::istringstream iss(mx_line);
		while (iss >> mx) {
			bxtomxs[bx].insert(strtoull(mx.c_str(), nullptr, 0));
		}
	}
	auto t = std::chrono::steady_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0);
	if (!silent) {
		std::cerr << "Time at readMxs (ms): " << diff.count() << '\n';
	}
	return bxtomxs;
}

static void
writeMxs(BxtoMxs bxtomxs, std::ostream& os, const std::string& opath, bool silent)
{
	for (const auto& item : bxtomxs) {
		const auto& bx = item.first;
		const auto& mxs = item.second;
		os << bx;
		char sep = '\t';
		for (const auto& mx : mxs) {
			os << sep << mx;
			sep = ' ';
		}
		os << '\n';
		assert_good(os, opath);
	}
	auto t = std::chrono::steady_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0);
	if (!silent) {
		std::cerr << "Time at writeMxs (ms): " << diff.count() << '\n';
	}
}

static MxtoCounts
countMxs(const BxtoMxs& bxtomxs, bool silent)
{
	MxtoCounts counts;
	for (const auto& item : bxtomxs) {
		const auto& mxs = item.second;
		for (const auto& mx : mxs) {
			if (counts.find(mx) == counts.end()) {
				counts[mx] = 1;
			} else {
				++counts[mx];
			}
		}
	}
	auto t = std::chrono::steady_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0);
	if (!silent) {
		std::cerr << "Time at countMxs (ms): " << diff.count() << '\n';
	}
	return counts;
}

static void
removeSingletonMxs(BxtoMxs& bxtomxs, bool silent)
{
	MxtoCounts counts = countMxs(bxtomxs, silent);
	std::cerr << "Counted " << counts.size() << " minimizers." << '\n';
	Mxs uniqueMxs;
	uint64_t singletons = 0;
	tsl::robin_map<Mx, bool> counted;
	for (const auto& item : bxtomxs) {
		Mxs not_singletons;
		not_singletons.reserve(item.second.size());
		for (const auto& mx : item.second) {
			if (counts[mx] >= 2) {
				not_singletons.insert(mx);
			} else {
				++singletons;
			}
		}
		uniqueMxs.insert(item.second.begin(), item.second.end()); // need to do before removing!
		bxtomxs[item.first] = std::move(not_singletons);
	}
	auto t = std::chrono::steady_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0);
	if (!silent) {
		std::cerr << "Time at removeSingletonMxs (ms): " << diff.count() << '\n';
		std::cerr << "Removed " << singletons << " minimizers that occur once of "
		          << uniqueMxs.size() << " (" << std::setprecision(1) << std::fixed
		          << 100.0 * singletons / uniqueMxs.size() << "%)\n";
	}
}

static void
physlr_filterbarcodes(
    std::istream& is,
    const std::string& ipath,
    std::ostream& os,
    const std::string& opath,
    const size_t n,
    const size_t N,
    bool silent)
{
	if (is.peek() == std::ifstream::traits_type::eof()) {
		std::cerr << "physlr-filterbarcodes: error: Empty input file: " << ipath << '\n';
		exit(EXIT_FAILURE);
	}
	BxtoMxs bxtomxs = readMxs(is, ipath, silent);
	unsigned initial_size = bxtomxs.size();
	removeSingletonMxs(bxtomxs, silent);
	unsigned too_few = 0, too_many = 0;
	std::cerr << "There are " << initial_size << " barcodes." << '\n';
	for (auto it = bxtomxs.begin(); it != bxtomxs.end();) {
		auto& mxs = it->second;
		if (mxs.size() < n) {
			++too_few;
			it = bxtomxs.erase(it);
		} else if (mxs.size() >= N) {
			++too_many;
			it = bxtomxs.erase(it);
		} else {
			++it;
		}
	}
	auto t = std::chrono::steady_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0);
	if (!silent) {
		std::cerr << "Time at filterbarcodes (ms): " << diff.count() << '\n';
		std::cerr << "Discarded " << too_few << " barcodes with too few minimizers of "
		          << initial_size << " (" << std::setprecision(1) << std::fixed
		          << 100.0 * too_few / initial_size << "%)\n";
		std::cerr << "Discarded " << too_many << " barcodes with too many minimizers of "
		          << initial_size << " (" << std::setprecision(1) << std::fixed
		          << 100.0 * too_many / initial_size << "%)\n";
		std::cerr << "Wrote " << initial_size - too_few - too_many << " barcodes\n";
	}
	writeMxs(bxtomxs, os, opath, silent);
}

int
main(int argc, char* argv[])
{
	t0 = std::chrono::steady_clock::now();
	auto progname = "physlr-filterbarcodes";
	int c;
	int optindex = 0;
	static int help = 0;
	unsigned n = 0;
	unsigned N = 0;
	bool silent = false;
	bool failed = false;
	bool n_set = false;
	bool N_set = false;
	char* end = nullptr;
	std::string outfile("/dev/stdout");
	static const struct option longopts[] = { { "help", no_argument, &help, 1 },
		                                      { nullptr, 0, nullptr, 0 } };
	while ((c = getopt_long(argc, argv, "o:s:n:N:", longopts, &optindex)) != -1) {
		switch (c) {
		case 0:
			break;
		case 'o':
			outfile.assign(optarg);
			break;
		case 's':
			silent = true;
			break;
		case 'n':
			n_set = true;
			n = strtoul(optarg, &end, 10);
			break;
		case 'N':
			N_set = true;
			N = strtoul(optarg, &end, 10);
			break;
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
	}
	if (!n_set) {
		n = 1;
	}
	if (!N_set) {
		N = INT_MAX;
	}
	if (n == 0) {
		printErrorMsg(progname, "option has incorrect argument -- 'n'");
		failed = true;
	}
	if (N == 0) {
		printErrorMsg(progname, "option has incorrect argument -- 'N'");
		failed = true;
	}
	if (infiles.empty()) {
		printErrorMsg(progname, "missing file operand");
		failed = true;
	}
	if (failed) {
		exit(EXIT_FAILURE);
	}
	std::ofstream ofs(outfile);
	assert_good(ofs, outfile);
	for (auto& infile : infiles) {
		if (infile == "-") {
			infile = "/dev/stdin";
		}
		std::ifstream ifs(infile);
		assert_good(ifs, infile);
		physlr_filterbarcodes(ifs, infile, ofs, outfile, n, N, silent);
	}
	ofs.flush();
	assert_good(ofs, outfile);
	return 0;
}
