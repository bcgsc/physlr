#include "robin_hood.h"
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
#include <boost/functional/hash.hpp>

static std::chrono::time_point<std::chrono::steady_clock> t0; // NOLINT(cert-err58-cpp)

struct PairHash {
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U>& p) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.first);
        boost::hash_combine(seed, p.second);
        return seed;
    }
};

// struct PairHash {
//     inline std::size_t operator()(const std::pair<int,int> & v) const {
//         return v.first*31+v.second;
//     }
// };

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
	          << "  -n n -N N [-s] [-p] [-o file] file...\n\n"
	             "  -s         silent; disable verbose output\n"
	             "  -o file    write output to file, default is stdout\n"
	             "  -n         minimum number of minimizers per barcode\n"
	             "  -N         maximum number of minimizers per barcode\n"
                 "  -p         indicate that the minimizers contain position\n"
	             "  --help     display this help and exit\n"
	             "  file       space separated list of FASTQ files\n";
}

using Mx = uint64_t;
using Mxs = robin_hood::unordered_set<Mx>;
using MxwithPos = std::pair<Mx, std::size_t>;
using MxswithPos = robin_hood::unordered_set<MxwithPos, PairHash>;
using Bx = std::string;
using BxtoMxs = robin_hood::unordered_map<Bx, MxswithPos>;
using MxtoCounts = robin_hood::unordered_map<Mx, unsigned>;

static BxtoMxs
readMxs(std::istream& is, const std::string& ipath, bool silent, bool positioned)
{
	if (is.peek() == std::ifstream::traits_type::eof()) {
		std::cerr << "physlr-filterbarcodes: error: Empty input file: " << ipath << '\n';
		exit(EXIT_FAILURE);
	}
	BxtoMxs bxtomxs;
	Bx bx;
	std::string mx;
    std::string mxpos;
	std::string mx_line;
    bool warned = false;
    bool hasPos = false;
	while ((is >> bx) && (getline(is, mx_line))) {
		std::istringstream iss(mx_line);
		if (positioned) {
            while (iss >> mx) {
                size_t mxsep = mx.find(":");
                if (mxsep == std::string::npos){
                    // print error message if the minimizer position is not given and exit
                    std::cerr << "physlr-filter-barcodes: error: minimizer position not provided: " << ipath << '\n';
                    std::cerr << " barcode: " << bx << " | minimizer: " << mx << '\n';
                    exit(EXIT_FAILURE);
                }
                mxpos = mx.substr(mxsep+1);
                mx = mx.substr(0, mxsep);
                bxtomxs[bx].insert(std::make_pair(strtoull(mx.c_str(), nullptr, 0), strtoull(mxpos.c_str(), nullptr, 0)));
            }
        } else { 
            while (iss >> mx) {
                if (!warned){
                    warned = true;
                    if (mx.find(":") != std::string::npos){
                        // print error message if the minimizer position is included in file and exit
                        std::cerr << "physlr-filter-barcodes: warning: skipping minimizer position." << '\n';
                        std::cerr << " Set -p to avoid skipping." << '\n';
                        hasPos = true;
                        mx = mx.substr(0, mx.find(":"));
                    }
                } else if (hasPos) {
                    mx = mx.substr(0, mx.find(":"));
                }
                bxtomxs[bx].insert(std::make_pair(strtoull(mx.c_str(), nullptr, 0), 0));
            }
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
writeMxs(BxtoMxs bxtomxs, std::ostream& os, const std::string& opath, bool silent, bool positioned)
{
	for (const auto& item : bxtomxs) {
		const auto& bx = item.first;
		const auto& mxs = item.second;
		os << bx;
		char sep = '\t';
		for (const auto& mx : mxs) {
			if (positioned)
                os << sep << mx.first << ":" << mx.second;
			else
                os << sep << mx.first;
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
            auto mxhvalue = mx.first;
			if (counts.find(mxhvalue) == counts.end()) {
				counts[mxhvalue] = 1;
			} else {
				++counts[mxhvalue];
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
	robin_hood::unordered_map<Mx, bool> counted;
	for (const auto& item : bxtomxs) {
		MxswithPos not_singletons;
        //std::unique_ptr<MxswithPos> not_singletons;
		//not_singletons->reserve(item.second.size());
		not_singletons.reserve(item.second.size());
		for (const auto& mx : item.second) {
			if (counts[mx.first] >= 2) {
				//not_singletons->insert({mx.first, mx.second});
				not_singletons.insert({mx.first, mx.second});
			} else {
				++singletons;
			}
		}
		
		uniqueMxs.insert(item.second.begin(), item.second.end()); // need to do before removing!
		//bxtomxs[item.first] = std::move(*not_singletons);
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
    bool silent,
    bool positioned)
{
	if (is.peek() == std::ifstream::traits_type::eof()) {
		std::cerr << "physlr-filterbarcodes: error: Empty input file: " << ipath << '\n';
		exit(EXIT_FAILURE);
	}
	BxtoMxs bxtomxs = readMxs(is, ipath, silent, positioned);
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
	writeMxs(bxtomxs, os, opath, silent, positioned);
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
    bool positioned = false;
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
		case 'p':
            positioned = true;
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
		physlr_filterbarcodes(ifs, infile, ofs, outfile, n, N, silent, positioned);
	}
	ofs.flush();
	assert_good(ofs, outfile);
	return 0;
}