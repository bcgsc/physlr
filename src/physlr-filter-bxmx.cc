// #include "tsl/robin_map.h"
// #include "tsl/robin_set.h"
#include "robin_hood.h"
#include <algorithm>
#include <cassert>
#include <cfenv>
#include <chrono>
#include <climits>
#include <cmath>
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
	          << "\nTry 'physlr-filter-barcodes-minimizers --help' for more information.\n";
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
	             "  -C         maximum minimizer multiplicity\n"
	             "  --help     display this help and exit\n"
	             "  file       space separated list of FASTQ files\n";
}

struct PairHash {
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U>& p) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.first);
        boost::hash_combine(seed, p.second);
        return seed;
    }
};

using Mx = uint64_t;
using Mxs = robin_hood::unordered_set<Mx>;
using MxwithPos = std::pair<Mx, std::size_t>;
using MxswithPos = robin_hood::unordered_set<MxwithPos, PairHash>;
using Bx = std::string;
using BxtoMxs = robin_hood::unordered_map<Bx, MxswithPos>;
using MxtoCounts = robin_hood::unordered_map<Mx, unsigned>;

static void
readMxs(std::istream& is, const std::string& ipath, bool silent, BxtoMxs& bxtomxs, bool positioned)
{
	if (is.peek() == std::ifstream::traits_type::eof()) {
		std::cerr << "physlr-filter-barcodes-minimizers: error: Empty input file: " << ipath
		          << '\n';
		exit(EXIT_FAILURE);
	}
	bxtomxs.clear();
	assert(bxtomxs.empty());
	Bx bx;
	// std::string mx_string;
	std::string mx;
    std::string mxpos;
	std::string mx_line;
	std::istringstream iss;
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
		
	// while ((is >> bx) && (getline(is, mx_line))) {
	// 	iss.clear();
	// 	iss.str(mx_line);
	// 	auto& mxs = bxtomxs[bx];
	// 	while (iss >> mx_string) {
	// 		Mx mx = strtoull(mx_string.c_str(), nullptr, 0);
	// 		mxs.insert(mx);
	// 	}
	// }
	if (!silent) {
		auto t = std::chrono::steady_clock::now();
		auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0);
		std::cerr << "Time after reading minimizers (ms): " << diff.count() << '\n';
		std::cerr << "Read " << bxtomxs.size() << " barcodes.\n\n";
	}
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
		std::cerr << "Time after writing minimizers (ms): " << diff.count() << '\n';
		std::cerr << "Wrote " << bxtomxs.size() << " barcodes.\n\n";
	}
}

static void
countMxs(const BxtoMxs& bxtomxs, MxtoCounts& counts, bool silent)
{
	counts.clear();
	assert(counts.empty());
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
		std::cerr << "Time after counting minimizers (ms): " << diff.count() << '\n';
		std::cerr << "Counted " << counts.size() << " minimizers.\n\n";
	}
}

static void
removeSingletonMxs(BxtoMxs& bxtomxs, MxtoCounts& counts, bool silent)
{
	countMxs(bxtomxs, counts, silent);
	assert(!counts.empty() && !bxtomxs.empty());
	Mxs singletons;
	for (auto it = bxtomxs.begin(); it != bxtomxs.end(); ++it) {
		const auto& bx = it->first;
		auto& mxs = it->second;
		MxswithPos not_singletons;
		not_singletons.reserve(mxs.size());
		for (const auto& mx : mxs) {
			if (counts[mx.first] >= 2) {
				//not_singletons->insert({mx.first, mx.second});
				not_singletons.insert({mx.first, mx.second});
			} else {
				singletons.insert({mx.first, mx.second});
			}
		}
		bxtomxs[bx] = std::move(not_singletons);
	}
	for (auto it = counts.begin(); it != counts.end();) {
		if (it->second < 2) {
			it = counts.erase(it);
		} else {
			++it;
		}
	}
	if (!silent) {
		auto t = std::chrono::steady_clock::now();
		auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0);
		std::cerr << "Time after removing singleton minimizers (ms): " << diff.count() << '\n';
		std::cerr << "Removed " << singletons.size() << " minimizers that occur once of "
		          << counts.size() << " (" << std::setprecision(1) << std::fixed
		          << 100.0 * singletons.size() / counts.size() << "%)\n\n";
	}
}

static void
filter_barcodes(const size_t n, const size_t N, bool silent, BxtoMxs& bxtomxs)
{
	unsigned initial_size = bxtomxs.size();
	unsigned too_few = 0, too_many = 0;
	for (auto it = bxtomxs.begin(); it != bxtomxs.end();) {
		const auto& mxs = it->second;
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
	if (!silent) {
		auto t = std::chrono::steady_clock::now();
		auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0);
		std::cerr << "Time after filtering barcodes (ms): " << diff.count() << '\n';
		std::cerr << "Discarded " << too_few << " barcodes with too few minimizers of "
		          << initial_size << " (" << std::setprecision(1) << std::fixed
		          << 100.0 * too_few / initial_size << "%)\n";
		std::cerr << "Discarded " << too_many << " barcodes with too many minimizers of "
		          << initial_size << " (" << std::setprecision(1) << std::fixed
		          << 100.0 * too_many / initial_size << "%)\n";
		std::cerr << "There are " << initial_size - too_few - too_many << " barcodes remaining\n\n";
	}
}

static std::vector<float>
quantile(std::vector<float> quantiles, std::vector<unsigned> values)
{
	std::sort(values.begin(), values.end());
	std::vector<float> qs;
	qs.reserve(quantiles.size());
	std::fesetround(FE_TONEAREST);
	for (const auto& p : quantiles) {
		qs.push_back(values[size_t(std::nearbyint(p * (values.size() - 1)))]);
	}
	return qs;
}

static void
filter_minimizers(bool silent, BxtoMxs& bxtomxs, MxtoCounts counts, unsigned C)
{
	std::vector<unsigned> values;
	values.reserve(counts.size());
	for (auto& item : counts) {
		values.push_back(item.second);
	}
	std::vector<float> q = { 0.25, 0.5, 0.75 };
	q = quantile(q, values);
	unsigned high_whisker = int(q[2] + 1.5 * (q[2] - q[0]));
	if (C == 0) {
		C = high_whisker;
	}
	Mxs repetitives;
	unsigned empty = 0;
	unsigned bx_size = bxtomxs.size();
	for (auto it = bxtomxs.begin(); it != bxtomxs.end();) {
		auto& bx = it->first;
		auto mxs = it->second;
		MxswithPos not_repetitives;
		not_repetitives.reserve(mxs.size());
		for (const auto& mx : mxs) {
			auto mxhvalue = mx.first;
			auto mxpos = mx.second;
			if (counts[mxhvalue] < C) {
				not_repetitives.insert({mxhvalue, mxpos});
			} else {
				repetitives.insert({mxhvalue, mxpos});
			}
		}
		bxtomxs[bx] = std::move(not_repetitives);
		mxs = it->second;
		if (mxs.empty()) {
			++empty;
			it = bxtomxs.erase(it);
		} else {
			++it;
		}
	}
	if (!silent) {
		auto t = std::chrono::steady_clock::now();
		auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0);
		std::cerr << "Time after removing repetitive minimizers (ms): " << diff.count() << '\n';
		std::cerr << "Minimizer frequency: Q1=" << q[0] << " Q2=" << q[1] << " Q3=" << q[2]
		          << " C= " << C << '\n';
		std::cerr << "Removed " << repetitives.size() << " most frequent minimizers of "
		          << counts.size() << '\n';
		std::cerr << "Removed " << empty << " empty barcodes of " << bx_size << '\n';
		std::cerr << "There are " << bxtomxs.size() << " barcodes remaining." << '\n';
	}
}

int
main(int argc, char* argv[])
{
	t0 = std::chrono::steady_clock::now();
	auto progname = "physlr-filter-barcodes-minimizers";
	int c;
	int optindex = 0;
	static int help = 0;
	unsigned n = 0;
	unsigned N = 0;
	unsigned C = 0;
	bool positioned = true; // temporarily set to true, set to false and modify physlr-make
	bool silent = false;
	bool failed = false;
	bool n_set = false;
	bool N_set = false;
	char* end = nullptr;
	std::string outfile("/dev/stdout");
	static const struct option longopts[] = { { "help", no_argument, &help, 1 },
		                                      { nullptr, 0, nullptr, 0 } };
	while ((c = getopt_long(argc, argv, "o:s:n:N:C:", longopts, &optindex)) != -1) {
		switch (c) {
		case 0:
			break;
		case 'o':
			outfile.assign(optarg);
			break;
		case 's':
			silent = true;
			break;
		case 'p':
			positioned = true;
			break;
		case 'n':
			n_set = true;
			n = strtoul(optarg, &end, 10);
			break;
		case 'N':
			N_set = true;
			N = strtoul(optarg, &end, 10);
			break;
		case 'C':
			C = strtoul(optarg, &end, 10);
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
		BxtoMxs bxtomxs;
		MxtoCounts counts;
		readMxs(ifs, infile, silent, bxtomxs, positioned);
		removeSingletonMxs(bxtomxs, counts, silent);
		filter_barcodes(n, N, silent, bxtomxs);
		removeSingletonMxs(bxtomxs, counts, silent);
		filter_minimizers(silent, bxtomxs, counts, C);
		writeMxs(bxtomxs, ofs, outfile, silent, positioned);
	}
	ofs.flush();
	assert_good(ofs, outfile);
	return 0;
}
