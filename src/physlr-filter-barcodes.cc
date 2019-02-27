#include <algorithm>
#include <cassert>
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

static inline void assert_good(const std::ios& stream, const std::string& path)
{
	if (!stream.good()) {
		std::cerr << "error: " << strerror(errno) << ": " << path << '\n';
		exit(EXIT_FAILURE);
	}
}

static void printErrorMsg(const std::string& progname, const std::string& msg)
{
	std::cerr << progname << ": " << msg << "\nTry 'physlr-indexlr --help' for more information.\n";
}

static void printUsage(const std::string& progname)
{
	std::cout << "Usage:  " << progname
		<< "  -k K -w W [-v] [-o file] file...\n\n"
		"  -k K       use K as k-mer size\n"
		"  -w W       use W as sliding-window size\n"
		"  -s         silent; disable verbose output\n"
		"  -o file    write output to file, default is stdout\n"
		"  -n         minimum number of minimizers per barcode\n"
		"  -N         maximum number of minimizers per barcode\n"
		"  --help     display this help and exit\n"
		"  file       space separated list of FASTQ files\n";
}

/*typedef uint64_t Mx;
typedef std::unordered_set<Mx> Mxs;
typedef std::string Bx;
typedef std::unordered_map<Bx, Mxs> BxtoMxs;
*/
using Mx = uint64_t;
using Mxs = std::unordered_set<Mx>;
using Bx = std::string;
using BxtoMxs = std::unordered_map<Bx, Mxs>;

static Mxs splitMinimizers(const std::string& mx_line) {
	Mxs mx_set;
	std::istringstream iss(mx_line);
	std::string mx;
	while (iss >> mx) {
		mx_set.insert(strtoull(mx.c_str(), nullptr, 0));
	}
	return mx_set;
}

static BxtoMxs read_minimizers(std::istream &is, const std::string &ipath) {
	if (is.peek() == std::ifstream::traits_type::eof()) {
		std::cerr << "physlr-filterbarcodes: error: Empty input file: " << ipath << '\n';
		exit(EXIT_FAILURE);
	}
	BxtoMxs bxtomxs;
	Bx bx;
	std::string mx_line;
	while ((is >> bx) && (getline(is, mx_line))) {
		auto mxs = splitMinimizers(mx_line);
		bxtomxs[bx].insert(mxs.begin(), mxs.end());
	}
	return bxtomxs;
}

static void write(BxtoMxs bxtomxs, std::ostream& os, const std::string& opath) {
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
}

static std::unordered_map<Mx, unsigned> counterForRemove(const BxtoMxs& bxtomxs) {
	std::unordered_map<Mx, unsigned> counts;
	for (const auto& item : bxtomxs)
	{
		const auto& mxs = item.second;
		for (const auto& mx : mxs) {
			if (counts.find(mx) == counts.end()) {
				counts[mx] = 1;
			} else {
				++counts[mx];
			}
		}
	}
	return counts;
}

static void remove_singleton_minimizers(BxtoMxs& bxtomxs, bool silent) {
	std::unordered_map<Mx, unsigned> counts = counterForRemove(bxtomxs); 
	Mxs singletons;
	for (const auto& item : counts) {
		if (item.second < 2) {
			singletons.insert(item.first);
		}
	}
	Mxs uniqueMxs;
	for (const auto& item : bxtomxs) {
		const auto& mxs = item.second;
		Mxs not_singletons;
		for (const auto& mx : mxs) {
			uniqueMxs.insert(mx);
			if (singletons.find(mx) == singletons.end()) {
				not_singletons.insert(mx);
			}
		}
		bxtomxs[item.first] = std::move(not_singletons);
	}
	if (!silent) {
		std::cerr << "Removed " << singletons.size() << " minimizers that occur once of " << uniqueMxs.size() << " (" << std::setprecision(1) << std::fixed << 100.0 * singletons.size() / uniqueMxs.size() << "%)\n";
	}
}

static void physlr_filterbarcodes(std::istream& is, const std::string& ipath, std::ostream& os, const std::string& opath, const size_t n, const size_t N, bool silent) {
	if (is.peek() == std::ifstream::traits_type::eof()) {
		std::cerr << "physlr-filterbarcodes: error: Empty input file: " << ipath << '\n';
		exit(EXIT_FAILURE);
	}
	BxtoMxs bxtomxs = read_minimizers(is, ipath);
	unsigned initial_size = bxtomxs.size();
	remove_singleton_minimizers(bxtomxs, silent);
	unsigned too_few = 0, too_many = 0;
	for (auto it = bxtomxs.begin(); it != bxtomxs.end(); ) {
		auto& mxs = it->second;
		if (mxs.size() < n) {
			++too_few;
			it = bxtomxs.erase(it);
		}
		else if (mxs.size() >= N) {
			++too_many;
			it = bxtomxs.erase(it);
		}
		else {
			++it;
		}
	}
	if (!silent) {
		std::cerr << "Discarded " << too_few << " barcodes with too few minimizers of " << initial_size << " (" << std::setprecision(1) << std::fixed << 100.0 * too_few / initial_size << "%)\n";
		std::cerr << "Discarded " << too_many << " barcodes with too many minimizers of " << initial_size << " (" << std::setprecision(1) << std::fixed << 100.0 * too_many / initial_size << "%)\n";
		std::cerr << "Wrote " << initial_size - too_few - too_many << " barcodes\n";
	}
	write(bxtomxs, os, opath);
}

int main(int argc, char *argv[])
{
	auto progname = "physlr-filterbarcodes";
	int c;
	int optindex = 0;
	static int help = 0;
	unsigned k = 0;
	unsigned w = 0;
	unsigned n = 0;
	unsigned N = 0;
	bool silent = false;
	bool failed = false;
	bool w_set = false;
	bool k_set = false;
	bool n_set = false;
	bool N_set = false;
	char* end = nullptr;
	std::string outfile("/dev/stdout");
	static const struct option longopts[] = { { "help", no_argument, &help, 1 },{ nullptr, 0, nullptr, 0 } };
	while ((c = getopt_long(argc, argv, "k:w:o:s:n:N:", longopts, &optindex)) != -1) {
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
	else if (!k_set) {
		printErrorMsg(progname, "missing option -- 'k'");
		failed = true;
	}
	else if (!w_set) {
		printErrorMsg(progname, "missing option -- 'w'");
		failed = true;
	}
	else if (!n_set) {
		n = 1;
	}
	else if (!N_set) {
		N = INT_MAX;
	}
	else if (k == 0) {
		printErrorMsg(progname, "option has incorrect argument -- 'k'");
		failed = true;
	}
	else if (w == 0) {
		printErrorMsg(progname, "option has incorrect argument -- 'w'");
		failed = true;
	}
	else if (n == 0) {
		printErrorMsg(progname, "option has incorrect argument -- 'n'");
		failed = true;
	}
	else if (N == 0) {
		printErrorMsg(progname, "option has incorrect argument -- 'N'");
		failed = true;
	}
	else if (infiles.empty()) {
		printErrorMsg(progname, "missing file operand");
		failed = true;
	}
	if (failed) {
		exit(EXIT_FAILURE);
	}
	std::ofstream ofs(outfile);
	assert_good(ofs, outfile);
	for (auto &infile : infiles) {
		if (infile == "-") {
			infile = "/dev/stdin";
		}
		std::ifstream ifs(infile);
		assert_good(ifs, infile);
		physlr_filterbarcodes(ifs, infile, ofs, outfile, n, N , silent);
	}
	ofs.flush();
	assert_good(ofs, outfile);
	return 0;
}
