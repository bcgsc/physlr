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
#include "tsl/robin_map.h"
#include "tsl/robin_set.h"
#include <boost/range/algorithm/count.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/assign.hpp>

using Mx = uint64_t;
using Mxs = std::vector<Mx>;
using Bx = std::string;
using BxtoMxs = tsl::robin_map<Bx, Mxs>;
using MxtoCount = std::unordered_map<Mx, unsigned>;

static std::chrono::time_point<std::chrono::steady_clock> t0; //NOLINT(cert-err58-cpp)

static inline void assert_good(const std::ios& stream, const std::string& path)
{
	if (!stream.good()) {
		std::cerr << "error: " << strerror(errno) << ": " << path << '\n';
		exit(EXIT_FAILURE);
	}
}

static void printErrorMsg(const std::string& progname, const std::string& msg)
{
	std::cerr << progname << ": " << msg << "\nTry 'physlr-filter-barcodes-minimizers --help' for more information.\n";
}

static void printUsage(const std::string& progname)
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

static void read_minimizers(std::istream &is, const std::string &ipath, bool silent, MxtoCount& counts, BxtoMxs& bxtomxs) {
	if (is.peek() == std::ifstream::traits_type::eof()) {
		std::cerr << "physlr-filter-barcodes-minimizers: error: Empty input file: " << ipath << '\n';
		exit(EXIT_FAILURE);
	}
	assert(bxtomxs.empty() && counts.empty()); 
	Bx bx;
	std::string mx;
	std::string mx_line;
	while ((is >> bx) && (getline(is, mx_line))) {
		bxtomxs[bx].reserve(boost::count(mx_line,' ')+1);
		std::istringstream iss(mx_line);
		while (iss >> mx) {
			Mx mx_num = strtoull(mx.c_str(), nullptr, 0);
			bxtomxs[bx].push_back(mx_num);
			if (counts.find(mx_num) == counts.end()) {
				counts[mx_num] = 1;
			} else {
				++counts[mx_num];
			}
		}
	}
	if (!silent) {
		auto t = std::chrono::steady_clock::now();
		auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0);
		std::cerr << "Time after reading (ms): " << diff.count() << '\n';
	}
}

static void filter_barcodes(const size_t n, const size_t N, bool silent, MxtoCount& counts, BxtoMxs& bxtomxs) {
	Mxs mxs_to_change; //to 'update' counts based on the mxs of the bxs that were removed	
	unsigned too_few = 0, too_many = 0;
	uint64_t singletons = 0;
    unsigned initial_size = bxtomxs.size();
	unsigned initial_count_size = counts.size();
	for (auto it = bxtomxs.begin(); it != bxtomxs.end(); ) {
		auto& bx = it->first;
		auto mxs = it->second; //this will be re-assigned so cannot use address(&)
		Mxs not_singletons;
        not_singletons.reserve(mxs.size());
		for (const auto& mx : mxs){
			if (counts[mx] >= 2) {
				not_singletons.push_back(mx);
			} else {
				++singletons;
			}
		}
        bxtomxs[bx] = std::move(not_singletons);
		mxs = it->second;
		if (mxs.size() < n) {
			++too_few;
			mxs_to_change.insert(mxs_to_change.end(), mxs.begin(), mxs.end());
			it = bxtomxs.erase(it);
		}
		else if (mxs.size() >= N) {
			++too_many;
			mxs_to_change.insert(mxs_to_change.end(), mxs.begin(), mxs.end());
			it = bxtomxs.erase(it);
		}
		else {
			++it;
		}
	}
	//bxtomxs has been updated, update counts by removing singleton counts AND UPDATING THE mxs_to_change COUNTS!!
	for(auto it = counts.begin(); it != counts.end(); ){
		auto& mx = it->first;
		if(counts[mx] == 1){
			it = counts.erase(it);
		}
		else{
			unsigned freq = std::count(mxs_to_change.begin(), mxs_to_change.end(), mx);
			if(freq != 0){
				counts[mx] = counts[mx] - freq;
			}
			++it;
		}
	}
	if (!silent) {
		auto t = std::chrono::steady_clock::now();
    	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0);
		std::cerr << "Time after filtering barcodes(ms): " << diff.count() << '\n';	
		std::cerr << "Removed " << singletons << " minimizers that occur once of " << initial_count_size << " (" << std::setprecision(1) << std::fixed << 100.0 * singletons / initial_count_size << "%)\n";
		std::cerr << "Discarded " << too_few << " barcodes with too few minimizers of " << initial_size << " (" << std::setprecision(1) << std::fixed << 100.0 * too_few / initial_size << "%)\n";
		std::cerr << "Discarded " << too_many << " barcodes with too many minimizers of " << initial_size << " (" << std::setprecision(1) << std::fixed << 100.0 * too_many / initial_size << "%)\n";
		std::cerr << "Wrote " << initial_size - too_few - too_many << " barcodes\n";
	}
}

static std::vector<float> quantile(std::vector<float> quantiles, std::vector<unsigned> values){
	std::sort(values.begin(), values.end());
	std::vector<float> qs;
	for (const auto& p : quantiles){
			qs.push_back(values[round(p*(values.size()-1))]);
	}
	return qs;
}

static void remove_singletons(BxtoMxs& bxtomxs, MxtoCount& counts){
	uint64_t singletons = 0;
	for (auto it = bxtomxs.begin(); it != bxtomxs.end(); ++it) {
		Mxs not_singletons;
		not_singletons.reserve(it->second.size());
		for (const auto& mx : it->second) {
				if (counts[mx] >= 2) {
						not_singletons.push_back(mx);
				} else {
						++singletons;
				}
		}
		bxtomxs[it->first] = std::move(not_singletons);
	}
	//update count
	for(auto it = counts.begin(); it != counts.end(); ){
		if(counts[it->first] == 1){
			it = counts.erase(it);
		}
		else{
			++it;
		}
	}
}

static void filter_minimizers(std::ostream& os, const std::string& opath, BxtoMxs bxtomxs, MxtoCount counts, unsigned C){	
	std::vector<unsigned> values;
	values.reserve(counts.size());
	boost::copy(counts | boost::adaptors::map_values, std::back_inserter(values));
	std::vector<float> q = {0.25,0.5,0.75};
	q = quantile(q, values);
	unsigned high_whisker = int(int(q[2]) + 1.5 * (int(q[2]) - int(q[0])));
	if(C == 0){ //means that it hasn't been set!!
    	C = high_whisker;
    }
	unsigned repetitives = 0;
	for(auto it = bxtomxs.begin(); it != bxtomxs.end(); ++it){
		auto& bx = it->first;
        auto mxs = it->second;
		Mxs not_repetitives;
        not_repetitives.reserve(mxs.size());
		for (const auto& mx : mxs){
			if (counts[mx] < C){
					not_repetitives.push_back(mx);
			}
			else{
					++repetitives;
			}
		}
		bxtomxs[bx] = std::move(not_repetitives);
        mxs = it->second;
		if(mxs.size() != 0){
			std::ostringstream oss;
	        oss << bx;
	        char sep = '\t';    
            for(const auto& mx : mxs){
         		oss << sep << mx;
            	sep = ' ';
            }
            oss << '\n';
            os << oss.str();
            assert_good(oss, opath);
		}
	}
	assert_good(os, opath);
	std::cerr << "Minimizer frequency: Q1=" << q[0] << " Q2=" << q[1] << " Q3=" << q[2] << "*(Q3-Q1)=" << high_whisker << "C= " << C << '\n';
	std::cerr << "Removed " << repetitives << " most frequent minimizers of" << counts.size() << '\n';
}

int main(int argc, char *argv[])
{
	t0 = std::chrono::steady_clock::now();
	auto progname = "physlr-filter-barcodes-minimizers";
	int c;
	int optindex = 0;
	static int help = 0;
	unsigned n = 0;
	unsigned N = 0;
	unsigned C = 0;
	bool silent = false;
	bool failed = false;
	bool n_set = false;
	bool N_set = false;
	char* end = nullptr;
	std::string outfile("/dev/stdout");
	static const struct option longopts[] = { { "help", no_argument, &help, 1 },{ nullptr, 0, nullptr, 0 } };
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
		case 'C':
			C = strtoul(optarg, &end, 10);
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
	for (auto &infile : infiles) {
		if (infile == "-") {
			infile = "/dev/stdin";
		}
		std::ifstream ifs(infile);
		assert_good(ifs, infile);
		MxtoCount counts; 
		BxtoMxs bxtomxs;
		read_minimizers(ifs, infile, silent, counts, bxtomxs);
		filter_barcodes(n, N , silent, counts, bxtomxs);
		remove_singletons(bxtomxs, counts);	
		filter_minimizers(ofs, outfile, bxtomxs, counts, C);
	}
	ofs.flush();
	assert_good(ofs, outfile);
	return 0;
}
