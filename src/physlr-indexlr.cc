// Convert linked-reads to minimizers using ntHash-2.0.0.
// Usage:  physlr-indexlr -k K -w W [-v] [-o file] file...
// Output: Each line of output is a barcode followed by a list of minimzers.

#include <algorithm>
#include <cassert>
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

// ntHash 2.0.0
#include "ntHashIterator.h"
#include "nthash.h"

// Return true if the second string is a prefix of the string s.
template <size_t N>
static bool startsWith(const std::string& s, const char (&prefix)[N])
{
	auto n = N - 1;
	return s.size() > n && equal(s.begin(), s.begin() + n, prefix);
}

// Hash the k-mers of a read using ntHash.
static std::vector<uint64_t> hashKmers(const std::string &readstr, const size_t k)
{
    std::vector<uint64_t> hashes;
    hashes.reserve(readstr.size() - k + 1);
    for (ntHashIterator iter(readstr, 1, k); iter != ntHashIterator::end(); ++iter) {
        hashes.push_back((*iter)[0]);
    }
    return hashes;
}

// Minimerize a sequence: Find the minimizers of a vector of hash values representing a sequence.
/* Algorithm
v is a vector of non-negative integers
w is the window size
Invariants
    0 <  w <= v.size() - 1
    0 <= l <= r <= v.size() - 1
Initial conditions
    M    = NIL       Final set of minimizers, empty initially
    min  = -1        Minimum element
    i    = -1        Index of minimum element
    prev = -1        Index of previous minimum element
    l    = 0         Index of left end of window
    r    = l + w - 1 Index of right end of window
Computation
At each window, if the previous minimum is out of scope, find the new, right-most, minimum
or else, check with only the right-most element to determine if that is the new minimum.
A minimizer is added to the final vector only if it's index has changed.
for each window of v bounded by [l, r]
    if (i < l)
        i = index of minimum element in [l, r], furthest from l.
    else if (v[r] <= v[i])
        i = r
    min = v[i]
    if (i != prev) {
        prev = i
        M <- M + m
    }
    l = l + 1        Move window's left bound by one element
    r = l + w - 1    Set window's right bound
}*/
static std::vector<uint64_t> getMinimizers(const std::vector<uint64_t> &hashes, const unsigned w)
{
    std::vector<uint64_t> minimizers;
    minimizers.reserve(hashes.size() / w);
    int i = -1, prev = -1;
    auto firstIt = hashes.begin();
    auto minIt   = hashes.end();
    for (auto leftIt = firstIt; leftIt < hashes.end() - w + 1; ++leftIt) {
        auto rightIt = leftIt + w;
        if (i < leftIt - firstIt) {
            // Use of operator '<=' returns the minimum that is furthest from left.
            minIt = std::min_element(leftIt, rightIt, std::less_equal<uint64_t>());
        }
        else if (*(rightIt - 1) <= *minIt) {
            minIt = rightIt - 1;
        }
        i = minIt - firstIt;
        if (i > prev) {
            prev = i;
            minimizers.push_back(*minIt);
        }
    }
    return minimizers;
}

// Test the condition of a I/O stream.
static inline void assert_good(const std::ios& stream, const std::string& path)
{
	if (!stream.good()) {
		std::cerr << "error: " << strerror(errno) <<  ": " << path << '\n';
		exit(EXIT_FAILURE);
	}
}

// Print the barcode and minimzers of a read.
static void writeMinimizedRead(std::ostream &os, const std::string &opath, const std::string &barcode, const std::vector<uint64_t> &minimizers)
{
    os << barcode;
    char sep = '\t';
    for (auto &m : minimizers) {
        os << sep << m;
        sep = ' ';
    }
    os << '\n';
    assert_good(os, opath);
}

// Read a FASTQ file and reduce each read to a set of minimizers
static void minimizeReads(std::istream &is, const std::string &ipath, std::ostream &os, const std::string &opath, const size_t k, const size_t w, bool verbose)
{
    // Check if input file is empty.
    if (is.peek() == std::ifstream::traits_type::eof()) {
            std::cerr << "physlr-indexlr: error: Empty input file: " << ipath << '\n';
            exit(EXIT_FAILURE);
    }
    size_t nread = 0, nline = 0;
    for (std::string id, barcode, sequence; is >> id;) {
        while (is.peek() == ' ') {
            is.ignore();
        }
        if (!getline(is, barcode)) {
            std::cerr << "physlr-indexlr: error: Failed to read header on line " << nline + 1 << '\n';
            exit(EXIT_FAILURE);
        }
        if (!getline(is, sequence)) {
            std::cerr << "physlr-indexlr: error: Failed to read sequence on line " << nline + 2 << '\n';
            exit(EXIT_FAILURE);
        }
        assert(!id.empty());
        if (id[0] == '@') {
            // Skip the FASTQ quality.
            if (is.peek() != '+') {
                std::cerr << "physlr-indexlr: error: " << nline + 3 << ": Expected + and saw: " << id << '\n';
                exit(EXIT_FAILURE);
            }
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            nline += 2;
        } else if (id[0] != '>') {
            std::cerr << "physlr-indexlr: error: " << nline + 1 << ": Expected > or @ and saw: " << id << '\n';
            exit(EXIT_FAILURE);
        }
        nline += 2;
        nread += 1;
        if (startsWith(barcode, "BX:Z:")) {
            auto pos = barcode.find(' ');
            if (pos != std::string::npos) {
                barcode.erase(pos);
            }
            barcode.erase(0, 5);
        } else {
            barcode = "NA";
        }
        // Validate parameters.
        if (k > sequence.size()) {
            if (verbose) {
                std::cerr << "physlr-indexlr: warning: Skip read " << nread << " on line " << nline - 2
                          << "; k > read length "
                          << "(k = " << k << ", read length = " << sequence.size() << ")\n";
            }
            continue;
        }
        // Hash the kmers.
        // NOTE: The predicate P(#kmers != #hashes) will be true when reads contains Ns, so check with the number
        // of hashes ntHash returns after hashing a sequence. ntHash's iterator takes care to skip Ns.
        auto hashes = hashKmers(sequence, k);
        if (w > hashes.size()) {
            if (verbose) {
                std::cerr << "physlr-indexlr: warning: Skip read " << nread << " on line " << nline - 2
                          << "; window size > #hashes (w = " << w << ", #hashes = " << hashes.size() << ")\n";
            }
            continue;
        }
        writeMinimizedRead(os, opath, barcode, getMinimizers(hashes, w));
    }
}

static void printErrorMsg(const std::string &progname, const std::string &msg)
{
    std::cerr << progname << ": " << msg << "\nTry 'physlr-indexlr --help' for more information.\n";
}

static void printUsage(const std::string &progname)
{
    std::cout << "Usage:  " << progname
         << "  -k K -w W [-v] [-o file] file...\n\n"
            "  -k K       use K as k-mer size\n"
            "  -w W       use W as sliding-window size\n"
            "  -v         enable verbose output\n"
            "  -o file    write output to file, default is stdout\n"
            "  --help     display this help and exit\n"
            "  file       space separated list of FASTQ files\n";
}

int main(int argc, char *argv[])
{
    auto progname = "physlr-indexlr";
    int      c;
    int      optindex = 0;
    static int help   = 0;
    unsigned k        = 0;
    unsigned w        = 0;
    bool     verbose  = false;
    bool     failed   = false;
    bool     w_set    = false;
    bool     k_set    = false;
    char     *end     = nullptr;
    std::string outfile("/dev/stdout");
    static const struct option longopts[] = {{"help", no_argument, &help, 1}, {nullptr, 0, nullptr, 0}};
    while ((c = getopt_long(argc, argv, "k:w:o:v", longopts, &optindex)) != -1) {
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
        minimizeReads(ifs, infile, ofs, outfile, k, w, verbose);
    }
    ofs.flush();
    assert_good(ofs, outfile);
    return 0;
}
