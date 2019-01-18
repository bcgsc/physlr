// Convert linked-reads to minimizers using ntHash-2.0.0.
// Usage:  physlr-indexlr -k K -w W [-v] file...
// Output: Each line of output is a barcode followed by a list of minimzers.

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
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
// NOTE: This snippet was copied over from abyss-dev/Common/StringUtil.h
template <size_t N>
static bool startsWith(const std::string& s, const char (&prefix)[N]) {
	size_t n = N - 1;
	return s.size() > n && equal(s.begin(), s.begin() + n, prefix);
}

// Hash the k-mers of a read using ntHash.
static std::vector<uint64_t> hashKmers(const std::string &readstr, const size_t k) {
    std::vector<uint64_t> hashes;
    hashes.reserve(readstr.size() - k + 1);
    for (ntHashIterator iter(readstr, 1, k); iter != ntHashIterator::end(); ++iter) {
        hashes.push_back((*iter)[0]);
    }
    return hashes;
}

// Function object implements '<=' comparison; passed to std::min_element().
template<typename T>
struct Less_or_equal {
    bool operator()(const T& a, const T& b) const { return a <= b; }
};

/*
Algorithm to find minimizers from a vector of hash values:
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
}
*/
// Find the minimizers of a vector of hash values using a sliding-window of size 'w'.
static std::vector<uint64_t> getMinimizers(const std::vector<uint64_t> &hashes, const unsigned w) {
    Less_or_equal<uint64_t> less_or_equal;
    std::vector<uint64_t> minimizers;
    minimizers.reserve(hashes.size() / w);
    int i = -1, prev = -1;
    auto firstIt = hashes.begin();
    auto minIt   = hashes.end();
    for (auto leftIt = firstIt; leftIt < hashes.end() - w + 1; ++leftIt) {
        auto rightIt = leftIt + w;
        if (i < leftIt - firstIt) {
            // Use of the comparison operator '<=' returns the right-most minimum.
            minIt = std::min_element(leftIt, rightIt, less_or_equal);
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
static inline void assert_good(const std::ios& stream) {
	if (!stream.good()) {
		std::cerr << "physlr-indexlr: error: " << std::strerror(errno) << '\n';
		exit(EXIT_FAILURE);
	}
}

// Print the barcode and minimzers of a read.
static void printMinimizedRead(const std::string &barcode, const std::vector<uint64_t> &minimizers) {
    std::cout << barcode;
    char sep = '\t';
    for (auto &m : minimizers) {
        std::cout << sep << m;
        sep = ' ';
    }
    std::cout << '\n';
    assert_good(std::cout);
}

// Read a FASTQ file and reduce each read to a set of minimizers
static void minimizeReads(std::istream& is, const size_t k, const size_t w, bool verbose) {
    // Check if input file is empty.
    if (is.peek() == std::ifstream::traits_type::eof()) {
            std::cerr << "physlr-indexlr: error: Empty input file\n";
            exit(EXIT_FAILURE);
    }
    size_t nread = 0, nline = 0;
    while (true) {
        std::string id, barcode, sequence;
        is >> id >> std::ws;
        if (is.eof()) {
            break;
        }
        if (!getline(is, barcode)) {
            std::cerr << "physlr-indexlr: error: Failed to read header on line " << nline + 1 << '\n';
            exit(EXIT_FAILURE);
        }
        if (!getline(is, sequence)) {
            std::cerr << "physlr-indexlr: error: Failed to read sequence on line " << nline + 2 << '\n';
            exit(EXIT_FAILURE);
        }
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // At this point four lines have been read, that is, one read and its associated information.
        nline += 4;
        nread += 1;
        if (!startsWith(barcode, "BX:Z:")) {
            std::cerr << "physlr-indexlr: error: Expected BX:Z:... and saw " << barcode << " at line "
                      << nline - 3 << "\n";
            exit(EXIT_FAILURE);
        }
        barcode = barcode.erase(0, 5);
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
        // of hashes ntHash returns after hashing a read. ntHash takes care to skip Ns.
        std::vector<uint64_t> hashes = hashKmers(sequence, k);
        if (w > hashes.size()) {
            if (verbose) {
                std::cerr << "physlr-indexlr: warning: Skip read " << nread << " on line " << nline - 2
                          << "; window size > #hashes (w = " << w << ", #hashes = " << hashes.size()
                          << ")\n";
            }
            continue;
        }
        // Minimerize the read, that is, pick the minimum hash values from the vector of hashes.
        std::vector<uint64_t> minimizers = getMinimizers(hashes, w);
        printMinimizedRead(barcode, minimizers);
    }
}

static void printErrorMsg(const std::string &progname, const std::string &msg) {
    std::cerr << progname << ": " << msg << "\nTry 'physlr-indexlr --help' for more information.\n";
}

static void printUsage(const std::string &progname) {
    std::cout << "Usage:  " << progname
         << "  -k K -w W [-v] file...\n\n"
            "  -k K     use K as k-mer size\n"
            "  -w W     use W as sliding-window size\n"
            "  -v       enable verbose output\n"
            "  --help   display this help and exit\n"
            "  file     space separated list of FASTQ files\n";
}

int main(int argc, char *argv[]) {
    std::string progname = "physlr-indexlr";
    int      c;
    int      optindex = 0;
    static int help   = 0;
    unsigned k        = 0;
    unsigned w        = 0;
    bool     verbose  = false;
    bool     failed   = false;
    bool     w_set    = false;
    bool     k_set    = false;
    char     *end;
    static const struct option longopts[] = {{"help", no_argument, &help, 1}, {nullptr, 0, nullptr, 0}};
    while ((c = getopt_long(argc, argv, "k:w:v", longopts, &optindex)) != -1) {
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
        case 'v':
            verbose = true;
            break;
        default:
            exit(EXIT_FAILURE);
        }
    }
    std::vector<std::string> infile(&argv[optind], &argv[argc]);
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
    } else if (infile.empty()) {
        printErrorMsg(progname, "missing file operand");
        failed = true;
    }
    if (failed) {
        exit(EXIT_FAILURE);
    }
    for (auto &f : infile) {
        std::ifstream ifs(f);
        if (!ifs) {
            std::cerr << "phslyr-indexlr: error: Failed to open: " << f << '\n';
            exit(EXIT_FAILURE);
        }
        minimizeReads(ifs, k, w, verbose);
    }
    return 0;
}
