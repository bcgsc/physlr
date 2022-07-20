#include "btllib/seq_reader.hpp"
#include "btllib/status.hpp"
#include "btllib/util.hpp"

#include <cassert>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <math.h>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

const static std::string PROGNAME = "add-bx-to-long";
const static std::string VERSION = "v0.1.0";
const static size_t MAX_THREADS = 6;

static void
print_error_msg(const std::string& msg)
{
	std::cerr << PROGNAME << ' ' << VERSION << ": " << msg << std::endl;
}

static void
print_usage()
{
	std::cerr
	    << "Usage: Add BX:Z: style barcodes to long reads. Note that each long read will be assigned a unique barcode." << PROGNAME
	    << " [--fasta -m M -t T -b B] READS "
	       "\n\n"
	       "  --fasta     Output in fasta format.\n"
	       "  -m M        M minimum read length for a read to be considered a molecule. [2000].\n"
	       "  -t T        Use T number of threads (max 6) per input file. [6]\n"
	       
	       "  -v          Show verbose output.\n"
	       "  --help      Display this help and exit.\n"
	       "  --version   Display version and exit.\n"
	       "  READS       Space separated list of long reads FASTA/Q files to be tagged with barcodes."
	    << std::endl;
}

int
main(int argc, char* argv[])
{
	int c;
	int optindex = 0;
	static int help = 0, version = 0;
	size_t t = 6, m = 2000;
	
	std::vector<size_t> read_lengths;
	static int with_fasta = 0;
	//std::string configFile("tigmint-long.params.tsv");
	//std::string bxMultiplicityFile("barcode_multiplicity.tsv");
	bool failed = false;
	static const struct option longopts[] = {
		{ "fasta", no_argument, &with_fasta, 1 },
		{ "help", no_argument, &help, 1 },
		{ "version", no_argument, &version, 1 },
		{ nullptr, 0, nullptr, 0 }
	};
	while ((c = getopt_long(argc, argv, "t:m:", longopts, &optindex)) != -1) {
		switch (c) {
		case 0:
			break;
		case 'm':
			m = std::stoul(optarg);
			break;
		case 't':
			t = std::stoul(optarg);
		default:
			std::exit(EXIT_FAILURE);
		}
	}

	std::vector<std::string> infiles(&argv[optind], &argv[argc]);
	if (argc < 2) {
		print_usage();
		std::exit(EXIT_FAILURE);
	}
	if (help != 0) {
		print_usage();
		std::exit(EXIT_SUCCESS);
	} else if (version != 0) {
		std::cerr << PROGNAME << ' ' << VERSION << std::endl;
		std::exit(EXIT_SUCCESS);
	}
	if (infiles.empty()) {
		print_error_msg("missing file operand");
		failed = true;
	}
	if (failed) {
		std::cerr << "Try '" << PROGNAME << " --help' for more information.\n";
		std::exit(EXIT_FAILURE);
	}

	if (t > MAX_THREADS) {
		t = MAX_THREADS;
		std::cerr << (PROGNAME + ' ' + VERSION + ": Using more than " +
		              std::to_string(MAX_THREADS) + " threads does not scale, reverting to " +
		              std::to_string(MAX_THREADS) + ".\n")
		          << std::flush;
	}

	char header_symbol = '@';

	if (with_fasta) {
		header_symbol = '>';
	}

	for (auto& infile : infiles) {
		unsigned flags = 0;
		flags |= btllib::SeqReader::Flag::LONG_MODE;
		btllib::SeqReader reader(infile, flags, t);
		btllib::SeqReader::Record record;
		while ((record = reader.read())) {
			std::string& seq = record.seq;
			size_t seq_size = seq.size();
			std::string& qual = record.qual;
      size_t qual_size = qual.size();

      std::cout << header_symbol << record.id << " BX:Z:" << record.num + 1 << '\n';
      std::cout << seq << '\n';
      if (!with_fasta) {
						std::cout << "+\n";
						if (qual_size == 0) {
							std::cout << std::string(seq_size, '#') << '\n';
						} else {
							std::cout << qual << '\n';
            }
      }
      std::cout << std::flush;
    }
	}
	return 0;
}
