/*
 * physlr-unitig-overlap.cc
 * Author: cjustin
 */

#include "tsl/robin_map.h"
#include "tsl/robin_set.h"
#include "include/IOUtil.h"
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#if _OPENMP
#include <omp.h>
#endif

/*
 * Protocol:
 * Read in fai file of contigs
 * Record lengths of each contig to their ID
 * Read in sam file
 * Record molecule positions relative to barcode
 * contigID_array
 * Barcode->contigID
 * Run overlap code
 *
 * New edge metrics:
 * Strength of edges based on number of reads in common with overlap
 * Strength of edges based on number of unitigs in overlap
 *
 * Alternate version:
 * contigID_array
 * Barcode_contigID_regionStartEnd
 * Overlap code with regionStartEnd Check
 *
 */

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt
{
static unsigned threads = 1;
static unsigned minN = 2;
static std::string unitigFilename = "";
static unsigned verbose = 1;
static int version = 0;
static unsigned window = 2500;
static int mapQThresh = 1;
static std::string barcodeFile = "";
} // namespace opt

#define PROGRAM "physlr-unitig-overlap"
#define PACKAGE_NAME "physlr"
#define GIT_REVISION "pre-autotools" //TEMP

using BarcodeID = uint32_t;
using UnitigID = uint32_t;
using RegionID = uint32_t;
using Count = uint16_t;

static void printVersion()
{
	const char VERSION_MESSAGE[] =
	PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	"Written by Justin Chu.\n"
	"\n"
	"Copyright 2019 Canada's Michael Smith Genome Science Centre\n";
	std::cerr << VERSION_MESSAGE << std::endl;
	exit(EXIT_SUCCESS);
}

static void printHelpDialog()
{
	static const char dialog[] =
			"Usage: physlr-overlap [OPTION]... -u [UNITIG.fa.fai] [READS.bam]\n"
					"Read a sketch of linked reads and find overlapping barcodes.\n"
					"  -t, --threads=INT Number of threads to use [1].\n"
					"  -u, --unitigs=STR Unitig fasta index file (.fai) [Required].\n"
					"  -n, --min-n=INT   Remove region associations with fewer than n reads [2].\n"
					"  -w, --window      Window size for each region [2500].\n"
					"  -b, --barcode     List of whitelisted barcodes.\n"
					"  -m, --mapQ        Mapping Quality filtering threshold [1].\n"
					"  -v, --verbose     Print more information for debugging.\n"
					"      --version     Print version\n"
					"Report bugs to <cjustin@bcgsc.ca>.";
	std::cerr << dialog << std::endl;
	exit(0);
}

// invertible hash function
struct fastHash
{
	std::size_t operator()(uint64_t key) const
	{
		key = (~key + (key << 21u)); // key = (key << 21) - key - 1;
		key = key ^ key >> 24u;
		key = ((key + (key << 3u)) + (key << 8u)); // key * 265
		key = key ^ key >> 14u;
		key = ((key + (key << 2u)) + (key << 4u)); // key * 21
		key = key ^ key >> 28u;
		key = (key + (key << 31u));
		return key;
	}
};

template<size_t N>
static inline std::string parseSAMTag(const std::string& s,
		const char (&tag)[N])
{
	size_t start = s.find(tag);
	if (start == std::string::npos)
		return std::string();
	start += N - 1;

	// Find the next whitespace or EOL after "BX:Z:".
	size_t end = s.find_first_of(" \t\r\n", start);
	if (end == std::string::npos)
		end = s.length();

	return s.substr(start, end - start);
}

/**
 * Extract the first Chromium barcode sequence from a given string
 * (FASTQ comment or list of SAM tags). The barcode is expected
 * to match the format BX:Z:<BARCODE>.
 */
static inline std::string parseBXTag(const std::string& s)
{
	return parseSAMTag(s, "BX:Z:");
}

// returns memory of program in kb
static int memory_usage()
{
	int mem = 0;
	std::ifstream proc("/proc/self/status");
	std::string s;
	while (getline(proc, s), !proc.fail())
	{
		if (s.substr(0, 6) == "VmSize")
		{
			std::stringstream convert(
					s.substr(s.find_last_of('\t'), s.find_last_of('k') - 1));
			if (!(convert >> mem))
			{
				return 0;
			}
			return mem;
		}
	}
	return mem;
}

int main(int argc, char* argv[])
{
	bool die = false;

	// switch statement variable
	int c;

	// long form arguments
	static struct option long_options[] =
	{
	{ "min-n", required_argument, nullptr, 'n' },
	{ "unitigs", required_argument, nullptr, 'u' },
	{ "threads", required_argument, nullptr, 't' },
	{ "window", required_argument, nullptr, 'w' },
	{ "barcode", required_argument, nullptr, 'b' },
	{ "mapQ", required_argument, nullptr, 'm' },
	{ "verbose", no_argument, nullptr, 'v' },
	{ "version", no_argument, &opt::version, true },
	{ nullptr, 0, nullptr, 0 } };

	int i = 0;
	while ((c = getopt_long(argc, argv, "u:t:vn:w:b:m:", long_options, &i))
			!= -1)
	{
		switch (c)
		{
		case 't':
		{
			std::stringstream convert(optarg);
			if (!(convert >> opt::threads))
			{
				std::cerr << "Error - Invalid parameters! t: " << optarg
						<< std::endl;
				return 0;
			}
			break;
		}
		case 'u':
		{
			std::stringstream convert(optarg);
			if (!(convert >> opt::unitigFilename))
			{
				std::cerr << "Error - Invalid parameters! u: " << optarg
						<< std::endl;
				return 0;
			}
			break;
		}
		case 'n':
		{
			std::stringstream convert(optarg);
			if (!(convert >> opt::minN))
			{
				std::cerr << "Error - Invalid parameters! n: " << optarg
						<< std::endl;
				return 0;
			}
			break;
		}
		case 'b':
		{
			std::stringstream convert(optarg);
			if (!(convert >> opt::barcodeFile))
			{
				std::cerr << "Error - Invalid parameters! b: " << optarg
						<< std::endl;
				return 0;
			}
			break;
		}
		case 'm':
		{
			std::stringstream convert(optarg);
			if (!(convert >> opt::mapQThresh))
			{
				std::cerr << "Error - Invalid parameters! m: " << optarg
						<< std::endl;
				return 0;
			}
			break;
		}
		case 'w':
		{
			std::stringstream convert(optarg);
			if (!(convert >> opt::window))
			{
				std::cerr << "Error - Invalid parameters! w: " << optarg
						<< std::endl;
				return 0;
			}
			break;
		}
		case 'v':
		{
			opt::verbose++;
			break;
		}
		default:
		{
			die = true;
			break;
		}
		}
	}

	if (opt::version)
	{
		printVersion();
	}

#if _OPENMP
	omp_set_num_threads(opt::threads);
#endif

	// Threads currently not supported.
	if (opt::threads > 1)
	{
		std::cerr << "Error: physlr-unitig-overlap: Threads not yet supported.\n";
		die = true;
	}

	// Stores alignment fasta names
	std::vector<std::string> inputFiles;

	while (optind < argc)
	{
		inputFiles.emplace_back(argv[optind]);
		optind++;
	}

	if (inputFiles.empty())
	{
		std::cerr << "Missing Input Files" << std::endl;
		die = true;
	}

	if (die)
	{
		printHelpDialog();
		exit(EXIT_FAILURE);
	}

#if _OPENMP
	double sTime = omp_get_wtime();
#endif

	typedef std::tuple<std::string, size_t> Region;
	struct regionHash: std::unary_function<Region, std::size_t>
	{
		size_t operator()(Region const& e) const
		{
			return std::hash<std::string>()(std::get<0>(e)) ^ std::get<1>(e);
		}
	};

	std::vector<Region> regions;
	tsl::robin_map<Region, RegionID, regionHash> regionIDs;

	// Read in fai file of contigs
	{
		std::string faiBuffer;
		size_t lengthBuffer;
		std::ifstream fh;
		fh.open(opt::unitigFilename);
		std::string line;
		if (opt::verbose > 1)
			std::cerr << "Loading file " << opt::unitigFilename << std::endl;
		while (getline(fh, line))
		{

			//Record lengths of each contig to their ID
			std::stringstream ss(line);
			ss >> faiBuffer;
			ss >> lengthBuffer;
			size_t currentPos = 0;
			while (currentPos < lengthBuffer)
			{
				regions.emplace_back(Region(faiBuffer, currentPos));
				regionIDs[regions.back()] = regions.size() - 1;
				currentPos += opt::window;
			}
		}
	}

	// barcode to ID table (index in vector)
	tsl::robin_map<std::string, BarcodeID> barcodes;

	// vector of barcodes
	tsl::robin_map<RegionID, tsl::robin_map<BarcodeID, Count>> regionToBarcode;
//	tsl::robin_map<UnitigID, tsl::robin_map<BarcodeID, std::tuple<unsigned, unsigned>>> unitigToBarcode;
	std::vector<std::string> barcodeToStr;

	tsl::robin_set<std::string> barcodeWhitelist;

	// load in barcode whitelist
	if (opt::barcodeFile != "")
	{
		std::ifstream fh;
		fh.open(opt::barcodeFile.c_str());
		std::string line;
		while (getline(fh, line))
		{
			barcodeWhitelist.insert(line);
		}
	}

	// Read in sam file
	// Barcode->contigID

	std::string barcodeBuffer;
//	A00228:66:H3C7LDMXX:2:1375:10438:12085	163	1593868	14319	40	109S19M22S	=	14319	19
//	TTTAAGCATTAAAATACAAAATTAAATAAATTTTTGATCGAGGCTAAAAGCTTTAGATTTGATAATACAAAATTAGTAAAAGGATATTCAGAACAATAGGAAACGTCCGAAAAAAAATTATCTGGTGTGACACCAATGCAACTTTTTCTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF8-FFFFFFFFFFFFFF8FFFFFFFFFFF-FFFFFFFFF
//	NM:i:0	MD:Z:19	MC:Z:79S19M29S	AS:i:19	XS:i:0	BX:Z:AAACACCAGAAACCAT-1
	for (std::vector<std::string>::const_iterator itr = inputFiles.begin();
			itr != inputFiles.end(); itr++)
	{
		std::ifstream fh;
		fh.open(itr->c_str());
		std::string line;
		if (opt::verbose > 1)
			std::cerr << "Loading file " << *itr << std::endl;

		BarcodeID currentBX = 0;
		RegionID currentRegion = 0;
		unsigned readCount = 0;

		while (getline(fh, line))
		{

			std::stringstream ss(line);
			std::string readName, scafName, cigar, rnext, seq, qual, tags;
			int flag, pos, mapq, pnext, tlen;

			ss >> readName >> flag >> scafName >> pos >> mapq >> cigar >> rnext
					>> pnext >> tlen >> seq >> qual >> std::ws;

			//filter mapq
			if (mapq >= opt::mapQThresh)
			{
				getline(ss, tags);
				barcodeBuffer = parseBXTag(tags); /* Parse the index from the readName */
				const auto& barcode = barcodes.find(barcodeBuffer);
				if (barcodeWhitelist.empty()
						|| barcodeWhitelist.find(barcodeBuffer)
								!= barcodeWhitelist.end())
				{
					if (barcode == barcodes.end())
					{
						barcodeToStr.emplace_back(barcodeBuffer);
						barcodes[barcodeBuffer] = barcodeToStr.size() - 1;
					}
					if (currentRegion
							!= regionIDs[Region(scafName,
									pos - pos % opt::window)]
							|| currentBX != barcodes[barcodeBuffer])
					{
						if (readCount > opt::minN)
						{
							regionToBarcode[currentRegion][currentBX] =
									readCount;
						}
						currentRegion = regionIDs[Region(scafName,
								pos - pos % opt::window)];
						currentBX = barcodes[barcodeBuffer];
						readCount = 0;
					}
					readCount++;
				}
			}
		}
	}

	//Streaming Algo proposition (position sorted):
	//Record current position, for each position updated update to current position
	//Keep list of current positions for each barcode, in addition to barcode counts
	//When barcode goes out of scope store barcode association with associated number of reads into edge graph
	//(optional) If edge is greater than threshold, print out edge (heurisitic)
	//Counts of number of unitigs associated, Counts of number of bases associated, Counts of number of reads associated?
	//Indicate number of maximum number of possible molecules associated

	//BX sorted bam file
	//Molecule extents -> overlap

	/*
	 * Proposed bin based overlap
	 * 1000bp regions as a unitig
	 */

// * New edge metrics:
// * Strength of edges based on number of reads in common with overlap
// * Strength of edges based on number of unitigs in overlap
// *
// * Alternate version:
// * contigID_array
// * Barcode_contigID_regionStartEnd
// * Overlap code with regionStartEnd Check
#if _OPENMP
//	if(opt::verbose) {
	std::cerr << "Finished constructing contigsToBarcodes in sec: " << omp_get_wtime() - sTime
	<< std::endl;
	sTime = omp_get_wtime();
//	}
#endif
	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576)
			<< "GB" << std::endl;

	// store into 2d matrix / hash table
	// todo revisit Counts? -> can be smaller
	using SimMat = tsl::robin_map<uint64_t, std::pair<Count, Count>, fastHash>;
	SimMat barcodeSimMat;

	// counts of vector
	std::vector<Count> barcodeCount(barcodes.size(), 0);

	std::cerr << "Populating Overlaps" << std::endl;
	std::cerr << "Total Unitigs: " << regionToBarcode.size() << std::endl;
	std::cerr << "Total Barcodes: " << barcodes.size() << std::endl;

	for (const auto& itr : regionToBarcode)
	{
		for (auto barcodeCount_i = itr.second.begin();
				barcodeCount_i != itr.second.end(); barcodeCount_i++)
		{
			barcodeCount[barcodeCount_i->first] += barcodeCount_i->second;
			auto barcodeCount_j =
					tsl::robin_map<BarcodeID, Count>::const_iterator(
							barcodeCount_i);
			for (barcodeCount_j++; barcodeCount_j != itr.second.end();
					barcodeCount_j++)
			{
				// assign "canonical edge"
				if (barcodeCount_i->first > barcodeCount_j->first)
				{
					uint64_t id = (static_cast<size_t>(barcodeCount_i->first)
							<< 32u | barcodeCount_j->first);
					barcodeSimMat[id].first += barcodeCount_i->second;
					barcodeSimMat[id].second += barcodeCount_j->second;
				}
				else
				{
					uint64_t id = (static_cast<size_t>(barcodeCount_j->first)
							<< 32u | barcodeCount_i->first);
					barcodeSimMat[id].first += barcodeCount_j->second;
					barcodeSimMat[id].second += barcodeCount_i->second;
				}
				std::cout << barcodeToStr[barcodeCount_j->first] << " "
						<< barcodeToStr[barcodeCount_j->first] << " "
						<< std::get<0>(regions[itr.first]) << " "
						<< std::get<1>(regions[itr.first]) << std::endl;
			}
		}
	}

#if _OPENMP
	std::cerr << "Finished computing overlaps in sec: " << omp_get_wtime() - sTime << std::endl;
#endif
	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576)
			<< "GB" << std::endl;
	std::cerr << "Total number of unfiltered edges: " << barcodeSimMat.size()
			<< std::endl;

	std::cout << "U\tn\n";
	std::string bufferString;
	// print out vertexes + counts
#pragma omp parallel private(bufferString)
	for (const auto& itr : barcodes)
	{
//		if (barcodeCount[itr.second] >= opt::minN)
//		{
			bufferString.clear();
			bufferString += itr.first;
			bufferString += "\t";
			bufferString += std::to_string(barcodeCount[itr.second]);
			bufferString += "\n";
			std::cout << bufferString;
//		}
	}

	size_t edgeCount = 0;
	std::cout << "\nU\tV\tn\n";
#pragma omp parallel private(bufferString)
	for (const auto& itr : barcodeSimMat)
	{
//		if (sqrt(itr.second.first * itr.second.second) >= opt::minN)
//		{
			bufferString.clear();
			bufferString += barcodeToStr[itr.first >> 32u];
			bufferString += "\t";
			bufferString += barcodeToStr[static_cast<uint32_t>(itr.first)];
			bufferString += "\t";
			bufferString += std::to_string(
					unsigned(sqrt(itr.second.first * itr.second.second) + 0.5));
			bufferString += "\n";
			std::cout << bufferString;
			++edgeCount;
//		}
	}
	std::cerr << "Total number of edges: " << edgeCount << std::endl;

	return 0;
}
