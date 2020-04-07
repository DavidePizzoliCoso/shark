#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <unordered_set>

#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#include "sdsl/int_vector.hpp"
#include "sdsl/util.hpp"

#include "common.hpp"
#include "argument_parser.hpp"
#include "simpleBF.h"
#include "bloomtree.h"
#include "BloomfilterFiller.hpp"
#include "KmerBuilder.hpp"
#include "FastaSplitter.hpp"
#include "FastqSplitter.hpp"
#include "ReadAnalyzer.hpp"
#include "ReadOutput.hpp"
#include "kmer_utils.hpp"

#include <fstream>

using namespace std;

auto start_t = chrono::high_resolution_clock::now();

void pelapsed(const string &s = "") {
	auto now_t = chrono::high_resolution_clock::now();
	cerr << "[shark/" << s << "] Time elapsed "<< chrono::duration_cast<chrono::milliseconds>(now_t - start_t).count()/1000<< endl;
}

/*****************************************
 * Main
 *****************************************/
int main(int argc, char *argv[]) {
	parse_arguments(argc, argv);

	// Transcripts
	gzFile ref_file = gzopen(opt::fasta_path.c_str(), "r");
	kseq_t *seq = kseq_init(ref_file);
	kseq_destroy(seq);
	gzclose(ref_file);
	
	SSBT tree(opt::bf_size);
	vector<string> legend_ID;
	
	{
	ref_file = gzopen(opt::fasta_path.c_str(), "r");
	kseq_t *refseq = kseq_init(ref_file);
	tbb::filter_t<void, vector<pair<string, string>>*> tr(tbb::filter::serial_in_order, FastaSplitter(refseq, 100));
		
	tbb::filter_t<vector<pair<string, string>>*, vector<pair<string,vector<uint64_t>>>*> kb(tbb::filter::parallel, KmerBuilder(opt::k));
		
	tbb::filter_t<vector<pair<string,vector<uint64_t>>>*, void> bff(tbb::filter::serial_out_of_order, BloomfilterFiller(&tree));

	tbb::filter_t<void, void> pipeline = tr & kb & bff;
	tbb::parallel_pipeline(opt::nThreads, pipeline);

	kseq_destroy(refseq);
	gzclose(ref_file);
	}
	
	return 0;
}
