/**
 * shark - Mapping-free filtering of useless RNA-Seq reads
 * Copyright (C) 2019 Tamara Ceccato, Luca Denti, Yuri Pirola, Marco Previtali
 *
 * This file is part of shark.
 *
 * shark is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * shark is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with shark; see the file LICENSE. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

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

	// Sample 1
	gzFile read1_file = gzopen(opt::sample1_path.c_str(), "r");
	seq = kseq_init(read1_file);
	kseq_destroy(seq);
	gzclose(read1_file);

	// Sample 2
	gzFile read2_file = nullptr;
	if(opt::paired_flag) 
	{
		read2_file = gzopen(opt::sample2_path.c_str(), "r");
		seq = kseq_init(read2_file);
		kseq_destroy(seq);
		gzclose(read2_file);
	}

	if(opt::verbose) 
	{
		cerr << "Reference texts: " << opt::fasta_path << endl;
		cerr << "Sample 1: " << opt::sample1_path << endl;
		if(opt::paired_flag)
			cerr << "Sample 2: " << opt::sample2_path << endl;
		cerr << "K-mer length: " << opt::k << endl;
		cerr << "Threshold value: " << opt::c << endl;
		cerr << "Only single associations: " << (opt::single ? "Yes" : "No") << endl;
		cerr << "Minimum base quality: " << static_cast<int>(opt::min_quality) << endl;
		cerr << endl;
	}
	
	/****************************************************************************/

	/*** 1. Second iteration over transcripts ***********************************/
  
	SSBT tree(opt::bf_size);
	vector<string> legend_ID;
	ref_file = gzopen(opt::fasta_path.c_str(), "r");
	seq = kseq_init(ref_file);
	int nidx = 0, seq_len;
	
	while ((seq_len = kseq_read(seq)) >= 0) 
	{
		string input_name = seq->name.s;
		legend_ID.push_back(input_name);
		++nidx;
	}
	kseq_destroy(seq);
	gzclose(ref_file);
	
	/****************************************************************************/
	
	SimpleBF* bloom;
	int i = 0, counter = 0;
	dequeue<SimpleBF*> coda;
	coda.clear();
	
	for(; i < nidx; i++) 
	{
		bloom = new SimpleBF(sbt->_size, i, nHash);
		coda.push_back(bloom);
	}
	vector<SimpleBF*> leaves(coda);
		
	if(diffSizes)
	{
		while (coda.size() > 1)
		{
			SimpleBF* sx = coda.front();
			coda.pop_front();
			SimpleBF* dx = coda.front();
			coda.pop_front();
			
			SimpleBF* node = new SimpleBF(max(sx->_size, dx->_size) << 1, -1, nHash);
			node->setSxChild(sx);
			node->setDxChild(dx);
			sx->parent = node;
			dx->parent = node;
			
			sx->support = (node->_size >> 1 != sx->_size);
			dx->support = (node->_size >> 1 != dx->_size);
			
			coda.push_back(make_pair(node, indexes));
		}
		tree.setRoot(coda.front().first);
	}
	else
	{
		while (coda.size() > 1)
		{
			SimpleBF* sx = coda.front();
			coda.pop_front();
			SimpleBF* dx = coda.front();
			coda.pop_front();
			
			SimpleBF* node = new SimpleBF(sbt->_size, -1, nHash);
			node->setSxChild(sx);
			node->setDxChild(dx);
			node->setBF(sx, dx);
			sx->parent = node;
			dx->parent = node;
			
			coda.push_back(node);
		}
		tree.setRoot(coda.front());
	}
	
	pelapsed("BF created from transcripts (" + to_string(nidx) + " genes)");
	
	/****************************************************************************/
	
	/*** 2. First iteration over transcripts ************************************/
	
	{
		ref_file = gzopen(opt::fasta_path.c_str(), "r");
		kseq_t *refseq = kseq_init(ref_file);
		
		tbb::filter_t<void, vector<pair<string, string>>*> 
			tr(tbb::filter::serial_in_order, FastaSplitter(refseq, 100));
		tbb::filter_t<vector<pair<string, string>>*, vector<pair<string,vector<size_t>>>*> 
			kb(tbb::filter::parallel, KmerBuilder(opt::k, opt::bf_size, opt::nHash));
		tbb::filter_t<vector<pair<string,vector<size_t>>>*, void> 
			bff(tbb::filter::serial_out_of_order, BloomfilterFiller(&tree, opt::nHash, opt::diff_sizes, &counter, &leaves));

		tbb::filter_t<void, void> pipeline = tr & kb & bff;
		tbb::parallel_pipeline(opt::nThreads, pipeline);
		
		kseq_destroy(refseq);
		gzclose(ref_file);
	}
	
	pelapsed("Transcript file processed");
	
	/****************************************************************************/

	/*** 3. Iteration over the sample *****************************************/
	// IF (FASE 1) COMMENT FROM HERE
	
	{
		kseq_t *sseq1 = nullptr, *sseq2 = nullptr;
		FILE *out1 = nullptr, *out2 = nullptr;
		read1_file = gzopen(opt::sample1_path.c_str(), "r");
		sseq1 = kseq_init(read1_file);
		if (opt::out1_path != "")
			out1 = fopen(opt::out1_path.c_str(), "w");
		if(opt::paired_flag) 
		{
			read2_file = gzopen(opt::sample2_path.c_str(), "r");
			sseq2 = kseq_init(read2_file);
			if (opt::out2_path != "")
				out2 = fopen(opt::out2_path.c_str(), "w");
		}
		
		tbb::filter_t<void, FastqSplitter::output_t*> 
			sr(tbb::filter::serial_in_order, FastqSplitter(sseq1, sseq2, 50000, opt::min_quality, out1 != nullptr));
		tbb::filter_t<FastqSplitter::output_t*, ReadAnalyzer::output_t*> 
			ra(tbb::filter::parallel, ReadAnalyzer(&tree, legend_ID, opt::k, opt::c, opt::single, opt::method, opt::nHash, opt::diff_sizes));
		tbb::filter_t<ReadAnalyzer::output_t*, void> 
			so(tbb::filter::serial_in_order, ReadOutput(out1, out2));

		tbb::filter_t<void, void> pipeline_reads = sr & ra & so;
		tbb::parallel_pipeline(opt::nThreads, pipeline_reads);

		kseq_destroy(sseq1);
		gzclose(read1_file);
		if(opt::paired_flag) 
		{
			kseq_destroy(sseq2);
			gzclose(read2_file);
		}
		if (out1 != nullptr) fclose(out1);
		if (out2 != nullptr) fclose(out2);
	}
	
	pelapsed("Sample completed");
	
	// IF (FASE 1) COMMENT UNTIL HERE
	/****************************************************************************/
	
	pelapsed("Association done");
	
	return 0;
}
