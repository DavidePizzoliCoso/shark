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

#ifndef READANALYZER_HPP
#define READANALYZER_HPP

#include "bloomtree.h"
#include "kmer_utils.hpp"
#include <vector>
#include <array>

using namespace std;

class ReadAnalyzer {
public:
	typedef vector<assoc_t> output_t;

	ReadAnalyzer(SSBT *tree, uint _k, double _c, bool _only_single = false) :
	_tree(tree), k(_k), c(_c), only_single(_only_single) {}

	output_t* operator()(vector<elem_t> *reads) const 
	{
		output_t* associations = new output_t();
		
		map<string, int> counter_occurence;
		for(const auto & p : *reads) 
		{
			counter_occurence.clear();
			const string& read_seq = p.first; // Get Sequence ID
			unsigned int len = 0; // Get read length
			for (unsigned int pos = 0; pos < read_seq.size(); ++pos) 
				len += to_int[read_seq[pos]] > 0 ? 1 : 0;
			
			if(len >= k) 
			{
				int pos = 0;
				uint64_t kmer = build_kmer(read_seq, pos, k);
				if(kmer == (uint64_t)-1) continue;
				uint64_t rckmer = revcompl(kmer, k);
				
				auto genes = _tree->get_genes(min(kmer, rckmer));
				int n = genes.size();
				for(int i = 0; i < n; i++)
					counter_occurence[genes[i]]++;
				
				for (; pos < (int)read_seq.size(); ++pos) 
				{
					uint8_t new_char = to_int[read_seq[pos]];
					if(new_char == 0)  // Found a char different from A, C, G, T
					{
						++pos; // we skip this character then we build a new kmer
						kmer = build_kmer(read_seq, pos, k);
						if(kmer == (uint64_t)-1) break;
						rckmer = revcompl(kmer, k);
						--pos; // p must point to the ending position of the kmer, it will be incremented by the for
					} 
					else 
					{
						--new_char; // A is 1 but it should be 0
						kmer = lsappend(kmer, new_char, k);
						rckmer = rsprepend(rckmer, reverse_char(new_char), k);
					}
					
				
					genes = _tree->get_genes(min(kmer, rckmer));
					n = genes.size();
					for(int i = 0; i < n; i++)
						counter_occurence[genes[i]]++;
				}
			}
			
			vector<string> best_genes;
			int max = -1;
			for(auto it=counter_occurence.cbegin(); it!=counter_occurence.cend(); ++it)
			{
				if(it->second == max)
				{
					best_genes.push_back(it->first);
				}
				else if(it->second > max)
				{
					best_genes.clear();
					max = it->second;
					best_genes.push_back(it->first);
				}
			}
			if(max >= c*len && (!only_single || best_genes.size() == 1)) 
				for(const auto idx : best_genes) 
					associations->push_back({ idx, std::move(get<1>(p)) });
		}
		delete reads;
	
		if(associations->size())
			return associations;
		else 
		{
			delete associations;
			return NULL;
		}
	}

private:
	SSBT * const _tree;
	const uint k;
	const double c;
	const bool only_single;

};

#endif
