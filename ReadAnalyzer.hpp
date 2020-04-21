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

	ReadAnalyzer(SSBT *tree, uint _k, double _c, bool _only_single = false, std::string _method = "base") :
	_tree(tree), k(_k), c(_c), only_single(_only_single), method(_method) {}

	output_t* operator()(vector<elem_t> *reads) const 
	{
		output_t* associations = new output_t();
		vector<string> best_genes;
		typedef pair<pair<unsigned int, unsigned int>, unsigned int> gene_cov_t;
		map<string, gene_cov_t> classification_id;
		
		for(const auto & p : *reads) 
		{
			classification_id.clear();
			const string& read_seq = p.first;
			unsigned int len = 0;
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
				for(int i=0; i < n; i++)
				{
					auto& gene_cov = classification_id[genes[i]];
					//cout<<genes[i]<<" "; // <<gene_cov.first.first<<" + min("<<k<<","<<pos-gene_cov.second<<") - ";
					gene_cov.first.first += min(k, pos - gene_cov.second);
					gene_cov.first.second = 1;
					gene_cov.second = pos - 1;
					//cout<<gene_cov.first.first<<" "<<gene_cov.first.second<<" "<<gene_cov.second<<endl;
				}
				
				for (; pos < (int)read_seq.size(); ++pos) 
				{
					uint8_t new_char = to_int[read_seq[pos]];
					if(new_char == 0)
					{
						++pos;
						kmer = build_kmer(read_seq, pos, k);
						if(kmer == (uint64_t)-1) break;
						rckmer = revcompl(kmer, k);
						--pos;
					} 
					else 
					{
						--new_char;
						kmer = lsappend(kmer, new_char, k);
						rckmer = rsprepend(rckmer, reverse_char(new_char), k);
					}
					
					genes = _tree->get_genes(min(kmer, rckmer));
					int n = genes.size();
					for(int i=0; i < n; i++)
					{
						auto& gene_cov = classification_id[genes[i]];
						//cout<<genes[i]<<" "; // <<gene_cov.first.first<<" + min("<<k<<","<<pos-gene_cov.second<<") /-> ";
						gene_cov.first.first += min(k, pos - gene_cov.second);
						gene_cov.first.second += 1;
						gene_cov.second = pos;
						//cout<<gene_cov.first.first<<" "<<gene_cov.first.second<<" "<<gene_cov.second<<endl;
					}
				}
			}
			
			//cout<<endl;
			unsigned int maxk = 0;
			unsigned int max = 0;
			best_genes.clear();
			//cout<<method<<endl;
			if(method == "kmer")
			{
				for(auto it=classification_id.cbegin(); it!=classification_id.cend(); ++it) 
				{
					//cout<<it->first<<" "<<it->second.first.first<<" "<<it->second.first.second<<" "<<it->second.second<<endl;
					if(it->second.first.second == maxk) 
					{
						best_genes.push_back(it->first);
					} 
					else if(it->second.first.second > maxk)
					{
						best_genes.clear();
						maxk = it->second.first.second;
						best_genes.push_back(it->first);
					}
				}
				//cout<<"maxk: "<<maxk<<endl;
				if(maxk >= c*(len-k+1) && (!only_single || best_genes.size() == 1)) 
					for(const auto idx : best_genes) 
						associations->push_back({ idx, std::move(get<1>(p)) });
			}
			else
			{
				for(auto it=classification_id.cbegin(); it!=classification_id.cend(); ++it) 
				{
					//cout<<it->first<<" "<<it->second.first.first<<" "<<it->second.first.second<<" "<<it->second.second<<endl;
					if(it->second.first.first == max && it->second.first.second == maxk) 
					{
						best_genes.push_back(it->first);
					} 
					else if(it->second.first.first > max || (it->second.first.first == max && it->second.first.second > maxk)) 
					{
						best_genes.clear();
						max = it->second.first.first;
						maxk = it->second.first.second;
						best_genes.push_back(it->first);
					}
				}
				//cout<<"maxk: "<<maxk<<" max: "<<max<<endl;
				if(max >= c*len && (!only_single || best_genes.size() == 1)) 
					for(const auto idx : best_genes) 
						associations->push_back({ idx, std::move(get<1>(p)) });
			}
			
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
	const std::string method;

};

#endif
