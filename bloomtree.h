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
 
#ifndef _BLOOM_TREE_HPP
#define _BLOOM_TREE_HPP

#include <algorithm>
#include <array>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <string>
#include <list>
#include "simpleBF.h"

#include <sys/mman.h>

#include "kmer_utils.hpp"

using namespace std;
using namespace sdsl;

class KmerBuilder;
class BloomfilterFiller;

class SSBT {
	friend class KmerBuilder;
	friend class BloomfilterFiller;

	public:
		typedef uint64_t kmer_t;
		typedef bit_vector bit_vector_t;
		
		SSBT(const size_t size) : _size(size) {}

		~SSBT()
		{
			delete root;
		}
		
		void setRoot(SimpleBF* node)
		{
			root = node;
		}
		
		vector<string> get_genes(const kmer_t &kmer) const 
		{
			list<SimpleBF*> coda;
			vector<string> genes;
			vector<uint64_t> hash = _get_hash(kmer);
			vector<size_t> bf_idx;
			
			bf_idx.clear();
			for(const auto h : hash)
				bf_idx.push_back(h % _size);
			
			coda.push_back(root);
			genes.clear();
			//int counter = 0;
			
			while(coda.size() > 0)
			{
				SimpleBF* node = coda.front();
				coda.pop_front();
				//counter++;
				bool flag = true;
				for(const auto index : bf_idx)
				{
					//cout<<node->_bf[bf_idx.back()]<<" "<<bf_idx.back()<<endl;
					flag = node->_bf[index] == 1;
					if(!flag) break;
				}
				
				if (flag)
				{
					if(node->sx == nullptr && node->dx == nullptr) // This is a leaf
						genes.push_back(node->_id);
					else // This is a node
					{
						if(node->sx != nullptr)
							coda.push_back(node->sx);
						
						if(node->dx != nullptr)
							coda.push_back(node->dx);
					}
				}
			}
			//if(counter > 1)
			//	cout<<kmer<<";"<<bf_idx<<";"<<counter<<endl;
			return genes;
		}
		
		void printTree()
		{
			list<SimpleBF*> coda;
			
			coda.push_back(root);
			
			while(coda.size() > 0)
			{
				SimpleBF* node = coda.front();
				coda.pop_front();
				node->printBF(cout);
				
				if(node->sx != nullptr)
					coda.push_back(node->sx);
				
				if(node->dx != nullptr)
					coda.push_back(node->dx);
			}
			
		}
		
		SSBT() = delete;
		const SSBT &operator=(const SSBT &) = delete;
		const SSBT &operator=(const SSBT &&) = delete;
		
		SimpleBF* root;
		size_t _size;
		
};

#endif
