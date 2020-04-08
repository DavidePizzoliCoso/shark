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

#ifndef _BLOOM_FILTER_HPP
#define _BLOOM_FILTER_HPP

#include <algorithm>
#include <array>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <string>

#include <sys/mman.h>

#include "kmer_utils.hpp"

using namespace std;
using namespace sdsl;

class KmerBuilder;
class BloomfilterFiller;
class SSBT;

#ifndef SHARK_HUGEPAGESIZE
#define SHARK_HUGEPAGESIZE (2 * 1024 * 1024)
#endif

class SimpleBF
{
	friend class KmerBuilder;
	friend class BloomfilterFiller;
	friend class SSBT;
	
	public:
		typedef uint64_t kmer_t;
		typedef uint64_t hash_t;
		typedef bit_vector bit_vector_t;
		typedef int_vector<16> index_kmer_t;
		
		// Costruttore
		SimpleBF(const size_t size, const string& id_gene): sx(nullptr), dx(nullptr), _size(size), _bf(size,0), _id(id_gene) {}
		
		SimpleBF(const SimpleBF& x): sx(x.sx), dx(x.dx), _size(x._size), _bf(x._bf), _id(x._id) {}
		
		// Distruttore
		~SimpleBF()
		{
			delete sx;
			delete dx;
		}
		
		void add_at(const uint64_t p)
		{
			_bf[p] = 1;
		}

		// Function to add a k-mer to the BF
		void add_kmer(const kmer_t kmer) 
		{
			uint64_t hash = _get_hash(kmer);
			_bf[hash % _size] = 1;
		}
		
		bool add_to_kmer(const uint64_t &kmer, const int &input_idx) 
		{
			uint64_t hash = _get_hash(kmer);
			size_t bf_idx = hash % _size;
			if (_bf[bf_idx]) 
			{
				return true;
			}
			return false;
		}
		
		pair<index_kmer_t::const_iterator, index_kmer_t::const_iterator> get_index(const kmer_t &kmer) const {
			
			return make_pair(_index_kmer.end(), _index_kmer.end());
		}
		
		void setSxChild(SimpleBF* _sx)
		{
			sx = _sx;
		}
		
		void setDxChild(SimpleBF* _dx)
		{
			dx = _dx;
		}
		
		void setBF(SimpleBF* A, SimpleBF* B)
		{
			for(size_t i=0; i < _size; i++)
				_bf[i] = A->_bf[i] || B->_bf[i];
			
		}
		
		void printBF(ostream& out)
		{
			out << _id << "\t";
			for(size_t i=0; i < _bf.size(); i++)
				out << (int)_bf[i];
			out << endl;
		}
		
	private:
		SimpleBF() = delete;
		
		SimpleBF* sx;
		SimpleBF* dx;
		size_t _size;
		bit_vector_t _bf;
		string _id;
		index_kmer_t _index_kmer;
};

#endif