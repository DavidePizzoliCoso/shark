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
		SimpleBF(const size_t size, const int id_gene, const int nHash): sx(nullptr), dx(nullptr), _size(size), _bf(size,0), _id(id_gene), _nHash(nHash) {}
		
		// Costruttore
		SimpleBF(const size_t size, const int nHash): sx(nullptr), dx(nullptr), _size(size), _bf(size,0), _nHash(nHash) {}
		
		// Costruttore di copia
		SimpleBF(const SimpleBF& x): sx(x.sx), dx(x.dx), _size(x._size), _bf(x._bf), _id(x._id), _nHash(x._nHash) {}
		
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
		void add_kmer(const kmer_t kmer, vector<size_t> &hash) 
		{
			hash = _get_hash(hash, kmer, _nHash, _size);
			for(const auto h : hash)
				_bf[h] = 1;
		}
		
		bool add_to_kmer(const uint64_t &kmer, const int &input_idx, vector<size_t> &hash) 
		{
			hash = _get_hash(hash, kmer, _nHash, _size);
			bool flag = true;
			
			for(const auto index : hash)
				flag = _bf[index] && flag;
			return flag;
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
		
		void changeBFSize(size_t size)
		{
			_size = size;
			_bf = bit_vector_t(size,0);
		}
		
		void set_id(int id)
		{
			_id = id;
		}
		
		void setBF(SimpleBF* A, SimpleBF* B)
		{
			for(size_t i=0; i < _size; i+=64)
				_bf.set_int(i, A->_bf.get_int(i,64) | B->_bf.get_int(i,64));
		}
		
		void printBF(ostream& out)
		{
			
			out << _id << " " << _size << " " << support << "\t";
			
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
		int _id;
		index_kmer_t _index_kmer;
		int _nHash;
		bool support;
};

#endif