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

#ifndef KMER_BUILDER_HPP
#define KMER_BUILDER_HPP

#include "kseq.h"
#include <zlib.h>
#include <string>
#include <vector>
#include <memory>
#include "simpleBF.h"
#include "kmer_utils.hpp"

using namespace std;

class KmerBuilder {

public:
  KmerBuilder(size_t _k, uint64_t _bf_size, int _nHash) : k(_k), bf_size(_bf_size), nHash(_nHash) {}

  vector<pair<string,vector<size_t>>>* operator()(vector<pair<string, string>> *texts) const {
    if(texts) {
	  vector<pair<string,vector<size_t>>>* ret = new vector<pair<string,vector<size_t>>>();
      vector<uint64_t> kmer_pos;
	  vector<size_t> hash;
      uint64_t kmer, rckmer, key;
	  
      for(const auto & p : *texts) {
		kmer_pos.clear();
		//cout<<">Number of elements: "<<p.second.size() - k + 1<<endl; // FASE 1
        if(p.second.size() >= k) {
          int _pos = 0;
          kmer = build_kmer(p.second, _pos, k);
          if(kmer == (uint64_t)-1) continue;
          rckmer = revcompl(kmer, k);
          key = min(kmer, rckmer);
		  
		  hash = _get_hash(hash, key, nHash, bf_size);
		  
		  kmer_pos.insert(kmer_pos.end(), hash.begin(), hash.end());
		  /*
		  cout<<key<<";"; // Inizio - FASE 1
		  for(const auto index : _get_hash(key, nHash))
			cout<<index % bf_size<<"|";
		  cout<<endl; // Fine - FASE 1
		  */
          for (int pos = _pos; pos < (int)p.second.size(); ++pos) {
            uint8_t new_char = to_int[p.second[pos]];
            if(new_char == 0) { // Found a char different from A, C, G, T
              ++pos; // we skip this character then we build a new kmer
              kmer = build_kmer(p.second, pos, k);
              if(kmer == (uint64_t)-1) break;
              rckmer = revcompl(kmer, k);
              --pos; // p must point to the ending position of the kmer, it will be incremented by the for
            } else {
              --new_char; // A is 1 but it should be 0
              kmer = lsappend(kmer, new_char, k);
              rckmer = rsprepend(rckmer, reverse_char(new_char), k);
            }
            key = min(kmer, rckmer);
			
			hash = _get_hash(hash, key, nHash, bf_size);
			  
			kmer_pos.insert(kmer_pos.end(), hash.begin(), hash.end());
			/*
			cout<<key<<";"; // Inizio - FASE 1
			for(const auto index : _get_hash(key, nHash))
				cout<<index % bf_size<<"|";
			cout<<endl; // Fine - FASE 1
			*/
          }
        }
		ret->push_back(make_pair(p.first,kmer_pos));
		//cout<<"="<<endl; // FASE 1
      }
      delete texts;
      return ret;
    }
    return NULL;
  }

private:
  size_t k;
  uint64_t bf_size;
  int nHash;
};

#endif
