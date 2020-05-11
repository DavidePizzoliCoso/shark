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
#include <string>
#include <vector>

#include "kmer_utils.hpp"

using namespace std;

class KmerBuilder;
class BloomfilterFiller;
class SSBT;


class SimpleBF {
  friend class KmerBuilder;
  friend class BloomfilterFiller;
  friend class SSBT;

public:
  typedef uint64_t kmer_t;
  typedef uint64_t hash_t;
  typedef vector<bool> bit_vector_t;

  // Costruttore
  SimpleBF(const size_t size, const int id_gene, const int nHash)
      : sx(nullptr), dx(nullptr), _size(size), _bf(size, 0), _id(id_gene),
        _nHash(nHash) {}

  // Costruttore
  SimpleBF(const size_t size, const int nHash)
      : sx(nullptr), dx(nullptr), _size(size), _bf(size, 0), _nHash(nHash) {}

  // Costruttore di copia
  SimpleBF(const SimpleBF &x)
      : sx(x.sx), dx(x.dx), _size(x._size), _bf(x._bf), _id(x._id),
        _nHash(x._nHash) {}

  // Distruttore
  ~SimpleBF() {
    delete sx;
    delete dx;
  }

  void add_at(const uint64_t p) { _bf[p] = 1; }

  void setSxChild(SimpleBF *_sx) { sx = _sx; }

  void setDxChild(SimpleBF *_dx) { dx = _dx; }

  void changeBFSize(size_t size) {
    _size = size;
    _bf = bit_vector_t(size, 0);
  }

  void set_id(int id) { _id = id; }

  SimpleBF() = delete;

private:
  SimpleBF *sx;
  SimpleBF *dx;
  size_t _size;
  bit_vector_t _bf;
  int _id;
  int _nHash;
  bool support;
};

#endif
