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

#include <vector>

#include "kmer_utils.hpp"

class SSBT;

class SimpleBF {
  friend class SSBT;

public:
  SimpleBF(const size_t size, const int id_gene = -1)
    : sx(nullptr), dx(nullptr), parent(nullptr), _bf(size, false), _id(id_gene), support(false) {}

  SimpleBF(SimpleBF *_sx, SimpleBF *_dx)
    : sx(_sx), dx(_dx), parent(nullptr), _bf(max(sx->_bf.size(), dx->_bf.size()) * 2, false),
        _id(-1), support(false) {
    sx->support = _bf.size() >> 1 != sx->_bf.size();
    dx->support = _bf.size() >> 1 != dx->_bf.size();
    sx->set_parent(this);
    dx->set_parent(this);
  }

  ~SimpleBF() {
    delete sx;
    delete dx;
  }

  void add_at(const uint64_t p) { _bf[p % _bf.size()] = true; }

  SimpleBF* get_parent() const { return parent; }
  void set_parent(SimpleBF *p) { parent = p; }

  bool get_support() const { return support; }

private:
  typedef std::vector<bool> bit_vector_t;

  SimpleBF *sx;
  SimpleBF *dx;
  SimpleBF *parent;
  bit_vector_t _bf;
  int _id;
  bool support;
};

#endif
