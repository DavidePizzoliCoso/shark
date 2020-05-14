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
  SimpleBF(const int id_gene = -1)
      : sx(nullptr), dx(nullptr), parent(nullptr), _bf(0, false), _id(id_gene) {
  }

  SimpleBF(SimpleBF *_sx, SimpleBF *_dx)
      : sx(_sx), dx(_dx), parent(nullptr), _bf(sx->size() * 2, false), _id(-1) {
    _sx->parent = this;
    _dx->parent = this;
  }

  ~SimpleBF() {
    delete sx;
    delete dx;
  }

  void add_at(const uint64_t p) { _bf[p & (_bf.size() - 1)] = true; }

  SimpleBF *get_parent() const { return parent; }

  size_t size() const { return _bf.size(); };
  void resize(const size_t size) {
    if (size == _bf.size()) return;
    _bf.clear();
    _bf.resize(size, false);
    if (sx != nullptr) sx->resize(size >> 1);
    if (dx != nullptr) dx->resize(size >> 1);
  }

private:
  typedef std::vector<bool> bit_vector_t;

  SimpleBF *const sx;
  SimpleBF *const dx;
  SimpleBF *parent;
  bit_vector_t _bf;
  const int _id;
};

#endif
