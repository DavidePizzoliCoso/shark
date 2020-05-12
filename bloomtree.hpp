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

#include "simpleBF.hpp"
#include <algorithm>
#include <array>
#include <deque>
#include <string>

#include "kmer_utils.hpp"

using namespace std;

class KmerBuilder;
class BloomfilterFiller;

class SSBT {
  friend class KmerBuilder;
  friend class BloomfilterFiller;

public:
  typedef uint64_t kmer_t;

  explicit
  SSBT(const size_t size) : _size(size) {}

  ~SSBT() { delete root; }

  void setRoot(SimpleBF *node) { root = node; }

  void get_genes(const kmer_t &kmer, vector<int> &genes, vector<size_t> &hash) const {
    genes.clear();
    _get_hash(hash, kmer, _size);
    inner_get_genes(root, _size - 1, hash, genes);
  }

  SSBT() = delete;
  const SSBT &operator=(const SSBT &) = delete;
  const SSBT &operator=(const SSBT &&) = delete;

private:
  void inner_get_genes(const SimpleBF *node, const size_t dynamic_mask,
                       const vector<size_t> &hash, vector<int> &genes) const {
    for (const auto index : hash) {
      if (!node->_bf[index & dynamic_mask])
        return;
    }

    if (node->sx == nullptr) {
      // This is a leaf
      genes.push_back(node->_id);
    } else {
      inner_get_genes(node->sx, dynamic_mask >> (1 + node->sx->support), hash,
                      genes);
      inner_get_genes(node->dx, dynamic_mask >> (1 + node->dx->support), hash,
                      genes);
    }
  }

  SimpleBF *root;
  size_t _size;
};

#endif
