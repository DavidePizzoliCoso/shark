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

#ifndef BF_FILLER_HPP
#define BF_FILLER_HPP

#include "bloomtree.hpp"
#include "simpleBF.hpp"
#include <cmath>
#include <deque>
#include <memory>
#include <string>
#include <vector>

using namespace std;

class BloomfilterFiller {
public:
  BloomfilterFiller(SSBT *_sbt, int _nHash, int *_cnt, vector<SimpleBF*> &_leaves)
      : sbt(_sbt), nHash(_nHash), counter(_cnt), leaves(_leaves) {}

  void operator()(vector<pair<string, vector<size_t>>> *genes) const {
	  
    SimpleBF *node;
    for (const auto &gene : *genes) {
      node = leaves[*counter];
	  
	  while(node)
	  {
		for (const auto position : gene.second)
          node->add_at(position);
		
		node = node->parent;
	  }
	  ++*counter;
    }
    delete genes;
  }

private:
  SSBT *sbt;
  const int nHash;
  int *counter;
  vector<SimpleBF*> &leaves;
};
#endif
