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

#include "kseq.h"
#include <zlib.h>
#include <string>
#include <vector>
#include <list>
#include <memory>
#include "bloomtree.h"
#include "simpleBF.h"

using namespace std;

class BloomfilterFiller {
public:
	BloomfilterFiller(SSBT *_sbt, int _nHash) : sbt(_sbt), nHash(_nHash) {}

	void operator()(vector<pair<string,vector<size_t>>> *genes) const 
	{
		list<SimpleBF*> coda;
		SimpleBF* bloom;
		int i = 0;
		
		for(const auto & gene : *genes) 
		{
			bloom = new SimpleBF(sbt->_size, i, nHash);
			
			for(const auto position : gene.second)
			{
				bloom->add_at(position);
			}
			coda.push_back(bloom);
			++i;
		}
		delete genes;
		
		while (coda.size() > 1)
		{
			SimpleBF* sx = coda.front();
			coda.pop_front();
			SimpleBF* dx = coda.front();
			coda.pop_front();
			
			SimpleBF* node = new SimpleBF(sbt->_size, nHash);
			node->setSxChild(sx);
			node->setDxChild(dx);
			node->setBF(sx, dx);
			coda.push_back(node);
		}
		sbt->setRoot(coda.front());
		//sbt->printTree();
		
	}

private:
	SSBT* sbt;
	int nHash;
};
#endif
