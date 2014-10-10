#ifndef _DJS_HYBRID_ALGORITHM_H_
#define _DJS_HYBRID_ALGORITHM_H_

/* 
    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by 
    the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version. 
 
    This program is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
    GNU General Public License for more details. 
 
    You should have received a copy of the GNU General Public License 
    along with this program.  If not, see <http://www.gnu.org/licenses/> 
*/

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "Tools.h"
#include <list>
#include "MemoryManager.h"
#include "DegeneracyTools.h"

/*! \file HybridAlgorithm.h

    \brief see HybridAlgorithm.cpp

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

void listAllMaximalCliquesHybridRecursive( long* cliqueCount,
                                           #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                           std::list<std::list<int>> &cliques,
                                           #endif
                                           std::list<int> &partialClique, 
                                           NeighborListArray** orderingArray,
                                           int* vertexSets, int* vertexLookup,
                                           int beginX, int beginP, int beginR,
                                           int* scratch);

long listAllMaximalCliquesHybrid( std::vector<std::list<int>> const &adjList, 
                                  int** adjacencyList, 
                                  #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                  std::list<std::list<int>> &cliques,
                                  #endif
                                  int* degree, 
                                  int size);

#endif
