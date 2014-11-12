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

/*! \file FasterDegeneracyAlgorithm.cpp

    \brief This file contains the algorithm for listing all cliques
           according to the algorithm of Eppstein et al. (ISAAC 2010/SEA 2011).

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly

    See the algorithm's description in http://dx.doi.org/10.1007/978-3-642-20662-7_31
    and http://dx.doi.org/10.1007/978-3-642-17517-6_36

    This algorithm first orders the vertices in a degeneracy order (vertices
    are removed from the graph in order by degree in the remaining subgraph 
    and placed in this order in the ordering).

    We then recursively call a modified version of the algorithm of Tomita et
    al. (2006) http://dx.doi.org/10.1016/j.tcs.2006.06.015, for each vertex
    v in the ordering, where R = {v}, P = v's neighbors that are after v
    in the degneracy order, and X = v's neighbors that are before v 
    in the degeneracy order.

    This is a recursive backtracking algorithm that maintains three 
    sets of vertices, R, a partial clique, P, the common neighbors
    of vertices in R that are candidates to add to the partial clique,
    and X, the set of common neighbors of R that have been listed 
    in a maximal clique with R already.

    The algorithm recursively adds vertices to R from P, then 
    updates the sets P and X to be the new common neighbors of R
    and recurses. When P and X are empty, R is a maximal clique,
    and is reported.

    Updating the sets P and X is done by iterating over the later
    neighbors of vertices in P and X to see if they contain the vertex
    v added to R as a neighbor, and then testing if v's later neighbors
    are in P or X. Neighbors of v remain in their respective sets. For sparse
    graphs (graphs with low degeneracy) the number of neighbors of any vertex
    that come later in the degeneracy order is expected to be small. Then
    testing for neighbors in this way should be fast.

    Besides the ordering, this algorithm includes more efficient 
    pivot computation for sparse graphs.
    To find a "good" pivot, we must find the vertex in P union X that
    has the most neighbors in P. We accompish this step by maintaining
    a data structure that stores, for each vertex in P union X, its
    neighbors in P.

*/

// local includes
#include "DegeneracyVertexSets.h"
#include "DegeneracyTools.h"
#include "Tools.h"

// system includes
#include <limits.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <list>
#include <vector>
#include <set>
#include <utility>
#include <tuple>
#include <algorithm>

using namespace std;

DegeneracyVertexSets::DegeneracyVertexSets(vector<vector<int>> &adjacencyList)
: VertexSets("degeneracy")
, beginX(0)
, beginP(0)
, beginR(adjacencyList.size())
, orderingArray()
, neighborsInP()
, numNeighbors()
, m_AdjacencyList(adjacencyList)
, vertexSets()
, vertexLookup()
, m_lDelineators()
, m_iCurrentTopLevelIndex(0)
{
}

DegeneracyVertexSets::~DegeneracyVertexSets()
{
}

void DegeneracyVertexSets::Initialize()
{
    m_lDelineators.reserve(10);

    vertexSets  .resize(m_AdjacencyList.size(), 0);
    vertexLookup.resize(m_AdjacencyList.size(), 0);
    numNeighbors.resize(m_AdjacencyList.size(), 0);
    neighborsInP.resize(m_AdjacencyList.size());

    orderingArray = std::move(computeDegeneracyOrderArray(m_AdjacencyList, m_AdjacencyList.size()));

    for (int i = 0; i < neighborsInP.size(); ++i) {
        neighborsInP[i].resize(orderingArray[i].laterDegree);
    }

    // indices indicating where each set P, X, R starts in 
    // vertexSets; initially, P contains all the vertices
    for (size_t i = 0; i < m_AdjacencyList.size(); ++i) {
        vertexSets[i] = i;
        vertexLookup[i] = i;
    }
}

void DegeneracyVertexSets::PrintSummary(int const line) const
{
    cout << line << ": X[size=" << beginP-beginX << "], P[size=" << beginR-beginP << "]" << endl; 
}

bool DegeneracyVertexSets::GetNextTopLevelPartition()
{
    if (m_bDoneWithTopLevelPartitions) return false;

    //TODO/DS: Compute next top level partition
    int const orderNumber(m_iCurrentTopLevelIndex++);
    int const vertex(orderingArray[orderNumber].vertex);

    beginX = 0;
    beginP = 0;

    // fill in X with earlier neighbors
    for (int const neighbor : orderingArray[vertex].earlier) {
        if (vertexLookup[vertexSets[beginP]] == beginP)
            vertexLookup[vertexSets[beginP]] = -1; // make sure no one else claims this position.
        vertexSets[beginP] = neighbor;
        vertexLookup[neighbor] = (beginP)++;
    }

    beginR = beginP;

    // fill in P with later neighbors
    for (int const neighbor : orderingArray[vertex].later) {
        if (vertexLookup[vertexSets[beginR]] == beginR)
            vertexLookup[vertexSets[beginR]] = -1; // make sure no one else claims this position.
        vertexSets[beginR] = neighbor;
        vertexLookup[neighbor] = (beginR)++;
    }

    if (vertexLookup[vertexSets[beginR]] == beginR)
        vertexLookup[vertexSets[beginR]] = -1; // make sure no one else claims this position.

    vertexSets[beginR] = vertex;
    vertexLookup[vertex] = beginR;

    for (int const neighbor : orderingArray[orderNumber].earlier) {

        int const numNeighborsNeeded(min(beginR-beginP, orderingArray[neighbor].laterDegree));

        if (neighborsInP[neighbor].size() < numNeighborsNeeded) {
            neighborsInP[neighbor].resize(2*numNeighborsNeeded);
        }

        numNeighbors[neighbor] = 0;

        // fill in NeighborsInP
        for (int const laterNeighbor : orderingArray[neighbor].later) {
            int laterNeighborLocation = vertexLookup[laterNeighbor];
            if (laterNeighborLocation >= beginP && laterNeighborLocation < beginR) {
                neighborsInP[neighbor][numNeighbors[neighbor]] = laterNeighbor;
                numNeighbors[neighbor]++;
            }
        }
    }

    // reset numNeighbors and neighborsInP for this vertex
    int j = beginP;
    while (j<beginR) {
        int vertexInP = vertexSets[j];

        numNeighbors[vertexInP] = 0;

        int const numNeighborsNeeded(min( beginR-beginP, 
                    orderingArray[vertexInP].laterDegree 
                    + orderingArray[vertexInP].earlierDegree));

        if (neighborsInP[vertexInP].size() < numNeighborsNeeded) {
            neighborsInP[vertexInP].resize(2*numNeighborsNeeded);
        }

        j++;
    }

    // count neighbors in P, and fill in array of neighbors in P
    j = beginP;
    while (j<beginR) {
        int vertexInP = vertexSets[j];

        int k = 0;
        while (k<orderingArray[vertexInP].laterDegree) {
            int laterNeighbor = orderingArray[vertexInP].later[k];
            int laterNeighborLocation = vertexLookup[laterNeighbor];

            if (laterNeighborLocation >= beginP && laterNeighborLocation < beginR) {
                neighborsInP[vertexInP][numNeighbors[vertexInP]] = laterNeighbor;
                numNeighbors[vertexInP]++;
                neighborsInP[laterNeighbor][numNeighbors[laterNeighbor]] = vertexInP;
                numNeighbors[laterNeighbor]++;
            }

            k++;
        }

        j++;
    }

    m_bDoneWithTopLevelPartitions = (m_iCurrentTopLevelIndex == m_AdjacencyList.size());

    return true;
}
