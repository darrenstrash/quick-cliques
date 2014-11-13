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
#include "CacheEfficientDegeneracyVertexSets.h"
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

////#define OLD_CREATE

using namespace std;

CacheEfficientDegeneracyVertexSets::CacheEfficientDegeneracyVertexSets(vector<vector<int>> &adjacencyList)
: VertexSets("degeneracy")
, beginX(0)
, beginP(0)
, beginR(adjacencyList.size())
, orderingArray()
, newNeighborsInP()
, neighborsIndex()
, m_AdjacencyList(adjacencyList)
, vertexSets()
, vertexLookup()
, m_lDelineators()
, m_iCurrentTopLevelIndex(0)
{
}

CacheEfficientDegeneracyVertexSets::~CacheEfficientDegeneracyVertexSets()
{
}

void CacheEfficientDegeneracyVertexSets::Initialize()
{
    m_lDelineators.reserve(10);

    vertexSets  .resize(m_AdjacencyList.size(), 0);
    vertexLookup.resize(m_AdjacencyList.size(), 0);

    orderingArray = std::move(computeDegeneracyOrderArray(m_AdjacencyList, m_AdjacencyList.size()));

    // indices indicating where each set P, X, R starts in 
    // vertexSets; initially, P contains all the vertices
    for (size_t i = 0; i < m_AdjacencyList.size(); ++i) {
        vertexSets[i] = i;
        vertexLookup[i] = i;
    }
}

void CacheEfficientDegeneracyVertexSets::PrintSummary(int const line) const
{
    cout << line << ": X[size=" << beginP-beginX << "], P[size=" << beginR-beginP << "]" << endl; 
}

int CacheEfficientDegeneracyVertexSets::GetIndexOfVertexInP(int const vertex) const
{
    for (int i = 0; i < neighborsIndex.size(); ++i) {
        if (std::get<0>(neighborsIndex[i]) == vertex)
            return i;
    }
    return -1;
};

bool CacheEfficientDegeneracyVertexSets::GetNextTopLevelPartition()
{
    if (m_bDoneWithTopLevelPartitions) return false;

    //TODO/DS: Compute next top level partition
    int const orderNumber(m_iCurrentTopLevelIndex++);
    int const vertex(orderingArray[orderNumber].vertex);

    beginX = 0;
    beginP = 0;

    int sizeOfNeighborsInP(0);

    // fill in X with earlier neighbors
    for (int const neighbor : orderingArray[vertex].earlier) {
        if (vertexLookup[vertexSets[beginP]] == beginP)
            vertexLookup[vertexSets[beginP]] = -1; // make sure no one else claims this position.
        vertexSets[beginP] = neighbor;
        vertexLookup[neighbor] = (beginP)++;
        sizeOfNeighborsInP += min(orderingArray[neighbor].laterDegree,orderingArray[vertex].laterDegree);
    }

////    cout << "Reserved    for X later neighbors: " << sizeOfNeighborsInP << endl << flush;
////
////    int const savedXLaterNeighbors(sizeOfNeighborsInP);

    beginR = beginP;

    int index(0);

    // fill in P with later neighbors
    for (int const neighbor : orderingArray[vertex].later) {
        if (vertexLookup[vertexSets[beginR]] == beginR)
            vertexLookup[vertexSets[beginR]] = -1; // make sure no one else claims this position.
        vertexSets[beginR] = neighbor;
        vertexLookup[neighbor] = (beginR)++;
////        sizeOfNeighborsInP += min(orderingArray[neighbor].laterDegree+index,orderingArray[vertex].laterDegree);
        sizeOfNeighborsInP += orderingArray[vertex].laterDegree;
        index++;
    }

////    cout << "Reserved    for P later neighbors: " << sizeOfNeighborsInP - savedXLaterNeighbors << endl << flush;

    newNeighborsInP.resize(sizeOfNeighborsInP, -1);

    int lastNeighborOfPIndex(0);
    int currentNeighborOfPIndex(0);

    if (vertexLookup[vertexSets[beginR]] == beginR)
        vertexLookup[vertexSets[beginR]] = -1; // make sure no one else claims this position.

    vertexSets[beginR] = vertex;
    vertexLookup[vertex] = beginR;

    neighborsIndex.clear();
    neighborsIndex.reserve(orderingArray[vertex].laterDegree + orderingArray[vertex].earlierDegree);

    for (int const neighbor : orderingArray[orderNumber].earlier) {
        lastNeighborOfPIndex = currentNeighborOfPIndex;

        // fill in newNeighborsInP
        for (int const laterNeighbor : orderingArray[neighbor].later) {
            int laterNeighborLocation = vertexLookup[laterNeighbor];
            if (laterNeighborLocation >= beginP && laterNeighborLocation < beginR) {
                newNeighborsInP[currentNeighborOfPIndex++] = laterNeighbor;
            }
        }

        neighborsIndex.push_back(make_tuple(neighbor, lastNeighborOfPIndex, currentNeighborOfPIndex));
    }

////    cout << "Actual used for X later neighbors: " << currentNeighborOfPIndex << endl << flush;
////
////    int const savedActualXLaterNeighbors(currentNeighborOfPIndex);

    //TODO/DS: Sort to make finding elements faster
    vector<vector<int>> vvTemp(orderingArray[vertex].laterDegree);

    int j = beginP;
    int const startOfPInNeighborsIndex(neighborsIndex.size());
    while (j<beginR) {
////        cout << "Size of neighborsIndex: " << neighborsIndex.size() << endl;
        int const vertexInP = vertexSets[j];
        neighborsIndex.push_back(make_tuple(vertexInP, j - beginP, -1));
////        cout << "Inserting <" << vertexInP << "," << j - beginP << "," << -1 << ">" << endl; 
        vvTemp[j-beginP].reserve(beginR-beginP);
        j++;
    }


    // count neighbors in P, and fill in array of neighbors in P
    j = beginP;
    while (j<beginR) {
        int const vertexInP = vertexSets[j];
        int const vertexIndex(j-beginP);

        int k = 0;
        while (k<orderingArray[vertexInP].laterDegree) {
            int const laterNeighbor = orderingArray[vertexInP].later[k];
            int const laterNeighborLocation = vertexLookup[laterNeighbor];

            if (laterNeighborLocation >= beginP && laterNeighborLocation < beginR) {
////                cout << __LINE__ << " : Index of " << vertexInP << "=" << vertexIndex << " in vvTemp" << endl;
                vvTemp[vertexIndex].push_back(laterNeighbor);
                vvTemp[GetIndexOfVertexInP(laterNeighbor) - startOfPInNeighborsIndex].push_back(vertexInP);
////                cout << __LINE__ << " : Index of " << laterNeighbor << "=" << GetIndexOfVertexInP(laterNeighbor) << " in vvTemp"  << endl;
////                cout << __LINE__ << " : Inserting edges (" << std::get<0>(neighborsIndex[vertexIndex + startOfPInNeighborsIndex]) << "," << std::get<0>(neighborsIndex[GetIndexOfVertexInP(laterNeighbor)]) << ")" << endl;
            }

            k++;
        }

        j++;
    }

    for (int i = 0; i < vvTemp.size(); ++i) {
        lastNeighborOfPIndex = currentNeighborOfPIndex;
        for (int const neighbor : vvTemp[i]) {
////            if (sizeOfNeighborsInP <= currentNeighborOfPIndex)
////                cout << "Current cnt for P later neighbors: " << currentNeighborOfPIndex - savedActualXLaterNeighbors << endl << flush;

            newNeighborsInP[currentNeighborOfPIndex++] = neighbor;
////            cout << __LINE__ << " : Inserting edge " << std::get<0>(neighborsIndex[startOfPInNeighborsIndex + i]) << "->" << neighbor << "" << endl;
        }

        std::get<1>(neighborsIndex[startOfPInNeighborsIndex + i]) = lastNeighborOfPIndex;
        std::get<2>(neighborsIndex[startOfPInNeighborsIndex + i]) = currentNeighborOfPIndex;
    }


////    cout << "Actual used for P later neighbors: " << currentNeighborOfPIndex - savedActualXLaterNeighbors << endl << flush;

    m_bDoneWithTopLevelPartitions = (m_iCurrentTopLevelIndex == m_AdjacencyList.size());

    return true;
}
