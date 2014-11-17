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
#include "DegeneracyIndependentSets.h"
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

DegeneracyIndependentSets::DegeneracyIndependentSets(vector<vector<int>> &adjacencyList)
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

DegeneracyIndependentSets::~DegeneracyIndependentSets()
{
}

void DegeneracyIndependentSets::Initialize()
{
    cout << "There are " << m_AdjacencyList.size() << " vertices" << endl;
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

void DegeneracyIndependentSets::PrintSummary(int const line) const
{
    cout << line << ": X[size=" << beginP-beginX << "], P[size=" << beginR-beginP << "]" << endl; 
}

// DONE, need to verify
bool DegeneracyIndependentSets::GetNextTopLevelPartition()
{
    if (m_bDoneWithTopLevelPartitions) return false;

    //TODO/DS: Compute next top level partition
    int const orderNumber(m_iCurrentTopLevelIndex++);
    int vertex(0);
    for (int i = 0; i < orderingArray.size(); ++i) {
        if (orderingArray[i].orderNumber == orderNumber) {
            vertex = orderingArray[i].vertex;
            break;
        }
    }

    cout << "Evaluating another toplevel vertex...only " << orderingArray.size() - orderNumber << " more to go." << endl;

    beginX = 0;
    beginP = 0;

    vector<bool> vMarkedVertices(m_AdjacencyList.size(), false);

    for (int const earlierNeighbor : orderingArray[vertex].earlier) {
        vMarkedVertices[earlierNeighbor] = true;
////        vertexLookup[earlierNeighbor] = -1;
    }

    for (int const laterNeighbor : orderingArray[vertex].later) {
        vMarkedVertices[laterNeighbor] = true;
////        vertexLookup[laterNeighbor] = -1;
    }

    int const vertexOrderNumber(orderingArray[vertex].orderNumber);

    // fill in earlier non-neighbors
    for (size_t i = 0; i < orderingArray.size(); ++i) {
        int const currentOrderNumber(orderingArray[i].orderNumber);
        if (currentOrderNumber < vertexOrderNumber && !vMarkedVertices[i]) {
            if (vertexLookup[vertexSets[beginP]] == beginP) {
                vertexLookup[vertexSets[beginP]] = -1; // make sure no one else claims this position.
////                cout << __LINE__ << " : invalidating position of vertex " << vertexSets[beginR] << endl;
            }
            vertexSets[beginP] = i;
            vertexLookup[i] = beginP++;
        }

////        cout << __LINE__ << " : vertex " << 75 << " in P?=" << InP(75) << endl;
    }

    beginR = beginP;

    // fill in later non-neighbors
    for (size_t i = 0; i < orderingArray.size(); ++i) {
        int const currentOrderNumber(orderingArray[i].orderNumber);
        if (currentOrderNumber > vertexOrderNumber && !vMarkedVertices[i]) {
            if (vertexLookup[vertexSets[beginR]] == beginR) {
                vertexLookup[vertexSets[beginR]] = -1; // make sure no one else claims this position.
////                cout << __LINE__ << " : invalidating position of vertex " << vertexSets[beginR] << endl;
            }
            vertexSets[beginR] = i;
            vertexLookup[i] = beginR++;
        }

////        cout << __LINE__ << " : vertex " << 75 << " in P?=" << InP(75) << endl;
    }

    if (vertexLookup[vertexSets[beginR]] == beginR)
        vertexLookup[vertexSets[beginR]] = -1; // make sure no one else claims this position.

    vertexSets[beginR] = vertex;
    vertexLookup[vertex] = beginR;

////    cout << __LINE__ << " : vertex " << 75 << " in P?=" << InP(75) << endl;

////    cout << __LINE__ << ": numNeighbors[11]=" << numNeighbors[11] << endl;

    int j = beginX;
    while (j < beginP) { // int const neighbor : orderingArray[orderNumber].earlier) {
        int const neighbor(vertexSets[j]);

////        cout << __LINE__ << " : evaluating earlier nonneighbor " << neighbor << endl;
        int const numNeighborsNeeded(min(beginR-beginP, orderingArray[neighbor].laterDegree));

        if (neighborsInP[neighbor].size() < numNeighborsNeeded) {
            neighborsInP[neighbor].resize(2*numNeighborsNeeded);
        }

        numNeighbors[neighbor] = 0;
////        cout << __LINE__ << ": numNeighbors[11]=" << numNeighbors[11] << endl;

        // fill in NeighborsInP
        for (int const laterNeighbor : orderingArray[neighbor].later) {
////            cout << __LINE__ << " : vertex " << 75 << " in neighborsOfP[11]?=" << (find(neighborsInP[11].begin(), neighborsInP[11].begin() + numNeighbors[11] + 1, 75) != neighborsInP[11].begin() + numNeighbors[11] + 1) << "InP?=" << InP(75) << endl;
////        if ((neighbor == 11 || laterNeighbor == 11) /*&& ( neighbor == 75 || laterNeighbor == 75)*/)
////            cout << __LINE__ << ": evaluating vertex 11" << endl;

            int const laterNeighborLocation = vertexLookup[laterNeighbor];
            if (laterNeighborLocation >= beginP && laterNeighborLocation < beginR) {
////        if ((neighbor == 11 || laterNeighbor == 11) /*&& ( neighbor == 75 || laterNeighbor == 75)*/)
////                    cout << __LINE__ << ": vertex " << neighbor << " has later neighbor " << laterNeighbor << " in P : InP?=" << InP(laterNeighbor) << endl;
                neighborsInP[neighbor][numNeighbors[neighbor]] = laterNeighbor;
                numNeighbors[neighbor]++;
            } else {
////        if ((neighbor == 11 || laterNeighbor == 11) /*&& ( neighbor == 75 || laterNeighbor == 75)*/)
////                    cout << __LINE__ << ": vertex " << neighbor << " doesn't have later neighbor " << laterNeighbor << " in P" << endl;
            }
        }

        j++;

////        cout << __LINE__ << " : vertex " << 75 << " in P?=" << InP(75) << endl;
    }

    // reset numNeighbors and neighborsInP for this vertex
    j = beginP;
    while (j<beginR) {
        int const vertexInP = vertexSets[j];
        numNeighbors[vertexInP] = 0;

        int const numNeighborsNeeded(min( beginR-beginP, 
                    orderingArray[vertexInP].laterDegree 
                    + orderingArray[vertexInP].earlierDegree));

////        cout << "vertex in P " << vertexInP << endl << flush;

        if (neighborsInP[vertexInP].size() < numNeighborsNeeded) {
            neighborsInP[vertexInP].resize(2*numNeighborsNeeded);
        }

        j++;

////        cout << __LINE__ << " : vertex " << 75 << " in P?=" << InP(75) << endl;
    }

////    cout << __LINE__ << ": numNeighbors[11]=" << numNeighbors[11] << endl;

    // count neighbors in P, and fill in array of neighbors in P
    j = beginP;
    while (j<beginR) {
        int const vertexInP = vertexSets[j];

        int k = 0;
        while (k<orderingArray[vertexInP].laterDegree) {
////            cout << __LINE__ << " : vertex " << 75 << " in neighborsOfP[11]?=" << (find(neighborsInP[11].begin(), neighborsInP[11].begin() + numNeighbors[11] + 1, 75) != neighborsInP[11].begin() + numNeighbors[11] + 1) << "InP?=" << InP(75) << endl;
            int const laterNeighbor = orderingArray[vertexInP].later[k];
            int const laterNeighborLocation = vertexLookup[laterNeighbor];
////        if (vertexInP == 11 || laterNeighbor == 11)
////            cout << __LINE__ << ": evaluating vertex 11" << endl;

            if (laterNeighborLocation >= beginP && laterNeighborLocation < beginR) {
////        if (vertexInP == 11 || laterNeighbor == 11)
////                    cout << __LINE__ << ": vertex " << vertexInP << " has later neighbor " << laterNeighbor << " in P" << endl;
                neighborsInP[vertexInP][numNeighbors[vertexInP]] = laterNeighbor;
                numNeighbors[vertexInP]++;
                neighborsInP[laterNeighbor][numNeighbors[laterNeighbor]] = vertexInP;
                numNeighbors[laterNeighbor]++;
            } else {
////        if (vertexInP == 11 || laterNeighbor == 11)
////                    cout << __LINE__ << ": vertex " << vertexInP << " doesn't have later neighbor " << laterNeighbor << " in P" << endl;
            }

            k++;
        }

        j++;

////        cout << __LINE__ << " : vertex " << 75 << " in P?=" << InP(75) << endl;
////        cout << __LINE__ << " : vertex " << 75 << " in neighborsOfP[11]?=" << (find(neighborsInP[11].begin(), neighborsInP[11].begin() + numNeighbors[11] + 1, 75) != neighborsInP[11].begin() + numNeighbors[11] + 1) << "InP?=" << InP(75) << endl;
    }

////    cout << __LINE__ << " : vertex " << 75 << " in P?=" << InP(75) << endl;

////    std::cout << "X and P: ";
////    std::cout << "X=";
////    for (int i = beginX; i < beginP; i++) {
////        std::cout << vertexSets[i];
////        if (orderingArray[vertexSets[i]].orderNumber >= vertexOrderNumber) cout << "*";
////        std::cout << " ";
////    }
////    cout << "P=";
////    for (int i = beginP; i < beginR; i++) {
////        std::cout << vertexSets[i];
////        if (orderingArray[vertexSets[i]].orderNumber <= vertexOrderNumber) cout << "*";
////        std::cout << " ";
////    }
////    std::cout << std::endl;
////
////    std::cout << "InX?";
////    for (int i = beginX; i < beginP; i++) {
////        int const location(vertexLookup[vertexSets[i]]);
////        std::cout << (location >= beginX && location < beginP && location == i) << " ";
////    }
////    std::cout << std::endl;
////
////    std::cout << "InP?";
////    for (int i = beginP; i < beginR; i++) {
////        int const location(vertexLookup[vertexSets[i]]);
////        std::cout << (location >= beginP && location < beginR && location == i) << " ";
////    }
////
////    std::cout << std::endl;

    m_bDoneWithTopLevelPartitions = (m_iCurrentTopLevelIndex == m_AdjacencyList.size());

////    return VerifyStartConfiguration();
    return true;
}

bool DegeneracyIndependentSets::VerifyStartConfiguration() const
{
    bool foundError(false);
    set<int> setVerticesInP;
    for (int u = beginP; u < beginR; ++u) {
        setVerticesInP.insert(vertexSets[u]);
    }

    int const vertexInR(vertexSets[beginR]);

    for (size_t u = beginX; u < beginR; ++u) {
        int const vertex = vertexSets[u];
        set<int> neighbors;
        neighbors.insert(m_AdjacencyList[vertex].begin(), m_AdjacencyList[vertex].end());

        // make sure non-neighbor of partial clique
        if (neighbors.find(vertex) != neighbors.end()) {
            cout << __LINE__ << " : ERROR: vertex " << vertexInR << " is a neighbor of : " << vertex << endl;
            foundError = true;
        }

        // check that neighbors in P only contains neighbors in P
        for (size_t v = 0; v < numNeighbors[vertex]; ++v) {
            int const neighbor(neighborsInP[vertex][v]);
            if (neighbors.find(neighbor) == neighbors.end()) {
                cout << __LINE__ << " : ERROR: nonneighbor of " << vertex << " is included in neighborsInP: " << neighbor << endl;
                foundError = true;
            }
            if (!InP(neighbor)) {
                cout << __LINE__ << " : ERROR: neighbor of " << vertex << " is not in P, but is included in neighborsInP: " << neighbor << endl;
                foundError = true;
            }
        }

        // check that all neighbors in P are in neighborsOfP
        for (int const neighbor : m_AdjacencyList[vertex]) {
            if (InP(neighbor) && setVerticesInP.find(neighbor) == setVerticesInP.end()) {
                cout << __LINE__ << " : ERROR: neighbor of " << vertex << " is missing from neighborsInP: " << neighbor << endl;
                foundError = true;
            }
            if (!InP(neighbor) && setVerticesInP.find(neighbor) != setVerticesInP.end()) {
                cout << __LINE__ << " : ERROR: neighbor of " << vertex << " is not in P, but is included in neighborsInP: " << neighbor << endl;
                foundError = true;
            }
        }
    }

    int const vertexInROrderNumber(orderingArray[vertexInR].orderNumber);

    for (size_t u = beginX; u < beginP; ++u) {
        int const currentVertex(vertexSets[u]);
        int const currentOrderNumber(orderingArray[currentVertex].orderNumber);
        if (currentOrderNumber >= vertexInROrderNumber) {
            cout << __LINE__ << " : ERROR: vertex " << currentVertex << " is in X, but its order number is not less than vertex " << vertexInR << endl;
            foundError = true;
        }
    }

    for (size_t u = beginP; u < beginR; ++u) {
        int const currentVertex(vertexSets[u]);
        int const currentOrderNumber(orderingArray[currentVertex].orderNumber);
        if (currentOrderNumber <= vertexInROrderNumber) {
            cout << __LINE__ << " : ERROR: vertex " << currentVertex << " is in P, but its order number is not greater than vertex " << vertexInR << endl;
            foundError = true;
        }
    }

    return !foundError;
}
