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
#include "CliqueTools.h"
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

    {
        int dumbGuess(0);
        cout << "Dumb guess: <Maximum independent set>" << flush;

        vector<bool> vMarkedVertices(m_AdjacencyList.size(), false);

        for (int orderNumber = 0; orderNumber < orderingArray.size(); ++orderNumber) {
            int vertex(-1);
            for (int i = 0; i < orderingArray.size(); ++i) {
                if (orderingArray[i].orderNumber == orderNumber) {
                    vertex = orderingArray[i].vertex;
                    break;
                }
            }

            if (!vMarkedVertices[vertex]) {
                vMarkedVertices[vertex] = true;
                dumbGuess++;
                for (int const neighbor : orderingArray[vertex].later) {
                    vMarkedVertices[neighbor] = true;
                }
            }
        }
        cout << dumbGuess << endl << flush;
    }

#if 0
    {
        vector<bool> vMarkedVertices(m_AdjacencyList.size(), false);

        int dumbGuess(0);
        cout << "Dumb guess: (upper bound) <Maximum independent set>" << flush;
        vector<vector<int>> cliques;
        CliqueTools::DecomposeIntoDisjointCliques(m_AdjacencyList, cliques);

        cout << cliques.size() << endl << flush;
    }
#endif


#if 0
    {
        vector<NeighborListArray> maximumOrderingArray = std::move(computeMaximumLaterOrderArray(m_AdjacencyList, m_AdjacencyList.size()));
        int dumbGuess(0);
        cout << "Dumb guess: <Minimum independent set>" << flush;

        vector<bool> vMarkedVertices(m_AdjacencyList.size(), false);

        for (int orderNumber = 0; orderNumber < maximumOrderingArray.size(); ++orderNumber) {
            int vertex(-1);
            for (int i = 0; i < orderingArray.size(); ++i) {
                if (maximumOrderingArray[i].orderNumber == orderNumber) {
                    vertex = maximumOrderingArray[i].vertex;
                    break;
                }
            }

            if (!vMarkedVertices[vertex]) {
                vMarkedVertices[vertex] = true;
                dumbGuess++;
                for (int const neighbor : maximumOrderingArray[vertex].later) {
                    vMarkedVertices[neighbor] = true;
                }
            }
        }
        cout << dumbGuess << endl << flush;
    }
#endif

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
    int vertex(-1);
    for (int i = 0; i < orderingArray.size(); ++i) {
        if (orderingArray[i].orderNumber == orderNumber) {
            vertex = orderingArray[i].vertex;
            break;
        }
    }

    if (vertex == -1) {
        cout << __LINE__ << ": ERROR: a vertex could not be selected with order number " << orderNumber << endl;
        m_bDoneWithTopLevelPartitions = true;
        return false;
    }

////    cout << "Evaluating another top level vertex: " << vertex << endl;

////    cout << "Evaluating another toplevel vertex...only " << orderingArray.size() - orderNumber << " more to go." << endl;

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

////    cout << "Filling in later neighbors of X" << endl;

    int j = beginX;
    while (j < beginP) { // int const neighbor : orderingArray[orderNumber].earlier) {
        int const neighbor(vertexSets[j]);
////        cout << "Evaluating vertex: " << neighbor << endl;

////        cout << __LINE__ << " : evaluating earlier nonneighbor " << neighbor << endl;
        int const numNeighborsNeeded(min(beginR-beginP, orderingArray[neighbor].laterDegree));

        if (neighborsInP[neighbor].size() < numNeighborsNeeded) {
            neighborsInP[neighbor].resize(2*numNeighborsNeeded);
        }

        numNeighbors[neighbor] = 0;
////        cout << __LINE__ << ": numNeighbors[11]=" << numNeighbors[11] << endl;

        // fill in NeighborsInP
        for (int const laterNeighbor : orderingArray[neighbor].later) {
////            cout << "     evaluating later neighbor: " << laterNeighbor << endl;
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

////        if (numNeighbors[neighbor] == 0) {
////            cout << "Could have pruined this recursion." << endl;
////        }

////        cout << __LINE__ << " : vertex " << 75 << " in P?=" << InP(75) << endl;
    }

////    cout << "Filling in later neighbors of P" << endl;
    // reset numNeighbors and neighborsInP for this vertex
    j = beginP;
    while (j<beginR) {
        int const vertexInP = vertexSets[j];
        numNeighbors[vertexInP] = 0;

        int const numNeighborsNeeded(min( beginR-beginP, 
                    orderingArray[vertexInP].laterDegree 
                    + orderingArray[vertexInP].earlierDegree));

////        if (vertexInP==11) {
////            cout << "Computing size for 11: " << numNeighborsNeeded << endl;
////            cout << "Size of P: " << beginR-beginP << endl;
////            cout << "Size of Neighbors of 11: " << (orderingArray[vertexInP].laterDegree + orderingArray[vertexInP].earlierDegree) << endl;
////        }


////        cout << "vertex in P " << vertexInP << endl << flush;

        if (neighborsInP[vertexInP].size() < numNeighborsNeeded) {
            neighborsInP[vertexInP].resize(2*numNeighborsNeeded);
        }

        j++;

////        cout << __LINE__ << " : vertex " << 75 << " in P?=" << InP(75) << endl;
    }

////    cout << "#slots reserved for 11: " << neighborsInP[11].size() << endl;

////    cout << __LINE__ << ": numNeighbors[11]=" << numNeighbors[11] << endl;

////    cout << __LINE__ << " : P:";
////    for (int i = beginP; i < beginR; i++) {
////        cout << vertexSets[i] << " ";
////    }
////    cout << endl;
////
////    cout << "Neighbors of 11 in P : " << endl;

    // count neighbors in P, and fill in array of neighbors in P
    j = beginP;
    while (j<beginR) {
        int const vertexInP = vertexSets[j];
////        cout << "Evaluating vertex: " << vertexInP << endl;

        int k = 0;
        while (k<orderingArray[vertexInP].laterDegree) {
////            cout << __LINE__ << " : vertex " << 75 << " in neighborsOfP[11]?=" << (find(neighborsInP[11].begin(), neighborsInP[11].begin() + numNeighbors[11] + 1, 75) != neighborsInP[11].begin() + numNeighbors[11] + 1) << "InP?=" << InP(75) << endl;
            int const laterNeighbor = orderingArray[vertexInP].later[k];
            int const laterNeighborLocation = vertexLookup[laterNeighbor];
////            cout << "     evaluating neighbor: " << laterNeighbor << endl;
////        if (vertexInP == 11 || laterNeighbor == 11)
////            cout << __LINE__ << ": evaluating vertex 11" << endl;

////            cout << "Evaluating later neighbor " << laterNeighbor << endl;

            if (laterNeighborLocation >= beginP && laterNeighborLocation < beginR) {
////        if (vertexInP == 11 || laterNeighbor == 11)
////                    cout << __LINE__ << ": vertex " << vertexInP << " has later neighbor " << laterNeighbor << " in P" << endl;
                neighborsInP[vertexInP][numNeighbors[vertexInP]] = laterNeighbor;
                numNeighbors[vertexInP]++;
                neighborsInP[laterNeighbor][numNeighbors[laterNeighbor]] = vertexInP;
                numNeighbors[laterNeighbor]++;
////                if (laterNeighbor == 11 || vertexInP == 11) {
////                    cout << "(" << vertexInP << "," << laterNeighbor << ") ";
////                }
            } else {
////        if (vertexInP == 11 || laterNeighbor == 11)
////                    cout << __LINE__ << ": vertex " << vertexInP << " doesn't have later neighbor " << laterNeighbor << " in P" << endl;
            }

            k++;
        }

        j++;
////        if (numNeighbors[vertexInP] == 0) {
////            cout << "Could have removed vertex from consideration" << endl;
////        }

////        cout << __LINE__ << " : vertex " << 75 << " in P?=" << InP(75) << endl;
////        cout << __LINE__ << " : vertex " << 75 << " in neighborsOfP[11]?=" << (find(neighborsInP[11].begin(), neighborsInP[11].begin() + numNeighbors[11] + 1, 75) != neighborsInP[11].begin() + numNeighbors[11] + 1) << "InP?=" << InP(75) << endl;
    }

////    cout << endl;

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

    vector<int> vDominatedVertices;
////    RemoveDominatedVertices(vDominatedVertices);

#if 0
    //sort (vDominatedVertices.begin(),  vDominatedVertices.end());
    cout << "D : ";
    for (int const dominatedVertex : vDominatedVertices) {
        cout << dominatedVertex << " ";
    }
    cout << endl;
#endif

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

void DegeneracyIndependentSets::RemoveDominatedVertices(vector<int> &removedVertices)
{
////        clock_t clockStart = clock();
        
////    cout << "X=" << beginD-beginX << ", P=" << beginR-beginP << endl;
    // for each vertex x in X

////    bool lookForDominatedNodes(true);
////    while (lookForDominatedNodes) {

    vector<bool> vMarkedNeighbors(m_AdjacencyList.size(), false);

    int const savedBeginR = beginR;
    for (int i = beginP; i < beginR; i++) {
////        lookForDominatedNodes = false;

        // mark neighbors in an array
        int const p(vertexSets[i]);
        int const potentialNeighborsOfP = numNeighbors[p]; // TODO/DS: optimize. Can we stop at the first one that isn't in P?
////        bool const debug(p == 68);
////        vMarkedNeighbors[p] = true;
////        if (debug) {
////            cout << "Evaluating " << p << "'s neighbors : ";
////        }
        for (int j = 0; j < potentialNeighborsOfP; ++j) {
            int const neighborOfP(neighborsInP[p][j]);
            int const neighborLocation(vertexLookup[neighborOfP]);
            if (InP(neighborOfP)) {
                vMarkedNeighbors[neighborOfP] = true;
////                if (debug)
////                cout << neighborOfP << " ";
            }
////            else
////                break;
        }

////        if (debug)
////        cout << endl;

        // for each vertex x in X: check that all of x's neighbors in P are neighbors of p
        for (int j = beginX; j < beginP; j++) {
            int const x(vertexSets[j]);
////            if (debug)
////                cout << "Does " << x << " dominate?" << endl;
            bool dominated(true);
            int const potentialNeighborsOfX = numNeighbors[x]; // TODO/DS: optimize. Can we stop at the first one that isn't in P?
////            if (debug)
////                cout << "Neighbors: ";
            for (int k = 0; k < potentialNeighborsOfX; ++k) { // TODO/DS: optimize this... Can we stop at the first one that isn't in P?
                int const neighborOfX(neighborsInP[x][k]);
                if (InP(neighborOfX)) { // in P
////                    if (debug)
////                        cout << neighborOfX << " ";
                    dominated = vMarkedNeighbors[neighborOfX];
                    if (!dominated) {
////                        if (debug)
////                        cout << endl << "    Does not dominate!" << endl;
                        break;
                    }
                }
////                else {
////                    break;
////                }
            } // for neighbors of p

            if (dominated) {
////                cout << "Dominated: " << p << " : ";
////                for (int k = 0; k < potentialNeighborsOfP; ++k) {
////                    int const neighborOfP(neighborsInP[p][k]);
////                    if (InP(neighborOfP))
////                        cout << neighborOfP << " ";
////                }
////                cout << endl;
////
////                cout << "Dominator: " << x << " : ";
////                for (int k = 0; k < potentialNeighborsOfX; ++k) {
////                    int const neighborOfX(neighborsInP[x][k]);
////                    if (InP(neighborOfX))
////                        cout << neighborOfX << " ";
////                }
////                cout << endl;
////                lookForDominatedNodes = true;
                // swap p from P to D
                beginR--;
                vertexSets[vertexLookup[p]] = vertexSets[beginR]; vertexLookup[vertexSets[beginR]] = vertexLookup[p]; // move vertex in beginning of P to p's position
                vertexSets[beginR] = p; vertexLookup[p] = beginR; // move p to beginning of P
                removedVertices.push_back(p);
                i--; // evaluate this position again
////                cout << "Found dominated non-neighbor" << endl;
                break;
            }

        } // for vertices x in X

        // unmark neighbors in the array
        for (int j = 0; j < numNeighbors[p]; ++j) {
            int const neighborOfP(neighborsInP[p][j]);
            ////int const neighborLocation(vertexLookup[neighborOfX]);
            ////if (beginP <= neighborLocation && neighborLocation <= beginR)
            vMarkedNeighbors[neighborOfP] = false;
        }
////        vMarkedNeighbors[p] = false;
    } // for p in P
////    } // while looking for dominated nodes

////    clock_t clockEnd = clock();
////
////    timeTestingDominancy += (clockEnd-clockStart);

#if 0 // TODO/DS: toggle, eventually to 0
    beginR = savedBeginR;
#endif
}

void DegeneracyIndependentSets::RemoveDominatedVerticesFromVector(vector<int> &vVerticesInP)
{
////        clock_t clockStart = clock();
        
////    cout << "X=" << beginD-beginX << ", P=" << beginR-beginP << endl;
    // for each vertex x in X

////    bool lookForDominatedNodes(true);
////    while (lookForDominatedNodes) {

    vector<bool> vMarkedNeighbors(m_AdjacencyList.size(), false);

    int numVertices(vVerticesInP.size());
    for (int i = 0; i < numVertices; i++) {
////        lookForDominatedNodes = false;

        // mark neighbors in an array
        int const p(vVerticesInP[i]);
        int const potentialNeighborsOfP = numNeighbors[p]; // TODO/DS: optimize. Can we stop at the first one that isn't in P?
        vMarkedNeighbors[p] = true;
        for (int j = 0; j < potentialNeighborsOfP; ++j) {
            int const neighborOfP(neighborsInP[p][j]);
            int const neighborLocation(vertexLookup[neighborOfP]);
            if (beginP <= neighborLocation && neighborLocation <= beginR)
                vMarkedNeighbors[neighborOfP] = true;
            else
                break;
        }

        // for each vertex x in X: check that all of x's neighbors in P are neighbors of p
        for (int j = beginX; j < beginP; j++) {
            int const x(vertexSets[j]);
            bool dominated(true);
            if (dominated) {
                int const potentialNeighborsOfX = numNeighbors[x]; // TODO/DS: optimize. Can we stop at the first one that isn't in P?
                for (int k = 0; k < potentialNeighborsOfX; ++k) { // TODO/DS: optimize this... Can we stop at the first one that isn't in P?
                    int const neighborOfX(neighborsInP[x][k]);
                    int const neighborLocation(vertexLookup[neighborOfX]);
                    if (neighborLocation >= beginP && neighborLocation < beginR) { // in P
                        dominated = !vMarkedNeighbors[neighborOfX];
                        if (!dominated) break;
                    } else {
                        break;
                    }
                } // for neighbors of p
            }

            if (dominated) {
////                lookForDominatedNodes = true;
                // swap p from P to D
                numVertices--;
                vVerticesInP[i] = vVerticesInP[numVertices];
                i--; // evaluate this position again
////                cout << "Found dominated non-neighbor" << endl;
                break;
            }

        } // for vertices p in P

        // unmark neighbors in the array
        for (int j = 0; j < numNeighbors[p]; ++j) {
            int const neighborOfP(neighborsInP[p][j]);
            ////int const neighborLocation(vertexLookup[neighborOfX]);
            ////if (beginP <= neighborLocation && neighborLocation <= beginR)
            vMarkedNeighbors[neighborOfP] = false;
        }
        vMarkedNeighbors[p] = false;
    } // for x in X
////    } // while looking for dominated nodes

////    clock_t clockEnd = clock();
////
////    timeTestingDominancy += (clockEnd-clockStart);

#if 0 // TODO/DS: toggle, eventually to 0
    beginR = savedBeginR;
#endif

    vVerticesInP.resize(numVertices);
}

// all dominated vertices are moved to space for R
// this method moves them back to P.
void DegeneracyIndependentSets::ReturnDominatedVertices(std::vector<int> const &vRemovedVertices)
{
    for (int const vertex : vRemovedVertices) {
        int const vertexLocation = vertexLookup[vertex];
        vertexSets[vertexLocation] = vertexSets[beginR];
        vertexSets[beginR] = vertex;
        vertexLookup[vertex] = beginR;
        vertexLookup[vertexSets[vertexLocation]] = vertexLocation;
        beginR++;
    }
}

size_t DegeneracyIndependentSets::RemainingSizeEstimate() const
{
#if 0
    vector<bool> vMarkedVertices(m_AdjacencyList.size(), false);

    int estimate(0);
    for (int i = beginP; i <  beginR; ++i) {
        int const vertex(vertexSets[i]);
        if (vMarkedVertices[vertex]) continue;
        vMarkedVertices[vertex] = true;
        estimate++;
        for (int j = 0; j < numNeighbors[vertex]; ++j) {
            int const neighbor(neighborsInP[vertex][j]);
            if (InP(neighbor) && !vMarkedVertices[neighbor]) {
                vMarkedVertices[neighbor] = true;
                break;
            }
        }
    }

////    cout << "Tighter: " << estimate << ", looser: " << SizeOfP() << endl;

    return estimate;
#else
    return SizeOfP();
#endif
}
