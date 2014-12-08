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
#include "PartialMatchDegeneracyVertexSets.h"
#include "DegeneracyTools.h"
#include "PartialMatchGraph.h"
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

PartialMatchDegeneracyVertexSets::PartialMatchDegeneracyVertexSets(vector<vector<int>> &adjacencyList)
: VertexSets("partial-match-degeneracy")
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

PartialMatchDegeneracyVertexSets::~PartialMatchDegeneracyVertexSets()
{
}

void PartialMatchDegeneracyVertexSets::Initialize()
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

void PartialMatchDegeneracyVertexSets::PrintSummary(int const line) const
{
    cout << line << ": X[size=" << beginP-beginX << "], P[size=" << beginR-beginP << "]" << endl; 
}

bool PartialMatchDegeneracyVertexSets::GetNextTopLevelPartition()
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

    vector<int> vDominatedVertices  = ComputeDominanceStandard();
    vector<int> vDominatedVertices2 = ComputeDominanceGraph();

#if 1
    sort (vDominatedVertices.begin(),  vDominatedVertices.end());
    sort (vDominatedVertices2.begin(), vDominatedVertices2.end());
    cout << "dominated(standard):";
    for (int const dominatedVertex : vDominatedVertices) {
        cout << dominatedVertex << " ";
    }
    cout << endl;

    cout << "dominated(graph   ):";
    for (int const dominatedVertex : vDominatedVertices2) {
        cout << dominatedVertex << " ";
    }
    cout << endl;

    return false;
#endif

    return true;

#ifdef PM_TEST

    std::vector<int> vValues;
    for (int i = 0; i < 5; ++i) {
        vValues.push_back(i);
    }

    std::vector<std::vector<int>> vOrderedValues(1, vValues);

    PartialMatchGraph graph(vValues, vOrderedValues);

    auto Test  = [&graph] (vector<int> const &vTestValues, string const &testName) {
        if (!graph.Contains(vTestValues)) {
            cout << testName << ": Failed" << endl;
        } else {
            cout << testName << ": Passed" << endl;
        }
    };

////    graph.Print();

    std::vector<int> vContiguousTestValues;
    vContiguousTestValues.push_back(2);
    vContiguousTestValues.push_back(3);
    vContiguousTestValues.push_back(4);
    
    std::vector<int> vNonContiguous;
    vNonContiguous.push_back(2);
    vNonContiguous.push_back(4);

    std::vector<int> vSingleStart;
    vSingleStart.push_back(1);

    std::vector<int> vSingleRegular;
    vSingleRegular.push_back(4);

    Test(vContiguousTestValues, "Contiguous     test");
    Test(vValues              , "All            test");
    Test(vNonContiguous       , "Non-Contiguous test");
    Test(vSingleStart         , "Single Start   test");
    Test(vSingleRegular       , "Single Regular test");

    return false;

#endif // RUN_ALGORITHM

}

vector<int> PartialMatchDegeneracyVertexSets::ComputeDominanceStandard()
{
    vector<bool> vMarkedNeighbors(m_AdjacencyList.size(), false);
    vector<int> dominatedVertices;

    int const savedBeginR = beginR;

////        clock_t clockStart = clock();
        
////    cout << "X=" << beginD-beginX << ", P=" << beginR-beginP << endl;
    // for each vertex x in X

////    bool lookForDominatedNodes(true);
////    while (lookForDominatedNodes) {
    for (int i = beginX; i < beginP; i++) {
////        lookForDominatedNodes = false;

        // mark neighbors in an array
        int const x(vertexSets[i]);
        int const potentialNeighborsOfX = min(beginR - beginP, numNeighbors[x]); // TODO/DS: optimize. Can we stop at the first one that isn't in P?
        for (int j = 0; j < potentialNeighborsOfX; ++j) {
            int const neighborOfX(neighborsInP[x][j]);
            int const neighborLocation(vertexLookup[neighborOfX]);
            if (beginP <= neighborLocation && neighborLocation <= beginR)
                vMarkedNeighbors[neighborOfX] = true;
            else
                break;
        }

        // for each vertex p in P: check that all of p's neighbors in P are neighbors of x
        for (int j = beginP; j < beginR; j++) {
            int const p(vertexSets[j]);
            bool dominated(vMarkedNeighbors[p]);
            if (dominated) {
                int const potentialNeighborsOfP = min(beginR - beginP, numNeighbors[p]); // TODO/DS: optimize. Can we stop at the first one that isn't in P?
                for (int k = 0; k < potentialNeighborsOfP; ++k) { // TODO/DS: optimize this... Can we stop at the first one that isn't in P?
                    int const neighborOfP(neighborsInP[p][k]);
                    int const neighborLocation(vertexLookup[neighborOfP]);
                    if (neighborLocation >= beginP && neighborLocation < beginR) { // in P
                        dominated = vMarkedNeighbors[neighborOfP];
                        if (!dominated) break;
                    } else {
                        break;
                    }
                } // for neighbors of p
            }

            // TODO/DS: finish computing dominated vertices for comparison with new partial match graph
            if (dominated) {
                dominatedVertices.push_back(p);
////                lookForDominatedNodes = true;
                // swap p from P to R
                beginR--;
                vertexSets[vertexLookup[p]] = vertexSets[beginR]; vertexLookup[vertexSets[beginR]] = vertexLookup[p]; // move vertex in beginning of P to p's position
                vertexSets[beginR] = p; vertexLookup[p] = beginR; // move p to beginning of P
////                beginP++; // move boundary, now D contains p
                j--; // evaluate this position again
////////                cout << "Found dominated non-neighbor" << endl;
            }

        } // for vertices p in P

        // unmark neighbors in the array
        for (int j = 0; j < numNeighbors[x]; ++j) {
            int const neighborOfX(neighborsInP[x][j]);
            ////int const neighborLocation(vertexLookup[neighborOfX]);
            ////if (beginP <= neighborLocation && neighborLocation <= beginR)
            vMarkedNeighbors[neighborOfX] = false;
        }
    } // for x in X
////    } // while looking for dominated nodes

////    clock_t clockEnd = clock();
////
////    timeTestingDominancy += (clockEnd-clockStart);

    beginR = savedBeginR;

    return dominatedVertices;
}

vector<int> PartialMatchDegeneracyVertexSets::ComputeDominanceGraph()
{
    vector<int> dominatedVertices;

    vector<int> vVertices;
    for (size_t u = beginP; u < beginR; ++u) {
        vVertices.push_back(vertexSets[u]);
    }

    vector<vector<int>> vSortedPaths;

    for (size_t u = beginX; u < beginP; ++u) {
        vector<int> vPath;
        int const vertex(vertexSets[u]);
        for (int const neighbor : orderingArray[vertex].later) {
            if (InP(neighbor)) vPath.push_back(neighbor);
        }
        sort(vPath.begin(), vPath.end());
        for (int const node : vPath) {
            cout << node << " ";
        }
        cout << endl;
        vSortedPaths.emplace_back(std::move(vPath));
    }

    PartialMatchGraph graph(vVertices, vSortedPaths);

    for (size_t u = beginP; u < beginR; ++u) {
        vector<int> vPath;
        int const vertex(vertexSets[u]);
        vPath.push_back(vertex);
        for (size_t v = 0; v < numNeighbors[vertex]; ++v) {
            int const neighbor(neighborsInP[vertex][v]);
            if (InP(neighbor)) vPath.push_back(neighbor);
        }
        sort(vPath.begin(), vPath.end());
        if (graph.Contains(vPath)) {
            cout << "Vertex and neighbors:" << endl;
            for (int const node : vPath) {
                cout << node << " ";
            }
            cout << endl;
            dominatedVertices.push_back(vertex);
            break;
        }

    }

    return dominatedVertices;
}
