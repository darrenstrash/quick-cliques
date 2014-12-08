#ifndef _DJS_DEGENERACY_INDEPENDENT_SETS_H_
#define _DJS_DEGENERACY_INDEPENDENT_SETS_H_

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

// local includes
#include "VertexSets.h"
#include "DegeneracyTools.h"

// system includes
#include <list>
#include <vector>
#include <cstring> // memcpy
#include <algorithm>
#include <iostream>

/*! \file FasterDegeneracyAlgorithm.h

    \brief see FasterDegeneracyAlgorithm.cpp

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

class DegeneracyIndependentSets : public VertexSets {
public:
    DegeneracyIndependentSets(std::vector<std::vector<int>> &adjacencyList);
    ~DegeneracyIndependentSets();

    DegeneracyIndependentSets           (DegeneracyIndependentSets const &sets) = delete;
    DegeneracyIndependentSets& operator=(DegeneracyIndependentSets const &sets) = delete;

    void MoveFromPToR(int const vertexInP) __attribute__((always_inline));
    void MoveFromRToX(int const vertexInP) __attribute__((always_inline));

    void ReturnVerticesToP(std::vector<int> const &vVertices) __attribute__((always_inline));

    std::vector<int> ChoosePivot() const __attribute__((always_inline));
    bool InP(int const vertex) const __attribute__((always_inline));

    bool PIsEmpty() const __attribute__((always_inline));
    bool XAndPAreEmpty() const __attribute__((always_inline));

    virtual size_t SizeOfX() const { return beginP - beginX; }
    virtual size_t SizeOfP() const { return beginR - beginP; }

    virtual size_t GetGraphSize() const { return m_AdjacencyList.size(); }

    virtual size_t RemainingSizeEstimate() const;

    void Initialize();

    virtual void PrintSummary(int const line) const;

    bool GetNextTopLevelPartition();

    void GetTopLevelPartialClique(std::list<int> &partialClique) const
    {
        if (m_iCurrentTopLevelIndex == 0 || m_iCurrentTopLevelIndex > GetGraphSize()) return;
        int const orderNumber(m_iCurrentTopLevelIndex-1);
        int vertex(-1);
        for (int i = 0; i < orderingArray.size(); ++i) {
            if (orderingArray[i].orderNumber == orderNumber) {
                vertex = orderingArray[i].vertex;
                partialClique.push_back(vertex);
                return;
            }
        }

        if (vertex == -1) {
            std::cout << __LINE__ << ": ERROR: could not find vertex with order number " << orderNumber << std::endl << std::flush;
            return;
        }

////        std::cout << "Putting vertex " << vertex << " in top level partial clique" << std::endl;
////        std::cout << "Partial clique contains vertex" << partialClique.front() << std::endl;
    }

protected: // methods
    bool VerifyStartConfiguration() const;
    void RemoveDominatedVerticesFromVector(std::vector<int> &vVerticesInP) const;
    void RemoveDominatedVertices(std::vector<int> &vRemovedVertices);
    void ReturnDominatedVertices(std::vector<int> const &vRemovedVertices);

protected: // members
    int beginX;
    int beginP;
    int beginR;
    std::vector<NeighborListArray> orderingArray;
    std::vector<std::vector<int>> neighborsInP;
    std::vector<int> numNeighbors;
    std::vector<std::vector<int>> &m_AdjacencyList;
    std::vector<int> vertexSets;
    std::vector<int> vertexLookup;
    std::vector<SetDelineator> m_lDelineators;
    int m_iCurrentTopLevelIndex;
};

inline void DegeneracyIndependentSets::ReturnVerticesToP(std::vector<int> const &vVertices)
{
    for (int const vertex : vVertices) {
        int const vertexLocation = vertexLookup[vertex];
        beginP--;
        vertexSets[vertexLocation] = vertexSets[beginP];
        vertexSets[beginP] = vertex;
        vertexLookup[vertex] = beginP;
        vertexLookup[vertexSets[vertexLocation]] = vertexLocation;
    }
}

// DONE, need to verify
inline void DegeneracyIndependentSets::MoveFromPToR(int const vertex)
{
////    bool const debug(vertex == 105 || vertex == 0);
////    if (debug) {
////        std::cout << "Adding " << vertex << " to clique: " << std::endl;
////        std::cout << "    0 has neighbors : ";
////        for (int i = 0; i < numNeighbors[0]; ++i) {
////            std::cout << neighborsInP[0][i] << " ";
////        }
////        std::cout << std::endl;
////    }
////    std::cout << "Adding vertex " << vertex << " to R" << std::endl;
    int const vertexLocation = vertexLookup[vertex];

    int const sizeOfP = beginR - beginP;

    //swap vertex into R and update beginR
    beginR--;
    vertexSets[vertexLocation] = vertexSets[beginR];
    vertexLookup[vertexSets[beginR]] = vertexLocation;
    vertexSets[beginR] = vertex;
    vertexLookup[vertex] = beginR;

    // Make new indices into vertexSets for recursive call
    // initially the new sets X, P, and R are empty.
    int newBeginX = beginP;
    int newBeginP = beginP;
    int newBeginR = beginP;

    // for each neighbor of vertex, ask if it is in X or P,
    // if it is, leave it there. Otherwise, swap it out.
    int j = beginX;
////    std::cout << "Iterate through vertices in X: ";
    while(j<newBeginX)
    {
        int const neighbor = vertexSets[j];
        int const neighborLocation = j;
////        std::cout << neighbor << " ";

        int incrementJ = 1;

////        int const numPotentialNeighbors = std::min(sizeOfP, numNeighbors[neighbor]);
        int const numPotentialNeighbors = numNeighbors[neighbor];
////        std::cout << "(later neighbors: " << numNeighbors[neighbor] << ") ";

        int k = 0;
        bool hasVertexAsNeighbor(false);
        while(k<numPotentialNeighbors)
        {
////            if (neighbor == 0 && vertex == 11) {
////                std::cout << "Does vertex 0 in X have 11 as a neighbor? : " << std::endl;
////            }

            if (neighborsInP[neighbor][k] == vertex) {
                hasVertexAsNeighbor = true;
                break;
////                if (neighbor == 0 && vertex == 11)
////                std::cout << "Yes" << std::endl;
            } else {
////                std::cout << "No" << std::endl;
            }

            k++;
        }


        // swap non-neighbors into X for recursion
        if (!hasVertexAsNeighbor) {
////            if (find(orderingArray[vertex].later.begin(), orderingArray[vertex].later.end(), neighbor) != orderingArray[vertex].later.end()  && InP(vertex) ||
////                find(orderingArray[vertex].earlier.begin(), orderingArray[vertex].earlier.end(), neighbor) != orderingArray[vertex].earlier.end() && InP(vertex)) {
////                std::cout << __LINE__ << ": ERROR - vertex " << vertex << " has neighbor " << neighbor << std::endl;
////            }
            newBeginX--;
            vertexSets[neighborLocation] = vertexSets[newBeginX];
            vertexLookup[vertexSets[newBeginX]] = neighborLocation;
            vertexSets[newBeginX] = neighbor;
            vertexLookup[neighbor] = newBeginX;
            incrementJ=0;
        }

        if (incrementJ) j++;
    }

    j = beginP;
    while(j<beginR) {
        int const neighbor = vertexSets[j];
        int const neighborLocation = j;

////        int const numPotentialNeighbors = std::min(sizeOfP, numNeighbors[neighbor]);
        int const numPotentialNeighbors = numNeighbors[neighbor];

        int k = 0;
        bool hasVertexAsNeighbor(false);
        while(k<numPotentialNeighbors) {
            if(neighborsInP[neighbor][k] == vertex) {
                hasVertexAsNeighbor = true;
                break;
            }

            k++;
        }

        // move non-neighbors into P.
        if (!hasVertexAsNeighbor) {
////            if (find(orderingArray[vertex].later.begin(), orderingArray[vertex].later.end(), neighbor) != orderingArray[vertex].later.end() ||
////                find(orderingArray[vertex].earlier.begin(), orderingArray[vertex].earlier.end(), neighbor) != orderingArray[vertex].earlier.end()) {
////                std::cout << __LINE__ << ": ERROR - vertex " << vertex << " has neighbor " << neighbor << std::endl;
////            }
            vertexSets[neighborLocation] = vertexSets[newBeginR];
            vertexLookup[vertexSets[newBeginR]] = neighborLocation;
            vertexSets[newBeginR] = neighbor;
            vertexLookup[neighbor] = newBeginR;
            newBeginR++;
        }

        j++;
    }

    j = newBeginX;

    while(j < newBeginR) {
        int const thisVertex = vertexSets[j];
////        int const numPotentialNeighbors = std::min(sizeOfP, numNeighbors[thisVertex]);
        int const numPotentialNeighbors = numNeighbors[thisVertex];
        int       numNeighborsInP = 0;

        int k = 0;
        while(k < numPotentialNeighbors) {
            int const neighbor = neighborsInP[thisVertex][k];
            int const neighborLocation = vertexLookup[neighbor];
            if(neighborLocation >= newBeginP && neighborLocation < newBeginR) {
                neighborsInP[thisVertex][k] = neighborsInP[thisVertex][numNeighborsInP];
                neighborsInP[thisVertex][numNeighborsInP] = neighbor;
                numNeighborsInP++;
            }
            k++;
        }

        j++;
////        if (numNeighborsInP == 0) {
////            std::cout << __LINE__ << ": could have pruned whole subtree" << std::endl; 
////        }
    }

    m_lDelineators.emplace_back(VertexSets::SetDelineator(beginX, beginP, beginR));

    beginX = newBeginX;
    beginP = newBeginP;
    beginR = newBeginR;
}

/*! \brief Move a vertex from the set R to the set X, and update all necessary pointers
           and arrays of neighbors in P

    \param vertex The vertex to move from R to X.
*/

inline void DegeneracyIndependentSets::MoveFromRToX(int const vertex)
{
    SetDelineator const &delineator = m_lDelineators.back();

    beginX = delineator.m_BeginX;
    beginP = delineator.m_BeginP;
    beginR = delineator.m_BeginR;

    m_lDelineators.pop_back();

    // the location of vertex may have changed
    int const vertexLocation = vertexLookup[vertex];

    //swap vertex into X and increment beginP and beginR
    vertexSets[vertexLocation] = vertexSets[beginP];
    vertexLookup[vertexSets[beginP]] = vertexLocation;
    vertexSets[beginP] = vertex;
    vertexLookup[vertex] = beginP;
    beginP = beginP + 1;
    beginR = beginR + 1;
}

/*! \brief Computes the vertex v in P union X that has the most neighbors in P,
  and places P \ {neighborhood of v} in an array. These are the 
  vertices to consider adding to the partial clique during the current
  recursive call of the algorithm.

  \param pivotNonNeighbors  An intially empty vector, which will contain the set 
  P \ {neighborhood of v} when this function completes.
 */

#define DJS_PIVOT_DONE

inline std::vector<int> DegeneracyIndependentSets::ChoosePivot() const
{
    clock_t clockStart = clock();
    int pivot = -1;
    int maxIntersectionSize = -1;

////    if (PIsEmpty()) return std::vector<int>();

////    std::cout << "beginP=" << beginP << ", beginR=" << beginR << std::endl;

////    set<tuple<int,int,int>> uniqueEdges;
////    int totalEdges(0);
////    for (int i = beginX; i < beginP; i++) {
////        for (int j = 0; j < neighborsInP[vertexSets[i]].size(); ++j) {
////            totalEdges++;
////            if (vertexLookup[neighborsInP[vertexSets[i]][j]] < beginP || vertexLookup[neighborsInP[vertexSets[i]][j]] >= beginR || (j == neighborsInP[vertexSets[i]].size() -1)) {
////                sort(neighborsInP[vertexSets[i]].begin(), neighborsInP[vertexSets[i]].begin() + j);
////                for (int k = 0; k < j-1; ++k) {
////                    if (k == 0)
////                        uniqueEdges.insert(make_tuple(k, neighborsInP[vertexSets[i]][k], neighborsInP[vertexSets[i]][k+1]));
////                    else 
////                        uniqueEdges.insert(make_tuple(1, neighborsInP[vertexSets[i]][k], neighborsInP[vertexSets[i]][k+1]));
////                }
////                if (vertexLookup[neighborsInP[vertexSets[i]][j-1]] >= beginP && vertexLookup[neighborsInP[vertexSets[i]][j]] < beginR) {
////                    if (j == 0)
////                        uniqueEdges.insert(make_tuple(0, neighborsInP[vertexSets[i]][j], -1));
////                    else if (j == 1)
////                        uniqueEdges.insert(make_tuple(0, neighborsInP[vertexSets[i]][j-1], neighborsInP[vertexSets[i]][j]));
////                    else
////                        uniqueEdges.insert(make_tuple(1, neighborsInP[vertexSets[i]][j-1], neighborsInP[vertexSets[i]][j]));
////                }
////                break;
////            }
////        }
////    }
////
////    std::cout << "total: " << totalEdges << "(from " << (beginP - beginX) << " vertices), unique: " << uniqueEdges.size() << std::endl;

////    std::cout << "P: ";
////    for (size_t u = beginP; u < beginR; ++u) {
////        std::cout << vertexSets[u] << " ";
////    }
////    std::cout << std::endl;

#ifdef DJS_PIVOT_DONE
    // iterate over each vertex in P union X 
    // to find the vertex with the most neighbors in P.
    int j = beginX;
    while(j<beginR)
    {
        int vertex = vertexSets[j];
////        int const numPotentialNeighbors = std::min(beginR - beginP, numNeighbors[vertex]);
        int const numPotentialNeighbors = numNeighbors[vertex];
////        std::cout << "evaluating neighbors of " << vertex << std::endl;

        int numNeighborsInP = 0;

        int k = 0;
        while(k<numPotentialNeighbors)
        {
            int neighbor = neighborsInP[vertex][k];
            int neighborLocation = vertexLookup[neighbor];

            if(neighborLocation >= beginP && neighborLocation < beginR)
            {
////                std::cout << neighbor << " is in P!" << std::endl;
                numNeighborsInP++;
            }
////             else {
////                break;
////            }

            k++;
        }

////        std::cout << "Comparing " << (beginR - beginP - numNeighborsInP) << " and " << maxIntersectionSize << std::endl;

        if((beginR - beginP - numNeighborsInP) > maxIntersectionSize) {
////            std::cout << "Updating pivot to " << pivot << std::endl;
            pivot = vertex;
            maxIntersectionSize = beginR - beginP - numNeighborsInP;
        }

        j++;
    }

    // compute non neighbors of pivot by marking its neighbors
    // and moving non-marked vertices into pivotNonNeighbors.
    // we must do this because this is an efficient way
    // to compute non-neighbors of a vertex in 
    // an adjacency list.

    // we initialize enough space for all of P; this is
    // slightly space inefficient, but it results in faster
    // computation of non-neighbors.
////    int numPivotNeighbors = std::min(beginR - beginP, numNeighbors[pivot]);
////    int numPivotNeighbors = beginR - beginP;
////    std::cout << "pivot: " << pivot << std::endl;
    int numPivotNeighbors = numNeighbors[pivot];
    std::vector<int> pivotNonNeighbors(numPivotNeighbors+1,0);

    // we will decrement numNonNeighbors as we find neighbors
    int numNonNeighbors = 0;

    if (InP(pivot)) pivotNonNeighbors[numNonNeighbors++] = pivot;

    // mark the neighbors of pivot that are in P.
    j = 0;
    while(j<numPivotNeighbors) {
        int const neighbor = neighborsInP[pivot][j];
        int const neighborLocation = vertexLookup[neighbor];

        if(neighborLocation >= beginP && neighborLocation < beginR) {
            pivotNonNeighbors[numNonNeighbors++] = neighbor;
        }
////         else {
////            break;
////        }

        j++;
    }

////    if (InP(pivot) && numNonNeighbors == 0) {
////        std::cout << "Could have pruned whole subtree..." << std::endl;
////    } else if (numNonNeighbors == 0) {
////        std::cout << "Did prune whole subtree..." << std::endl;
////    }

    pivotNonNeighbors.resize(numNonNeighbors);

    ////RemoveDominatedVerticesFromVector(pivotNonNeighbors);

    return pivotNonNeighbors;
#else
    std::vector<int> vVerticesInP;
    for (int i = beginP; i < beginR; ++i)
        vVerticesInP.push_back(vertexSets[i]);

    return vVerticesInP;
#endif
}

inline bool DegeneracyIndependentSets::InP(int const vertex) const
{
    int const vertexLocation(vertexLookup[vertex]);
    return (vertexLocation >= beginP && vertexLocation < beginR);
}

inline bool DegeneracyIndependentSets::PIsEmpty() const
{
    return (beginP >= beginR);
}

inline bool DegeneracyIndependentSets::XAndPAreEmpty() const
{
    return (beginX >= beginP && beginP >= beginR);
}

#endif // _DJS_DEGENERACY_INDEPENDENT_SETS_H_
