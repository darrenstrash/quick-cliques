#ifndef _DJS_DEGENERACY_VERTEX_SETS_H_
#define _DJS_DEGENERACY_VERTEX_SETS_H_

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
#include <set>
#include <vector>
#include <cstring> // memcpy
#include <algorithm>
#include <iostream>

////#define ORDERED // TODO/DS: take out?

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

class DegeneracyVertexSets : public VertexSets {
public:
    DegeneracyVertexSets(std::vector<std::vector<int>> &adjacencyList);
    DegeneracyVertexSets(std::vector<std::vector<int>> &adjacencyList, std::set<int> const &verticesToExclude);
    ~DegeneracyVertexSets();

    DegeneracyVertexSets           (DegeneracyVertexSets const &sets) = delete;
    DegeneracyVertexSets& operator=(DegeneracyVertexSets const &sets) = delete;

    void MoveFromPToR(int const vertexInP) __attribute__((always_inline));
    void MoveFromRToX(int const vertexInP) __attribute__((always_inline));

    void ReturnVerticesToP(std::vector<int> const &vVertices) __attribute__((always_inline));

    std::vector<int> ChoosePivot() const __attribute__((always_inline));
    bool InP(int const vertex) const __attribute__((always_inline));

    bool PIsEmpty() const __attribute__((always_inline));
    bool XAndPAreEmpty() const __attribute__((always_inline));

    virtual size_t SizeOfX() const { return beginP - beginX; }
    virtual size_t SizeOfP() const { return beginR - beginP; }

    virtual void Initialize();

    virtual void PrintSummary(int const line) const;

    virtual bool GetNextTopLevelPartition();

    virtual void GetTopLevelPartialClique(std::list<int> &partialClique)  const {
        if (m_iCurrentTopLevelIndex == 0 || m_iCurrentTopLevelIndex > m_AdjacencyList.size()) return;

#ifdef ORDERED
        for (int i = 0; i < orderingArray.size(); ++i) {
            if (orderingArray[i].orderNumber == m_iCurrentTopLevelIndex-1) {
                partialClique.push_back(orderingArray[i].vertex);
                return;
            }
        }
#else
        partialClique.push_back(orderingArray[m_iCurrentTopLevelIndex-1].vertex);
        return ;
#endif
    }

    virtual size_t GetGraphSize() const { return m_AdjacencyList.size(); }

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
    std::set<int> m_EmptySet;
    std::set<int> const &m_VerticesToExclude;
};

inline void DegeneracyVertexSets::ReturnVerticesToP(std::vector<int> const &vVertices)
{
    for (int const vertex : vVertices) {
////        std::cout << "Returning " << vertex << " to P" << std::endl;
        int const vertexLocation = vertexLookup[vertex];
        beginP--;
        vertexSets[vertexLocation] = vertexSets[beginP];
        vertexSets[beginP] = vertex;
        vertexLookup[vertex] = beginP;
        vertexLookup[vertexSets[vertexLocation]] = vertexLocation;
    }
}

inline void DegeneracyVertexSets::MoveFromPToR(int const vertex)
{
////    std::cout << "Moving " << vertex << " to R" << std::endl;
    int const vertexLocation = vertexLookup[vertex];

    int sizeOfP = beginR - beginP;

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
        int neighbor = vertexSets[j];
        int neighborLocation = j;
////        std::cout << neighbor << " ";

        int incrementJ = 1;

        int numPotentialNeighbors = std::min(sizeOfP, numNeighbors[neighbor]);
////        int numPotentialNeighbors = numNeighbors[neighbor];
////        std::cout << "(later neighbors: " << numNeighbors[neighbor] << ") ";

        int k = 0;
        while(k<numPotentialNeighbors)
        {
////            if (neighbor == 0 && vertex == 11) {
////                std::cout << "Does vertex 0 in X have 11 as a neighbor? : " << std::endl;
////            }

            if (neighborsInP[neighbor][k] == vertex)
            {
////                if (neighbor == 0 && vertex == 11)
////                std::cout << "Yes" << std::endl;
                (newBeginX)--;
                vertexSets[neighborLocation] = vertexSets[(newBeginX)];
                vertexLookup[vertexSets[(newBeginX)]] = neighborLocation;
                vertexSets[(newBeginX)] = neighbor;
                vertexLookup[neighbor] = (newBeginX);
                incrementJ=0;
            } else {
////                std::cout << "No" << std::endl;
            }

            k++;
        }

        if(incrementJ) j++;
    }

    j = beginP;
    while(j<beginR)
    {
        int neighbor = vertexSets[j];
        int neighborLocation = j;

        int numPotentialNeighbors = std::min(sizeOfP, numNeighbors[neighbor]);
////       int numPotentialNeighbors = numNeighbors[neighbor];

        int k = 0;
        while(k<numPotentialNeighbors)
        {
            if(neighborsInP[neighbor][k] == vertex)
            {
                vertexSets[neighborLocation] = vertexSets[(newBeginR)];
                vertexLookup[vertexSets[(newBeginR)]] = neighborLocation;
                vertexSets[(newBeginR)] = neighbor;
                vertexLookup[neighbor] = (newBeginR);
                (newBeginR)++;
            }

            k++;
        }

        j++;
    }

    j = (newBeginX);

    while(j < newBeginR) {
        int thisVertex = vertexSets[j];

        int numPotentialNeighbors = std::min(sizeOfP, numNeighbors[thisVertex]);
////       int numPotentialNeighbors = numNeighbors[thisVertex];

        int numNeighborsInP = 0;

        int k = 0;
        while(k < numPotentialNeighbors)
        {
            int neighbor = neighborsInP[thisVertex][k];
            int neighborLocation = vertexLookup[neighbor];
            if(neighborLocation >= newBeginP && neighborLocation < newBeginR)
            {
                neighborsInP[thisVertex][k] = neighborsInP[thisVertex][numNeighborsInP];
                neighborsInP[thisVertex][numNeighborsInP] = neighbor;
                numNeighborsInP++;
            }
            k++;
        }

        j++;
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

inline void DegeneracyVertexSets::MoveFromRToX(int const vertex)
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

inline std::vector<int> DegeneracyVertexSets::ChoosePivot() const
{
    clock_t clockStart = clock();
    int pivot = -1;
    int maxIntersectionSize = -1;

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

    // iterate over each vertex in P union X 
    // to find the vertex with the most neighbors in P.
    int j = beginX;
    while(j<beginR)
    {
        int vertex = vertexSets[j];
        int numPotentialNeighbors = std::min(beginR - beginP, numNeighbors[vertex]);

        int numNeighborsInP = 0;

        int k = 0;
        while(k<numPotentialNeighbors)
        {
            int neighbor = neighborsInP[vertex][k];
            int neighborLocation = vertexLookup[neighbor];

            if(neighborLocation >= beginP && neighborLocation < beginR)
            {
                numNeighborsInP++;
            } else {
                break;
            }

            k++;
        }

        if(numNeighborsInP > maxIntersectionSize) {
            pivot = vertex;
            maxIntersectionSize = numNeighborsInP;
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
    std::vector<int> pivotNonNeighbors(beginR-beginP,0);
    memcpy(&pivotNonNeighbors[0], &vertexSets[beginP], (beginR-beginP)*sizeof(int));

    // we will decrement numNonNeighbors as we find neighbors
    int numNonNeighbors = beginR-beginP;

    int numPivotNeighbors = std::min(beginR - beginP, numNeighbors[pivot]);

    // mark the neighbors of pivot that are in P.
    j = 0;
    while(j<numPivotNeighbors)
    {
        int neighbor = neighborsInP[pivot][j];
        int neighborLocation = vertexLookup[neighbor];

        if(neighborLocation >= beginP && neighborLocation < beginR)
        {
            pivotNonNeighbors[neighborLocation-beginP] = -1;
        }
        else
        {
            break;
        }

        j++;
    }

    // move non-neighbors of pivot in P to the beginning of
    // pivotNonNeighbors and set numNonNeighbors appriopriately.

    // if a vertex is marked as a neighbor, the we move it
    // to the end of pivotNonNeighbors and decrement numNonNeighbors.
    j = 0;
    while (j<numNonNeighbors) {
        int vertex = pivotNonNeighbors[j];

        if (vertex == -1) {
            numNonNeighbors--;
            pivotNonNeighbors[j] = pivotNonNeighbors[numNonNeighbors];
            continue;
        }

        j++;
    }

    pivotNonNeighbors.resize(numNonNeighbors);

    return pivotNonNeighbors; 
}

inline bool DegeneracyVertexSets::InP(int const vertex) const
{
    int const vertexLocation(vertexLookup[vertex]);
    return (vertexLocation >= beginP && vertexLocation < beginR);
}

inline bool DegeneracyVertexSets::PIsEmpty() const
{
    return (beginP >= beginR);
}

inline bool DegeneracyVertexSets::XAndPAreEmpty() const
{
    return (beginX >= beginP && beginP >= beginR);
}

#endif //_DJS_DEGENERACY_ALGORITHM_H_
