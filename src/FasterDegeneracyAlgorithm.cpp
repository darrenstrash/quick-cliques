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

#include <limits.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>

#include "Tools.h"
#include <list>
#include <vector>
#include <set>
#include <utility>
#include <tuple>
#include <algorithm>
#include "MemoryManager.h"
#include "DegeneracyTools.h"

#include "FasterDegeneracyAlgorithm.h"

using namespace std;

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

static clock_t timeComputingPivot(0);
static clock_t timeMovingFromRtoX(0);
static clock_t timeMovingToR(0);
static clock_t timeMovingXToP(0);
static clock_t timeFillInPX(0);

FasterDegeneracyAlgorithm::FasterDegeneracyAlgorithm(vector<vector<int>> &adjacencyArray)
 : Algorithm("faster-degeneracy")
 , m_AdjacencyArray(adjacencyArray)
{
}

FasterDegeneracyAlgorithm::~FasterDegeneracyAlgorithm()
{
}

long FasterDegeneracyAlgorithm::Run(list<list<int>> &cliques)
{
    return listAllMaximalCliquesFasterDegeneracy(
                m_AdjacencyArray,
#ifdef RETURN_CLIQUES_ONE_BY_ONE
                cliques,
#endif
                m_AdjacencyArray.size());
}


/*! \brief Computes the vertex v in P union X that has the most neighbors in P,
           and places P \ {neighborhood of v} in an array. These are the 
           vertices to consider adding to the partial clique during the current
           recursive call of the algorithm.

    \param pivotNonNeighbors  An intially unallocated pointer, which will contain the set 
                              P \ {neighborhood of v} when this function completes.

    \param numNonNeighbors A pointer to a single integer, which has been preallocated,
                           which will contain the number of elements in pivotNonNeighbors.
 
    \param vertexSets An array containing sets of vertices divided into sets X, P, R, and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param neighborsInP Maps vertices to arrays of neighbors such that 
                        neighbors in P fill the first cells

    \param numNeighbors An the neighbor of neighbors a vertex had in P,
                        the first time this function is called, this bound is 
                        used to keep us from allocating more than linear space.
 
    \param beginX The index where set X begins in vertexSets.
 
    \param beginP The index where set P begins in vertexSets.

    \param beginR The index where set R begins in vertexSets.

*/

inline int findBestPivotNonNeighborsFasterDegeneracy( int** pivotNonNeighbors, int* numNonNeighbors,
                                                int* vertexSets, int* vertexLookup,
                                                vector<vector<int>> &neighborsInP, int* numNeighbors,
                                                int beginX, int beginP, int beginR)
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
////    cout << "total: " << totalEdges << "(from " << (beginP - beginX) << " vertices), unique: " << uniqueEdges.size() << endl;

    // iterate over each vertex in P union X 
    // to find the vertex with the most neighbors in P.
    int j = beginX;
    while(j<beginR)
    {
        int vertex = vertexSets[j];
        int numPotentialNeighbors = min(beginR - beginP, numNeighbors[vertex]);

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
    *pivotNonNeighbors = (int*)Calloc(beginR-beginP, sizeof(int));
    memcpy(*pivotNonNeighbors, &vertexSets[beginP], (beginR-beginP)*sizeof(int));

    // we will decrement numNonNeighbors as we find neighbors
    *numNonNeighbors = beginR-beginP;

    int numPivotNeighbors = min(beginR - beginP, numNeighbors[pivot]);

    // mark the neighbors of pivot that are in P.
    j = 0;
    while(j<numPivotNeighbors)
    {
        int neighbor = neighborsInP[pivot][j];
        int neighborLocation = vertexLookup[neighbor];

        if(neighborLocation >= beginP && neighborLocation < beginR)
        {
            (*pivotNonNeighbors)[neighborLocation-beginP] = -1;
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
    while (j<*numNonNeighbors) {
        int vertex = (*pivotNonNeighbors)[j];

        if (vertex == -1) {
            (*numNonNeighbors)--;
            (*pivotNonNeighbors)[j] = (*pivotNonNeighbors)[*numNonNeighbors];
            continue;
        }

        j++;
    }

    clock_t clockEnd = clock();

    timeComputingPivot += (clockEnd - clockStart);

    return pivot; 
}

#if 0
void fillInPXFast()
{
    *pNewBeginX = 0;
    *pNewBeginP = 0;
    for (int const neighbor : orderingArray[vertex].earlier) {
        vertexSets[*pNewBeginP] = neighbor;
        vertexLookup[neighbor] = (*pNewBeginP)++;
    }

    *pNewBeginR = *pNewBeginP;

    for (int const neighbor : orderingArray[vertex].later) {
        vertexSets[*pNewBeginR] = neighbor;
        vertexLookup[neighbor] = (*pNewBeginR)++;
    }

    vertexSets[*pNewBeginR] = vertex;
    vertexLookup[vertex] = *pNewBeginR;

    // TODO/DS: Finish filling in VertexSets, VertexLookup and NeighborsInP (FAST!).
}
#endif //0

/*! \brief Move vertex to R, set P to vertex's later neighbors and
           set X to vertex's earlier neighbors.

    \param vertex The vertex to move to R.

    \param orderNumber The position of vertex in the ordering.

    \param vertexSets An array containing sets of vertices divided into sets X, P, and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param orderingArray A degeneracy order of the input graph.

    \param neighborsInP Maps vertices to arrays of neighbors such that 
                        neighbors in P fill the first cells

    \param numNeighbors An the neighbor of neighbors a vertex had in P,
                        the first time this function is called, this bound is 
                        used to keep us from allocating more than linear space.
 
    \param pBeginX The index where set X begins in vertexSets.
 
    \param pBeginP The index where set P begins in vertexSets.

    \param pBeginR The index where set R begins in vertexSets.

    \param pNewBeginX After function, contains the new index where set X begins
                      in vertexSets after adding vertex to R.
 
    \param pNewBeginP After function, contains the new index where set P begins
                      in vertexSets after adding vertex to P.

    \param pNewBeginR After function, contains the new index where set R begins
                      in vertexSets after adding vertex to R.
*/

#define FAST

inline void fillInPandXForRecursiveCallFasterDegeneracy( int vertex, int orderNumber,
                                                   int* vertexSets, int* vertexLookup, 
                                                   vector<NeighborListArray> &orderingArray,
                                                   vector<vector<int>> &neighborsInP, int* numNeighbors,
                                                   int* pBeginX, int *pBeginP, int *pBeginR, 
                                                   int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{
    clock_t startClock = clock();

    int j = 0;
#ifdef FAST

    *pNewBeginX = 0;
    *pNewBeginP = 0;
    for (int const neighbor : orderingArray[vertex].earlier) {
        if (vertexLookup[vertexSets[*pNewBeginP]] == *pNewBeginP)
            vertexLookup[vertexSets[*pNewBeginP]] = -1; // make sure no one else claims this position.
        vertexSets[*pNewBeginP] = neighbor;
        vertexLookup[neighbor] = (*pNewBeginP)++;
    }

    *pNewBeginR = *pNewBeginP;

    for (int const neighbor : orderingArray[vertex].later) {
        if (vertexLookup[vertexSets[*pNewBeginR]] == *pNewBeginR)
            vertexLookup[vertexSets[*pNewBeginR]] = -1; // make sure no one else claims this position.
        vertexSets[*pNewBeginR] = neighbor;
        vertexLookup[neighbor] = (*pNewBeginR)++;

////        if (neighbor == 1420)
////            cout << "Adding vertex " << neighbor << " to P" << endl;
    }

    if (vertexLookup[vertexSets[*pNewBeginR]] == *pNewBeginR)
        vertexLookup[vertexSets[*pNewBeginR]] = -1; // make sure no one else claims this position.
    vertexSets[*pNewBeginR] = vertex;
    vertexLookup[vertex] = *pNewBeginR;

////    cout << __LINE__ << ": 1420 in position " << vertexLookup[1420] << endl << flush;
    for (int const neighbor : orderingArray[orderNumber].earlier) {

        int const numNeighborsNeeded(min(*pNewBeginR-*pNewBeginP, orderingArray[neighbor].laterDegree));

        if (neighborsInP[neighbor].size() < numNeighborsNeeded) {
            neighborsInP[neighbor].resize(2*numNeighborsNeeded);
        }

        numNeighbors[neighbor] = 0;

        // fill in NeighborsInP
        for (int const laterNeighbor : orderingArray[neighbor].later) {
            int laterNeighborLocation = vertexLookup[laterNeighbor];
            if (laterNeighborLocation >= *pNewBeginP && laterNeighborLocation < *pNewBeginR) {
                neighborsInP[neighbor][numNeighbors[neighbor]] = laterNeighbor;
                numNeighbors[neighbor]++;
            }
        }
    }

#else
    int const vertexLocation = vertexLookup[vertex];

    (*pBeginR)--;
    vertexSets[vertexLocation] = vertexSets[*pBeginR];
    vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
    vertexSets[*pBeginR] = vertex;
    vertexLookup[vertex] = *pBeginR;

    *pNewBeginR = *pBeginR;
    *pNewBeginP = *pBeginR;

    // swap later neighbors of vertex into P section of vertexSets
    while (j < orderingArray[orderNumber].laterDegree) {
        int neighbor = orderingArray[orderNumber].later[j];
        int neighborLocation = vertexLookup[neighbor];

        (*pNewBeginP)--;

        vertexSets[neighborLocation] = vertexSets[*pNewBeginP];
        vertexLookup[vertexSets[*pNewBeginP]] = neighborLocation;
        vertexSets[*pNewBeginP] = neighbor;
        vertexLookup[neighbor] = *pNewBeginP;

        j++; 
    }

    *pNewBeginX = *pNewBeginP;

    // swap earlier neighbors of vertex into X section of vertexSets
    j = 0;
    while (j<orderingArray[orderNumber].earlierDegree) {
        int neighbor = orderingArray[orderNumber].earlier[j];
        int neighborLocation = vertexLookup[neighbor];

        (*pNewBeginX)--;
        vertexSets[neighborLocation] = vertexSets[*pNewBeginX];
        vertexLookup[vertexSets[*pNewBeginX]] = neighborLocation;
        vertexSets[*pNewBeginX] = neighbor;
        vertexLookup[neighbor] = *pNewBeginX;

        int const numNeighborsNeeded(min(*pNewBeginR-*pNewBeginP,orderingArray[neighbor].laterDegree));

        if (neighborsInP[neighbor].size() < numNeighborsNeeded) {
            neighborsInP[neighbor].resize(2*numNeighborsNeeded);
        }

        numNeighbors[neighbor] = 0;

        // fill in NeighborsInP
        int k = 0;
        while (k<orderingArray[neighbor].laterDegree) {
            int laterNeighbor = orderingArray[neighbor].later[k];
            int laterNeighborLocation = vertexLookup[laterNeighbor];
            if (laterNeighborLocation >= *pNewBeginP && laterNeighborLocation < *pNewBeginR) {
                neighborsInP[neighbor][numNeighbors[neighbor]] = laterNeighbor;
                numNeighbors[neighbor]++;
            }

            k++;
        }

        j++; 
    }
#endif

    // reset numNeighbors and neighborsInP for this vertex
    j = *pNewBeginP;
    while (j<*pNewBeginR) {
        int vertexInP = vertexSets[j];

        numNeighbors[vertexInP] = 0;

        int const numNeighborsNeeded(min( *pNewBeginR-*pNewBeginP, 
                    orderingArray[vertexInP].laterDegree 
                    + orderingArray[vertexInP].earlierDegree));

        if (neighborsInP[vertexInP].size() < numNeighborsNeeded) {
            neighborsInP[vertexInP].resize(2*numNeighborsNeeded);
        }

        j++;
    }

    // count neighbors in P, and fill in array of neighbors
    // in P
    j = *pNewBeginP;
    while (j<*pNewBeginR) {
        int vertexInP = vertexSets[j];

        int k = 0;
        while (k<orderingArray[vertexInP].laterDegree) {
            int laterNeighbor = orderingArray[vertexInP].later[k];
            int laterNeighborLocation = vertexLookup[laterNeighbor];

            if (laterNeighborLocation >= *pNewBeginP && laterNeighborLocation < *pNewBeginR) {
                neighborsInP[vertexInP][numNeighbors[vertexInP]] = laterNeighbor;
                numNeighbors[vertexInP]++;
                neighborsInP[laterNeighbor][numNeighbors[laterNeighbor]] = vertexInP;
                numNeighbors[laterNeighbor]++;
            }

            k++;
        }

        j++;
    }

    clock_t endClock = clock();
    timeFillInPX += (endClock - startClock);

////    cout << __LINE__ << ": 1420 in position " << vertexLookup[1420] << endl << flush;
}

/*! \brief List all maximal cliques in a given graph using the algorithm
  by Eppstein et al. (ISAAC 2010/SEA 2011).

  \param adjList An array of linked lists, representing the input graph in the
  "typical" adjacency list format.

  \param cliques A linked list of cliques to return. <b>(only available when compiled 
  with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

  \param degree An array, indexed by vertex, containing the degree of that vertex. (not currently used)

  \param size The number of vertices in the graph.

  \return the number of maximal cliques of the input graph.
 */

static unsigned long largestDifference(0);
static unsigned long numLargeJumps;
static unsigned long stepsSinceLastReportedClique(0);

long listAllMaximalCliquesFasterDegeneracy( vector<vector<int>> &adjArray, 
#ifdef RETURN_CLIQUES_ONE_BY_ONE
        list<list<int>> &cliques,
#endif
        int size)
{
    // vertex sets are stored in an array like this:
    // |--X--|--P--|
    int* vertexSets = (int*)Calloc(size, sizeof(int));

    // vertex i is stored in vertexSets[vertexLookup[i]]
    int* vertexLookup = (int*)Calloc(size, sizeof(int));

    int* numNeighbors = (int*)Calloc(size, sizeof(int));
    
    // compute the degeneracy order
    clock_t clockStart = clock();
    vector<NeighborListArray> orderingArray = std::move(computeDegeneracyOrderArray(adjArray, size));
    clock_t clockEnd = clock();
    clock_t timeDegeneracyOrder = clockEnd - clockStart;

    vector<vector<int>> neighborsInP(size);

    for (int i = 0; i < neighborsInP.size(); ++i) {
        neighborsInP[i].resize(orderingArray[i].laterDegree);
    }

    int i = 0;
    long cliqueCount = 0;
#if 1
    while (i<size) {
        vertexLookup[i] = i;
        vertexSets[i] = i;
        numNeighbors[i] = 1;
        i++;
    }

    int beginX = 0;
    int beginP = 0;
    int beginR = size;

    list<int> partialClique;

    // for each vertex
    for(i=0;i<size;i++)
    {
        int const vertex = orderingArray[i].vertex;

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("%d ", vertex);
        #endif

        // add vertex to partial clique R
        partialClique.push_back(vertex);

        int newBeginX, newBeginP, newBeginR;

        // set P to be later neighbors and X to be be earlier neighbors
        // of vertex
        fillInPandXForRecursiveCallFasterDegeneracy( i, vertex, 
                                               vertexSets, vertexLookup, 
                                               orderingArray,
                                               neighborsInP, numNeighbors,
                                               &beginX, &beginP, &beginR, 
                                               &newBeginX, &newBeginP, &newBeginR);

        // recursively compute maximal cliques containing vertex, some of its
        // later neighbors, and avoiding earlier neighbors
        listAllMaximalCliquesFasterDegeneracyRecursive( &cliqueCount,
                                                  #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                                  cliques,
                                                  #endif
                                                  partialClique, 
                                                  vertexSets, vertexLookup,
                                                  neighborsInP, numNeighbors,
                                                  newBeginX, newBeginP, newBeginR); 

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        beginR = beginR + 1;

        partialClique.pop_back();
    }
#endif

    cerr << endl;
    cerr << "Largest Difference  : " << largestDifference << endl;
    cerr << "Num     Differences : " << numLargeJumps << endl;
    cerr << "Time Computing Pivot: " << ((double)(timeComputingPivot)/(double)(CLOCKS_PER_SEC)) << endl;
    cerr << "Time Moving R to X  : " << ((double)(timeMovingFromRtoX)/(double)(CLOCKS_PER_SEC)) << endl;
    cerr << "Time Moving   to R  : " << ((double)(timeMovingToR)/(double)(CLOCKS_PER_SEC)) << endl;
    cerr << "Time Moving X to P  : " << ((double)(timeMovingXToP)/(double)(CLOCKS_PER_SEC)) << endl;
    cerr << "Time Making X and P : " << ((double)(timeFillInPX)/(double)(CLOCKS_PER_SEC)) << endl;
    cerr << "Time Degeneracy Ordr: " << ((double)(timeDegeneracyOrder)/(double)(CLOCKS_PER_SEC)) << endl;

    Free(vertexSets);
    Free(vertexLookup);

    Free(numNeighbors);

    return cliqueCount;
}

/*! \brief Move a vertex to the set R, and update sets P and X
           and the arrays of neighbors in P

    \param vertex The vertex to move to R.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R, and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param neighborsInP Maps vertices to arrays of neighbors such that 
                        neighbors in P fill the first cells

    \param numNeighbors An the neighbor of neighbors a vertex had in P,
                        the first time this function is called, this bound is 
                        used to keep us from allocating more than linear space.

    \param pBeginX The index where set X begins in vertexSets.
 
    \param pBeginP The index where set P begins in vertexSets.

    \param pBeginR The index where set R begins in vertexSets.

    \param pNewBeginX After function, contains the new index where set X begins
                      in vertexSets after adding vertex to R.
 
    \param pNewBeginP After function, contains the new index where set P begins
                      in vertexSets after adding vertex to P.

    \param pNewBeginR After function, contains the new index where set R begins
                      in vertexSets after adding vertex to R.
*/

inline void moveToRFasterDegeneracy( int vertex, 
                               int* vertexSets, int* vertexLookup, 
                               vector<vector<int>> &neighborsInP, int* numNeighbors,
                               int* pBeginX, int *pBeginP, int *pBeginR, 
                               int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{

    clock_t clockStart = clock();
    int vertexLocation = vertexLookup[vertex];
////    if (vertex == 11)
////        cout << "Adding vertex " << vertex << " to partial clique." << endl << flush;

    (*pBeginR)--;
    vertexSets[vertexLocation] = vertexSets[*pBeginR];
    vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
    vertexSets[*pBeginR] = vertex;
    vertexLookup[vertex] = *pBeginR;

    // this is not a typo, initially newX is empty
    *pNewBeginX = *pBeginP;
    *pNewBeginP = *pBeginP;
    *pNewBeginR = *pBeginP;

    int sizeOfP = *pBeginR - *pBeginP;

    int j = *pBeginX;
////    cout << "Iterate through vertices in X: ";
    while(j<*pNewBeginX)
    {
        int neighbor = vertexSets[j];
        int neighborLocation = j;
////        cout << neighbor << " ";

        int incrementJ = 1;

        int numPotentialNeighbors = min(sizeOfP, numNeighbors[neighbor]);
////        cout << "(later neighbors: " << numNeighbors[neighbor] << ") ";

        int k = 0;
        while(k<numPotentialNeighbors)
        {
////            if (neighbor == 0 && vertex == 11) {
////                cout << "Does vertex 0 in X have 11 as a neighbor? : " << endl;
////            }

            if (neighborsInP[neighbor][k] == vertex)
            {
////                if (neighbor == 0 && vertex == 11)
////                cout << "Yes" << endl;
                (*pNewBeginX)--;
                vertexSets[neighborLocation] = vertexSets[(*pNewBeginX)];
                vertexLookup[vertexSets[(*pNewBeginX)]] = neighborLocation;
                vertexSets[(*pNewBeginX)] = neighbor;
                vertexLookup[neighbor] = (*pNewBeginX);
                incrementJ=0;
            } else {
////                if (neighbor == 0 && vertex == 11)
////                cout << "No" << endl;
            }

            k++;
        }

        if(incrementJ) j++;
    }

    j = (*pBeginP);
    while(j<(*pBeginR))
    {
        int neighbor = vertexSets[j];
        int neighborLocation = j;

        int numPotentialNeighbors = min(sizeOfP, numNeighbors[neighbor]);

        int k = 0;
        while(k<numPotentialNeighbors)
        {
            if(neighborsInP[neighbor][k] == vertex)
            {
                vertexSets[neighborLocation] = vertexSets[(*pNewBeginR)];
                vertexLookup[vertexSets[(*pNewBeginR)]] = neighborLocation;
                vertexSets[(*pNewBeginR)] = neighbor;
                vertexLookup[neighbor] = (*pNewBeginR);
                (*pNewBeginR)++;
            }

            k++;
        }

        j++;
    }

    j = (*pNewBeginX);

    while(j < *pNewBeginR) {
        int thisVertex = vertexSets[j];

        int numPotentialNeighbors = min(sizeOfP, numNeighbors[thisVertex]);

        int numNeighborsInP = 0;

        int k = 0;
        while(k < numPotentialNeighbors)
        {
            int neighbor = neighborsInP[thisVertex][k];
            int neighborLocation = vertexLookup[neighbor];
            if(neighborLocation >= *pNewBeginP && neighborLocation < *pNewBeginR)
            {
                neighborsInP[thisVertex][k] = neighborsInP[thisVertex][numNeighborsInP];
                neighborsInP[thisVertex][numNeighborsInP] = neighbor;
                numNeighborsInP++;
            }
            k++;
        }

        j++;
    }

    clock_t clockEnd = clock();

    timeMovingToR += (clockEnd - clockStart);
}

/*! \brief Move a vertex from the set R to the set X, and update all necessary pointers
           and arrays of neighbors in P

    \param vertex The vertex to move from R to X.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R, and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param pBeginX The index where set X begins in vertexSets.
 
    \param pBeginP The index where set P begins in vertexSets.

    \param pBeginR The index where set R begins in vertexSets.

*/

inline void moveFromRToXFasterDegeneracy( int vertex, 
                                    int* vertexSets, int* vertexLookup, 
                                    int* pBeginX, int* pBeginP, int* pBeginR )
{
    clock_t clockStart = clock();
    int vertexLocation = vertexLookup[vertex];

    //swap vertex into X and increment beginP and beginR
    vertexSets[vertexLocation] = vertexSets[*pBeginP];
    vertexLookup[vertexSets[*pBeginP]] = vertexLocation;
    vertexSets[*pBeginP] = vertex;
    vertexLookup[vertex] = *pBeginP;

    *pBeginP = *pBeginP + 1;
    *pBeginR = *pBeginR + 1;

    clock_t clockEnd = clock();

    timeMovingFromRtoX += (clockEnd - clockStart);
}

/*! \brief Recursively list all maximal cliques containing all of
           all vertices in R, some vertices in P and no vertices in X.

    \param cliqueCount A pointer to the number of maximal cliques computed 
                       thus far.

    \param cliques A linked list of cliques to return. <b>(only available when compiled 
                   with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

    \param partialClique A linked list storing R, the partial clique for this
                         recursive call. 

    \param vertexSets An array containing sets of vertices divided into sets X, P, R and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param neighborsInP Maps vertices to arrays of neighbors such that 
                        neighbors in P fill the first cells

    \param numNeighbors An the neighbor of neighbors a vertex had in P,
                        the first time this function is called, this bound is 
                        used to keep us from allocating more than linear space.

    \param beginX The index where set X begins in vertexSets.
 
    \param beginP The index where set P begins in vertexSets.

    \param beginR The index where set R begins in vertexSets.

*/

static unsigned long recursionNode(0);

void listAllMaximalCliquesFasterDegeneracyRecursive( long* cliqueCount,
                                               #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                               list<list<int>> &cliques,
                                               #endif
                                               list<int> &partialClique, 
                                               int* vertexSets, int* vertexLookup,
                                               vector<vector<int>> &neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR)
{
    int const currentRecursionNode(recursionNode++);

////    for (int i = beginX; i < beginP; ++i) {
////        cout << vertexSets[i] << " ";
////    }
////    cout << "], ";
////
////    cout << "P[";
////    for (int i = beginP; i < beginR; ++i) {
////        cout << vertexSets[i] << " ";
////    }
////    cout << "], ";
////
////    cout << "R[";
////    for (int const vertex : partialClique) {
////        cout << vertex << " ";
////    }
////    cout << "]" << endl;

////    cout << __LINE__ << ": 1420 in position " << vertexLookup[1420] << endl << flush;
    stepsSinceLastReportedClique++;

    // if X is empty and P is empty, process partial clique as maximal
    if(beginX >= beginP && beginP >= beginR)
    {
////        cout << currentRecursionNode << ": Clique: ";
////        for (int const vertex : partialClique) {
////            cout << vertex << " ";
////        }
////        cout << endl;

        (*cliqueCount)++;

        if (stepsSinceLastReportedClique > partialClique.size()) {
            numLargeJumps++;
            //cout << "steps: " << stepsSinceLastReportedClique << ">" << partialClique.size() << endl;
            if (largestDifference < (stepsSinceLastReportedClique - partialClique.size())) {
                largestDifference = stepsSinceLastReportedClique - partialClique.size();
            }
        }

        stepsSinceLastReportedClique = 0;

        processClique( 
                       #ifdef RETURN_CLIQUES_ONE_BY_ONE
                       cliques,
                       #endif
                       partialClique );
////    cout << currentRecursionNode << ": " << endl; // <<  "X[";
        return;
    }

    // avoid work if P is empty.
    if(beginP >= beginR) {
////    cout << currentRecursionNode << ": " << endl; // <<  "X[";
        return;
    }

////    cout << currentRecursionNode << ": "; // <<  "X[";

    int* myCandidatesToIterateThrough;
    int numCandidatesToIterateThrough;

    // get the candidates to add to R to make a maximal clique
    findBestPivotNonNeighborsFasterDegeneracy( &myCandidatesToIterateThrough,
                                         &numCandidatesToIterateThrough,
                                         vertexSets, vertexLookup,
                                         neighborsInP, numNeighbors,
                                         beginX, beginP, beginR);

    // add candiate vertices to the partial clique one at a time and 
    // search for maximal cliques
    if(numCandidatesToIterateThrough != 0)
    {
    int iterator = 0;
    while(iterator < numCandidatesToIterateThrough)
    {
        // vertex to be added to the partial clique
        int vertex = myCandidatesToIterateThrough[iterator];

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("%d ", vertex);
        #endif

        int newBeginX, newBeginP, newBeginR;

        // add vertex into partialClique, representing R.
////        if (vertex == 11) {
////            cout << "Adding 11 to R" << endl;
////        }

        partialClique.push_back(vertex);

////        cout << __LINE__ << ": 1420 in position " << vertexLookup[1420] << endl << flush;
        // swap vertex into R and update all data structures 
        moveToRFasterDegeneracy( vertex, 
                           vertexSets, vertexLookup, 
                           neighborsInP, numNeighbors,
                           &beginX, &beginP, &beginR, 
                           &newBeginX, &newBeginP, &newBeginR);
////        cout << __LINE__ << ": 1420 in position " << vertexLookup[1420] << endl << flush;

        // recursively compute maximal cliques with new sets R, P and X
        listAllMaximalCliquesFasterDegeneracyRecursive( cliqueCount,
                                                  #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                                  cliques,
                                                  #endif
                                                  partialClique, 
                                                  vertexSets, vertexLookup,
                                                  neighborsInP, numNeighbors,
                                                  newBeginX, newBeginP, newBeginR);

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        // remove vertex from partialClique
        partialClique.pop_back();

////        cout << __LINE__ << ": 1420 in position " << vertexLookup[1420] << endl << flush;
        moveFromRToXFasterDegeneracy( vertex, 
                                vertexSets, vertexLookup,
                                &beginX, &beginP, &beginR );

////        cout << __LINE__ << ": 1420 in position " << vertexLookup[1420] << endl << flush;
        iterator++;
    }

    // swap vertices that were moved to X back into P, for higher recursive calls.
    iterator = 0;

////    cout << __LINE__ << ": 1420 in position " << vertexLookup[1420] << endl << flush;
    clock_t clockStart = clock();
    while(iterator < numCandidatesToIterateThrough)
    {
        int vertex = myCandidatesToIterateThrough[iterator];
        int vertexLocation = vertexLookup[vertex];

        beginP--;
        vertexSets[vertexLocation] = vertexSets[beginP];
        vertexSets[beginP] = vertex;
        vertexLookup[vertex] = beginP;
        vertexLookup[vertexSets[vertexLocation]] = vertexLocation;

        iterator++;
    }

////    cout << __LINE__ << ": 1420 in position " << vertexLookup[1420] << endl << flush;
    clock_t clockEnd = clock();
    timeMovingXToP += (clockEnd - clockStart);
    }

    // don't need to check for emptiness before freeing, since
    // something will always be there (we allocated enough memory
    // for all of P, which is nonempty)
    Free(myCandidatesToIterateThrough);
}
