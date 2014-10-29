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
#include <time.h>
#include <iostream>
#include <algorithm>

#include "Tools.h"
#include <list>
#include <vector>
#include "MemoryManager.h"
#include "TimeDelayMaxDegreeAlgorithm.h"
#include "MaximalCliquesAlgorithm.h"

using namespace std;

/*! \file TimeDelayMaxDegreeAlgorithm.cpp

    \brief This file contains the main algorithm for listing all cliques
           according to the algorithm of Tomita et al. (TCS 2006) with one
           crucial change: the input graph is represented as an adjacency list.

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly

    See the main algorithm's description in 
    http://dx.doi.org/10.1016/j.tcs.2006.06.015, and a description of this
    variant of the algorithm in http://dx.doi.org/10.1007/978-3-642-20662-7_31

    This is a recursive backtracking algorithm that maintains three 
    sets of vertices, R, a partial clique, P, the common neighbors
    of vertices in R that are candidates to add to the partial clique,
    and X, the set of common neighbors of R that have been listed 
    in a maximal clique with R already.

    The algorithm recursively adds vertices to R from P, then 
    updates the sets P and X to be the new common neighbors of R
    and recurses. When P and X are empty, R is a maximal clique,
    and is reported.

    Updating the sets P and X is done by iterating over the
    neighbors of the new vertex v added to R and testing
    if they are in P or X. Neighbors of v remain in their
    respective sets.

*/

TimeDelayMaxDegreeAlgorithm::TimeDelayMaxDegreeAlgorithm(vector<vector<int>> const &adjacencyList)
 : MaximalCliquesAlgorithm("timedelay-maxdegree")
 , m_AdjacencyList(adjacencyList)
{
}

TimeDelayMaxDegreeAlgorithm::~TimeDelayMaxDegreeAlgorithm()
{
}

long TimeDelayMaxDegreeAlgorithm::Run(list<list<int>> &cliques)
{
    return listAllMaximalCliquesTimeDelayMaxDegree(
                m_AdjacencyList,
#ifdef RETURN_CLIQUES_ONE_BY_ONE
                cliques,
#endif
                m_AdjacencyList.size());
}

/*! \brief List all maximal cliques in a given graph using the algorithm
           by Tomita et al. (TCS 2006), modified to use an adjacency list
           representation of the graph instead of an adjacency matrix. 
 
    \param adjList An array of linked lists, representing the input graph in the
                   "typical" adjacency list format. (not currently used)

    \param adjacencyList an array of arrays, representing the input graph in a more
                         compact and cache-friendly adjacency list format.

    \param cliques A linked list of cliques to return. <b>(only available when compiled 
                   with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

    \param size The number of vertices in the graph.

    \return The number of maximal cliques of the input graph.
*/

static unsigned long largestDifference(0);
static unsigned long numLargeJumps;
static unsigned long stepsSinceLastReportedClique(0);

vector<bool> vMarkedNeighbors;

vector<vector<int>> ComputeInducedSubgraph(vector<vector<int>> const &adjacencyList, int const vertex, int &subgraphSize)
{
    vector<bool> inGraph(adjacencyList.size(), false);
    vector<vector<int>> subgraphInducedByVertex(adjacencyList.size());
    subgraphInducedByVertex[vertex] = adjacencyList[vertex];
    inGraph[vertex] = true;
    for (int const neighbor : subgraphInducedByVertex[vertex]) {
        inGraph[neighbor] = true;
    }

    subgraphSize = subgraphInducedByVertex[vertex].size() + 1;

    // for every other vertex in the subgraph, fill in its neighbors in the induced subgraph
    for (int const neighbor : subgraphInducedByVertex[vertex]) {
        for (int const nNeighbor : adjacencyList[neighbor]) {
            if (inGraph[nNeighbor]) {
                subgraphInducedByVertex[neighbor].push_back(nNeighbor);
            }
        }
    }

    return subgraphInducedByVertex;
}

//Pre-conditions:
// subraph with all undirected edges, containing vertex v, and all of v's neighbors.
// subgraph size ia the number of vertices in this graph
// size of inX is such that any valid vertex can index into it.

//Post-conditions: sets X,D,P,R are propagated correctly with vertex v and its neighbors
// X will contain all neighbors of v that are "true" in inX.
// P will contain all other neighbors.
// D will be empty
// R will contain v

void FillInSetsAndLookupTables(vector<vector<int>> const &subgraphInducedByVertex, int const inducedSubgraphSize, vector<bool> const &inX, int *vertexSets, int *vertexLookup, int &beginX, int &beginD, int &beginP, int &beginR, int const vertex)
{
    // fill in array representing sets X,D,P,R

    // v goes into R:
    vertexSets[inducedSubgraphSize-1] = vertex;
    vertexLookup[vertex] = inducedSubgraphSize-1;
    beginR = inducedSubgraphSize-1;

    beginP = beginR;
    beginD = beginP;

    // prepare X to receive vertices
    beginX = 0;
    int endX(beginX);

    // v's neighbors go into X or P.
    for (int const neighbor : subgraphInducedByVertex[vertex]) {
        if (inX[neighbor])  {
            vertexSets[endX] = neighbor;
            vertexLookup[neighbor] = endX;
            ++endX;
        } else {
            --beginP;
            --beginD;
            vertexSets[beginP] = neighbor;
            vertexLookup[neighbor] = beginP;
        }
    }
}

long listAllMaximalCliquesTimeDelayMaxDegree( vector<vector<int>> const &adjacencyList, 
                                         #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                         list<list<int> &cliques,
                                         #endif
                                         int size)
{

    // vertex sets are stored in an array like this:
    // |--X--|--P--|
    int* vertexSets = (int*)Calloc(size, sizeof(int));

    // vertex i is stored in vertexSets[vertexLookup[i]]
    int* vertexLookup = (int*)Calloc(size,sizeof(int));

    list<int> partialClique;

    for (int i=size-1; i > -1; i--) {
        vertexSets[i] = i;
        vertexLookup[i] = i;
    }

    vMarkedNeighbors.resize(size, false);

    // indices indicating where each set P, X, R starts in 
    // vertexSets; initially, P contains all the vertices
    int beginX = 0; 
    int beginP = 0;
    int beginD = 0;
    int beginR = size;

    long cliqueCount = 0;

    vector<bool> inX(size, false);

    for (int vertex = 0; vertex < size; vertex++) {
        // create a subgraph out of vertex's neighborhood
        int inducedSubgraphSize(0);
        vector<vector<int>> subgraphInducedByVertex = std::move(ComputeInducedSubgraph(adjacencyList, vertex, inducedSubgraphSize));

        FillInSetsAndLookupTables(subgraphInducedByVertex, inducedSubgraphSize, inX, vertexSets, vertexLookup, beginX, beginD, beginP, beginR, vertex);

        partialClique.push_back(vertex);

        listAllMaximalCliquesTimeDelayMaxDegreeRecursive(&cliqueCount,
                                                     #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                                     cliques,
                                                     #endif
                                                     partialClique, 
                                                     subgraphInducedByVertex,
                                                     vertexSets, vertexLookup, inducedSubgraphSize,
                                                     beginX, beginD, beginP, beginR);

        partialClique.pop_back();
        inX[vertex] = true;
   }

    cout << "Largest Difference : " << largestDifference << endl;
    cout << "Num     Differences: " << numLargeJumps << endl;

    Free(vertexSets);
    Free(vertexLookup);

    return cliqueCount;
}

void moveDominatedVerticesFromPtoDMaxDegree(std::vector<std::vector<int>> const &adjacencyList, int* vertexSets,
                                   int* vertexLookup, int size, int const beginX, int const beginD, int &beginP, int &beginR, list<int> &newlyDominatedVertices)
{
////    bool lookForDominated(true);
////    while (lookForDominated) {
////        lookForDominated = false;
    // for each vertex x in X
    for (int i = beginX; i < beginD; i++) {

        // mark neighbors in an array
        int const x(vertexSets[i]);
        for (int const neighborOfX : adjacencyList[x]) {
            vMarkedNeighbors[neighborOfX] = true;
        }

        // for each vertex p in P: check that all of p's neighbors in P are neighbors of x
        for (int j = beginP; j < beginR; j++) {
            int const p(vertexSets[j]);
            bool dominated(vMarkedNeighbors[p]);
            if (dominated) {
                for (int const neighborOfP : adjacencyList[p]) {
                    if (vertexLookup[neighborOfP] >= beginP && vertexLookup[neighborOfP] < beginR) { // in P
                        dominated = vMarkedNeighbors[neighborOfP];
                        if (!dominated) break;
                    }
                } // for neighbors of p
            }

            if (dominated) {
////                lookForDominated = true; // should do another round of domination checking
                // swap p from P to D
                vertexSets[j] = vertexSets[beginP]; vertexLookup[vertexSets[beginP]] = j; // move vertex in beginning of P to p's position
                vertexSets[beginP] = p;             vertexLookup[p] = beginP; // move p to beginning of P
                beginP++; // move boundary, now D contains p
                newlyDominatedVertices.push_back(p);
            }

        } // for vertices p in P

        // unmark neighbors in the array
        for (int const neighborOfX : adjacencyList[x]) {
            vMarkedNeighbors[neighborOfX] = false;
        }
        
    } // for x in X
////    } // while looking for dominated vertices
}

void moveDominatedVerticesFromNonNeighborstoDMaxDegree(std::vector<std::vector<int>> const &adjacencyList, int* vertexSets,
                                   int* vertexLookup, int* nonNeighbors, int &numNonNeighbors, int const start, int const size, int const beginX, int const beginD, int &beginP, int &beginR, list<int> &newlyDominatedVertices)
{
////    cout << "beginX=" << beginX << endl << flush;
    // for each vertex x in X
    for (int i = beginX; i < beginD; i++) {

        // mark neighbors in an array
        int const x(vertexSets[i]);
        for (int const neighborOfX : adjacencyList[x]) {
            vMarkedNeighbors[neighborOfX] = true;
        }

        // for each vertex p in P: check that all of p's neighbors in P are neighbors of x
        for (int j = start; j < numNonNeighbors; j++) {
            int const p(nonNeighbors[j]);
            bool dominated(vMarkedNeighbors[p]);
            if (dominated) {
                for (int const neighborOfP : adjacencyList[p]) {
                    if (vertexLookup[neighborOfP] >= beginP && vertexLookup[neighborOfP] < beginR) { // in P
                        dominated = vMarkedNeighbors[neighborOfP];
                        if (!dominated) break;
                    }
                } // for neighbors of p
            }

            if (dominated) {
                // swap p from P to D
                vertexSets[vertexLookup[p]] = vertexSets[beginP]; vertexLookup[vertexSets[beginP]] = vertexLookup[p]; // move vertex in beginning of P to p's position
                vertexSets[beginP] = p; vertexLookup[p] = beginP; // move p to beginning of P
                beginP++; // move boundary, now D contains p
                newlyDominatedVertices.push_back(p);
                numNonNeighbors--;
                nonNeighbors[j] = nonNeighbors[numNonNeighbors];
                j--; // evaluate this position again
            }

        } // for vertices p in P

        // unmark neighbors in the array
        for (int const neighborOfX : adjacencyList[x]) {
            vMarkedNeighbors[neighborOfX] = false;
        }
        
    } // for x in X
}

/*! \brief Computes the vertex v in P union X that has the most neighbors in P,
           and places P \ {neighborhood of v} in an array. These are the 
           vertices to consider adding to the partial clique during the current
           recursive call of the algorithm.

    \param pivotNonNeighbors  An intially unallocated pointer, which will contain the set 
                              P \ {neighborhood of v} when this function completes.

    \param numNonNeighbors A pointer to a single integer, which has been preallocated,
                           which will contain the number of elements in pivotNonNeighbors.
 
    \param adjacencyList An array of arrays, representing the input graph in a more
                         compact and cache-friendly adjacency list format.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R, and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param size The number of vertices in the graph. 
 
    \param beginX The index where set X begins in vertexSets.
 
    \param beginP The index where set P begins in vertexSets.

    \param beginR The index where set R begins in vertexSets.

 */

int findBestPivotNonNeighborsTimeDelayMaxDegree( int** pivotNonNeighbors, int* numNonNeighbors,
                                            vector<vector<int>> const &adjacencyList,
                                            int* vertexSets, int* vertexLookup, int size,
                                            int beginX, int beginD, int beginP, int beginR )
{


    int pivot = -1;
    int maxIntersectionSize = -1;
    int i = beginX;

    // loop through all vertices in P union X
    while (i < beginR) {
        if (i >= beginD && i < beginP) { // skip nodes in D
            i = beginP;
            continue;
        }

        int vertex = vertexSets[i];
        int neighborCount = 0;
        int numNeighbors = 0;

        // count the number of neighbors vertex has in P.
        // only count them if the degree of the vertex
        // is greater than the the best count so far.
        int j = 0;
        int const degreeOfVertex(adjacencyList[vertex].size());
        if (degreeOfVertex > maxIntersectionSize) {
            // count the number of neighbors vertex has in P.
            while (j < degreeOfVertex) {
                int const neighbor = adjacencyList[vertex][j];

                int neighborLocation = vertexLookup[neighbor];

                if (neighborLocation >= beginP && neighborLocation < beginR)
                    neighborCount++;

                j++;
            }

            // if vertex has more neighbors in P, then update the pivot
            if (neighborCount > maxIntersectionSize) {
                maxIntersectionSize = neighborCount;
                pivot = vertexSets[i];
            }
        }

        i++;
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

    // mark neighbors of pivot that are in P.
    int j = 0;
    int const degreeOfPivot(adjacencyList[pivot].size());
    while (j < degreeOfPivot) {
        int neighbor = adjacencyList[pivot][j];

        int neighborLocation = vertexLookup[neighbor];

        // if the neighbor is in P, mark it as -1
        if (neighborLocation >= beginP && neighborLocation < beginR) {
            (*pivotNonNeighbors)[neighborLocation-beginP] = -1;
        }
 
        j++;
    }

    // put the non-neighbors at the beginning of the array
    // and update numNonNeighbors appropriately
    i = 0; 
    while (i < *numNonNeighbors) {
        if ((*pivotNonNeighbors)[i] == -1) {
            (*numNonNeighbors)--;
            (*pivotNonNeighbors)[i] = (*pivotNonNeighbors)[*numNonNeighbors];
        } else
            i++;
    }

    return pivot; 
}

/*! \brief Recursively list all maximal cliques containing all of
           all vertices in R, some vertices in P and no vertices in X.

    \param cliqueCount A pointer to the number of maximal cliques computed 
                       thus far.

    \param cliques A linked list of cliques to return. <b>(only available when compiled 
                   with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

    \param partialClique A linked list storing R, the partial clique for this
                         recursive call. 

    \param adjacencyList An array of arrays, representing the input graph in a more
                         compact and cache-friendly adjacency list format.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param size The number of vertices in the graph. 
 
    \param beginX The index where set X begins in vertexSets.
 
    \param beginP The index where set P begins in vertexSets.

    \param beginR The index where set R begins in vertexSets.

*/

#define SIZE 1

static bool printDebug(false); 
static unsigned long recursionNode(1);
static unsigned long recursionId(0);

void listAllMaximalCliquesTimeDelayMaxDegreeRecursive( long* cliqueCount,
                                                  #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                                  list<list<int>> &cliques,
                                                  #endif
                                                  list<int> &partialClique, 
                                                  vector<vector<int>> const &adjacencyList,
                                                  int* vertexSets, int* vertexLookup, int size,
                                                  int beginX, int beginD, int beginP, int beginR )
{
    int const currentRecursionNode(recursionNode++);

////    cout << "Starting node " << currentRecursionNode << endl;

    stepsSinceLastReportedClique++;

    // if X is empty and P is empty, return partial clique as maximal
    if (beginX >= beginD && beginP >= beginR) {
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

////        cout << __LINE__  << ": Done in node " << currentRecursionNode << endl;
        return;
    }

    // avoid work if P is empty.
    if (beginP >= beginR) {
////        cout << __LINE__  << ": Done in node " << currentRecursionNode << endl;
        return;
    }

    int* myCandidatesToIterateThrough;
    int numCandidatesToIterateThrough;

    list<int> newlyDominatedVertices;

    //moveDominatedVerticesFromPtoDMaxDegree(adjacencyList, vertexSets, vertexLookup, size, beginX, beginD, beginP, beginR, newlyDominatedVertices);
    //if (beginP >= beginR) return;

    // get the candidates to add to R to make a maximal clique
    findBestPivotNonNeighborsTimeDelayMaxDegree( &myCandidatesToIterateThrough,
                                            &numCandidatesToIterateThrough,
                                            adjacencyList, vertexSets, vertexLookup, size,
                                            beginX, beginD, beginP, beginR);

    moveDominatedVerticesFromNonNeighborstoDMaxDegree(adjacencyList, vertexSets, vertexLookup, myCandidatesToIterateThrough, numCandidatesToIterateThrough, 0, size, beginX, beginD, beginP, beginR, newlyDominatedVertices);

    if (beginP >= beginR) {
        //cout << __LINE__  << ": Done in node " << currentRecursionNode << endl;
        return;
    }

    // add candiate vertices to the partial clique one at a time and 
    // search for maximal cliques
    if (numCandidatesToIterateThrough != 0) {
    int iterator = 0;

    while (iterator < numCandidatesToIterateThrough) {

        // vertex is to be added to the partial clique
        int const vertex         = myCandidatesToIterateThrough[iterator];
        int const vertexLocation = vertexLookup[vertex];

        //swap vertex into R and update beginR
        beginR--;
        vertexSets[vertexLocation] = vertexSets[beginR];
        vertexLookup[vertexSets[beginR]] = vertexLocation;
        vertexSets[beginR] = vertex;
        vertexLookup[vertex] = beginR;

        partialClique.push_back(vertex);

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("%d ", vertex);
        #endif

        // Make new indices into vertexSets for recursive call
        // initially the new sets X, P, and R are empty.
        int newBeginX = beginD;
        int newBeginR = beginP;

        // for each neighbor of vertex, ask if it is in X or P,
        // if it is, leave it there. Otherwise, swap it out.
        int j = 0;
        int const degreeOfVertex(adjacencyList[vertex].size());
        while (j < degreeOfVertex) {

            int const neighbor = adjacencyList[vertex][j];
            int const neighborLocation = vertexLookup[neighbor];

            // if in X
            if (neighborLocation >= beginX && neighborLocation < beginD) {
                // swap into new X territory
                newBeginX--;
                vertexSets[neighborLocation] = vertexSets[newBeginX];
                vertexLookup[vertexSets[newBeginX]] = neighborLocation;
                vertexSets[newBeginX] = neighbor;
                vertexLookup[neighbor] = newBeginX;
            }


            //TODO/DS: need to understand this, I might be missing something with the new dominated set D.
            //if in P
            else if (neighborLocation >= beginP && neighborLocation < beginR) {
                // swap into new P territory
                vertexSets[neighborLocation] = vertexSets[newBeginR];
                vertexLookup[vertexSets[newBeginR]] = neighborLocation;
                vertexSets[newBeginR] = neighbor;
                vertexLookup[neighbor] = newBeginR;
                newBeginR++;
            }

            j++;

        }

        // recursively compute maximal cliques with new sets R, P and X
        listAllMaximalCliquesTimeDelayMaxDegreeRecursive( cliqueCount,
                                                     #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                                     cliques,
                                                     #endif
                                                     partialClique,
                                                     adjacencyList,
                                                     vertexSets, vertexLookup, size,
                                                     newBeginX, beginD, beginP, newBeginR );

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        // remove vertex from partialCliques
        partialClique.pop_back();

        // the location of vertex may have changed
        int const newVertexLocation = vertexLookup[vertex];

        /* ----swap vertex into X and increment beginD, beginP and beginR---- */

        int const firstVertexFromR(vertexSets[beginR]);
        int const firstVertexFromD(vertexSets[beginD]);
        int const firstVertexFromP(vertexSets[beginP]);
        int const dIsEmpty(beginD == beginP);
        int const pIsEmpty(beginP == beginR);

        // move first element in R to fill space of vertex moving to X.

        vertexSets[newVertexLocation] = firstVertexFromR; vertexLookup[vertexSets[newVertexLocation]] = newVertexLocation;

        // move first element in P to end, to fill hole in R.
        if (!pIsEmpty) {
            vertexSets[beginR] = firstVertexFromP; vertexLookup[firstVertexFromP] = beginR;
        }

        beginR++;

        // next, move first vertex from D to end, to fill hole in P, then increase D so that it contains this vertex.
        if (!dIsEmpty) {
            vertexSets[beginP] = firstVertexFromD; vertexLookup[firstVertexFromD] = beginP;
        }

        beginP++;

        // move vertex into X
        vertexSets[beginD] = vertex; vertexLookup[vertex] = beginD;

        beginD++;

        iterator++;
        moveDominatedVerticesFromNonNeighborstoDMaxDegree(adjacencyList, vertexSets, vertexLookup, myCandidatesToIterateThrough, numCandidatesToIterateThrough, iterator, size, beginX, beginD, beginP, beginR, newlyDominatedVertices);

    } // for non-neighbors
    } // if there are non-neighbors

    /*---- swap vertices that were moved to X and D back into P, for higher recursive calls. ---- */

    // swap each (no longer) dominated vertex to end of D, and decrement beginP so it contains those values.
    for (int const dominatedVertex : newlyDominatedVertices) {
        int const dominatedVertexLocation(vertexLookup[dominatedVertex]);
        beginP--;
        // swap last element in D to fill in hole in D left by removing dominated vertex.
        vertexSets[dominatedVertexLocation] = vertexSets[beginP]; vertexLookup[vertexSets[dominatedVertexLocation]] = dominatedVertexLocation;
        vertexSets[beginP] = dominatedVertex; vertexLookup[dominatedVertex] = beginP; // swap previously dominated vertex into P
    }

    // swap vertices added to X into P. Must swap around vertices in D too... otherwise they get overwritten.
    for (int i=0; i < numCandidatesToIterateThrough; i++) {

        int const vertexInX = myCandidatesToIterateThrough[i];
        beginD--;
        beginP--;

        int const lastElementOfD(vertexSets[beginP]);
        int const lastElementofX(vertexSets[beginD]);
        int const vertexLocationInX = vertexLookup[vertexInX];
        bool const dIsEmpty(beginD == beginP);

        vertexSets[vertexLocationInX] = lastElementofX; vertexLookup[lastElementofX] = vertexLocationInX; // move last element of X into the hole in X...

        // move last element of D into the beginning of D.
        if (!dIsEmpty) {
            vertexSets[beginD] = lastElementOfD;
            vertexLookup[lastElementOfD] = beginD;
        }

        vertexSets[beginP] = vertexInX; vertexLookup[vertexInX] = beginP; // move vertex in X into P.
    }

    // don't need to check for emptiness before freeing, since
    // something will always be there (we allocated enough memory
    // for all of P, which is nonempty)
    Free(myCandidatesToIterateThrough);

    stepsSinceLastReportedClique++;
}
