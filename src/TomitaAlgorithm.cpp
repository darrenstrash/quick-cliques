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
#include "Tools.h"
#include "MemoryManager.h"
#include "TomitaAlgorithm.h"
#include "Algorithm.h"

// system includes
#include <list>
#include <vector>
#include <climits>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>

using namespace std;

/*! \file TomitaAlgorithm.cpp

    \brief This file contains the algorithm for listing all cliques
           according to the algorithm of Tomita et al. (TCS 2006).

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly

    See the algorithm's description in 
    http://dx.doi.org/10.1016/j.tcs.2006.06.015  

    This is a recursive backtracking algorithm that maintains three 
    sets of vertices, R, a partial clique, P, the common neighbors
    of vertices in R that are candidates to add to the partial clique,
    and X, the set of common neighbors of R that have been listed 
    in a maximal clique with R already.

    The algorithm recursively adds vertices to R from P, then 
    updates the sets P and X to be the new common neighbors of R
    and recurses. When P and X are empty, R is a maximal clique,
    and is reported.

    Updating the sets P and X is done by testing if the vertices
    in P (X) have thee new vertex v added to R as a neighbor, by
    looking up the edge in the adjacency matrix. 
    Neighbors of v remain in their respective sets.

*/

TomitaAlgorithm::TomitaAlgorithm(char **ppAdjacencyMatrix, int const numVertices)
 : Algorithm("tomita")
 , m_ppAdjacencyMatrix(ppAdjacencyMatrix)
 , m_iNumVertices(numVertices)
{
}

TomitaAlgorithm::~TomitaAlgorithm()
{
}

long TomitaAlgorithm::Run(list<list<int>> &cliques)
{
    return listAllMaximalCliquesMatrix(
                m_ppAdjacencyMatrix,
                m_iNumVertices);
}

/*! \brief List all maximal cliques in a given graph using the algorithm
           by Tomita et al. (TCS 2006). 
 
    \param adjacencyMatrix An input graph in the adjacency matrix format.

    \param numVertices The number of vertices in the graph.

    \return the number of maximal cliques of the input graph.
*/

////static unsigned long largestDifference(0);
////static unsigned long numLargeJumps(0);

long TomitaAlgorithm::listAllMaximalCliquesMatrix( char** adjacencyMatrix,
                                  int    numVertices )
{

    int i;

    // vertex sets are stored in an array like this:
    // |--X--|--P--|
    int* vertexSets = (int*)Calloc(numVertices, sizeof(int));

    // vertex i is stored in vertexSets[vertexLookup[i]]
    int* vertexLookup = (int*)Calloc(numVertices, sizeof(int));

    list<int> partialClique;

    for(i=numVertices-1;i>-1;i--)
    {
        vertexSets[i] = i;
        vertexLookup[i] = i;
    }

    // indices indicating where each set P, X, R starts in 
    // vertexSets; initially, P contains all the vertices
    int beginX = 0;
    int beginP = 0;
    int beginR = numVertices;

    long cliqueCount = 0;

    long stepsSinceLastReportedClique = 0;

    listAllMaximalCliquesMatrixRecursive( &cliqueCount,
                                          partialClique, 
                                          adjacencyMatrix,
                                          vertexSets, vertexLookup, numVertices,
                                          beginX, beginP, beginR, stepsSinceLastReportedClique);

////    cout << "Largest Difference : " << largestDifference << endl;
////    cout << "Num     Differences: " << numLargeJumps << endl;

    Free(vertexSets);
    Free(vertexLookup);

    partialClique.clear();

    return cliqueCount;
}

/*! \brief Computes the vertex v in P union X that has the most neighbors in P,
           and places P \ {neighborhood of v} in an array. These are the 
           vertices to consider adding to the partial clique during the current
           recursive call of the algorithm.

    \param pivotNonNeighbors  An intially unallocated pointer, which will contain the set 
                              P \ {neighborhood of v} when this function completes.

    \param numNonNeighbors A pointer to a single integer, which has been preallocated,
                           which will contain the number of elements in pivotNonNeighbors.
 
    \param adjacencyMatrix The input graph in adjacency matrix format.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R, and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param size The number of vertices in the graph. 
 
    \param beginX The index where set X begins in vertexSets.
 
    \param beginP The index where set P begins in vertexSets.

    \param beginR The index where set R begins in vertexSets.

*/

int TomitaAlgorithm::findBestPivotNonNeighborsMatrix( int** pivotNonNeighbors, int* numNonNeighbors,
                                     char** adjacencyMatrix, 
                                     int* vertexSets, int* vertexLookup, int size,
                                     int beginX, int beginP, int beginR )
{
    int pivot = -1;
    int maxIntersectionSize = -1;

    // loop through all vertices in P union X
    int i = beginX;
    while(i < beginR)
    {
        int vertex = vertexSets[i];
        int neighborCount = 0;

        // count the number of neighbors vertex has in P.
        int j = beginP;
        while(j < beginR) 
        {
            if(adjacencyMatrix[vertexSets[j]][vertex])
                neighborCount++;

            j++;
        }

        // if vertex has more neighbors in P, then update the pivot
        if(neighborCount > maxIntersectionSize)
        {
            maxIntersectionSize = neighborCount;
            pivot = vertex;
        }

        i++;
    }

    // now pivot is the vertex with the most neighbors in P
    *numNonNeighbors = 0;

    // make an array with the chosen pivot's non-neighbors
    if(beginR-beginP-maxIntersectionSize > 0)
    {
        *pivotNonNeighbors = (int*)Calloc(beginR-beginP-maxIntersectionSize, sizeof(int));

        int j = beginP;

        while(j < beginR)
        {
            if(!adjacencyMatrix[pivot][vertexSets[j]])
            {
                (*pivotNonNeighbors)[*numNonNeighbors] = vertexSets[j];
                (*numNonNeighbors)++;
            }

            j++;
        }
    }

    return pivot; 
}

/*! \brief Move a vertex to the set R, and update the sets P and X

    \param vertex The vertex to move to R.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R, and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param adjacencyMatrix An input graph in the adjacency matrix format.
 
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

inline void moveToRMatrix( int vertex, 
                           int* vertexSets, int* vertexLookup, 
                           char** adjacencyMatrix,
                           int* pBeginX, int *pBeginP, int *pBeginR, 
                           int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{
    int vertexLocation = vertexLookup[vertex];

    // swap vertex into R and update beginR
    (*pBeginR)--;
    vertexSets[vertexLocation] = vertexSets[*pBeginR];
    vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
    vertexSets[*pBeginR] = vertex;
    vertexLookup[vertex] = *pBeginR;

    *pNewBeginX = *pBeginP;
    *pNewBeginP = *pBeginP;
    *pNewBeginR = *pBeginP;

    // for each vertex in X, ask if it has vertex as a neighbor,
    // if it does, keep it in X. Otherwise, swap it out.
    int j = *pBeginX;
    while(j<*pNewBeginX)
    {
        int neighbor = vertexSets[j];
        int neighborLocation = j;

        if(adjacencyMatrix[vertex][neighbor])
        {
            // swap into new X territory
            (*pNewBeginX)--;
            vertexSets[neighborLocation] = vertexSets[*pNewBeginX];
            vertexLookup[vertexSets[*pNewBeginX]] = neighborLocation;
            vertexSets[*pNewBeginX] = neighbor;
            vertexLookup[neighbor] = *pNewBeginX;
        }
        else
        {
            j++;
        }

    }

    // for each vertex in P, ask if it has vertex as a neighbor,
    // if it does, keep it in P. Otherwise, swap it out.
    j = *pBeginP;
    while(j<*pBeginR)
    {
        int neighbor = vertexSets[j];
        int neighborLocation = j;

        if(adjacencyMatrix[vertex][neighbor])
        {
            // swap into new P territory
            vertexSets[neighborLocation] = vertexSets[*pNewBeginR];
            vertexLookup[vertexSets[*pNewBeginR]] = neighborLocation;
            vertexSets[*pNewBeginR] = neighbor;
            vertexLookup[neighbor] = *pNewBeginR;

            (*pNewBeginR)++;
        }

        j++;
    }
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

inline void moveFromRToXMatrix( int vertex, 
                                int* vertexSets, int* vertexLookup, 
                                int* pBeginX, int *pBeginP, int *pBeginR )
{
   int vertexLocation = vertexLookup[vertex];

   //swap vertex into X and increment beginP and beginR
   vertexSets[vertexLocation] = vertexSets[*pBeginP];
   vertexLookup[vertexSets[*pBeginP]] = vertexLocation;
   vertexSets[*pBeginP] = vertex;
   vertexLookup[vertex] = *pBeginP;

   *pBeginP = *pBeginP + 1;
   *pBeginR = *pBeginR + 1;
}

/*! \brief Recursively list all maximal cliques containing all of
           all vertices in R, some vertices in P and no vertices in X.

    \param cliqueCount A pointer to the number of maximal cliques computed 
                       thus far.

    \param partialClique A linked list storing R, the partial clique for this
                         recursive call. 

    \param adjacencyMatrix The input graph in adjacency matrix format.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param size The number of vertices in the graph. 
 
    \param beginX The index where set X begins in vertexSets.
 
    \param beginP The index where set P begins in vertexSets.

    \param beginR The index where set R begins in vertexSets.

*/

void TomitaAlgorithm::listAllMaximalCliquesMatrixRecursive( long* cliqueCount,
                                           list<int> &partialClique, 
                                           char** adjacencyMatrix,
                                           int* vertexSets, int* vertexLookup, int size,
                                           int beginX, int beginP, int beginR , long &stepsSinceLastReportedClique)
{
////    stepsSinceLastReportedClique++;
    // if X is empty and P is empty, return partial clique as maximal
    if(beginX >= beginP && beginP >= beginR)
    {
        (*cliqueCount)++;

////        if (stepsSinceLastReportedClique > partialClique.size()) {
////            numLargeJumps++;
////            //cout << "steps: " << stepsSinceLastReportedClique << ">" << partialClique.size() << endl;
////            if (largestDifference < (stepsSinceLastReportedClique - partialClique.size())) {
////                largestDifference = stepsSinceLastReportedClique - partialClique.size();
////            }
////        }
////
////        stepsSinceLastReportedClique = 0;
        ExecuteCallBacks(partialClique);
        processClique(partialClique);

        return;
    }

    // avoid work if P is empty.
    if(beginP >= beginR)
        return;

    int* myCandidatesToIterateThrough;
    int numCandidatesToIterateThrough;

    // get the candidates to add to R to make a maximal clique
    findBestPivotNonNeighborsMatrix( &myCandidatesToIterateThrough,
                                     &numCandidatesToIterateThrough,
                                     adjacencyMatrix,
                                     vertexSets, vertexLookup, size,
                                     beginX, beginP, beginR);

    // add candiate vertices to the partial clique one at a time and 
    // search for maximal cliques
    if(numCandidatesToIterateThrough != 0)
    {
    int iterator = 0;
    while(iterator < numCandidatesToIterateThrough)
    {
        // vertex is to be added to the partial clique
        int vertex = myCandidatesToIterateThrough[iterator];

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("%d ", vertex);
        #endif

        int newBeginX, newBeginP, newBeginR;

        // add vertex into partialClique, representing R.
        partialClique.push_back(vertex);
        list<int>::iterator vertexLink = partialClique.end();
        --vertexLink;

        // swap vertex into R and update all data structures 
        moveToRMatrix( vertex, 
                       vertexSets, vertexLookup, 
                       adjacencyMatrix,
                       &beginX, &beginP, &beginR, 
                       &newBeginX, &newBeginP, &newBeginR);

        // recursively compute maximal cliques with new sets R, P and X
        listAllMaximalCliquesMatrixRecursive( cliqueCount, 
                                              partialClique,
                                              adjacencyMatrix,
                                              vertexSets, vertexLookup, size,
                                              newBeginX, newBeginP, newBeginR, stepsSinceLastReportedClique );

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        // remove vertex from partialCliques
        partialClique.erase(vertexLink);

        moveFromRToXMatrix( vertex, 
                            vertexSets, vertexLookup, 
                            &beginX, &beginP, &beginR );

        iterator++;
    }

    // swap vertices that were moved to X back into P, for higher recursive calls.
    iterator = 0;
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

    Free(myCandidatesToIterateThrough);
    }

    stepsSinceLastReportedClique++;
}
