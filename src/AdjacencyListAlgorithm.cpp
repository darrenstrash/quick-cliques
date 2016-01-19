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

#include "Tools.h"
#include <list>
#include <vector>
#include "MemoryManager.h"
#include "AdjacencyListAlgorithm.h"
#include "Algorithm.h"

using namespace std;

/*! \file AdjacencyListAlgorithm.cpp

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

AdjacencyListAlgorithm::AdjacencyListAlgorithm(vector<vector<int>> const &adjacencyList)
 : Algorithm("adjlist")
 , m_AdjacencyList(adjacencyList)
 , m_pDegree(nullptr)
{
    m_pDegree = new int[m_AdjacencyList.size()];
    for (size_t i = 0; i < m_AdjacencyList.size(); ++i) {
        m_pDegree[i] = m_AdjacencyList[i].size();
    }
}

AdjacencyListAlgorithm::~AdjacencyListAlgorithm()
{
    delete[] m_pDegree;
}

long AdjacencyListAlgorithm::Run(std::list<std::list<int>> &cliques)
{
    return listAllMaximalCliquesAdjacencyList(
                m_AdjacencyList,
                m_pDegree,
                m_AdjacencyList.size());
}

/*! \brief List all maximal cliques in a given graph using the algorithm
           by Tomita et al. (TCS 2006), modified to use an adjacency list
           representation of the graph instead of an adjacency matrix. 
 
    \param adjList An array of linked lists, representing the input graph in the
                   "typical" adjacency list format. (not currently used)

    \param adjacencyList an array of arrays, representing the input graph in a more
                         compact and cache-friendly adjacency list format.

    \param degree An array, indexed by vertex, containing the degree of that vertex.

    \param size The number of vertices in the graph.

    \return The number of maximal cliques of the input graph.
*/

////static unsigned long largestDifference(0);
////static unsigned long numLargeJumps;
////static unsigned long stepsSinceLastReportedClique(0);

long AdjacencyListAlgorithm::listAllMaximalCliquesAdjacencyList(vector<vector<int>> const &adjacencyList, 
                                         int* degree, 
                                         int size)
{

    int i;

    // vertex sets are stored in an array like this:
    // |--X--|--P--|
    int* vertexSets = (int*)Calloc(size, sizeof(int));

    // vertex i is stored in vertexSets[vertexLookup[i]]
    int* vertexLookup = (int*)Calloc(size,sizeof(int));

    list<int> partialClique;

    for(i=size-1;i>-1;i--)
    {
        vertexSets[i] = i;
        vertexLookup[i] = i;
    }

    // indices indicating where each set P, X, R starts in 
    // vertexSets; initially, P contains all the vertices
    int beginX = 0;
    int beginP = 0;
    int beginR = size;

    long cliqueCount = 0;

    listAllMaximalCliquesAdjacencyListRecursive(&cliqueCount,
                                                partialClique, 
                                                adjacencyList, degree,
                                                vertexSets, vertexLookup, size,
                                                beginX, beginP, beginR);

////    cerr << "Largest Difference : " << largestDifference << endl;
////    cerr << "Num     Differences: " << numLargeJumps << endl;

    Free(vertexSets);
    Free(vertexLookup);

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
 
    \param adjacencyList An array of arrays, representing the input graph in a more
                         compact and cache-friendly adjacency list format.

    \param degree An array, indexed by vertex, containing the degree of that vertex.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R, and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param size The number of vertices in the graph. 
 
    \param beginX The index where set X begins in vertexSets.
 
    \param beginP The index where set P begins in vertexSets.

    \param beginR The index where set R begins in vertexSets.

 */

int AdjacencyListAlgorithm::findBestPivotNonNeighborsAdjacencyList(int** pivotNonNeighbors, int* numNonNeighbors,
                                            vector<vector<int>> const &adjacencyList, int* degree,
                                            int* vertexSets, int* vertexLookup, int size,
                                            int beginX, int beginP, int beginR)
{

    int pivot = -1;
    int maxIntersectionSize = -1;
    int i = beginX;

    // loop through all vertices in P union X
    while(i < beginR)
    {

        int vertex = vertexSets[i];
        int neighborCount = 0;
        int numNeighbors = 0;

        // count the number of neighbors vertex has in P.
        // only count them if the degree of the vertex
        // is greater than the the best count so far.
        int j = 0;
        if(degree[vertex] > maxIntersectionSize)
        {
            // count the number of neighbors vertex has in P.
            while(j < degree[vertex])
            {
                int const neighbor = adjacencyList[vertex][j];

                int neighborLocation = vertexLookup[neighbor];

                if(neighborLocation >= beginP && neighborLocation < beginR)
                    neighborCount++;

                j++;
            }

            // if vertex has more neighbors in P, then update the pivot
            if(neighborCount > maxIntersectionSize)
            {
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
    while(j < degree[pivot])
    {
        int neighbor = adjacencyList[pivot][j];

        int neighborLocation = vertexLookup[neighbor];

        // if the neighbor is in P, mark it as -1
        if(neighborLocation >= beginP && neighborLocation < beginR)
        {
            (*pivotNonNeighbors)[neighborLocation-beginP] = -1;
        }
 
        j++;
    }

    // put the non-neighbors at the beginning of the array
    // and update numNonNeighbors appropriately
    i = 0; 
    while(i < *numNonNeighbors)
    {
        if((*pivotNonNeighbors)[i] == -1)
        {
            (*numNonNeighbors)--;
            (*pivotNonNeighbors)[i] = (*pivotNonNeighbors)[*numNonNeighbors];
        }
        else
            i++;
    }

    return pivot; 
}

/*! \brief Recursively list all maximal cliques containing all of
           all vertices in R, some vertices in P and no vertices in X.

    \param cliqueCount A pointer to the number of maximal cliques computed 
                       thus far.

    \param partialClique A linked list storing R, the partial clique for this
                         recursive call. 

    \param adjacencyList An array of arrays, representing the input graph in a more
                         compact and cache-friendly adjacency list format.

    \param degree An array, indexed by vertex, containing the degree of that vertex.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param size The number of vertices in the graph. 
 
    \param beginX The index where set X begins in vertexSets.
 
    \param beginP The index where set P begins in vertexSets.

    \param beginR The index where set R begins in vertexSets.

*/

void AdjacencyListAlgorithm::listAllMaximalCliquesAdjacencyListRecursive(long* cliqueCount,
                                                  list<int> &partialClique, 
                                                  vector<vector<int>> const &adjacencyList, int* degree,
                                                  int* vertexSets, int* vertexLookup, int size,
                                                  int beginX, int beginP, int beginR)
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
    findBestPivotNonNeighborsAdjacencyList( &myCandidatesToIterateThrough,
                                            &numCandidatesToIterateThrough,
                                            adjacencyList, degree, vertexSets, vertexLookup, size,
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

        int vertexLocation = vertexLookup[vertex];

        //swap vertex into R and update beginR
        beginR--;
        vertexSets[vertexLocation] = vertexSets[beginR];
        vertexLookup[vertexSets[beginR]] = vertexLocation;
        vertexSets[beginR] = vertex;
        vertexLookup[vertex] = beginR;

        // add vertex into partialClique, representing R.
        partialClique.push_back(vertex);
        list<int>::iterator vertexLink(partialClique.end());
        --vertexLink;

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("%d ", vertex);
        #endif

        // Make new indices into vertexSets for recursive call
        // initially the new sets X, P, and R are empty.
        int newBeginX = beginP;
        int newBeginP = beginP;
        int newBeginR = beginP;

        // for each neighbor of vertex, ask if it is in X or P,
        // if it is, leave it there. Otherwise, swap it out.
        int j = 0;
        while(j<degree[vertex])
        {
            int const neighbor = adjacencyList[vertex][j];
            int const neighborLocation = vertexLookup[neighbor];

            // if in X
            if(neighborLocation >= beginX && neighborLocation < beginP)
            {
                // swap into new X territory
                newBeginX--;
                vertexSets[neighborLocation] = vertexSets[newBeginX];
                vertexLookup[vertexSets[newBeginX]] = neighborLocation;
                vertexSets[newBeginX] = neighbor;
                vertexLookup[neighbor] = newBeginX;
            }

            //if in P
            else if(neighborLocation >= beginP && neighborLocation < beginR)
            {
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
        listAllMaximalCliquesAdjacencyListRecursive(cliqueCount,
                                                    partialClique,
                                                    adjacencyList, degree,
                                                    vertexSets, vertexLookup, size,
                                                    newBeginX, newBeginP, newBeginR);


        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        // remove vertex from partialCliques
        partialClique.erase(vertexLink);

        // the location of vertex may have changed
        vertexLocation = vertexLookup[vertex];

        //swap vertex into X and increment beginP and beginR
        vertexSets[vertexLocation] = vertexSets[beginP];
        vertexLookup[vertexSets[beginP]] = vertexLocation;
        vertexSets[beginP] = vertex;
        vertexLookup[vertex] = beginP;
        beginP = beginP + 1;
        beginR = beginR + 1;

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
    }

    // don't need to check for emptiness before freeing, since
    // something will always be there (we allocated enough memory
    // for all of P, which is nonempty)
    Free(myCandidatesToIterateThrough);

////    stepsSinceLastReportedClique++;
}
