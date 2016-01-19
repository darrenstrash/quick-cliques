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

#include <climits>
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "Tools.h"
#include <list>
#include <vector>
#include "MemoryManager.h"
#include "DegeneracyTools.h"
#include "HybridAlgorithm.h"

using namespace std;

/*! \file HybridAlgorithm.cpp

    \brief This file contains an algorithm for listing all cliques
           according to the "hybrid" algorithm of Eppstein and Strash (SEA 2011).

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly

    See the algorithm's description in http://dx.doi.org/10.1007/978-3-642-17517-6_36

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
    has the most neighbors in P. We accompish this step iterating through each
    vertex in X and P, and marking each vertex with the number of later neighbors
    that are in P. (if the vertex is in P, we add one to the count for the 
    vertex itself as well) and then we iterate over the vertices has the highest count.

*/

HybridAlgorithm::HybridAlgorithm(vector<list<int>> const &adjacencyList)
 : Algorithm("hybrid")
 , m_AdjacencyList(adjacencyList)
 , m_pDegree(nullptr)
{
    m_pDegree = new int[m_AdjacencyList.size()];
    for (size_t i = 0; i < m_AdjacencyList.size(); ++i) {
        m_pDegree[i] = m_AdjacencyList[i].size();
    }
}

HybridAlgorithm::~HybridAlgorithm()
{
    delete[] m_pDegree;
}

long HybridAlgorithm::Run(list<list<int>> &cliques)
{
    return listAllMaximalCliquesHybrid(
                m_AdjacencyList,
                m_pDegree,
                m_AdjacencyList.size());
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

    \param beginX The index where set X begins in vertexSets.
 
    \param beginP The index where set P begins in vertexSets.

    \param beginR The index where set R begins in vertexSets.

    \param orderingArray A array storing the degeneracy order, with each vertex
                         containing arrays of its later and earlier neighbrs in the order.
    
    \param scratch An array that we use for scratch space to avoid allocating and
                   deallocating memory often. We use this space to count how many
                   neighbors each vertex has in P.

    \return the pivot that is chosen from P union X

*/

int findBestPivotNonNeighborsHybrid( int** pivotNonNeighbors, int* numNonNeighbors,
                                     int* vertexSets, int* vertexLookup,
                                     int beginX, int beginP, int beginR,
                                     NeighborListArray** orderingArray, 
                                     int* scratch)
{

    int pivot = -1;
    int maxIntersectionSize = -1;

    // count the number of later neighbors each 
    // vertex in X has in P
    int j = beginX;
    while(j<beginP)
    {
        int vertex = vertexSets[j];

        int k = 0;
        while(k < orderingArray[vertex]->laterDegree)
        {
            int neighbor = orderingArray[vertex]->later[k];
            int neighborLocation = vertexLookup[neighbor];

            if(neighborLocation >= beginP && neighborLocation < beginR)
                scratch[vertex]++;

            k++;
        }

        j++;
    }

    // count the number of neighbors each vertex in P has in P
    // and also increment counts for earlier neighbors in X
    // missed by the previous loop
    int candidateToEvaluate = beginP;
    while(candidateToEvaluate < beginR)
    {
        int vertex = vertexSets[candidateToEvaluate];

        int j = 0;
        while(j < orderingArray[vertex]->laterDegree)
        {
            int neighbor = (orderingArray[vertex]->later)[j];
            int neighborLocation = vertexLookup[neighbor];
           
            if(neighborLocation >= beginP && neighborLocation < beginR)
            {
                scratch[vertex]++; 
                scratch[neighbor]++; 
            }
            else if(neighborLocation >= beginX && neighborLocation < beginP)
            {
                scratch[neighbor]++; 
            }
            
            j++;
        }

        candidateToEvaluate++;
    }

    // loop through all vertices in P union X
    candidateToEvaluate = beginX;
    while(candidateToEvaluate < beginR)
    {
        int vertex = vertexSets[candidateToEvaluate];

        // if vertex has more neighbors in P, then update the pivot
        if(scratch[vertex] > maxIntersectionSize)
        {
            pivot = vertex;
            maxIntersectionSize = scratch[vertex];
        }

        // reset scratch to 0
        scratch[vertex] = 0;

        candidateToEvaluate++;
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

    // mark later neighbors of pivot that are in P.
    j = 0;
    while(j<orderingArray[pivot]->laterDegree)
    {
        int neighbor = orderingArray[pivot]->later[j];
        int neighborLocation = vertexLookup[neighbor];

        if(neighborLocation >= beginP && neighborLocation < beginR)
        {
            (*pivotNonNeighbors)[neighborLocation-beginP] = -1;
        }

        j++;
    }

    // move non-neighbors of pivot in P to the beginning of
    // pivotNonNeighbors and set numNonNeighbors appriopriately.

    // if a vertex is marked as a neighbor, the we move it
    // to the end of pivotNonNeighbors.

    // for other vertices, we iterate over their later neighbors
    // to determine if they have the pivot as a later neighbor.
    j = 0;
    while(j<*numNonNeighbors)
    {
        int vertex = (*pivotNonNeighbors)[j];
        int incrementJ = 1;

        if(vertex == -1)
        {
            (*numNonNeighbors)--;
            (*pivotNonNeighbors)[j] = (*pivotNonNeighbors)[*numNonNeighbors];
            continue;
        }

        if(orderingArray[vertex]->orderNumber < orderingArray[pivot]->orderNumber)
        {
            int k = 0;
            while(k<orderingArray[vertex]->laterDegree)
            {
                if(orderingArray[vertex]->later[k]==pivot)
                {
                    (*numNonNeighbors)--;
                    (*pivotNonNeighbors)[j] = (*pivotNonNeighbors)[*numNonNeighbors];
                    incrementJ = 0;
                }

                k++;
            }
        }

        if(incrementJ)
            j++;
    }

    return pivot; 
}

/*! \brief Move vertex to R, set P to vertex's later neighbors and
           set X to vertex's earlier neighbors.

    \param vertex The vertex to move to R.

    \param orderNumber The position of vertex in the ordering.

    \param vertexSets An array containing sets of vertices divided into sets X, P, and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param orderingArray A degeneracy order of the input graph.
 
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

inline void fillInPandXForRecursiveCallHybrid( int vertex, int orderNumber,
                                               int* vertexSets, int* vertexLookup, 
                                               NeighborListArray** orderingArray,
                                               int* pBeginX, int *pBeginP, int *pBeginR, 
                                               int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{
    int vertexLocation = vertexLookup[vertex];

    // move vertex into R section of vertexSets
    (*pBeginR)--;
    vertexSets[vertexLocation] = vertexSets[*pBeginR];
    vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
    vertexSets[*pBeginR] = vertex;
    vertexLookup[vertex] = *pBeginR;

    *pNewBeginX = 0;
    *pNewBeginP = *pNewBeginX;

    // swap earlier neighbors of vertex into X section of vertexSets
    int j = 0;
    while(j<orderingArray[orderNumber]->earlierDegree)
    {
        int neighbor = orderingArray[orderNumber]->earlier[j];
        int neighborLocation = vertexLookup[neighbor];

        vertexSets[neighborLocation] = vertexSets[*pNewBeginP];
        vertexLookup[vertexSets[*pNewBeginP]] = neighborLocation;
        vertexSets[*pNewBeginP] = neighbor;
        vertexLookup[neighbor] = *pNewBeginP;
        (*pNewBeginP)++;

        j++; 

    }

    *pNewBeginR = *pNewBeginP;

    // swap later neighbors of vertex into P section of vertexSets
    j = 0;
    while(j<orderingArray[orderNumber]->laterDegree)
    {
        int neighbor = orderingArray[orderNumber]->later[j];
        int neighborLocation = vertexLookup[neighbor];

        vertexSets[neighborLocation] = vertexSets[*pNewBeginR];
        vertexLookup[vertexSets[*pNewBeginR]] = neighborLocation;
        vertexSets[*pNewBeginR] = neighbor;
        vertexLookup[neighbor] = *pNewBeginR;

        (*pNewBeginR)++;

        j++; 
    }
}

/*! \brief List all maximal cliques in a given graph using the "hybrid" algorithm
           by Eppstein et al. (SEA 2011).

    \param adjList An array of linked lists, representing the input graph in the
                   "typical" adjacency list format.
 
    \param degree An array, indexed by vertex, containing the degree of that vertex. (not currently used)

    \param size The number of vertices in the graph.

    \return the number of maximal cliques of the input graph.

*/

long HybridAlgorithm::listAllMaximalCliquesHybrid( vector<list<int>> const &adjList, 
                                  int* degree, 
                                  int size)
{
    // vertex sets are stored in an array like this:
    // |--X--|--P--|
    int* vertexSets = (int*)Calloc(size, sizeof(int));

    // vertex i is stored in vertexSets[vertexLookup[i]]
    int* vertexLookup = (int*)Calloc(size, sizeof(int));

    // scratch array to make pivot computation faster
    int* scratch = (int*)Calloc(sizeof(int), size);

    // compute the degeneracy order
    NeighborListArray** orderingArray = computeDegeneracyOrderArray(adjList, size);

    int i = 0;

    while(i<size)
    {
        vertexLookup[i] = i;
        vertexSets[i] = i;
        i++;
    }

    int beginX = 0;
    int beginP = 0;
    int beginR = size;

    long cliqueCount = 0;

    list<int> partialClique;

    // for each vertex in order
    for(i=0;i<size;i++)
    {
        int const vertex = orderingArray[i]->vertex;

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("%d ", vertex);
        #endif

        // add vertex to partial clique R
        partialClique.push_back(vertex);

        int newBeginX, newBeginP, newBeginR;

        // set P to be later neighbors and X to be be earlier neighbors
        // of vertex
        fillInPandXForRecursiveCallHybrid( i, vertex, 
                                           vertexSets, vertexLookup, 
                                           orderingArray,
                                           &beginX, &beginP, &beginR, 
                                           &newBeginX, &newBeginP, &newBeginR);

        // recursively compute maximal cliques containing vertex, some of its
        // later neighbors, and avoiding earlier neighbors
        listAllMaximalCliquesHybridRecursive( &cliqueCount,
                                              partialClique, 
                                              orderingArray,
                                              vertexSets, vertexLookup,
                                              newBeginX, newBeginP, newBeginR,
                                              scratch);

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        beginR = beginR + 1;

        partialClique.pop_back();
    }

    partialClique.clear();

    Free(vertexSets);
    Free(vertexLookup);
    Free(scratch);

    for(i = 0; i<size; i++)
    {
        delete orderingArray[i];
    }

    Free(orderingArray);

    return cliqueCount;
}

/*! \brief Move a vertex to the set R, and update the sets P and X

    \param vertex The vertex to move to R.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R, and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param orderingArray A degeneracy order of the input graph.
 
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

inline void moveToRHybrid( int vertex, 
                           int* vertexSets, int* vertexLookup, 
                           NeighborListArray** orderingArray,
                           int* pBeginX, int *pBeginP, int *pBeginR, 
                           int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{
    int vertexLocation = vertexLookup[vertex];

    // initially new sets X, P and R to be empty
    *pNewBeginX = *pBeginP;
    *pNewBeginP = *pBeginP;
    *pNewBeginR = *pBeginP;

    // swap vertex into R and update beginR
    (*pBeginR)--;
    vertexSets[vertexLocation] = vertexSets[*pBeginR];
    vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
    vertexSets[*pBeginR] = vertex;
    vertexLookup[vertex] = *pBeginR;


    // for each vertex in X determine it if has vertex
    // as a later neighbor, if so, move it to new X
    int j = *pBeginX;
    while(j<*pNewBeginX)
    {
        int neighbor = vertexSets[j];
        int neighborLocation = j;

        int incrementJ = 1;
        int k = 0;
        while(k<orderingArray[neighbor]->laterDegree)
        {
            //fprintf(stderr, "got here\n");
            if(orderingArray[neighbor]->later[k] == vertex)
            {
                (*pNewBeginX)--;
                vertexSets[neighborLocation] = vertexSets[*pNewBeginX];
                vertexLookup[vertexSets[*pNewBeginX]] = neighborLocation;
                vertexSets[*pNewBeginX] = neighbor;
                vertexLookup[neighbor] = *pNewBeginX;
                incrementJ=0;
            }
                
            k++;
        }

        if(incrementJ) j++;
    }

    // for each vertex in P determine it if has vertex
    // as a later neighbor, if so, move it to new P
    j = *pBeginP;
    while(j<*pBeginR)
    {
        int neighbor = vertexSets[j];
        int neighborLocation = j;

        int k = 0;
        while(k<orderingArray[neighbor]->laterDegree)
        {
            if(orderingArray[neighbor]->later[k] == vertex)
            {
                vertexSets[neighborLocation] = vertexSets[*pNewBeginR];
                vertexLookup[vertexSets[*pNewBeginR]] = neighborLocation;
                vertexSets[*pNewBeginR] = neighbor;
                vertexLookup[neighbor] = *pNewBeginR;
                (*pNewBeginR)++;
            }

            k++;
        }

        j++;
    }

    // vertex, determine if it has later neighbors in P (X)
    // if so, move them into new P (new X)
    j = 0;
    while(j<orderingArray[vertex]->laterDegree)
    {
        int neighbor = orderingArray[vertex]->later[j];
        int neighborLocation = vertexLookup[neighbor];

        // if later neghbor in P
        if(neighborLocation >= *pBeginP && neighborLocation < *pBeginR)
        {
            // swap neighbor into new P territory
            vertexSets[neighborLocation] = vertexSets[*pNewBeginR];
            vertexLookup[vertexSets[*pNewBeginR]] = neighborLocation;
            vertexSets[*pNewBeginR] = neighbor;
            vertexLookup[neighbor] = *pNewBeginR;
            (*pNewBeginR)++;
        }

        // if later neghbor in X
        if(neighborLocation >= *pBeginX && neighborLocation < *pBeginP)
        {
            // swap neighbor into new X territory
            (*pNewBeginX)--;
            vertexSets[neighborLocation] = vertexSets[*pNewBeginX];
            vertexLookup[vertexSets[*pNewBeginX]] = neighborLocation;
            vertexSets[*pNewBeginX] = neighbor;
            vertexLookup[neighbor] = *pNewBeginX;
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

inline void moveFromRToXHybrid( int vertex, 
                                int* vertexSets, int* vertexLookup, 
                                int* pBeginX, int *pBeginP, int *pBeginR )
{
    int vertexLocation = vertexLookup[vertex];

    //swap vertex into X and increment *pBeginP and *pBeginR
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

    \param orderingArray A array storing the degeneracy order, with each vertex
                         containing arrays of its later and earlier neighbrs in the order.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param beginX The index where set X begins in vertexSets.
 
    \param beginP The index where set P begins in vertexSets.

    \param beginR The index where set R begins in vertexSets.

    \param scratch An array that we use for scratch space to avoid allocating and
                   deallocating memory often. We use this space to count how many
                   neighbors each vertex has in P.

 */

void HybridAlgorithm::listAllMaximalCliquesHybridRecursive( long* cliqueCount,
                                           list<int> &partialClique, 
                                           NeighborListArray** orderingArray,
                                           int* vertexSets, int* vertexLookup,
                                           int beginX, int beginP, int beginR,
                                           int* scratch )
{

    // if X is empty and P is empty, return partial clique as maximal
    if(beginX >= beginP && beginP >= beginR)
    {
        (*cliqueCount)++;

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
    findBestPivotNonNeighborsHybrid(&myCandidatesToIterateThrough,
                                    &numCandidatesToIterateThrough,
                                    vertexSets, vertexLookup,
                                    beginX, beginP, beginR,
                                    orderingArray, 
                                    scratch);

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
        partialClique.push_back(vertex);
        list<int>::iterator vertexLink = partialClique.end();
        --vertexLink;

        // swap vertex into R and update all data structures 
        moveToRHybrid( vertex, 
                       vertexSets, vertexLookup, 
                       orderingArray,
                       &beginX, &beginP, &beginR, 
                       &newBeginX, &newBeginP, &newBeginR);

        // recursively compute maximal cliques with new sets R, P and X
        listAllMaximalCliquesHybridRecursive(cliqueCount,
                                              partialClique, 
                                              orderingArray,
                                              vertexSets, vertexLookup,
                                              newBeginX, newBeginP, newBeginR,
                                              scratch);

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        // remove vertex from partialCliques
        partialClique.erase(vertexLink);

        moveFromRToXHybrid( vertex, 
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
    }

    // don't need to check for emptiness before freeing, since
    // something will always be there (we allocated enough memory
    // for all of P, which is nonempty)
    Free(myCandidatesToIterateThrough);
}
