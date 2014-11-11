#include "AdjacencyListVertexSets.h"

#include <vector>
#include <iostream>
#include <cstring> // memcpy

using namespace std;

////    std::vector<int> m_VertexSets;
////    std::vector<int> m_VertexLocation;
////    std::list<SetDelineator> m_vDelineators;

AdjacencyListVertexSets::AdjacencyListVertexSets(vector<vector<int>> const &adjacencyList)
: beginX(0)
, beginP(0)
, beginR(adjacencyList.size())
, m_AdjacencyList(adjacencyList)
, vertexSets(adjacencyList.size(), 0)
, vertexLookup(adjacencyList.size(), 0)
, degree(m_AdjacencyList.size(), 0)
, m_lDelineators()
{
    m_lDelineators.reserve(10);
    for (size_t i = 0; i < m_AdjacencyList.size(); ++i) {
        degree[i] = m_AdjacencyList[i].size();
    }

    // indices indicating where each set P, X, R starts in 
    // vertexSets; initially, P contains all the vertices
    for (size_t i = 0; i < m_AdjacencyList.size(); ++i) {
        vertexSets[i] = i;
        vertexLookup[i] = i;
    }
}

void AdjacencyListVertexSets::PrintSummary(int const line) const
{
    cout << line << ": X[size=" << beginP-beginX << "], P[size=" << beginR-beginP << "]" << endl; 
}

AdjacencyListVertexSets::~AdjacencyListVertexSets()
{
}

////void AdjacencyListVertexSets::ReturnVerticesToP(std::vector<int> const &vVertices)
////{
////    for (int const vertex : vVertices) {
////        int const vertexLocation = vertexLookup[vertex];
////        beginP--;
////        vertexSets[vertexLocation] = vertexSets[beginP];
////        vertexSets[beginP] = vertex;
////        vertexLookup[vertex] = beginP;
////        vertexLookup[vertexSets[vertexLocation]] = vertexLocation;
////    }
////}
////
////void AdjacencyListVertexSets::MoveFromPToR(int const vertex)
////{
////    int const vertexLocation = vertexLookup[vertex];
////
////    //swap vertex into R and update beginR
////    beginR--;
////    vertexSets[vertexLocation] = vertexSets[beginR];
////    vertexLookup[vertexSets[beginR]] = vertexLocation;
////    vertexSets[beginR] = vertex;
////    vertexLookup[vertex] = beginR;
////
////    // Make new indices into vertexSets for recursive call
////    // initially the new sets X, P, and R are empty.
////    int newBeginX = beginP;
////    int newBeginP = beginP;
////    int newBeginR = beginP;
////
////    // for each neighbor of vertex, ask if it is in X or P,
////    // if it is, leave it there. Otherwise, swap it out.
////    int j = 0;
////    while (j<degree[vertex]) {
////        int const neighbor = m_AdjacencyList[vertex][j];
////        int const neighborLocation = vertexLookup[neighbor];
////
////        // if in X
////        if (neighborLocation >= beginX && neighborLocation < beginP) {
////            // swap into new X territory
////            newBeginX--;
////            vertexSets[neighborLocation] = vertexSets[newBeginX];
////            vertexLookup[vertexSets[newBeginX]] = neighborLocation;
////            vertexSets[newBeginX] = neighbor;
////            vertexLookup[neighbor] = newBeginX;
////        }
////
////        //if in P
////        else if (neighborLocation >= beginP && neighborLocation < beginR) {
////            // swap into new P territory
////            vertexSets[neighborLocation] = vertexSets[newBeginR];
////            vertexLookup[vertexSets[newBeginR]] = neighborLocation;
////            vertexSets[newBeginR] = neighbor;
////            vertexLookup[neighbor] = newBeginR;
////            newBeginR++;
////        }
////
////        j++;
////    }
////
////    m_lDelineators.emplace_back(VertexSets::SetDelineator(beginX, beginP, beginR));
////
////    beginX = newBeginX;
////    beginP = newBeginP;
////    beginR = newBeginR;
////}
////
////void AdjacencyListVertexSets::MoveFromRToX(int const vertex)
////{
////    SetDelineator const &delineator = m_lDelineators.back();
////
////    beginX = delineator.m_BeginX;
////    beginP = delineator.m_BeginP;
////    beginR = delineator.m_BeginR;
////
////    m_lDelineators.pop_back();
////
////    // the location of vertex may have changed
////    int const vertexLocation = vertexLookup[vertex];
////
////    //swap vertex into X and increment beginP and beginR
////    vertexSets[vertexLocation] = vertexSets[beginP];
////    vertexLookup[vertexSets[beginP]] = vertexLocation;
////    vertexSets[beginP] = vertex;
////    vertexLookup[vertex] = beginP;
////    beginP = beginP + 1;
////    beginR = beginR + 1;
////
////}
////
/////*! \brief Computes the vertex v in P union X that has the most neighbors in P,
////  and places P \ {neighborhood of v} in an array. These are the 
////  vertices to consider adding to the partial clique during the current
////  recursive call of the algorithm.
////
////  \param pivotNonNeighbors  An intially empty vector, which will contain the set 
////  P \ {neighborhood of v} when this function completes.
//// */
////
////int AdjacencyListVertexSets::ChoosePivot(std::vector<int> &pivotNonNeighbors) const
////{
////    int pivot = -1;
////    int maxIntersectionSize = -1;
////    int i = beginX;
////
////    // loop through all vertices in P union X
////    while (i < beginR) {
////        int vertex = vertexSets[i];
////        int neighborCount = 0;
////        int numNeighbors = 0;
////
////        // count the number of neighbors vertex has in P.
////        // only count them if the degree of the vertex
////        // is greater than the the best count so far.
////        int j = 0;
////        if (degree[vertex] > maxIntersectionSize) {
////            // count the number of neighbors vertex has in P.
////            while(j < degree[vertex]) {
////                int const neighbor = m_AdjacencyList[vertex][j];
////
////                int const neighborLocation = vertexLookup[neighbor];
////
////                if (neighborLocation >= beginP && neighborLocation < beginR)
////                    neighborCount++;
////
////                j++;
////            }
////
////            // if vertex has more neighbors in P, then update the pivot
////            if (neighborCount > maxIntersectionSize) {
////                maxIntersectionSize = neighborCount;
////                pivot = vertexSets[i];
////            }
////        }
////
////        i++;
////    }
////
////    // compute non neighbors of pivot by marking its neighbors
////    // and moving non-marked vertices into pivotNonNeighbors.
////    // we must do this because this is an efficient way
////    // to compute non-neighbors of a vertex in 
////    // an adjacency list.
////
////    // we initialize enough space for all of P; this is
////    // slightly space inefficient, but it results in faster
////    // computation of non-neighbors.
////    pivotNonNeighbors.resize(beginR-beginP, 0);
////    memcpy(&pivotNonNeighbors[0], &vertexSets[beginP], (beginR-beginP)*sizeof(int));
////
////    // we will decrement numNonNeighbors as we find neighbors
////    size_t numNonNeighbors = beginR-beginP;
////
////    // mark neighbors of pivot that are in P.
////    int j = 0;
////    while (j < degree[pivot]) {
////        int const neighbor = m_AdjacencyList[pivot][j];
////        int const neighborLocation = vertexLookup[neighbor];
////
////        // if the neighbor is in P, mark it as -1
////        if (neighborLocation >= beginP && neighborLocation < beginR) {
////            pivotNonNeighbors[neighborLocation-beginP] = -1;
////        }
//// 
////        j++;
////    }
////
////    // put the non-neighbors at the beginning of the array
////    // and update numNonNeighbors appropriately
////    i = 0; 
////    while (i < numNonNeighbors) {
////        if (pivotNonNeighbors[i] == -1) {
////            numNonNeighbors--;
////            pivotNonNeighbors[i] = pivotNonNeighbors[numNonNeighbors];
////        } else
////            i++;
////    }
////
////    pivotNonNeighbors.resize(numNonNeighbors);
////
////    return pivot;
////}
////
////bool AdjacencyListVertexSets::InP(int const vertex) const
////{
////    int const vertexLocation(vertexLookup[vertex]);
////    return (vertexLocation >= beginP && vertexLocation < beginR);
////}
////
////bool AdjacencyListVertexSets::PIsEmpty() const
////{
////    return (beginP >= beginR);
////}
////
////bool AdjacencyListVertexSets::XAndPAreEmpty() const
////{
////    return (beginX >= beginP && beginP >= beginR);
////}
