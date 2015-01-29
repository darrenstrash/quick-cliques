#ifndef ADJACENCY_MATRIX_VERTEX_SETS_MAX_H
#define ADJACENCY_MATRIX_VERTEX_SETS_MAX_H

// local includes
#include "VertexSets.h"
#include "CliqueColoringStrategy.h"

// system includes
#include <vector>
#include <cstring> // memcopy

class AdjacencyMatrixVertexSetsMax : public VertexSets {
public:
    AdjacencyMatrixVertexSetsMax(std::vector<std::vector<char>> const &adjacencyMatrix);
    ~AdjacencyMatrixVertexSetsMax();

    AdjacencyMatrixVertexSetsMax           (AdjacencyMatrixVertexSetsMax const &sets) = delete;
    AdjacencyMatrixVertexSetsMax& operator=(AdjacencyMatrixVertexSetsMax const &sets) = delete;

    void MoveFromPToR(int const vertexInP) __attribute__((always_inline));
    void MoveFromRToX(int const vertexInP) __attribute__((always_inline));

    void ReturnVerticesToP(std::vector<int> const &vVertices) __attribute__((always_inline));

    std::vector<int> ChoosePivotNonConst() __attribute__((always_inline));
    bool InP(int const vertex) const __attribute__((always_inline));

    bool PIsEmpty() const __attribute__((always_inline));
    bool XAndPAreEmpty() const __attribute__((always_inline));

    virtual size_t SizeOfX() const { return beginP - beginX; }
    virtual size_t SizeOfP() const { return beginR - beginP; }

    virtual size_t RemainingSizeEstimate() const
    {
        if (stackColor.empty() || stackColor.back().empty()) return 0;
        return stackColor.back().back();
    }

    virtual size_t GetGraphSize() const { return m_AdjacencyMatrix.size(); }

    virtual int  GetNextVertexToEvaluate(std::vector<int> &/*vVertices*/)
    {
        if (stackP.empty() || stackP.back().empty()) return -1;

        int const nextVertex(stackP.back().back());
        stackP.back().pop_back();
        stackColor.back().pop_back();

        return nextVertex;
    }

    void Initialize();

    void GetTopLevelPartialClique(std::list<int> &/*partialClique*/) const { }

    void PrintSummary(int const line) const;

    bool GetNextTopLevelPartition();

private: // members
    int beginX;
    int beginP;
    int beginR;
    std::vector<std::vector<char>> const &m_AdjacencyMatrix;
    std::vector<int> vertexSets;
    std::vector<int> vertexLookup;
    std::vector<int> degree;
    std::vector<SetDelineator> m_lDelineators;
    std::vector<std::vector<int>> stackP;
    std::vector<std::vector<int>> stackColor;
    CliqueColoringStrategy coloringStrategy;
////    SparseCliqueColoringStrategy coloringStrategy;
};

inline void AdjacencyMatrixVertexSetsMax::ReturnVerticesToP(std::vector<int> const &vVertices)
{
    for (int const vertex : vVertices) {
        int const vertexLocation = vertexLookup[vertex];
        beginP--;
        vertexSets[vertexLocation] = vertexSets[beginP];
        vertexSets[beginP] = vertex;
        vertexLookup[vertex] = beginP;
        vertexLookup[vertexSets[vertexLocation]] = vertexLocation;
    }

    stackP.pop_back();
    stackColor.pop_back();
}

inline void AdjacencyMatrixVertexSetsMax::MoveFromPToR(int const vertex)
{
    int const vertexLocation = vertexLookup[vertex];

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
    int j = 0;
    for (int index = beginP-1; index >= beginX; index--) {
        int const vertexInX = vertexSets[index];
        if (m_AdjacencyMatrix[vertex][vertexInX]) {
            //swap into X
            newBeginX--;
            vertexSets[index] = vertexSets[newBeginX];
            vertexLookup[vertexSets[newBeginX]] = index;
            vertexSets[newBeginX] = vertexInX;
            vertexLookup[vertexInX] = newBeginX;
        }
    }

    for (int index = beginP; index < beginR; index++) {
        int const vertexInP = vertexSets[index];
        if (m_AdjacencyMatrix[vertex][vertexInP]) {
            //swap into P
            vertexSets[index] = vertexSets[newBeginR];
            vertexLookup[vertexSets[newBeginR]] = index;
            vertexSets[newBeginR] = vertexInP;
            vertexLookup[vertexInP] = newBeginR;
            newBeginR++;
        }
    }

#if 0
    while (j<degree[vertex]) {
        int const neighbor = m_AdjacencyList[vertex][j];
        int const neighborLocation = vertexLookup[neighbor];

        // if in X
        if (neighborLocation >= beginX && neighborLocation < beginP) {
            // swap into new X territory
            newBeginX--;
            vertexSets[neighborLocation] = vertexSets[newBeginX];
            vertexLookup[vertexSets[newBeginX]] = neighborLocation;
            vertexSets[newBeginX] = neighbor;
            vertexLookup[neighbor] = newBeginX;
        }

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
#endif // 0

    m_lDelineators.emplace_back(VertexSets::SetDelineator(beginX, beginP, beginR));

    beginX = newBeginX;
    beginP = newBeginP;
    beginR = newBeginR;

    // maintain order of items in P.
    std::vector<int> &currentP    (stackP.back());
    std::vector<int> &currentColor(stackColor.back());
    std::vector<int> newP(currentP.size());
    std::vector<int> newColor(currentP.size());

    int index(0);
    for (size_t vertexIndex = 0; vertexIndex < currentP.size(); ++vertexIndex) {
        int const vertex(currentP[vertexIndex]);
        if (InP(vertex)) {
            newP[index] = vertex;
            newColor[index] = currentColor[vertexIndex];
            index++;
        }
    }
    newP.resize(index);
    newColor.resize(index);

    coloringStrategy.Color(m_AdjacencyMatrix, currentP /* evaluation order */, currentP /* color order */, currentColor);

    stackP.emplace_back(std::move(newP));
    stackColor.emplace_back(std::move(newColor));
}

inline void AdjacencyMatrixVertexSetsMax::MoveFromRToX(int const vertex)
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

    stackP.pop_back();
    stackColor.pop_back();
}

/*! \brief Computes the vertex v in P union X that has the most neighbors in P,
  and places P \ {neighborhood of v} in an array. These are the 
  vertices to consider adding to the partial clique during the current
  recursive call of the algorithm.

  \param pivotNonNeighbors  An intially empty vector, which will contain the set 
  P \ {neighborhood of v} when this function completes.
 */

inline std::vector<int> AdjacencyMatrixVertexSetsMax::ChoosePivotNonConst()
{
#if 0
    int pivot = -1;
    int maxIntersectionSize = -1;

    // loop through all vertices in P union X
    int i = beginX;
    while(i < beginR) {
        int vertex = vertexSets[i];
        int neighborCount = 0;

        // count the number of neighbors vertex has in P.
        int j = beginP;
        while(j < beginR) {
            if (m_AdjacencyMatrix[vertexSets[j]][vertex])
                neighborCount++;

            j++;
        }

        // if vertex has more neighbors in P, then update the pivot
        if (neighborCount > maxIntersectionSize) {
            maxIntersectionSize = neighborCount;
            pivot = vertex;
        }

        i++;
    }

    // now pivot is the vertex with the most neighbors in P
    size_t numNonNeighbors = 0;
    std::vector<int> pivotNonNeighbors;

    // make an array with the chosen pivot's non-neighbors
    if(beginR-beginP-maxIntersectionSize > 0) {
        pivotNonNeighbors.resize(beginR-beginP-maxIntersectionSize);

        int j = beginP;
        while (j < beginR) {
            if (!m_AdjacencyMatrix[pivot][vertexSets[j]]) {
                pivotNonNeighbors[numNonNeighbors] = vertexSets[j];
                numNonNeighbors++;
            }

            j++;
        }
    }

    return pivotNonNeighbors; 
#endif // 0

#if 0
    int pivot = -1;
    int maxIntersectionSize = -1;
    int i = beginX;

    // loop through all vertices in P union X
    while (i < beginR) {
        int vertex = vertexSets[i];
        int neighborCount = 0;
        int numNeighbors = 0;

        // count the number of neighbors vertex has in P.
        // only count them if the degree of the vertex
        // is greater than the the best count so far.
        int j = 0;
        if (degree[vertex] > maxIntersectionSize) {
            // count the number of neighbors vertex has in P.
            while(j < degree[vertex]) {
                int const neighbor = m_AdjacencyList[vertex][j];

                int const neighborLocation = vertexLookup[neighbor];

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
    std::vector<int> pivotNonNeighbors(beginR-beginP);
    memcpy(&pivotNonNeighbors[0], &vertexSets[beginP], (beginR-beginP)*sizeof(int));

    // we will decrement numNonNeighbors as we find neighbors
    std::size_t numNonNeighbors = beginR-beginP;

    // mark neighbors of pivot that are in P.
    int j = 0;
    while (j < degree[pivot]) {
        int const neighbor = m_AdjacencyList[pivot][j];
        int const neighborLocation = vertexLookup[neighbor];

        // if the neighbor is in P, mark it as -1
        if (neighborLocation >= beginP && neighborLocation < beginR) {
            pivotNonNeighbors[neighborLocation-beginP] = -1;
        }
 
        j++;
    }

    // put the non-neighbors at the beginning of the array
    // and update numNonNeighbors appropriately
    i = 0; 
    while (i < numNonNeighbors) {
        if (pivotNonNeighbors[i] == -1) {
            numNonNeighbors--;
            pivotNonNeighbors[i] = pivotNonNeighbors[numNonNeighbors];
        } else
            i++;
    }

    pivotNonNeighbors.resize(numNonNeighbors);

    return pivotNonNeighbors;
#endif // 0

    // snapshot current P
////    std::cout << "stackP-back-size (Before)=" << stackP.back().size() << std::endl;
    std::vector<int> currentP(stackP.back());
    stackP.emplace_back(std::move(currentP));
////    std::cout << "stackP-back-size (After )=" << stackP.back().size() << std::endl;

    std::vector<int> currentColor(stackColor.back());
    stackColor.emplace_back(std::move(currentColor));

    return stackP.back();
}

inline bool AdjacencyMatrixVertexSetsMax::InP(int const vertex) const
{
    int const vertexLocation(vertexLookup[vertex]);
    return (vertexLocation >= beginP && vertexLocation < beginR);
}

inline bool AdjacencyMatrixVertexSetsMax::PIsEmpty() const
{
    return (beginP >= beginR);
}

inline bool AdjacencyMatrixVertexSetsMax::XAndPAreEmpty() const
{
    return (beginX >= beginP && beginP >= beginR);
}

#endif //ADJACENCY_MATRIX_VERTEX_SETS_MAX_H
