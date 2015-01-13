#ifndef ISOLATES2_H
#define ISOLATES2_H

#include "Set.h"
#include "ArraySet.h"

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <ctime>

class Isolates2
{
public:
    Isolates2(std::vector<std::vector<int>> &adjacencyArray);
    ~Isolates2();

    void RemoveVertexAndNeighbors(int const vertex, std::vector<int> &vRemoved);
    void RemoveVertex(int const vertex);

    void RemoveAllIsolates(int const independentSetSIze, std::vector<int> &vIsolateVertices, std::vector<int> &vOtherRemovedVertices, std::vector<std::pair<int,int>> &vAddedEdges);
    void ReplaceAllRemoved(std::vector<int> const &vRemoved);

    int NextVertexToRemove(std::vector<int> &vVertices);
    int NextVertexToRemove();

    size_t size() const { return isolates.Size(); }

    ArraySet const& GetIsolates() const { return isolates; }
    ArraySet const& GetInGraph()  const { return inGraph;  }
    std::vector<std::vector<int>> const& Neighbors()  const { return neighbors;  }

////    void RemoveEdges(std::vector<std::pair<int,int>> const &vEdges);

    int GetAlternativeVertex(int const vertex) const;

protected: // methods
    bool RemoveIsolatedClique    (int const vertex, std::vector<int> &vIsolateVertices, std::vector<int> &vOtherRemovedVertices);
////    bool RemoveIsolatedPath      (int const vertex, std::vector<int> &vIsolateVertices,  std::vector<int> &vOtherRemovedVertices, std::vector<std::pair<int,int>> &vAddedEdges);
    void SwapOutNeighbor(int const vertex, int const neighborToMove);
    void SwapOutNonGraphNeighbors(int const vertex);
    void SwapInGraphNeighbors(int const vertex);

protected: // members
    std::vector<std::vector<int>> &m_AdjacencyArray;
    std::vector<int>               degree;
    std::vector<std::vector<int>>  neighbors;
    ArraySet inGraph;
    ArraySet isolates;
    ArraySet remaining;
    std::vector<bool> vMarkedVertices;
    std::map<int,int> m_AlternativeVertices;
    clock_t timer;
    clock_t removeTimer;
    clock_t replaceTimer;
};

#endif //ISOLATES2_H
