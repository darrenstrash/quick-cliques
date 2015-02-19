#ifndef ISOLATES_4_H
#define ISOLATES_4_H

#include "Set.h"
#include "ArraySet.h"
#include "SparseArraySet.h"

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <ctime>

////#define TIMERS
////#define SPARSE

template <typename NeighborSet> class Isolates4
{
public:
    Isolates4(std::vector<std::vector<int>> const &adjacencyArray);
    ~Isolates4();

    void RemoveVertexAndNeighbors(int const vertex, std::vector<int> &vRemoved);
    void RemoveVertex(int const vertex);

    void RemoveAllIsolates(int const independentSetSIze, std::vector<int> &vIsolateVertices, std::vector<int> &vOtherRemovedVertices, std::vector<std::pair<int,int>> &vAddedEdges, bool const bConsiderAllVertices);
    void ReplaceAllRemoved(std::vector<int> const &vRemoved);

    int NextVertexToRemove(std::vector<int> &vVertices);
    int NextVertexToRemove();

    size_t size() const { return isolates.Size(); }

    ArraySet const& GetIsolates() const { return isolates; }
    ArraySet const& GetInGraph()  const { return inGraph;  }
#ifdef SPARSE
    std::vector<SparseArraySet> const& Neighbors()  const { return neighbors;  }
#else
    std::vector<NeighborSet> const& Neighbors()  const { return neighbors;  }
#endif //SPARSE

////    void RemoveEdges(std::vector<std::pair<int,int>> const &vEdges);

    int GetAlternativeVertex(int const vertex) const;

    void SetConnectedComponent(std::vector<int> const &vVertices);

protected: // methods
    bool RemoveIsolatedClique    (int const vertex, std::vector<int> &vIsolateVertices, std::vector<int> &vOtherRemovedVertices);
////    bool RemoveIsolatedPath      (int const vertex, std::vector<int> &vIsolateVertices,  std::vector<int> &vOtherRemovedVertices, std::vector<std::pair<int,int>> &vAddedEdges);

protected: // members
    std::vector<std::vector<int>> const &m_AdjacencyArray;
#ifdef SPARSE
    std::vector<SparseArraySet>     neighbors;
#else
    std::vector<NeighborSet>     neighbors;
#endif // SPARSE
    ArraySet inGraph;
    ArraySet isolates;
    ArraySet remaining;
    std::vector<bool> vMarkedVertices;
    std::map<int,int> m_AlternativeVertices;
#ifdef TIMERS
    clock_t timer;
    clock_t removeTimer;
    clock_t replaceTimer;
    clock_t sortDuringNextTimer;
    clock_t removeOneDuringNextTimer;
    clock_t removeDuringNextTimer;
    clock_t replaceDuringNextTimer;
    #endif // TIMERS
    bool m_bConnectedComponentMode;
};

#endif //ISOLATES_4_H
