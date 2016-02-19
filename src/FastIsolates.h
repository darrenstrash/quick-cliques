#ifndef FAST_ISOLATES_H
#define FAST_ISOLATES_H

#include "Set.h"
#include "ArraySet.h"
#include "SparseArraySet.h"
#include "Reduction.h"

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <ctime>

////#define TIMERS
////#define SPARSE

template <typename NeighborSet> class FastIsolates
{
public:
    FastIsolates(std::vector<std::vector<int>> const &adjacencyArray);
    ~FastIsolates();

    void RemoveVertex(int const vertex, std::vector<Reduction> &vReductions);

    void RemoveAllIsolates(int const independentSetSIze, std::vector<int> &vIsolateVertices, std::vector<int> &vOtherRemovedVertices, std::vector<Reduction> &vReductions, bool const bConsiderAllVertices);

#ifdef SPARSE
    std::vector<SparseArraySet> const& Neighbors()  const { return neighbors;  }
#else
    std::vector<NeighborSet> const& Neighbors()  const { return neighbors;  }
#endif //SPARSE

////    void RemoveEdges(std::vector<std::pair<int,int>> const &vEdges);

    int GetAlternativeVertex(int const vertex) const;

    void SetConnectedComponent(std::vector<int> const &vVertices);

    size_t GetFoldedVertexCount() const { return foldedVertexCount; }

    void SetAllowVertexFolds(bool const allow) { m_bAllowVertexFolds = allow; }

    size_t GetReductionCount() const { return m_uReductionCount; }

    size_t GetRemainingGraphSize() const { return m_uRemainingGraphSize; }

protected: // methods
    bool RemoveIsolatedClique    (int const vertex, std::vector<int> &vIsolateVertices, std::vector<int> &vOtherRemovedVertices, std::vector<Reduction> &vReductions);
    bool FoldVertex(int const vertex, std::vector<int> &vIsolateVertices,  std::vector<int> &vOtherRemovedVertices, std::vector<Reduction> &vReductions);

protected: // members
    std::vector<std::vector<int>> const &m_AdjacencyArray;
#ifdef SPARSE
    std::vector<SparseArraySet>     neighbors;
#else
    std::vector<NeighborSet>     neighbors;
#endif // SPARSE
    std::vector<bool> m_vbInGraph;
    ArraySet remaining;
    std::vector<bool> vMarkedVertices;
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
    size_t foldedVertexCount;
    bool m_bAllowVertexFolds;
    size_t m_uReductionCount;
    size_t m_uRemainingGraphSize;
};

#endif //FAST_ISOLATES_H
