#ifndef LIGHTWEIGHT_REDUCTION_SPARSE_MISQ_H
#define LIGHTWEIGHT_REDUCTION_SPARSE_MISQ_H

#include "Algorithm.h"
#include "SparseIndependentSetColoringStrategy2.h"
#include "Isolates2.h"

#include <vector>
#include <list>

class LightWeightReductionSparseMISQ : public Algorithm
{
public:
    LightWeightReductionSparseMISQ(std::vector<std::vector<char>> const &vAdjacencyMatrix, std::vector<std::vector<int>> const &vAdjacencyArray);
    virtual ~LightWeightReductionSparseMISQ();

    virtual long Run(std::list<std::list<int>> &cliques);

    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);
    virtual void GetNewOrder(std::vector<int> &vNewVertexOrder, std::vector<int> &vVertexOrder, std::vector<int> const &P, int const chosenVertex, std::vector<int> &vCliqueVertices, std::vector<int> &vRemoved);
    virtual void PostProcessOrder(std::vector<int> &vVertexOrder, int const chosenVertex) {}

    virtual void RunRecursive(std::vector<int> &P, std::vector<int> &vVertexOrder, std::list<std::list<int>> &cliques, std::vector<int> &vColors);

////    void SetInvert(bool const invert);

protected:
    std::vector<std::vector<char>> const &m_AdjacencyMatrix;
    std::vector<std::vector<int>>  const &m_AdjacencyArray;
    size_t m_uMaximumCliqueSize;
    std::vector<int> R;
    std::vector<std::vector<int>> stackP;
    std::vector<std::vector<int>> stackColors;
    std::vector<std::vector<int>> stackOrder;
    int nodeCount;
    Isolates2 isolates;
    SparseIndependentSetColoringStrategy2 coloringStrategy;
////    bool m_bInvert;
};
#endif //LIGHTWEIGHT_REDUCTION_SPARSE_MISQ_H
