#ifndef LIGHTWEIGHT_REDUCTION_SPARSE_MISQ_H
#define LIGHTWEIGHT_REDUCTION_SPARSE_MISQ_H

#include "MaxSubgraphAlgorithm.h"
#include "SparseIndependentSetColoringStrategy.h"
#include "Isolates2.h"

#include <vector>
#include <list>
#include <ctime>

class LightWeightReductionSparseMISQ : public MaxSubgraphAlgorithm
{
public:
    LightWeightReductionSparseMISQ(std::vector<std::vector<int>> const &vAdjacencyArray);

    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);
    virtual void GetNewOrder(std::vector<int> &vNewVertexOrder, std::vector<int> &vVertexOrder, std::vector<int> const &P, int const chosenVertex);
    virtual void ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex);
    virtual void ProcessOrderBeforeReturn(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors);

protected:
    std::vector<std::vector<int>>  const &m_AdjacencyArray;
    SparseIndependentSetColoringStrategy coloringStrategy;
    std::vector<std::vector<int>> stackClique;
    std::vector<std::vector<int>> stackOther;
    std::vector<std::vector<int>> stackPersistentClique;
    std::vector<std::vector<int>> stackPersistentOther;
    Isolates2 isolates;
};
#endif //LIGHTWEIGHT_REDUCTION_SPARSE_MISQ_H
