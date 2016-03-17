#ifndef LIGHTWEIGHT_SPARSE_STATIC_ORDER_MCS_H
#define LIGHTWEIGHT_SPARSE_STATIC_ORDER_MCS_H

#include "LightWeightSparseMCQ.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightSparseStaticOrderMCS : public LightWeightSparseMCQ
{
public:
    LightWeightSparseStaticOrderMCS(std::vector<std::vector<int>> const &vAdjacencyArray);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);

    virtual void GetNewOrder(std::vector<int> &vNewVertexOrder, std::vector<int> &vVertexOrder, std::vector<int> const &P, int const chosenVertex);

    virtual void ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex);

protected:
    std::vector<std::vector<int>> stackOrder;
};

#endif //LIGHTWEIGHT_SPARSE_STATIC_ORDER_MCS_H
