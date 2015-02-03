#ifndef LIGHTWEIGHT_REDUCTION_SPARSE_STATIC_ORDER_MISS_H
#define LIGHTWEIGHT_REDUCTION_SPARSE_STATIC_ORDER_MISS_H

#include "LightWeightReductionSparseMISQ.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightReductionSparseStaticOrderMISS : public LightWeightReductionSparseMISQ
{
public:
    LightWeightReductionSparseStaticOrderMISS(std::vector<std::vector<int>> const &vAdjacencyArray);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);

    virtual void GetNewOrder(std::vector<int> &vNewVertexOrder, std::vector<int> &vVertexOrder, std::vector<int> const &P, int const chosenVertex);

////    virtual void ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex);

protected:
    std::vector<std::vector<int>> stackOrder;
};

#endif //LIGHTWEIGHT_REDUCTION_SPARSE_STATIC_ORDER_MISS_H
