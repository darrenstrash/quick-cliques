#ifndef LIGHTWEIGHT_STATIC_ORDER_MCS_H
#define LIGHTWEIGHT_STATIC_ORDER_MCS_H

#include "LightWeightMCQ.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightStaticOrderMCS : public LightWeightMCQ
{
public:
    LightWeightStaticOrderMCS(std::vector<std::vector<char>> const &vAdjacencyMatrix);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);

    virtual void GetNewOrder(std::vector<int> &vNewVertexOrder, std::vector<int> &vVertexOrder, std::vector<int> const &P, int const chosenVertex);

    virtual void ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex);

protected:
    std::vector<std::vector<int>> stackOrder;
};

#endif //LIGHTWEIGHT_STATIC_ORDER_MCS_H
