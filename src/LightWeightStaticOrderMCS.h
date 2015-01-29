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
////    virtual ~LightWeightStaticOrderMCS();

    virtual long Run(std::list<std::list<int>> &cliques);

    virtual void RunRecursive(std::vector<int> &P, std::vector<int> &vVertexOrder, std::list<std::list<int>> &cliques, std::vector<int> &vColors);
protected:
    std::vector<std::vector<int>> stackOrder;
};

#endif //LIGHTWEIGHT_STATIC_ORDER_MCS_H
