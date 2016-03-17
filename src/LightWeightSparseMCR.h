#ifndef LIGHTWEIGHT_SPARSE_MCR_H
#define LIGHTWEIGHT_SPARSE_MCR_H

#include "LightWeightSparseMCQ.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightSparseMCR : public LightWeightSparseMCQ
{
public:
    LightWeightSparseMCR(std::vector<std::vector<int>> const &vAdjacencyArray);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);
};
#endif //LIGHTWEIGHT_SPARSE_MCR_H
