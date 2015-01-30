#ifndef LIGHTWEIGHT_MCR_H
#define LIGHTWEIGHT_MCR_H

#include "LightWeightMCQ.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightMCR : public LightWeightMCQ
{
public:
    LightWeightMCR(std::vector<std::vector<char>> const &vAdjacencyMatrix);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);
};
#endif //LIGHTWEIGHT_MCR_H
