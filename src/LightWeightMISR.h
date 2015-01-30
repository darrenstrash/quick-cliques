#ifndef LIGHTWEIGHT_MISR_H
#define LIGHTWEIGHT_MISR_H

#include "LightWeightMISQ.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightMISR : public LightWeightMISQ
{
public:
    LightWeightMISR(std::vector<std::vector<char>> const &vAdjacencyMatrix);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);
};
#endif //LIGHTWEIGHT_MISR_H
