#ifndef LIGHTWEIGHT_FULL_MCS_H
#define LIGHTWEIGHT_FULL_MCS_H

#include "LightWeightStaticOrderMCS.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightFullMCS : public LightWeightStaticOrderMCS
{
public:
    LightWeightFullMCS(std::vector<std::vector<char>> const &vAdjacencyMatrix);

    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);
};

#endif //LIGHTWEIGHT_FULL_MCS_H
