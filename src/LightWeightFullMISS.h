#ifndef LIGHTWEIGHT_FULL_MISS_H
#define LIGHTWEIGHT_FULL_MISS_H

#include "LightWeightStaticOrderMISS.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightFullMISS : public LightWeightStaticOrderMISS
{
public:
    LightWeightFullMISS(std::vector<std::vector<char>> const &vAdjacencyMatrix);

    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);
};

#endif //LIGHTWEIGHT_FULL_MISS_H
