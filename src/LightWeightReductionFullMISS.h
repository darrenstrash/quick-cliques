#ifndef LIGHTWEIGHT_REDUCTION_FULL_MISS_H
#define LIGHTWEIGHT_REDUCTION_FULL_MISS_H

#include "LightWeightReductionStaticOrderMISS.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightReductionFullMISS : public LightWeightReductionStaticOrderMISS
{
public:
    LightWeightReductionFullMISS(std::vector<std::vector<char>> const &vAdjacencyMatrix, std::vector<std::vector<int>> const &vAdjacencyArray);

    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);
};

#endif //LIGHTWEIGHT_REDUCTION_FULL_MISS_H
