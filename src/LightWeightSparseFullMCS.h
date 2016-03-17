#ifndef LIGHTWEIGHT_SPARSE_FULL_MCS_H
#define LIGHTWEIGHT_SPARSE_FULL_MCS_H

#include "LightWeightSparseStaticOrderMCS.h"
#include "SparseCliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightSparseFullMCS : public LightWeightSparseStaticOrderMCS
{
public:
    LightWeightSparseFullMCS(std::vector<std::vector<int>> const &vAdjacencyArray);

    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);
};

#endif //LIGHTWEIGHT_SPARSE_FULL_MCS_H
