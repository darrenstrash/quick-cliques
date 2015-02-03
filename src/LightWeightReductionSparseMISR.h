#ifndef LIGHTWEIGHT_REDUCTION_SPARSE_MISR_H
#define LIGHTWEIGHT_REDUCTION_SPARSE_MISR_H

#include "LightWeightReductionSparseMISQ.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightReductionSparseMISR : public LightWeightReductionSparseMISQ
{
public:
    LightWeightReductionSparseMISR(std::vector<std::vector<int>> const &adjacencyArray);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);
};
#endif //LIGHTWEIGHT_REDUCTION_SPARSE_MISR_H
