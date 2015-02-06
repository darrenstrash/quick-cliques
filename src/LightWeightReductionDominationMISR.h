#ifndef LIGHTWEIGHT_REDUCTION_DOMINATION_MISR_H
#define LIGHTWEIGHT_REDUCTION_DOMINATION_MISR_H

#include "LightWeightReductionDominationMISQ.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightReductionDominationMISR : public LightWeightReductionDominationMISQ
{
public:
    LightWeightReductionDominationMISR(std::vector<std::vector<char>> const &vAdjacencyMatrix, std::vector<std::vector<int>> const &adjacencyArray);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);
};
#endif //LIGHTWEIGHT_REDUCTION_DOMINATION_MISR_H
