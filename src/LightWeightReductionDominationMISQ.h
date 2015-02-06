#ifndef LIGHTWEIGHT_REDUCTION_DOMINATION_MISQ_H
#define LIGHTWEIGHT_REDUCTION_DOMINATION_MISQ_H

#include "LightWeightMISQ.h"
#include "IndependentSetColoringStrategy.h"
#include "ArraySet.h"
#include "Reducer.h"

#include <ctime>
#include <vector>
#include <list>

class LightWeightReductionDominationMISQ : public LightWeightMISQ 
{
public:
    LightWeightReductionDominationMISQ(std::vector<std::vector<char>> const &vAdjacencyMatrix, std::vector<std::vector<int>> const &vAdjacencyArray);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);
    virtual void GetNewOrder(std::vector<int> &vNewVertexOrder, std::vector<int> &vVertexOrder, std::vector<int> const &P, int const chosenVertex);
    virtual void ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex);
    virtual void ProcessOrderBeforeReturn(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors);

////    virtual void PrintState() const;

protected:
    std::vector<std::vector<int>>  const &m_AdjacencyArray;
////    std::vector<std::vector<int>> stackOrder;
    std::vector<std::vector<int>> stackX;
    std::vector<std::vector<int>> stackClique;
    std::vector<std::vector<int>> stackOther;
    std::vector<std::vector<int>> stackPersistentClique;
    std::vector<std::vector<int>> stackPersistentOther;
    IsolateDominationReducer reducer;
////    bool m_bInvert;
};
#endif //LIGHTWEIGHT_REDUCTION_DOMINATION_MISQ_H
