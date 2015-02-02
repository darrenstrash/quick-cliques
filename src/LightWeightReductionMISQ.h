#ifndef LIGHTWEIGHT_REDUCTION_MISQ_H
#define LIGHTWEIGHT_REDUCTION_MISQ_H

#include "LightWeightMISQ.h"
#include "IndependentSetColoringStrategy.h"
#include "Isolates2.h"

#include <ctime>

#include <vector>
#include <list>

class LightWeightReductionMISQ : public LightWeightMISQ 
{
public:
    LightWeightReductionMISQ(std::vector<std::vector<char>> const &vAdjacencyMatrix, std::vector<std::vector<int>> const &vAdjacencyArray);
////    virtual ~LightWeightReductionMISQ();

////    virtual long Run(std::list<std::list<int>> &cliques);

////    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);
    virtual void GetNewOrder(std::vector<int> &vNewVertexOrder, std::vector<int> &vVertexOrder, std::vector<int> const &P, int const chosenVertex);

    virtual void ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex);
    virtual void ProcessOrderBeforeReturn(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors);

////    virtual void RunRecursive(std::vector<int> &P, std::vector<int> &vVertexOrder, std::list<std::list<int>> &cliques, std::vector<int> &vColors);

////    void SetInvert(bool const invert);

protected:
////    std::vector<std::vector<char>> const &m_AdjacencyMatrix;
    std::vector<std::vector<int>>  const &m_AdjacencyArray;
////    IndependentSetColoringStrategy coloringStrategy;
////    size_t m_uMaximumCliqueSize;
////    std::vector<int> R;
////    std::vector<std::vector<int>> stackP;
////    std::vector<std::vector<int>> stackColors;
////    std::vector<std::vector<int>> stackOrder;
    std::vector<std::vector<int>> stackClique;
    std::vector<std::vector<int>> stackOther;
    std::vector<std::vector<int>> stackPersistentClique;
    std::vector<std::vector<int>> stackPersistentOther;
////    int nodeCount;
////    int depth;
    Isolates2 isolates;
////    clock_t startTime;
////    bool m_bInvert;
};
#endif //LIGHTWEIGHT_REDUCTION_MISQ_H
