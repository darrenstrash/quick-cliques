#ifndef LIGHTWEIGHT_MISQ_H
#define LIGHTWEIGHT_MISQ_H

#include "MaxSubgraphAlgorithm.h"
#include "IndependentSetColoringStrategy.h"

#include <vector>
#include <list>
#include <ctime>

class LightWeightMISQ : public MaxSubgraphAlgorithm
{
public:
    LightWeightMISQ(std::vector<std::vector<char>> const &vAdjacencyMatrix);
////    virtual ~LightWeightMISQ();

////    virtual long Run(std::list<std::list<int>> &cliques);

    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);
    virtual void GetNewOrder(std::vector<int> &vNewVertexOrder, std::vector<int> &vVertexOrder, std::vector<int> const &P, int const chosenVertex);
    virtual void ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex);
    virtual void ProcessOrderBeforeReturn(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors) {}

////    void SetInvert(bool const invert);

protected:
    std::vector<std::vector<char>> const &m_AdjacencyMatrix;
    IndependentSetColoringStrategy coloringStrategy;
////    bool m_bInvert;
};
#endif
