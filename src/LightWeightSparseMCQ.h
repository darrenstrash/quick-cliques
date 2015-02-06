#ifndef LIGHTWEIGHT_SPARSE_MCQ_H
#define LIGHTWEIGHT_SPARSE_MCQ_H

#include "MaxSubgraphAlgorithm.h"
#include "SparseCliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightSparseMCQ : public MaxSubgraphAlgorithm
{
public:
    LightWeightSparseMCQ(std::vector<std::vector<int>> const &vAdjacencyArray);

    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);
    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);
    virtual void GetNewOrder(std::vector<int> &vNewVertexOrder, std::vector<int> &vVertexOrder, std::vector<int> const &P, int const chosenVertex);
    virtual void ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex);
    virtual void ProcessOrderBeforeReturn(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors) {}

////    void SetInvert(bool const invert);

protected:
    std::vector<std::vector<int>> const &m_AdjacencyArray;
    SparseCliqueColoringStrategy coloringStrategy;
    std::vector<bool> vMarkedVertices;
////    bool m_bInvert;
};
#endif //LIGHTWEIGHT_SPARSE_MCQ_H
