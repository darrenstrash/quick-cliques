#ifndef TESTER_MISQ_H
#define TESTER_MISQ_H

#include "LightWeightMISQ.h"
#include "IndependentSetColoringStrategy.h"
#include "ArraySet.h"
#include "Isolates2.h"
#include "Isolates3.h"
#include "Isolates4.h"
#include "IsolatesWithMatrix.h"

#include <ctime>
#include <vector>
#include <list>

class TesterMISQ : public LightWeightMISQ 
{
public:
    TesterMISQ(std::vector<std::vector<char>> const &vAdjacencyMatrix, std::vector<std::vector<int>> const &vAdjacencyArray);

////    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);

    size_t ComputeConnectedComponents(std::vector<std::vector<int>> &vComponents);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);
    virtual void GetNewOrder(std::vector<int> &vNewVertexOrder, std::vector<int> &vVertexOrder, std::vector<int> const &P, int const chosenVertex);
    virtual void ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex);
    virtual void ProcessOrderBeforeReturn(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors);

////    void SetIsolates(Isolates3<ArraySet> const &newIsolates)
////    {
////        isolates = newIsolates;
////    }
////
////    virtual void PrintState() const;

protected:
    std::vector<std::vector<int>>  const &m_AdjacencyArray;
////    std::vector<std::vector<int>> stackOrder;
////    std::vector<std::vector<int>> stackX;
    std::vector<std::vector<int>> stackClique;
    std::vector<std::vector<int>> stackOther;
    std::vector<std::vector<Reduction>> stackReductions;
    std::vector<std::vector<int>> stackPersistentClique;
    std::vector<std::vector<int>> stackPersistentOther;
    std::vector<std::vector<Reduction>> stackPersistentReductions;
////    Isolates4<SparseArraySet> isolates;
    Isolates4<ArraySet> isolates;
////    IsolatesWithMatrix<ArraySet> isolates;
    std::vector<bool> vRemoveIsolates;
////    std::vector<int>  vFoldedVertexCount;
////    bool m_bInvert;
};
#endif //TESTER_MISQ_H
