#ifndef CONNECTED_COMPONENT_MISS_2_H
#define CONNECTED_COMPONENT_MISS_2_H

#include "LightWeightReductionMISQ.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class ConnectedComponentMISS2 : public LightWeightReductionMISQ
{
public:
    ConnectedComponentMISS2(std::vector<std::vector<char>> const &vAdjacencyMatrix, std::vector<std::vector<int>> const &vAdjacencyArray);
    ~ConnectedComponentMISS2();

    void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors);

    virtual void GetNewOrder(std::vector<int> &vNewVertexOrder, std::vector<int> &vVertexOrder, std::vector<int> const &P, int const chosenVertex);

#if 1
    virtual void RunRecursive(std::vector<int> &P, std::vector<int> &vVertexOrder, std::list<std::list<int>> &cliques, std::vector<int> &vColors);
    void RunRecursive(std::vector<int> &P, std::vector<int> &vVertexOrder, std::list<std::list<int>> &cliques, std::vector<int> &vColors, bool const checkForConnectedComponents);
#endif // 0

    void SetIsolates(Isolates3<ArraySet> const &newIsolates)
    {
        isolates = newIsolates;
    }

////    virtual void PostProcessOrder(std::vector<int> &vVertexOrder, int const chosenVertex);
protected:
    size_t m_uSubgraphDepth;
    std::vector<int> m_vSubgraphClique;
    std::vector<int> m_vStackDelta;
};

#endif // CONNECTED_COMPONENT_MISS_2_H
