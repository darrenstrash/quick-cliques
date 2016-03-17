#ifndef SPARSE_CLIQUE_COLORING_STRATEGY
#define SPARSE_CLIQUE_COLORING_STRATEGY

#include "ColoringStrategy.h"

#include <vector>
#include <list>

class SparseCliqueColoringStrategy : public ColoringStrategy
{
public:
    SparseCliqueColoringStrategy(std::vector<std::vector<int>> const &adjacencyList);
    virtual void Color(std::vector<std::vector<int>> const &adjacencyList, std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);
    virtual void Color(std::vector<std::vector<int>> const &adjacencyList, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);

    virtual void Recolor(std::vector<std::vector<int>> const &adjacencyArray, std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors, int const currentBestCliqueSize, int const currentCliqueSize);
    bool HasConflict(int const vertex, std::vector<int> const &vVerticesWithColor);
    int  GetConflictingVertex(int const vertex, std::vector<int> const &vVerticesWithColor);
    bool Repair(int const vertex, int const color, int const iBestCliqueDelta);

protected:
    std::vector<std::vector<int>> const &m_AdjacencyArray;
    std::vector<std::list<int>> m_vlColorToVertices;
    std::vector<int> m_vVertexToColor;
    std::vector<bool> m_vNeighborHasColor;
    // From SparseIndependentSetColoringStrategy
    std::vector<std::vector<int>> m_vvVerticesWithColor;
    std::vector<int>              m_vNeighborColorCount;
    std::vector<bool>             m_vbNeighbors;
    std::vector<bool>             m_vbConflictNeighbors;
};

#endif //SPARSE_CLIQUE_COLORING_STRATEGY
