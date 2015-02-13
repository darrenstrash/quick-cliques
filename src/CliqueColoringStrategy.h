#ifndef CLIQUE_COLORING_STRATEGY
#define CLIQUE_COLORING_STRATEGY

#include "ColoringStrategy.h"

class CliqueColoringStrategy : public ColoringStrategy
{
public:
    CliqueColoringStrategy(std::vector<std::vector<char>> const &adjacencyMatrix);
    virtual void Color(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);
    virtual void Recolor(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors, int const currentBestCliqueSize, int const currentCliqueSize);
    bool HasConflict(int const vertex, std::vector<int> const &vVerticesWithColor);
    int  GetConflictingVertex(int const vertex, std::vector<int> const &vVerticesWithColor);
    bool Repair(int const vertex, int const color, int const iBestCliqueDelta);
////    virtual void Recolor();
////    virtual void RemoveVertex(int const vertex);
////    virtual void PeekAtNextVertexAndColor(int &vertex, int &color);
protected:
    std::vector<std::vector<char>> const &m_AdjacencyMatrix;
    std::vector<std::vector<int>> m_vvVerticesWithColor;
////    std::vector<int> m_Colors;
////    std::vector<int> m_VertexOrder;
};

#endif //CLIQUE_COLORING_STRATEGY
