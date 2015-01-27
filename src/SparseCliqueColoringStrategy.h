#ifndef SPARSE_CLIQUE_COLORING_STRATEGY
#define SPARSE_CLIQUE_COLORING_STRATEGY

#include "ColoringStrategy.h"

#include <vector>
#include <list>

class SparseCliqueColoringStrategy : public ColoringStrategy
{
public:
    SparseCliqueColoringStrategy(std::vector<std::vector<int>> const &adjacencyList);
    virtual void Color(std::vector<std::vector<int>> const &adjacencyList, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);
////    virtual void Recolor();
////    virtual void RemoveVertex(int const vertex);
////    virtual void PeekAtNextVertexAndColor(int &vertex, int &color);
protected:
    std::vector<std::list<int>> m_vlColorToVertices;
    std::vector<int> m_vVertexToColor;
    std::vector<bool> m_vNeighborHasColor;
};

#endif //SPARSE_CLIQUE_COLORING_STRATEGY
