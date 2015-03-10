#ifndef ISOLATES_IS_COLORING_STRATEGY
#define ISOLATES_IS_COLORING_STRATEGY

#include "ColoringStrategy.h"

#include "Isolates3.h"
#include "ArraySet.h"
#include <ctime>

template <typename IsolatesType>
class IsolatesIndependentSetColoringStrategy : public ColoringStrategy
{
public:
    IsolatesIndependentSetColoringStrategy(IsolatesType const &isolates, size_t const numVerticesInOriginalGraph);
    ~IsolatesIndependentSetColoringStrategy();
    virtual void Color(IsolatesType const &isolates, std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);
    virtual void Recolor(IsolatesType const &isolates, std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors, int const currentBestCliqueSize, int const currentCliqueSize);
    bool HasConflict(int const vertex, std::vector<int> const &vVerticesWithColor);
    int  GetConflictingVertex(int const vertex, std::vector<int> const &vVerticesWithColor);
    bool Repair(int const vertex, int const color, int const iBestCliqueDelta);
////    virtual void Recolor();
////    virtual void RemoveVertex(int const vertex);
////    virtual void PeekAtNextVertexAndColor(int &vertex, int &color);
protected:
    std::vector<std::vector<int>> m_vvVerticesWithColor;
    std::vector<int>              m_vVertexToColor;
    std::vector<int>              m_vNeighborColorCount;
    std::vector<bool>             m_vbNeighbors;
    std::vector<bool>             m_vbConflictNeighbors;
    IsolatesType                  m_Isolates;
    clock_t                       m_ColorTimer;
////    std::vector<int> m_Colors;
////    std::vector<int> m_VertexOrder;
};

#endif //ISOLATES_IS_COLORING_STRATEGY
