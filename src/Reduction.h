#ifndef REDUCTION_H
#define REDUCTION_H

#include <vector>
#include <utility>

enum ReductionType {ISOLATED_VERTEX, FOLDED_VERTEX, DOMINATED_VERTEX, REMOVED_VERTEX, REMOVED_VERTEX_AND_NEIGHBORS};

class Reduction
{

public:
    Reduction(ReductionType const type)
    : m_iVertex(-1)
    , m_vNeighbors()
    , m_vRemovedEdges()
    , reductionType(type) {
    }

    void SetVertex(int const vertex)
    {
        m_iVertex = vertex;
    }

    int GetVertex() const { return m_iVertex; }

    void AddNeighbor(int const neighbor) {
        m_vNeighbors.push_back(neighbor);
    }

    std::vector<int> const &GetNeighbors() const
    {
        return m_vNeighbors;
    }

    void AddRemovedEdge(int const v1, int const v2)
    {
        m_vRemovedEdges.push_back(std::make_pair(v1, v2));
    }

    std::vector<std::pair<int,int>> const &GetRemovedEdges() const
    {
        return m_vRemovedEdges;
    }

    ReductionType GetType() const { return reductionType; }

private:
    int m_iVertex;
    std::vector<int>                m_vNeighbors;
    std::vector<std::pair<int,int>> m_vRemovedEdges;
    ReductionType                   reductionType;
};

#endif // REDUCTION_H
