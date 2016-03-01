#ifndef _BI_DOUBLE_GRAPH_H_
#define _BI_DOUBLE_GRAPH_H_

#include <vector>
#include <set>

class BiDoubleGraph
{
public:
    BiDoubleGraph(std::vector<std::vector<int>> adjacencyList);
    ~BiDoubleGraph();

    std::vector<int> const &Neighbors(int const vertex) const;

    bool InLeftSide(int const vertex) const;

    bool ComputeResidualPath(std::vector<int> const &vMatching, std::vector<int> &vPath) const;

    void ComputeMaximumMatching(std::vector<int> &vMatching) const;

private:

    std::vector<std::vector<int>> m_AdjacencyList;
};

#endif // _BI_DOUBLE_GRAPH_H_
