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

    bool ComputeResidualPath(std::vector<int> const &vMatching, std::vector<int> &vPath);

    void ComputeMaximumMatching(std::vector<int> &vMatching);

private:
    void PushOnStack(int const vertex);
    int  PopOffStack();
    bool IsEvaluated(int const vertex) const;

    std::vector<std::vector<int>> m_AdjacencyList;
    std::vector<int> m_Stack;
    std::vector<int> m_Evaluated;
    int m_iCurrentRound;
};

#endif // _BI_DOUBLE_GRAPH_H_
