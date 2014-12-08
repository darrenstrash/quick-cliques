#ifndef PARTIAL_MATCH_GRAPH_H
#define PARTIAL_MATCH_GRAPH_H

#include <vector>
#include <set>

class NodeAndIndices {
public:
    NodeAndIndices(int const _node, int const _start, int const _end) : node(_node), startIndex(_start), endIndex(_end) { }

    int node;
    int startIndex;
    int endIndex;
};

class PartialMatchGraph
{
public:
    PartialMatchGraph(std::vector<int> const &values,
                      std::vector<std::vector<int>> orderedValues);

    ~PartialMatchGraph();

    bool IsInPartialMatchGraph(std::vector<int> const &values) const;
    bool IsEmpty() const;
    void Print() const;

    bool Contains(std::vector<int> const &values) const;

private: // methods
    void CreateStartNode(int const start);
    void CreateRegularNode(int const regular);

    NodeAndIndices GetStartNodeAndIndices(int const start) const;
    NodeAndIndices GetRegularNodeAndIndices(int const regular) const;

    void Connect(int const start, int const end, bool const isStartNode);

    bool IsConnected(int const start, int const end, bool const isStartNode) const;

    bool PathExists(int const start, int const end, bool const checkStart) const;

    int GetIndexInRegularNodeAndIndices(int const regular) const;
    int GetSuccessorIndexInRegularNodeAndIndices(int const regular) const;

    int GetIndexInStartNodeAndIndices(int const regular) const;

private: // members
    std::vector<std::set<int>> m_vTempGraph;
    std::vector<int> m_vGraph; // stored contiguously for speed.
    std::vector<NodeAndIndices> m_vStartNodeIndices;   // position in vGraph
    std::vector<NodeAndIndices> m_vRegularNodeIndices; // position in vGraph
};

#endif //PARTIAL_MATCH_GRAPH_H
