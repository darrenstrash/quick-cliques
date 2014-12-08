#include "PartialMatchGraph.h"

// system includes
#include <vector>
#include <list>
#include <iostream>
#include <algorithm>

using namespace std;

PartialMatchGraph::PartialMatchGraph(std::vector<int> const &vValues,
                      std::vector<std::vector<int>> vOrderedValues)
: m_vGraph()
, m_vStartNodeIndices()
, m_vRegularNodeIndices()
{
    for (vector<int> const &values : vOrderedValues) {
        if (values.size() == 1) {
            CreateStartNode(values[0]);
        }
        for (int i = 1; i < values.size(); ++i) {
            Connect(values[i-1], values[i], i==1 /* start node */);
        }
    }

    cout << "#StartNodes  =" << m_vStartNodeIndices.size() << endl;
    cout << "#RegularNodes=" << m_vRegularNodeIndices.size() << endl;

    auto sortNodesAndIndices = [] (NodeAndIndices const &left, NodeAndIndices const &right)
    {
        return left.node < right.node;
    };

    sort(m_vStartNodeIndices.begin(),   m_vStartNodeIndices.end(),   sortNodesAndIndices);
    sort(m_vRegularNodeIndices.begin(), m_vRegularNodeIndices.end(), sortNodesAndIndices);

    for (NodeAndIndices &start : m_vStartNodeIndices) {
        set<int> const &vTemp(m_vTempGraph[start.startIndex]);

        start = NodeAndIndices(start.node, m_vGraph.size(), m_vGraph.size() + vTemp.size());
        m_vGraph.insert(m_vGraph.end(), vTemp.begin(), vTemp.end());
    }

    for (NodeAndIndices &regular : m_vRegularNodeIndices) {
        set<int> const &vTemp(m_vTempGraph[regular.startIndex]);

        regular = NodeAndIndices(regular.node, m_vGraph.size(), m_vGraph.size() + vTemp.size());
        m_vGraph.insert(m_vGraph.end(), vTemp.begin(), vTemp.end());
    }

    Print();

    m_vTempGraph.clear();
}

PartialMatchGraph::~PartialMatchGraph()
{
}

// TODO/DS: make more efficient
NodeAndIndices PartialMatchGraph::GetStartNodeAndIndices(int const start) const
{
    for (NodeAndIndices const &startNodeIndex : m_vStartNodeIndices)
        if (startNodeIndex.node == start) return startNodeIndex;

    return NodeAndIndices(-1,-1,-1);
}

int PartialMatchGraph::GetIndexInRegularNodeAndIndices(int const regular) const
{
    for (int i = 0; i < m_vRegularNodeIndices.size(); ++i) {
        NodeAndIndices const &regularNodeIndex(m_vRegularNodeIndices[i]);
        if (regularNodeIndex.node == regular) return i;
    }

    return -1;
}

int PartialMatchGraph::GetSuccessorIndexInRegularNodeAndIndices(int const regular) const
{
    for (int i = 0; i < m_vRegularNodeIndices.size(); ++i) {
        NodeAndIndices const &regularNodeIndex(m_vRegularNodeIndices[i]);
        if (regularNodeIndex.node >= regular) return i;
    }

    return -1;
}

int PartialMatchGraph::GetIndexInStartNodeAndIndices(int const regular) const
{
    for (int i = 0; i < m_vStartNodeIndices.size(); ++i) {
        NodeAndIndices const &startNodeIndex(m_vStartNodeIndices[i]);
        if (startNodeIndex.node == regular) return i;
    }

    return -1;
}

// TODO/DS: make more efficient
NodeAndIndices PartialMatchGraph::GetRegularNodeAndIndices(int const regular) const
{
    for (NodeAndIndices const &regularNodeIndex : m_vRegularNodeIndices)
        if (regularNodeIndex.node == regular) return regularNodeIndex;

    return NodeAndIndices(-1,-1,-1);
}

// TODO/DS: make more efficient
void PartialMatchGraph::CreateStartNode(int const node)
{
    for (NodeAndIndices const &nodeAndIndex : m_vStartNodeIndices) {
        if (nodeAndIndex.node == node) return;
    }

    cout << "Creating start node: " << node << endl;

    m_vStartNodeIndices.push_back(NodeAndIndices(node, m_vTempGraph.size(), -1));
    m_vTempGraph.push_back(set<int>());
}

// TODO/DS: make more efficient
void PartialMatchGraph::CreateRegularNode(int const node)
{
    for (NodeAndIndices const &nodeAndIndex : m_vRegularNodeIndices) {
        if (nodeAndIndex.node == node) return;
    }

    cout << "Creating regular node: " << node << endl;

    m_vRegularNodeIndices.push_back(NodeAndIndices(node, m_vTempGraph.size(), -1));
    m_vTempGraph.push_back(set<int>());
}

// TODO/DS: make more efficient
void PartialMatchGraph::Connect(int const start, int const end, bool const isStartNode)
{
    cout << "Connecting " << start << (isStartNode? "(start)": "(normal)") << " and " << end << endl;
    int index(-1);
    if (isStartNode) {
        CreateStartNode(start);
        for (NodeAndIndices const &nodeAndIndex : m_vStartNodeIndices) {
            if (nodeAndIndex.node == start) index = nodeAndIndex.startIndex;
        }
    } else {
        CreateRegularNode(start);
        for (NodeAndIndices const &nodeAndIndex : m_vRegularNodeIndices) {
            if (nodeAndIndex.node == start) index = nodeAndIndex.startIndex;
        }
    }

    CreateRegularNode(end);
    m_vTempGraph[index].insert(end);
}

bool PartialMatchGraph::IsConnected(int const start, int const _end, bool const isStartNode) const
{
    NodeAndIndices const &nodeAndIndices = (isStartNode ? GetStartNodeAndIndices(start) : GetRegularNodeAndIndices(start));
    if (nodeAndIndices.node == -1) return false;

    if (!m_vTempGraph.empty()) {
        int const index(nodeAndIndices.startIndex);
        return find(m_vTempGraph[index].begin(), m_vTempGraph[index].end(), _end) != m_vTempGraph[index].end();
    }

    return find(m_vGraph.begin() + nodeAndIndices.startIndex, m_vGraph.begin() + nodeAndIndices.endIndex + 1, _end) != (m_vGraph.begin() + nodeAndIndices.endIndex + 1);
}

// TODO/DS: make more efficient
bool PartialMatchGraph::PathExists(int const start, int const end, bool const checkStart) const
{
    cout << "Checking if " << end << " is reachable from " << start << endl;
    int const indexOfStartNode(GetIndexInRegularNodeAndIndices(start));
    cout << "Index of start node: " << indexOfStartNode << endl;
    if (indexOfStartNode != -1) {
        int currentIndex(indexOfStartNode);
        std::vector<bool> reachable(m_vRegularNodeIndices.size(), 0);
        reachable[indexOfStartNode] = true;
        while (currentIndex <= reachable.size() && m_vRegularNodeIndices[currentIndex].node < end) {
            if (!reachable[currentIndex]) {
                currentIndex++; 
                continue;
            }

            NodeAndIndices const &currentNode = m_vRegularNodeIndices[currentIndex];
            int const currentStartIndex(currentNode.startIndex);
            int const currentEndIndex  (currentNode.endIndex);
            for (int i = currentStartIndex; i < currentEndIndex; ++i) {
                int const neighbor(m_vGraph[i]);
                cout << "Vertex " << neighbor << " is reachable" << endl;
                if (neighbor == end) return true;
                reachable[GetIndexInRegularNodeAndIndices(neighbor)] = true;
            }
            currentIndex++;
        }
    }

    if (checkStart) {
        int const indexOfStartNode(GetIndexInStartNodeAndIndices(start));
        if (indexOfStartNode == -1)
            return false;

        int currentIndex(indexOfStartNode);
        std::vector<bool> reachable(m_vRegularNodeIndices.size(), 0);

        // start nodes are not in the regular nodes array
        // first, evaluate all neighbors, then start incrementing 
        NodeAndIndices const &currentNode = m_vStartNodeIndices[currentIndex];
        int const currentStartIndex1(currentNode.startIndex);
        int const currentEndIndex1  (currentNode.endIndex);
        for (int i = currentStartIndex1; i < currentEndIndex1; ++i) {
            int const neighbor(m_vGraph[i]);
            if (neighbor == end) return true;
            reachable[GetIndexInRegularNodeAndIndices(neighbor)] = true;
        }

        if (currentStartIndex1 >= currentEndIndex1) return false;

        currentIndex = GetSuccessorIndexInRegularNodeAndIndices(start);
        while (currentIndex <= reachable.size() && m_vRegularNodeIndices[currentIndex].node < end) {
            if (!reachable[currentIndex]) {
                currentIndex++; 
                continue;
            }

            NodeAndIndices const &currentNode = m_vRegularNodeIndices[currentIndex];
            int const currentStartIndex(currentNode.startIndex);
            int const currentEndIndex  (currentNode.endIndex);
            for (int i = currentStartIndex; i < currentEndIndex; ++i) {
                int const neighbor(m_vGraph[i]);
                if (neighbor == end) return true;
                reachable[GetIndexInRegularNodeAndIndices(neighbor)] = true;
            }
            currentIndex++;
        }
    }

    return false;
}

bool PartialMatchGraph::Contains(std::vector<int> const &values) const
{
    if (values.empty()) return false;
    if (values.size() == 1) return GetRegularNodeAndIndices(values[0]).node != -1 || GetStartNodeAndIndices(values[0]).node != -1;
    int currentVertexIndex(1);
    int previousVertex(values[currentVertexIndex-1]);
    int thisVertex(values[currentVertexIndex]);

    while (PathExists(previousVertex, thisVertex, currentVertexIndex==1 /* check start */)) {
        if (currentVertexIndex == values.size()-1) return true;
        previousVertex = thisVertex;
        thisVertex = values[++currentVertexIndex];
    }

    return false;
}

bool PartialMatchGraph::IsEmpty() const
{
    return m_vStartNodeIndices.empty() && m_vRegularNodeIndices.empty();
}

void PartialMatchGraph::Print() const
{
    for (NodeAndIndices const &startNode : m_vStartNodeIndices) {
        int const startIndex(startNode.startIndex);
        int const endIndex(startNode.endIndex);
        cout << startNode.node << " [start] : ";
        for (int u = startIndex; u < endIndex; ++u) {
            cout << m_vGraph[u] << " ";
        }
        if (endIndex <= startIndex) {
            cout << "(empty)";
        }
        cout << endl;
    }

    for (NodeAndIndices const &normalNode : m_vRegularNodeIndices) {
        int const startIndex(normalNode.startIndex);
        int const endIndex(normalNode.endIndex);
        cout << normalNode.node << " [normal] : ";
        for (int u = startIndex; u < endIndex; ++u) {
            cout << m_vGraph[u] << " ";
        }
        if (endIndex <= startIndex) {
            cout << "(empty)";
        }
        cout << endl;
    }
}
