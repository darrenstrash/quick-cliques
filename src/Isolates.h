#ifndef ISOLATES_H
#define ISOLATES_H

#include <vector>
#include <map>
#include <set>
#include <utility>

class Isolates
{
public:
    Isolates(std::vector<std::vector<int>> &adjacencyArray);
    ~Isolates();

    void RemoveVertexAndNeighbors(int const vertex, std::vector<int> &vRemoved);
    void RemoveVertex(int const vertex);

    void RemoveAllIsolates(std::vector<int> &vIsolateVertices, std::vector<int> &vOtherRemovedVertices, std::vector<std::pair<int,int>> &vAddedEdges);
    void ReplaceAllRemoved(std::vector<int> const &vRemoved);

    int NextVertexToRemove();

    size_t size() const { return isolates.size(); }

    std::set<int> const& GetIsolates() const { return isolates; }
    std::set<int> const& GetInGraph()  const { return inGraph;  }
    std::vector<std::set<int>> const& Neighbors()  const { return neighbors;  }

    void RemoveEdges(std::vector<std::pair<int,int>> const &vEdges);

    int GetAlternativeVertex(int const vertex) const;

protected: // methods
    bool RemoveIsolatedClique    (int const vertex, std::vector<int> &vIsolateVertices, std::vector<int> &vOtherRemovedVertices);
    bool RemoveIsolatedPath      (int const vertex, std::vector<int> &vIsolateVertices,  std::vector<int> &vOtherRemovedVertices, std::vector<std::pair<int,int>> &vAddedEdges);

protected: // members
    std::vector<std::vector<int>> &m_AdjacencyArray;
    std::vector<std::set<int>>     neighbors;
    std::set<int> inGraph;
    std::set<int> isolates;
    std::set<int> removed;
    std::set<int> remaining;
    std::vector<bool> vMarkedVertices;
    std::vector<std::vector<int>> vvRemovedVertices;
    std::map<int,int> m_AlternativeVertices;
};

#endif //ISOLATES_H
