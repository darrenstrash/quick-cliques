#ifndef GRAPH_TOOLS_H
#define GRAPH_TOOLS_H

#include <vector>
#include <set>
#include <map>

namespace GraphTools
{
    void ComputeInducedSubgraph(std::vector<std::vector<int>> &adjacencyList, std::set<int> const &vertices, std::vector<std::vector<int>> &subraph, std::map<int,int> &remapping);
    std::vector<int> OrderVerticesByDegree(std::vector<std::vector<int>> const &adjacencyList, bool const ascending);
////    void RemoveVertices(vector<vector<int>> &adjacencyList, vector<int> const &vVertices);
};

#endif //GRAPH_TOOLS_H
