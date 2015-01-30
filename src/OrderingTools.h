#ifndef ORDERING_TOOLS_H
#define ORDERING_TOOLS_H

#include <vector>

namespace OrderingTools
{
    void InitialOrderingMCQ(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring);
    void InitialOrderingMISQ(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring);
    void InitialOrderingMCR(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring, size_t &cliqueSize);
    void InitialOrderingMISR(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring, size_t &cliqueSize);
};

#endif //ORDERING_TOOLS_H
