#ifndef ORDERING_TOOLS_H
#define ORDERING_TOOLS_H

#include <vector>

namespace OrderingTools
{
    void InitialOrderingMCQ(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring);
    void InitialOrderingMCR(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring, size_t &cliqueSize);

    void InitialOrderingMCQ(std::vector<std::vector<int>> const &adjacencyArray, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring);

    void InitialOrderingMISQ(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring);

    void InitialOrderingMISR(std::vector<std::vector<int>> const &adjacencyArray, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring, size_t &cliqueSize);
    void InitialOrderingMISR(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring, size_t &cliqueSize);

    void InitialOrderingConnectedComponent(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring);
};

#endif //ORDERING_TOOLS_H
