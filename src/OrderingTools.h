#ifndef ORDERING_TOOLS_H
#define ORDERING_TOOLS_H

#include <vector>
#include <cstddef>
#include <string>

#include "Isolates4.h"
#include "Isolates3.h"
#include "SparseArraySet.h"
#include "ArraySet.h"


namespace OrderingTools
{
    void InitialOrderingMCQ(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring);
    void InitialOrderingMCR(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring, size_t &cliqueSize);

    void InitialOrderingMCQ(std::vector<std::vector<int>> const &adjacencyArray, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring);

    void InitialOrderingMISQ(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring);

    void InitialOrderingMISR(std::vector<std::vector<int>> const &adjacencyArray, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring, size_t &cliqueSize);

    template <typename IsolatesType>
    void InitialOrderingMISR(std::vector<std::vector<int>> const &adjacencyArray, IsolatesType const &isolates, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring, size_t &cliqueSize);

    void InitialOrderingMISR(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring, size_t &cliqueSize);

    void InitialOrderingConnectedComponent(std::string const &filename, std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vOrderedVertices, std::vector<int> &vColoring);

    template <typename IsolatesType>
    void InitialOrderingReduction(IsolatesType &isolates, std::vector<int> &vOrderedVertices);
};

#endif //ORDERING_TOOLS_H
