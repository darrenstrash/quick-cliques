#ifndef ORDERING_TOOLS_H
#define ORDERING_TOOLS_H

#include <vector>

namespace OrderingTools
{
    std::vector<int> InitialOrderingMCQ(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> const &degree);

    std::vector<int> InitialOrderingMCR(std::vector<std::vector<char>> const &adjacencyMatrix);
};

#endif //ORDERING_TOOLS_H
