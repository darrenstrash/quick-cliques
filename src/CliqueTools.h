#ifndef CLIQUE_TOOLS_H
#define CLIQUE_TOOLS_H

#include <vector>
#include <list>
#include <set>

namespace CliqueTools
{
    bool IsMaximalClique(std::vector<std::vector<int>> &adjacencyArray, std::list<int> const&clique, bool const verbose);
    bool IsClique(std::vector<std::vector<char>> &adjacencyMatrix, std::list<int> const &clique, bool const verbose);
};

#endif //CLIQUE_TOOLS_H
