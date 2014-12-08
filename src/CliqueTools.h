#ifndef CLIQUE_TOOLS_H
#define CLIQUE_TOOLS_H

#include <vector>

namespace CliqueTools
{
    void ComputeCliqueGraph(std::vector<std::vector<int>> &adjacencyList, std::vector<std::vector<int>> &cliqueGraphAdjacencyList, std::vector<std::vector<int>> &vertexToClique);

    void DecomposeIntoDisjointCliques(std::vector<std::vector<int>> &adjacencyList, std::vector<std::vector<int>> &cliquesToReturn);

    void FindMaximalIndependentSetInCliqueGraph(std::vector<std::vector<int>> &adjacencyList);
};

#endif //CLIQUE_TOOLS_H
