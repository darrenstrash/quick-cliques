#ifndef _MATCHING_TOOLS_H_
#define _MATCHING_TOOLS_H_

#include <set>
#include <vector>

namespace MatchingTools
{

    std::set<int> ComputeLeftMIS(std::vector<std::vector<int>> const &biDoubleGraph);
    std::set<int> ComputeCriticalSet(std::vector<std::vector<int>> const &adjacencyList);
    std::set<int> ComputeBiDoubleMIS(std::vector<std::vector<int>> const &biDoubleGraph, std::vector<bool> const &vInGraph, std::set<int> const &setInGraph);
    std::set<int> ComputeLeftMIS(std::vector<std::vector<int>> const &biDoubleGraph, std::vector<bool> const &vInGraph, std::set<int> const &setInGraph);
};

#endif // _MATCHING_TOOLS_H_
