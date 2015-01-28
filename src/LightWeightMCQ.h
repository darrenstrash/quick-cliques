#ifndef LIGHTWEIGHT_MCQ_H

#include "Algorithm.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightMCQ : public Algorithm
{
public:
    LightWeightMCQ(std::vector<std::vector<char>> const &vAdjacencyMatrix);
    virtual ~LightWeightMCQ();

    virtual long Run(std::list<std::list<int>> &cliques);

    void RunRecursive(std::vector<int> &P, std::list<std::list<int>> &cliques, std::vector<int> &vColors);

private:
    std::vector<std::vector<char>> const &m_AdjacencyMatrix;
    CliqueColoringStrategy coloringStrategy;
    size_t m_uMaximumCliqueSize;
    std::vector<int> R;
    std::vector<std::vector<int>> stackP;
    std::vector<std::vector<int>> stackColors;
};
#endif
