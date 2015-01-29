#ifndef LIGHTWEIGHT_MCR_H
#define LIGHTWEIGHT_MCR_H

#include "LightWeightMCQ.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class LightWeightMCR : public LightWeightMCQ
{
public:
    LightWeightMCR(std::vector<std::vector<char>> const &vAdjacencyMatrix);
////    virtual ~LightWeightMCR();

    virtual long Run(std::list<std::list<int>> &cliques);

////    virtual void RunRecursive(std::vector<int> &P, std::list<std::list<int>> &cliques, std::vector<int> &vColors);
};
#endif //LIGHTWEIGHT_MCR_H
