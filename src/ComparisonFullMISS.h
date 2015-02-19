#ifndef COMPARISON_FULL_MISS_H
#define COMPARISON_FULL_MISS_H

#include "ComparisonStaticOrderMISS.h"
#include "CliqueColoringStrategy.h"

#include <vector>
#include <list>

class ComparisonFullMISS : public ComparisonStaticOrderMISS
{
public:
    ComparisonFullMISS(std::vector<std::vector<char>> const &vAdjacencyMatrix, std::vector<std::vector<int>> const &vAdjacencyArray);

    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);
};

#endif //COMPARISON_FULL_MISS_H
