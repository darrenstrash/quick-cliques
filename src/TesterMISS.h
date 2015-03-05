#ifndef TESTER_MISS_H
#define TESTER_MISS_H

#include "TesterStaticOrderMISS.h"
#include "CliqueColoringStrategy.h"
#include "IsolatesIndependentSetColoringStrategy.h"

#include <vector>
#include <list>

class TesterMISS : public TesterStaticOrderMISS
{
public:
    TesterMISS(std::vector<std::vector<char>> const &vAdjacencyMatrix, std::vector<std::vector<int>> const &vAdjacencyArray);

    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);

protected:
////    SparseIndependentSetColoringStrategy sparseColoringStrategy;
    IsolatesIndependentSetColoringStrategy isolatesColoringStrategy;
};

#endif //TESTER_MISS_H
