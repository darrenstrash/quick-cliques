#ifndef FORWARD_SEARCH_MISS_H
#define FORWARD_SEARCH_MISS_H

#include "ForwardSearchStaticOrderMISS.h"
#include "CliqueColoringStrategy.h"
#include "IsolatesIndependentSetColoringStrategy.h"
#include "Isolates3.h"
#include "Isolates4.h"
#include "IsolatesWithMatrix.h"
#include "ArraySet.h"

#include <vector>
#include <list>

class ForwardSearchMISS : public ForwardSearchStaticOrderMISS
{
public:
    ForwardSearchMISS(std::vector<std::vector<char>> const &vAdjacencyMatrix, std::vector<std::vector<int>> const &vAdjacencyArray);

    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors);

protected:
////    SparseIndependentSetColoringStrategy sparseColoringStrategy;
////    IsolatesIndependentSetColoringStrategy<Isolates4<SparseArraySet>> isolatesColoringStrategy;
    IsolatesIndependentSetColoringStrategy<Isolates4<ArraySet>> isolatesColoringStrategy;
};

#endif //FORWARD_SEARCH_MISS_H
