#include "LightWeightReductionSparseMISR.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <iostream>

using namespace std;

LightWeightReductionSparseMISR::LightWeightReductionSparseMISR(vector<vector<int>> const &vAdjacencyArray)
: LightWeightReductionSparseMISQ(vAdjacencyArray)
{
    SetName("reduction-sparse-misr");
}

void LightWeightReductionSparseMISR::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMISR(m_AdjacencyArray, P, vColors, m_uMaximumCliqueSize);
    vVertexOrder = P;
}
