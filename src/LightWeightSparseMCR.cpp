#include "LightWeightSparseMCR.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <iostream>

using namespace std;

LightWeightSparseMCR::LightWeightSparseMCR(vector<vector<int>> const &vAdjacencyArray)
: LightWeightSparseMCQ(vAdjacencyArray)
{
    SetName("sparse-mcr");
}

void LightWeightSparseMCR::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMCR(m_AdjacencyArray, P, vColors, m_uMaximumCliqueSize);
    vVertexOrder = P;
}
