#include "LightWeightReductionDominationMISR.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <iostream>

using namespace std;

LightWeightReductionDominationMISR::LightWeightReductionDominationMISR(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
: LightWeightReductionDominationMISQ(vAdjacencyMatrix, vAdjacencyArray)
{
    SetName("reduction-domination-misr");
}

void LightWeightReductionDominationMISR::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMISR(m_AdjacencyMatrix, P, vColors, m_uMaximumCliqueSize);
    vVertexOrder = P;
}
