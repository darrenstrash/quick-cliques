#include "LightWeightMCR.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <iostream>

using namespace std;

LightWeightMCR::LightWeightMCR(vector<vector<char>> const &vAdjacencyMatrix)
: LightWeightMCQ(vAdjacencyMatrix)
{
    SetName("mcr");
}

void LightWeightMCR::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMCR(m_AdjacencyMatrix, P, vColors, m_uMaximumCliqueSize);
    vVertexOrder = P;
}
