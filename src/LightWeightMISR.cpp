#include "LightWeightMISR.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <iostream>

using namespace std;

LightWeightMISR::LightWeightMISR(vector<vector<char>> const &vAdjacencyMatrix)
: LightWeightMISQ(vAdjacencyMatrix)
{
    SetName("mcr");
}

void LightWeightMISR::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMISR(m_AdjacencyMatrix, P, vColors, m_uMaximumCliqueSize);
    vVertexOrder = P;
}
