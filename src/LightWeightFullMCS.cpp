#include "LightWeightFullMCS.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

LightWeightFullMCS::LightWeightFullMCS(vector<vector<char>> const &vAdjacencyMatrix)
: LightWeightStaticOrderMCS(vAdjacencyMatrix)
{
    SetName("mcs");
}

void LightWeightFullMCS::Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors)
{
    coloringStrategy.Recolor(m_AdjacencyMatrix, vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors, static_cast<int>(m_uMaximumCliqueSize), static_cast<int>(R.size()));
}
