#include "LightWeightSparseFullMCS.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

LightWeightSparseFullMCS::LightWeightSparseFullMCS(vector<vector<int>> const &vAdjacencyArray)
: LightWeightSparseStaticOrderMCS(vAdjacencyArray)
{
    SetName("sparse-mcs");
}

void LightWeightSparseFullMCS::Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors)
{
    coloringStrategy.Recolor(m_AdjacencyArray, vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors, static_cast<int>(m_uMaximumCliqueSize), static_cast<int>(R.size()));
}
