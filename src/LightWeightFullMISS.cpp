#include "LightWeightFullMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

LightWeightFullMISS::LightWeightFullMISS(vector<vector<char>> const &vAdjacencyMatrix)
: LightWeightStaticOrderMISS(vAdjacencyMatrix)
{
    SetName("mcs");
}

void LightWeightFullMISS::Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors)
{
    coloringStrategy.Recolor(m_AdjacencyMatrix, vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors, static_cast<int>(m_uMaximumCliqueSize), static_cast<int>(R.size()));
}
