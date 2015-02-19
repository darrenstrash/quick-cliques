#include "ComparisonFullMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

ComparisonFullMISS::ComparisonFullMISS(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
: ComparisonStaticOrderMISS(vAdjacencyMatrix, vAdjacencyArray)
{
    SetName("reduction-miss");
}

void ComparisonFullMISS::Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors)
{
    coloringStrategy.Recolor(m_AdjacencyMatrix, vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors, static_cast<int>(m_uMaximumCliqueSize), static_cast<int>(R.size()));
}
