#include "ForwardSearchMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

ForwardSearchMISS::ForwardSearchMISS(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
: ForwardSearchStaticOrderMISS(vAdjacencyMatrix, vAdjacencyArray)
////, sparseColoringStrategy(vAdjacencyArray)
, isolatesColoringStrategy(isolates, vAdjacencyArray.size())
{
    SetName("tester-miss");
}

void ForwardSearchMISS::Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors)
{
////    sparseColoringStrategy.Recolor(m_AdjacencyArray, vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors, static_cast<int>(m_uMaximumCliqueSize), static_cast<int>(R.size()));
    isolatesColoringStrategy.Recolor(isolates, vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors, static_cast<int>(m_uMaximumCliqueSize), static_cast<int>(R.size() + isolates.GetFoldedVertexCount()));
////    coloringStrategy.Recolor(m_AdjacencyMatrix, vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors, static_cast<int>(m_uMaximumCliqueSize), static_cast<int>(R.size()));
}
