#include "LightWeightSparseMCQ.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <iostream>

using namespace std;

LightWeightSparseMCQ::LightWeightSparseMCQ(vector<vector<int>> const &vAdjacencyArray)
: MaxSubgraphAlgorithm("sparse-mcq")
, m_AdjacencyArray(vAdjacencyArray)
, coloringStrategy(vAdjacencyArray)
, vMarkedVertices(vAdjacencyArray.size(), false)
////, m_bInvert(0)
{
    R.reserve(m_AdjacencyArray.size());

    stackP.resize(m_AdjacencyArray.size() + 1);
    stackColors.resize(m_AdjacencyArray.size() + 1);
    stackOrder.resize(m_AdjacencyArray.size() + 1);

    // don't reserve for 0-th vectors, they get std::move'd by InitialOrdering
    for (int index = 0; index < stackP.size(); ++index) {
        stackP[index].reserve(m_AdjacencyArray.size());
        stackColors[index].reserve(m_AdjacencyArray.size());
        stackOrder[index].reserve(m_AdjacencyArray.size());
    }
}

////void LightWeightSparseMCQ::SetInvert(bool const invert)
////{
////    m_bInvert = invert;
////}

void LightWeightSparseMCQ::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMCQ(m_AdjacencyArray, P, vColors);
    vVertexOrder = P;
}

void LightWeightSparseMCQ::Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors)
{
    coloringStrategy.Color(m_AdjacencyArray, vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors);

////    cout << "Colors:";
////    for (int const color : vColors) {
////        cout << color << " ";
////    }
////    cout << endl;
}

void LightWeightSparseMCQ::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
////    cout << "OldP    :";
////    for (int const vertex : P) {
////        cout << vertex << " ";
////    }
////    cout << endl;

    vNewVertexOrder.resize(P.size());

    for (int const neighbor : m_AdjacencyArray[chosenVertex]) {
        vMarkedVertices[neighbor] = true;
    }

    size_t uNewIndex(0);
    for (int const candidate : P) {
        ////            if (!m_bInvert && m_AdjacencyMatrix[chosenVertex][candidate]) vNewVertexOrder[uNewIndex++] = candidate;
        if (vMarkedVertices[candidate]) vNewVertexOrder[uNewIndex++] = candidate;
    }
    vNewVertexOrder.resize(uNewIndex);

    for (int const neighbor : m_AdjacencyArray[chosenVertex]) {
        vMarkedVertices[neighbor] = false;
    }

////    cout << "NewOrder:";
////    for (int const vertex : vNewVertexOrder) {
////        cout << vertex << " ";
////    }
////    cout << endl;

    R.push_back(chosenVertex);
}

void LightWeightSparseMCQ::ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex)
{
    if (chosenVertex != -1) R.pop_back();
}
