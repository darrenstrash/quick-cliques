#include "LightWeightMISQ.h"
#include "OrderingTools.h"
#include "GraphTools.h"
#include "Tools.h"

#include <cmath>
#include <iostream>

using namespace std;

LightWeightMISQ::LightWeightMISQ(vector<vector<char>> const &vAdjacencyMatrix)
: MaxSubgraphAlgorithm("misq")
, m_AdjacencyMatrix(vAdjacencyMatrix)
, coloringStrategy(m_AdjacencyMatrix)
////, m_uMaximumCliqueSize(0)
////, stackP(vAdjacencyMatrix.size() + 1)
////, stackColors(vAdjacencyMatrix.size() + 1)
////, stackOrder(vAdjacencyMatrix.size() + 1)
////, nodeCount(0)
////, depth(-1)
////, startTime(clock())
////, m_bInvert(0)
{
    R.reserve(m_AdjacencyMatrix.size());

    stackP.resize(m_AdjacencyMatrix.size() + 1);
    stackColors.resize(m_AdjacencyMatrix.size() + 1);
    stackOrder.resize(m_AdjacencyMatrix.size() + 1);

    // don't reserve for 0-th vectors, they get std::move'd by InitialOrdering
    for (int index = 0; index < stackP.size(); ++index) {
        stackP[index].reserve(m_AdjacencyMatrix.size());
        stackColors[index].reserve(m_AdjacencyMatrix.size());
        stackOrder[index].reserve(m_AdjacencyMatrix.size());
    }
}

////void LightWeightMISQ::SetInvert(bool const invert)
////{
////    m_bInvert = invert;
////}

void LightWeightMISQ::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMISQ(m_AdjacencyMatrix, P, vColors);
    vVertexOrder = P;
}

void LightWeightMISQ::Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors)
{
    coloringStrategy.Color(m_AdjacencyMatrix, vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors);
}

void LightWeightMISQ::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
    vNewVertexOrder.resize(P.size());
    {
        size_t uNewIndex(0);
        for (int const candidate : P) {
////            if (!m_bInvert && m_AdjacencyMatrix[chosenVertex][candidate]) vNewVertexOrder[uNewIndex++] = candidate;
////            cout << depth << ": accessing " << chosenVertex << "," << candidate << " in " << m_AdjacencyMatrix.size() << "," << m_AdjacencyMatrix[chosenVertex].size() << endl << flush;
////            cout << depth << ": put in index " << uNewIndex << "/" << vNewVertexOrder.size() << endl << flush;
////            cout << depth << ": size of stackP : " << stackP.size() << endl;
////            cout << depth << ": size of stackOrder: " << stackOrder.size() << endl;
            if (chosenVertex == candidate) continue;
            if (!m_AdjacencyMatrix[chosenVertex][candidate]) vNewVertexOrder[uNewIndex++] = candidate;
        }
        vNewVertexOrder.resize(uNewIndex);
    }

    R.push_back(chosenVertex);
}

void LightWeightMISQ::ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex)
{
////    cout << "LightWeightMISQ::ProcessOrderAfterRecursion" << endl;
////    Color(vVertexOrder, P, vColors);
    if (chosenVertex != -1) R.pop_back();
}

