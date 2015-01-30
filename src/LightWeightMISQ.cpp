#include "LightWeightMISQ.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <iostream>

using namespace std;

LightWeightMISQ::LightWeightMISQ(vector<vector<char>> const &vAdjacencyMatrix)
: Algorithm("MISQ")
, m_AdjacencyMatrix(vAdjacencyMatrix)
, coloringStrategy(m_AdjacencyMatrix)
, m_uMaximumCliqueSize(0)
, stackP(vAdjacencyMatrix.size())
, stackColors(vAdjacencyMatrix.size())
, stackOrder(vAdjacencyMatrix.size())
, nodeCount(0)
////, m_bInvert(0)
{
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

long LightWeightMISQ::Run(list<std::list<int>> &cliques)
{
    R.reserve(m_AdjacencyMatrix.size());

    // don't reserve for 0-th vectors, they get std::move'd by InitialOrdering
    for (int index = 1; index < stackP.size(); ++index) {
        stackP.reserve(m_AdjacencyMatrix.size());
        stackColors.reserve(m_AdjacencyMatrix.size());
        stackOrder.reserve(m_AdjacencyMatrix.size());
    }

    // Initial coloring should be 1 to maxDegree, then the color the rest maxDegree+1.
    vector<int> &P(stackP[0]);
    vector<int> &vColors(stackColors[0]);
    vector<int> &vVertexOrder(stackOrder[0]);

    InitializeOrder(P, vVertexOrder, vColors);

    cliques.push_back(list<int>());

    RunRecursive(P, vVertexOrder, cliques, vColors);
    return cliques.size();
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
            if (!m_AdjacencyMatrix[chosenVertex][candidate]) vNewVertexOrder[uNewIndex++] = candidate;
        }
        vNewVertexOrder.resize(uNewIndex);
    }
}

void LightWeightMISQ::RunRecursive(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors)
{
    nodeCount++;
    vector<int> &vNewP(stackP[R.size() + 1]);
    vector<int> &vNewColors(stackColors[R.size() + 1]);
    while (!P.empty()) {
        int const largestColor(vColors.back());
        if (R.size() + largestColor <= m_uMaximumCliqueSize) {
            return;
        }

        vColors.pop_back();
        int const nextVertex(P.back()); P.pop_back();
        R.push_back(nextVertex);

        vector<int> &vNewVertexOrder(stackOrder[R.size()]);
        GetNewOrder(vNewVertexOrder, vVertexOrder, P, nextVertex);

        if (!vNewVertexOrder.empty()) {
            vNewP.resize(vNewVertexOrder.size());
            vNewColors.resize(vNewVertexOrder.size());
            Color(vNewVertexOrder/* evaluation order */, vNewP /* color order */, vNewColors);
            RunRecursive(vNewP, vNewVertexOrder, cliques, vNewColors);
        } else if (R.size() > m_uMaximumCliqueSize) {
            cliques.back().clear();
            cliques.back().insert(cliques.back().end(), R.begin(), R.end());
            ExecuteCallBacks(cliques.back());
            m_uMaximumCliqueSize = R.size();
        }

        PostProcessOrder(vVertexOrder, nextVertex);
        R.pop_back();
    }

    vNewColors.clear();
    vNewP.clear();
}


LightWeightMISQ::~LightWeightMISQ()
{
    cerr << "Largest Clique     : " << m_uMaximumCliqueSize << endl;
    cerr << "Node    Count      : " << nodeCount << endl;
}
