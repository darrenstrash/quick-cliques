#include "LightWeightMISQ.h"
#include "OrderingTools.h"
#include "GraphTools.h"
#include "Tools.h"

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
, depth(-1)
, startTime(clock())
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

    if (R.size() < m_uMaximumCliqueSize) {
        cliques.back().clear();
        cliques.back().insert(cliques.back().end(), P.begin(), P.begin() + m_uMaximumCliqueSize);
        ExecuteCallBacks(cliques.back());
    }

    ProcessOrderAfterRecursion(vVertexOrder, P, vColors, -1 /* no vertex chosen for removal */);

    if (R.size() > m_uMaximumCliqueSize) {
        cliques.back().clear();
        cliques.back().insert(cliques.back().end(), R.begin(), R.end());
        ExecuteCallBacks(cliques.back());
    }

    depth++;
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

    R.push_back(chosenVertex);
}

void LightWeightMISQ::ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex)
{
    if (chosenVertex != -1) R.pop_back();
}

void LightWeightMISQ::RunRecursive(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors)
{
    nodeCount++;
    vector<int> &vNewP(stackP[R.size() + 1]);
    vector<int> &vNewColors(stackColors[R.size() + 1]);

    if (nodeCount%10000 == 0) {
        cout << "Evaluated " << nodeCount << " nodes. " << GetTimeInSeconds(clock() - startTime) << endl;
    }

    while (!P.empty()) {

        if (depth == 0) {
            cout << "Only " << P.size() << " more vertices to go! " << GetTimeInSeconds(clock() - startTime) << endl;
        }
        int const largestColor(vColors.back());
        if (R.size() + largestColor <= m_uMaximumCliqueSize) {
            ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
            return;
        }

        vColors.pop_back();
        int const nextVertex(P.back()); P.pop_back();

        vector<int> &vNewVertexOrder(stackOrder[R.size()]);
        GetNewOrder(vNewVertexOrder, vVertexOrder, P, nextVertex);

        if (!vNewVertexOrder.empty()) {
            vNewP.resize(vNewVertexOrder.size());
            vNewColors.resize(vNewVertexOrder.size());
            Color(vNewVertexOrder/* evaluation order */, vNewP /* color order */, vNewColors);
            depth++;
            RunRecursive(vNewP, vNewVertexOrder, cliques, vNewColors);
            depth--;
        } else if (R.size() > m_uMaximumCliqueSize) {
            cliques.back().clear();
            cliques.back().insert(cliques.back().end(), R.begin(), R.end());
            ExecuteCallBacks(cliques.back());
            m_uMaximumCliqueSize = R.size();
        }

        bool bPIsEmpty(P.empty());
        ProcessOrderAfterRecursion(vVertexOrder, P, vColors, nextVertex);

        if (!bPIsEmpty && P.empty()) {
            if (R.size() > m_uMaximumCliqueSize) {
                cliques.back().clear();
                cliques.back().insert(cliques.back().end(), R.begin(), R.end());
                ExecuteCallBacks(cliques.back());
                m_uMaximumCliqueSize = R.size();
            }
        }
    }

    ProcessOrderBeforeReturn(vVertexOrder, P, vColors);

    vNewColors.clear();
    vNewP.clear();
}


LightWeightMISQ::~LightWeightMISQ()
{
    cerr << "Largest Clique     : " << m_uMaximumCliqueSize << endl;
    cerr << "Node    Count      : " << nodeCount << endl;
}
