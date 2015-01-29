#include "LightWeightStaticOrderMCS.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

LightWeightStaticOrderMCS::LightWeightStaticOrderMCS(vector<vector<char>> const &vAdjacencyMatrix)
: LightWeightMCQ(vAdjacencyMatrix)
, stackOrder(vAdjacencyMatrix.size())
{
    SetName("static-order-mcs");
}

long LightWeightStaticOrderMCS::Run(list<std::list<int>> &cliques)
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

    OrderingTools::InitialOrderingMCR(m_AdjacencyMatrix, P, vColors, m_uMaximumCliqueSize);

    cliques.push_back(list<int>());

    // if we found a clique in the ordering phase, it's sitting in the first slots of the ordering array.
    if (m_uMaximumCliqueSize > 0) {
        cliques.back().clear();
        cliques.back().insert(cliques.back().end(), P.begin(), P.begin() + m_uMaximumCliqueSize);
        ExecuteCallBacks(cliques.back());
    }

    vVertexOrder = std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyMatrix, false /* non-increasing */));

    RunRecursive(P, vVertexOrder, cliques, vColors);
    return cliques.size();
}

void LightWeightStaticOrderMCS::RunRecursive(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors)
{
    vector<int> &vNewP(stackP[R.size() + 1]);
    vector<int> &vNewColors(stackColors[R.size() + 1]);
    while (!P.empty()) {
        vNewP.resize(P.size());
        vNewColors.resize(vColors.size());
        int const largestColor(vColors.back());
        if (R.size() + largestColor <= m_uMaximumCliqueSize) {
            vNewColors.clear();
            vNewP.clear();
            return;
        }

        vColors.pop_back();
        int const nextVertex(P.back()); P.pop_back();
        R.push_back(nextVertex);

        vector<int> &vNewVertexOrder(stackOrder[R.size()]);
        vNewVertexOrder.resize(vVertexOrder.size());
        {
            size_t uNewIndex(0);
            for (int const candidate : vVertexOrder) {
                if (m_AdjacencyMatrix[nextVertex][candidate]) vNewVertexOrder[uNewIndex++] = candidate;
            }
            vNewVertexOrder.resize(uNewIndex);
            vNewP.resize(uNewIndex);
            vNewColors.resize(uNewIndex);
        }

        if (vNewP.empty() && R.size() > m_uMaximumCliqueSize) {
            cliques.back().clear();
            cliques.back().insert(cliques.back().end(), R.begin(), R.end());
            ExecuteCallBacks(cliques.back());
            m_uMaximumCliqueSize = R.size();
        }

        if (!vNewP.empty()) {
            coloringStrategy.Color(m_AdjacencyMatrix, vNewVertexOrder/* evaluation order */, vNewP /* color order */, vNewColors);
            RunRecursive(vNewP, vNewVertexOrder, cliques, vNewColors);
            vNewVertexOrder.clear();
        }

        vVertexOrder.erase(find(vVertexOrder.begin(), vVertexOrder.end(), nextVertex));
        // try searching from end, might be faster in general...
////        size_t indexAfterVertex(0);
////        for (indexAfterVertex = vVertexOrder.size(); indexAfterVertex >= 1; indexAfterVertex--) {
////            if (vVertexOrder[indexAfterVertex-1] == nextVertex) {
////                break;
////            }
////        }
////
////        for (; indexAfterVertex < vVertexOrder.size(); indexAfterVertex++) {
////            vVertexOrder[indexAfterVertex-1] = vVertexOrder[indexAfterVertex];
////        }
////        vVertexOrder.pop_back();
        R.pop_back();
    }

    vNewColors.clear();
    vNewP.clear();
}

////LightWeightStaticOrderMCS::~LightWeightStaticOrderMCS()
////{
////    cerr << "Largest Clique     : " << m_uMaximumCliqueSize << endl;
////}
