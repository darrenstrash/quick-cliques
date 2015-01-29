#include "LightWeightMCQ.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <iostream>

using namespace std;

LightWeightMCQ::LightWeightMCQ(vector<vector<char>> const &vAdjacencyMatrix)
: Algorithm("MCQ")
, m_AdjacencyMatrix(vAdjacencyMatrix)
, coloringStrategy(m_AdjacencyMatrix)
, m_uMaximumCliqueSize(0)
, stackP(vAdjacencyMatrix.size())
, stackColors(vAdjacencyMatrix.size())
{
}

long LightWeightMCQ::Run(list<std::list<int>> &cliques)
{
    R.reserve(m_AdjacencyMatrix.size());

    // don't reserve for 0-th vectors, they get std::move'd by InitialOrdering
    for (int index = 1; index < stackP.size(); ++index) {
        stackP.reserve(m_AdjacencyMatrix.size());
        stackColors.reserve(m_AdjacencyMatrix.size());
    }

    // Initial coloring should be 1 to maxDegree, then the color the rest maxDegree+1.
    vector<int> &P(stackP[0]);
    vector<int> &vColors(stackColors[0]);

    OrderingTools::InitialOrderingMCQ(m_AdjacencyMatrix, P, vColors);

    cliques.push_back(list<int>());

    RunRecursive(P, cliques, vColors);
    return cliques.size();
}

void LightWeightMCQ::RunRecursive(vector<int> &P, list<list<int>> &cliques, vector<int> &vColors)
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

        size_t      uNewIndex(0);
        for (size_t index = 0; index < P.size(); index++) {
            int const vertex(P[index]);
            if (m_AdjacencyMatrix[nextVertex][vertex]) {
                vNewP[uNewIndex] = vertex;
                vNewColors[uNewIndex] = vColors[index];
                uNewIndex++;
            }
        }
        vNewP.resize(uNewIndex);
        vNewColors.resize(uNewIndex);

        if (vNewP.empty() && R.size() > m_uMaximumCliqueSize) {
            cliques.back().clear();
            cliques.back().insert(cliques.back().end(), R.begin(), R.end());
            ExecuteCallBacks(cliques.back());
            m_uMaximumCliqueSize = R.size();
        }

        if (!vNewP.empty()) {
            coloringStrategy.Color(m_AdjacencyMatrix, vNewP /* evaluation order */, vNewP /* color order */, vNewColors);
            RunRecursive(vNewP, cliques, vNewColors);
        }

        R.pop_back();
    }

    vNewColors.clear();
    vNewP.clear();
}

LightWeightMCQ::~LightWeightMCQ()
{
    cerr << "Largest Clique     : " << m_uMaximumCliqueSize << endl;
}
