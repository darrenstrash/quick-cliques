#include "LightWeightMCQ.h"
#include "GraphTools.h"

#include <iostream>

using namespace std;

LightWeightMCQ::LightWeightMCQ(vector<vector<char>> const &vAdjacencyMatrix)
: Algorithm("MCQ")
, m_AdjacencyMatrix(vAdjacencyMatrix)
, coloringStrategy(m_AdjacencyMatrix)
, m_uMaximumCliqueSize(0)
{
}

long LightWeightMCQ::Run(list<std::list<int>> &cliques)
{
    vector<int> P(std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyMatrix, false /* descending */)));
    R.reserve(m_AdjacencyMatrix.size());

    vector<int> vColors(P.size(), 0);
    coloringStrategy.Color(m_AdjacencyMatrix, P, vColors);

    cliques.push_back(list<int>());

    RunRecursive(P, cliques, vColors);
    return cliques.size();
}

void LightWeightMCQ::RunRecursive(vector<int> &P, list<list<int>> &cliques, vector<int> &vColors)
{
    coloringStrategy.Color(m_AdjacencyMatrix, P, vColors);
    while (!P.empty()) {
        int const largestColor(vColors.back());
        if (R.size() + largestColor < m_uMaximumCliqueSize) return;

        vColors.pop_back();
        int const nextVertex(P.back()); P.pop_back();
        R.push_back(nextVertex);

        vector<int> vNewP(P);
        vector<int> vNewColors(vColors);
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
            RunRecursive(vNewP, cliques, vNewColors);
        }

        R.pop_back();
    }
}

LightWeightMCQ::~LightWeightMCQ()
{
    cerr << "Largest Clique     : " << m_uMaximumCliqueSize << endl;
}
