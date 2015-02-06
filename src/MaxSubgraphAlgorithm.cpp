#include "MaxSubgraphAlgorithm.h"
#include "Tools.h"

#include <iostream>

using namespace std;

MaxSubgraphAlgorithm::MaxSubgraphAlgorithm(string const &name)
: Algorithm(name)
, m_uMaximumCliqueSize(0)
, stackP()
, stackColors()
, stackOrder()
, nodeCount(0)
, depth(-1)
, startTime(clock())
, m_bQuiet(false)
////, m_bInvert(0)
{
}

////void MaxSubgraphAlgorithm::SetInvert(bool const invert)
////{
////    m_bInvert = invert;
////}

long MaxSubgraphAlgorithm::Run(list<std::list<int>> &cliques)
{
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
        m_uMaximumCliqueSize = R.size();
    }

    depth++;
    RunRecursive(P, vVertexOrder, cliques, vColors);
    return cliques.size();
}

////long MaxSubgraphAlgorithm::Run(vector<int> const &startingVertices, list<std::list<int>> &cliques)
////{
////    vector<int> &P(stackP[0]);
////    vector<int> &vColors(stackColors[0]);
////    vector<int> &vVertexOrder(stackOrder[0]);
////
////    InitializeOrder(P, vVertexOrder, vColors);
////
////    cliques.push_back(list<int>());
////
////    if (R.size() < m_uMaximumCliqueSize) {
////        cliques.back().clear();
////        cliques.back().insert(cliques.back().end(), P.begin(), P.begin() + m_uMaximumCliqueSize);
////        ExecuteCallBacks(cliques.back());
////    }
////
////    ProcessOrderAfterRecursion(vVertexOrder, P, vColors, -1 /* no vertex chosen for removal */);
////
////    if (R.size() > m_uMaximumCliqueSize) {
////        cliques.back().clear();
////        cliques.back().insert(cliques.back().end(), R.begin(), R.end());
////        ExecuteCallBacks(cliques.back());
////        m_uMaximumCliqueSize = R.size();
////    }
////
////    depth++;
////    RunRecursive(P, vVertexOrder, cliques, vColors);
////    return cliques.size();
////}

void MaxSubgraphAlgorithm::RunRecursive(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors)
{
    nodeCount++;
    vector<int> &vNewP(stackP[R.size()+1]);
    vector<int> &vNewColors(stackColors[R.size()+1]);
    vector<int> &vNewVertexOrder(stackOrder[R.size()+1]);

////    cout << depth << ": P: ";
////    for (int const p : P) {
////        cout << p << " ";
////    }
////    cout << endl;

    if (nodeCount%10000 == 0) {
        if (!m_bQuiet) {
            cout << "Evaluated " << nodeCount << " nodes. " << GetTimeInSeconds(clock() - startTime) << endl;
        }
////        PrintState();
    }

    while (!P.empty()) {

        if (depth == 0) {
            if (!m_bQuiet) {
                cout << "Only " << P.size() << " more vertices to go! " << GetTimeInSeconds(clock() - startTime) << endl;
            }
        }

        int const largestColor(vColors.back());
        if (R.size() + largestColor <= m_uMaximumCliqueSize) {
            ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
            return;
        }

        vColors.pop_back();
        int const nextVertex(P.back()); P.pop_back();

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

////        if (R.size() > m_uMaximumCliqueSize && bPIsEmpty && P.empty()) {
////            cout << "ERROR!" << endl << flush;
////        }

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


MaxSubgraphAlgorithm::~MaxSubgraphAlgorithm()
{
    if (!m_bQuiet) {
        cerr << "Largest Clique     : " << m_uMaximumCliqueSize << endl;
        cerr << "Node    Count      : " << nodeCount << endl;
    }
}

void MaxSubgraphAlgorithm::PrintState() const
{
    cout << "(";
    for (size_t index = 0; index <= R.size(); ++index) {
        cout << stackP[index].size();
        if (index != R.size()) cout << ", ";
    }
    cout << ")" << endl << flush;
}
