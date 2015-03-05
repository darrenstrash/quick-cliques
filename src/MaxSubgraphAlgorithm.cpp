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
, timeToLargestClique(0)
, m_bQuiet(false)
, stackEvaluatedHalfVertices()
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

    bool bMaximumCliqueSizeIsSet(m_uMaximumCliqueSize != 0);

    InitializeOrder(P, vVertexOrder, vColors);

    cliques.push_back(list<int>());

    if (!bMaximumCliqueSizeIsSet && R.size() < m_uMaximumCliqueSize) {
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
        timeToLargestClique = clock() - startTime;
    }

    depth++;
    RunRecursive(P, vVertexOrder, cliques, vColors);
    return cliques.size();
}

void MaxSubgraphAlgorithm::RunRecursive(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors)
{
    nodeCount++;
    vector<int> &vNewP(stackP[R.size()+1]);
    vector<int> &vNewColors(stackColors[R.size()+1]);
    vector<int> &vNewVertexOrder(stackOrder[R.size()+1]);

////    bool &bEvaluatedHalfVertices = stackEvaluatedHalfVertices[R.size() + 1];
////    stackEvaluatedHalfVertices[depth + 1] = (rand()%(depth+1) == depth);
////    stackEvaluatedHalfVertices[depth + 1] = true;
////    stackEvaluatedHalfVertices[depth + 1] = false;
////    stackEvaluatedHalfVertices[depth + 1] = (depth<=2);

////    size_t halfVertices(0);
////    size_t index = P.size()+1;
////    for (; index > 0; --index) {
////        if (vColors[index-1] + R.size() <= m_uMaximumCliqueSize) {
////            halfVertices = (index - 1) + 0.9*(P.size() - (index - 1));
////            break;
////        }
////    }

////    stackEvaluatedHalfVertices[depth + 1] = ((P.size() - index) > 50); //(depth<=2);
    stackEvaluatedHalfVertices[depth + 1] = true;

    size_t const uOriginalPSize(P.size());

    if (nodeCount%10000 == 0) {
        if (!m_bQuiet) {
            cout << "Evaluated " << nodeCount << " nodes. " << Tools::GetTimeInSeconds(clock() - startTime) << endl;
            PrintState();
        }
    }

    while (!P.empty()) {
////    if (!stackEvaluatedHalfVertices[depth + 1]) {
////        stackEvaluatedHalfVertices[depth + 1] = (rand()%(depth+1) == depth);
////        stackEvaluatedHalfVertices[depth + 1] = (rand()%2 == 1);

////        stackEvaluatedHalfVertices[depth+1] = (rand() % 20 == 1); ////!stackEvaluatedHalfVertices[depth+1];

////        if (!stackEvaluatedHalfVertices[depth + 1]) {
////            if (uOriginalPSize >= 100 && P.size() > 1 && vColors[P.size()-5] + R.size() <= m_uMaximumCliqueSize) {
////                stackEvaluatedHalfVertices[depth + 1] = true;
////            }
////            size_t index = P.size();
////            for (; index > 0; --index) {
////                if (vColors[index-1] + R.size() <= m_uMaximumCliqueSize) {
////                    halfVertices = (index - 1 + uOriginalPSize)/2;
////                    break;
////                }
////            }
////
////            if (P.size() - index == 1) {
////                stackEvaluatedHalfVertices[depth + 1] = true;
////            }

////            if (P.size() <= halfVertices) {
////                stackEvaluatedHalfVertices[depth + 1] = true;
////            }
////////            else if (rand() % 2 == 1) { ////(P.size() - halfVertices) == 1) {
////////                stackEvaluatedHalfVertices[R.size() + 1] = true;
////////            }
////            else if (rand() % P.size() <= index) {
////                stackEvaluatedHalfVertices[depth + 1] = true;
////            }
////            else if (rand() % 2 == 1) {
////                stackEvaluatedHalfVertices[depth+1] = stackEvaluatedHalfVertices[depth];
////            }
////        }

////        cout << depth << ": P: ";
////        for (int const p : P) {
////            cout << p << " ";
////        }
////        cout << endl;

        if (depth == 0) {
            if (!m_bQuiet) {
                cout << "Only " << P.size() << " more vertices to go! " << Tools::GetTimeInSeconds(clock() - startTime) << endl;
            }
        }

        int const largestColor(vColors.back());
        if (R.size() + largestColor <= m_uMaximumCliqueSize) {
            ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
            P.clear();
            return;
        }

        vColors.pop_back();
        int const nextVertex(P.back()); P.pop_back();

////        cout << depth << ": Choosing next vertex: " << nextVertex << endl;

        GetNewOrder(vNewVertexOrder, vVertexOrder, P, nextVertex);

        if (!vNewVertexOrder.empty()) {
            vNewP.resize(vNewVertexOrder.size());
            vNewColors.resize(vNewVertexOrder.size());
            Color(vNewVertexOrder/* evaluation order */, vNewP /* color order */, vNewColors);
#ifdef PREPRUNE
            if (R.size() + vNewColors.back() > m_uMaximumCliqueSize) {
                depth++;
                RunRecursive(vNewP, vNewVertexOrder, cliques, vNewColors);
                depth--;
            }
#else
            depth++;
            RunRecursive(vNewP, vNewVertexOrder, cliques, vNewColors);
            depth--;
#endif // PREPRUNE
        } else if (R.size() > m_uMaximumCliqueSize) {
            cliques.back().clear();
            cliques.back().insert(cliques.back().end(), R.begin(), R.end());
            ExecuteCallBacks(cliques.back());
            m_uMaximumCliqueSize = R.size();
            timeToLargestClique = clock() - startTime;
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
                timeToLargestClique = clock() - startTime;
            }
        }
    }

    ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
    P.clear();

    vNewColors.clear();
    vNewP.clear();
}


MaxSubgraphAlgorithm::~MaxSubgraphAlgorithm()
{
    if (!m_bQuiet) {
        cerr << "Largest Clique     : " << m_uMaximumCliqueSize << endl;
        cerr << "Time to Clique     : " << Tools::GetTimeInSeconds(timeToLargestClique) << endl;
        cerr << "Search Nodes       : " << nodeCount << endl;
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
