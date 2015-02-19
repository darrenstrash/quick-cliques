#include "TesterStaticOrderMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"
#include "Tools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

TesterStaticOrderMISS::TesterStaticOrderMISS(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
: TesterMISQ(vAdjacencyMatrix, vAdjacencyArray)
, onlyConsider(vAdjacencyMatrix.size())
, vMarkedVertices(vAdjacencyMatrix.size(), false)
{
    SetName("tester-static-order-miss");
}

void TesterStaticOrderMISS::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMISR(m_AdjacencyMatrix, P, vColors, m_uMaximumCliqueSize);
    vVertexOrder = P; ////std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyMatrix, true /* non-decreasing */)); //// = P; //?
}

void TesterStaticOrderMISS::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
    vector<int> &vCliqueVertices(stackClique[depth+1]); vCliqueVertices.clear(); vCliqueVertices.push_back(chosenVertex);
    vector<int> &vRemoved(stackOther[depth+1]); vRemoved.clear();
    vector<pair<int,int>> vAddedEdgesUnused;
    isolates.RemoveVertexAndNeighbors(chosenVertex, vRemoved);

#ifdef NO_ISOLATES_P_LEFT_10
    vector<int> const &vColors(stackColors[depth]);
    bool bRemoveIsolates((vColors.size() < 10) || (vColors[vColors.size()-10] + R.size() <= m_uMaximumCliqueSize));
#else
    bool bRemoveIsolates(true);
#endif //NO_ISOLATES_P_LEFT_10
////    size_t index = vColors.size();
////    for (; index > 0; --index) {
////        if (R.size() + vColors[index] <= m_uMaximumCliqueSize) { bRemoveIsolates = !(P.size() - index < 10); break; }
////    }

////    double const density(isolates.GetDensity());
////    size_t const maxDegree(isolates.GetMaxDegree());
////    bool const &bRemoveIsolates(stackEvaluatedHalfVertices[depth + 1]);
    if (bRemoveIsolates)
////        if (onlyConsider.Size() > 50) {
////        if (depth <= 5)
            isolates.RemoveAllIsolates(0/*unused*/, vCliqueVertices, vRemoved, vAddedEdgesUnused /* unused */, false /* only consider updated vertices */);
////        } else {
////            isolates.RemoveAllIsolates(0/*unused*/, vCliqueVertices, vRemoved, vAddedEdgesUnused /* unused */, false /* only consider updated vertices */, onlyConsider);
////        }
////    }

////    onlyConsider.Clear();

////    cout << __LINE__ << ": density=" << density << ", max-degree=" << maxDegree << ", clique-vertices=" << vCliqueVertices.size() << ", other-removed=" << vRemoved.size()  << ", percent-total-removed=" << (vCliqueVertices.size() + vRemoved.size())/static_cast<double>(P.size())*100 << "%" << endl;
////
    vNewVertexOrder.resize(vVertexOrder.size());
    size_t uNewIndex(0);
    for (int const candidate : vVertexOrder) {
////    if (!m_bInvert && m_AdjacencyMatrix[chosenVertex][candidate]) vNewVertexOrder[uNewIndex++] = candidate;
        if (isolates.GetInGraph().Contains(candidate)) vNewVertexOrder[uNewIndex++] = candidate;
    }
    vNewVertexOrder.resize(uNewIndex);

    R.insert(R.end(), vCliqueVertices.begin(), vCliqueVertices.end());

////    vector<vector<int>> vComponents;
////    GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());
////
////    if (vComponents.size() > 1) {
////        cerr << "# connected components       : " << vComponents.size() << endl << flush;
////        cerr << "size of connected components : ";
////        cout << "[ ";
////        for (vector<int> const& vComponent : vComponents) {
////            cout << vComponent.size() << " ";
////        }
////        cout << "]" << endl << flush;
////    }
}

void TesterStaticOrderMISS::GetNewOrderNoIsolates(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
////    cout << "Old order: ";
////    for (int const vertex : vVertexOrder) {
////        cout << vertex << " ";
////    }
////    cout << endl;
    vNewVertexOrder.resize(vVertexOrder.size());
    {
        size_t uNewIndex(0);
        for (int const candidate : vVertexOrder) {
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

////    cout << "New order: ";
////    for (int const vertex : vNewVertexOrder) {
////        cout << vertex << " ";
////    }
////    cout << endl;
}

void TesterStaticOrderMISS::ProcessOrderAfterRecursionNoIsolates(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex)
{
////    cout << "LightWeightStaticOrderMISS::ProcessOrderAfterRecursion" << endl;
    if (chosenVertex == -1) return;
////    cout << "    order: ";
////    for (int const vertex : vVertexOrder) {
////        cout << vertex << " ";
////    }
////    cout << endl;
////    cout << "# vertices=" << vVertexOrder.size() << endl << flush;
    ////        vVertexOrder.erase(find(vVertexOrder.begin(), vVertexOrder.end(), chosenVertex));
    // try searching from end, might be faster in general...
////    if (find(vVertexOrder.begin(), vVertexOrder.end(), chosenVertex) != vVertexOrder.end()) {
////        cout << "vertex " << chosenVertex << " is in ordering..." << endl << flush;
////    } else {
////        cout << "vertex " << chosenVertex << " is not in ordering..." << endl << flush;
////    }
    size_t indexAfterVertex(0);
    for (indexAfterVertex = vVertexOrder.size(); indexAfterVertex >= 1; indexAfterVertex--) {
        if (vVertexOrder[indexAfterVertex-1] == chosenVertex) {
            break;
        }
    }

    for (; indexAfterVertex < vVertexOrder.size(); indexAfterVertex++) {
        vVertexOrder[indexAfterVertex-1] = vVertexOrder[indexAfterVertex];
    }

    vVertexOrder.pop_back();
    R.pop_back();
}

void TesterStaticOrderMISS::RunRecursive(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors)
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
            cout << "Evaluated " << nodeCount << " nodes. " << GetTimeInSeconds(clock() - startTime) << endl;
            PrintState();
        }
    }

    size_t numLeft = P.size();
    for (; numLeft > 0; --numLeft) {
        if (R.size() + vColors[numLeft-1] <= m_uMaximumCliqueSize) { break; }
    }

    numLeft = P.size() - numLeft;

    vector<int> evaluateP;

    bool pivot(false);

    if (!P.empty() && depth <= 3) {
        int const nextVertexToEvaluate(P.back());
        if (numLeft > isolates.Neighbors()[nextVertexToEvaluate].Size()) {
            pivot = true;
            vMarkedVertices[nextVertexToEvaluate] = true;
            for (int const neighbor : isolates.Neighbors()[nextVertexToEvaluate]) {
                vMarkedVertices[neighbor] = true;
            }

            for (size_t index = 0; index < P.size(); ++index) {
                int const vertex(P[index]);
                if (vMarkedVertices[vertex]) {
                    evaluateP.push_back(vertex);
                }
            }

            vMarkedVertices[nextVertexToEvaluate] = false;
            for (int const neighbor : isolates.Neighbors()[nextVertexToEvaluate]) {
                vMarkedVertices[neighbor] = false;
            }
        }
    }

    while (!evaluateP.empty()) {
        if (depth == 0) {
            if (!m_bQuiet) {
                cout << "Only " << evaluateP.size() << " more vertices to go! " << GetTimeInSeconds(clock() - startTime) << endl;
            }
        }

        int const nextVertex(evaluateP.back()); evaluateP.pop_back();

        if (!isolates.GetInGraph().Contains(nextVertex)) {
            continue;
        }

        size_t const numLeft = evaluateP.size();
        if (numLeft > 10) {
            cout << "depth = " << depth << ", P.size = " << P.size() << ", P.left = " << numLeft << ", neighbors=" << isolates.Neighbors()[nextVertex].Size() << endl;
        }

        ////        cout << depth << ": Choosing next vertex: " << nextVertex << endl;

        GetNewOrder(vNewVertexOrder, vVertexOrder, P, nextVertex);

        if (!vNewVertexOrder.empty()) {
            vNewP.resize(vNewVertexOrder.size());
            vNewColors.resize(vNewVertexOrder.size());
            Color(vNewVertexOrder/* evaluation order */, vNewP /* color order */, vNewColors);
            bool const bSwitchToNoIsolatesAlgorithm((vNewColors.size() < 5) || (vNewColors[vNewColors.size()-5] + R.size() <= m_uMaximumCliqueSize));
            if (bSwitchToNoIsolatesAlgorithm) {
                depth++;
                RunRecursiveNoIsolates(vNewP, vNewVertexOrder, cliques, vNewColors);
                depth--;
            } else {
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
            }
        } else if (R.size() > m_uMaximumCliqueSize) {
            cliques.back().clear();
            cliques.back().insert(cliques.back().end(), R.begin(), R.end());
            ExecuteCallBacks(cliques.back());
            m_uMaximumCliqueSize = R.size();
            timeToLargestClique = clock() - startTime;
        }

        bool bPIsEmpty(evaluateP.empty());
        ProcessOrderAfterRecursion(vVertexOrder, P, vColors, nextVertex);

        ////        if (R.size() > m_uMaximumCliqueSize && bPIsEmpty && P.empty()) {
        ////            cout << "ERROR!" << endl << flush;
        ////        }

        if (!bPIsEmpty && evaluateP.empty()) {
            if (R.size() > m_uMaximumCliqueSize) {
                cliques.back().clear();
                cliques.back().insert(cliques.back().end(), R.begin(), R.end());
                ExecuteCallBacks(cliques.back());
                m_uMaximumCliqueSize = R.size();
                timeToLargestClique = clock() - startTime;
            }
        }
    }

    while (!P.empty() && !pivot) {
        if (depth == 0) {
            if (!m_bQuiet) {
                cout << "Only " << P.size() << " more vertices to go! " << GetTimeInSeconds(clock() - startTime) << endl;
            }
        }

        int const largestColor(vColors.back());
        if (R.size() + largestColor <= m_uMaximumCliqueSize) {
            TesterMISQ::ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
            P.clear();
            return;
        }

        vColors.pop_back();
        int const nextVertex(P.back()); P.pop_back();


    numLeft = P.size();
    for (; numLeft > 0; --numLeft) {
        if (R.size() + vColors[numLeft-1] <= m_uMaximumCliqueSize) { break; }
    }

    numLeft = P.size() - numLeft;

////    if (numLeft > 10) {
////        cout << "depth = " << depth << ", P.size = " << P.size() << ", P.left = " << numLeft << ", neighbors=" << isolates.Neighbors()[nextVertex].Size() << endl;
////    }

////        cout << depth << ": Choosing next vertex: " << nextVertex << endl;

        GetNewOrder(vNewVertexOrder, vVertexOrder, P, nextVertex);

        if (!vNewVertexOrder.empty()) {
            vNewP.resize(vNewVertexOrder.size());
            vNewColors.resize(vNewVertexOrder.size());
            Color(vNewVertexOrder/* evaluation order */, vNewP /* color order */, vNewColors);
            bool const bSwitchToNoIsolatesAlgorithm((vNewColors.size() < 5) || (vNewColors[vNewColors.size()-5] + R.size() <= m_uMaximumCliqueSize));
            if (bSwitchToNoIsolatesAlgorithm) {
                depth++;
                RunRecursiveNoIsolates(vNewP, vNewVertexOrder, cliques, vNewColors);
                depth--;
            } else {
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
            }
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

    TesterMISQ::ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
    P.clear();

    vNewColors.clear();
    vNewP.clear();
}

void TesterStaticOrderMISS::RunRecursiveNoIsolates(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors)
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
            cout << "Evaluated " << nodeCount << " nodes. " << GetTimeInSeconds(clock() - startTime) << endl;
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
                cout << "Only " << P.size() << " more vertices to go! " << GetTimeInSeconds(clock() - startTime) << endl;
            }
        }

        int const largestColor(vColors.back());
        if (R.size() + largestColor <= m_uMaximumCliqueSize) {
            LightWeightMISQ::ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
            P.clear();
            return;
        }

        vColors.pop_back();
        int const nextVertex(P.back()); P.pop_back();

////        cout << depth << ": Choosing next vertex: " << nextVertex << endl;

        GetNewOrderNoIsolates(vNewVertexOrder, vVertexOrder, P, nextVertex);

        if (!vNewVertexOrder.empty()) {
            vNewP.resize(vNewVertexOrder.size());
            vNewColors.resize(vNewVertexOrder.size());
            Color(vNewVertexOrder/* evaluation order */, vNewP /* color order */, vNewColors);
#ifdef PREPRUNE
            if (R.size() + vNewColors.back() > m_uMaximumCliqueSize) {
                depth++;
                RunRecursiveNoIsolates(vNewP, vNewVertexOrder, cliques, vNewColors);
                depth--;
            }
#else
            depth++;
            RunRecursiveNoIsolates(vNewP, vNewVertexOrder, cliques, vNewColors);
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
        ProcessOrderAfterRecursionNoIsolates(vVertexOrder, P, vColors, nextVertex);

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

    LightWeightMISQ::ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
    P.clear();

    vNewColors.clear();
    vNewP.clear();
}

