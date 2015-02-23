#include "TesterStaticOrderMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"
#include "Tools.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>

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

    vector<int> X;

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

////    for (size_t index = 0; /*P.size() - numLeft;*/ index < P.size(); ++index) {
////        vector<vector<int>> vComponents;
////////        int const vertexToRemove(P[index]);
////////        isolates.RemoveVertex(vertexToRemove);
////        GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());
////////        vector<int> const vToReplace(1, vertexToRemove);
////////        isolates.ReplaceAllRemoved(vToReplace);
////        if (vComponents.size() > 1) {
////        bool printMessage(false);
////        ostringstream strm;
////        strm << depth << ": components : ";
////        for (size_t index = 0; /*P.size() - numLeft;*/ index < P.size(); ++index) {
////            size_t firstComponent(string::npos);
////            for (size_t componentIndex = 0; componentIndex < vComponents.size(); ++componentIndex) {
////                if (find(vComponents[componentIndex].begin(), vComponents[componentIndex].end(), P[index]) != vComponents[componentIndex].end()) {
////                    if (firstComponent == string::npos) {
////                        firstComponent = componentIndex;
////                    } else if (firstComponent != componentIndex) {
////                        printMessage = true;
////                    }
////                    strm << componentIndex << " ";
////                }
////            }
////        }
////////        if (printMessage) {
////        if (true) {
////            cout << strm.str();
////            cout << endl;
////            cerr << "# connected components       : " << vComponents.size() << endl << flush;
////            cerr << "size of connected components : ";
////            cout << "[ ";
////            for (vector<int> const& vComponent : vComponents) {
////                cout << vComponent.size() << "(mindegree=";
////                size_t minDegree(string::npos);
////                for (int const vertex : vComponent) {
////                    if (isolates.Neighbors()[vertex].Size() < minDegree) {
////                        minDegree = isolates.Neighbors()[vertex].Size();
////                    }
////                }
////                cout << minDegree << ") ";
////            }
////            cout << "]" << endl << flush;
////        }
////        }
////    }

    vector<int> evaluateP;
////    vector<int> evaluateColors;

    bool pivot(false);

    if (false ) {////!P.empty() /*&& depth <= 2*/) {
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
////                    evaluateColors.push_back(vertex);
                }
            }

            vMarkedVertices[nextVertexToEvaluate] = false;
            for (int const neighbor : isolates.Neighbors()[nextVertexToEvaluate]) {
                vMarkedVertices[neighbor] = false;
            }
        }
    }

    if (!evaluateP.empty()) {
        if (!m_bQuiet) cout << depth << ": Switching to pivot mode." << endl;
    }

    while (!evaluateP.empty()) {
        if (depth == 0) {
            if (!m_bQuiet) {
                cout << "Only " << evaluateP.size() << " more vertices to go! " << GetTimeInSeconds(clock() - startTime) << endl;
            }
        }

        int const nextVertex(evaluateP.back()); evaluateP.pop_back();

////        // check if we should return to non-pivot version of the algorithm
////        P.erase(find(P.begin(), P.end(), nextVertex));
////        vVertexOrder.erase(find(vVertexOrder.begin(), vVertexOrder.end(), nextVertex));
        Color(vVertexOrder, P, vColors);
        numLeft = P.size();
        for (; numLeft > 0; --numLeft) {
            if (R.size() + vColors[numLeft-1] <= m_uMaximumCliqueSize) { break; }
        }
        numLeft = P.size() - numLeft;
        if (numLeft < 3.0*evaluateP.size()/4.0) {
            if (!m_bQuiet) cout << depth << ": Switching out of pivot mode" << endl;
            pivot = false;
            break;
        }

        if (!isolates.GetInGraph().Contains(nextVertex)) { // was filtered by reducer...
            X.push_back(nextVertex);
            continue;
        }

        for (int const neighbor : isolates.Neighbors()[nextVertex]) {
            vMarkedVertices[neighbor] = true;
        }

        bool dominates(false);
        for (int const x : X) {
            dominates = true;
            for (int const neighborX : m_AdjacencyArray[x]) {
                if (isolates.GetInGraph().Contains(neighborX) && !vMarkedVertices[neighborX]) {
                    dominates = false;
                    break;
                }
            }

            if (dominates) {
                break;
            }
        }

        for (int const neighbor : isolates.Neighbors()[nextVertex]) {
            vMarkedVertices[neighbor] = false;
        }

        if (dominates) {
            if (!m_bQuiet)
                cout << depth << ": domination check removed " << nextVertex << endl;
            isolates.RemoveVertex(nextVertex);
            vector<int> &vRemovedVerticesToReplace(stackPersistentOther[depth+1]);
            vRemovedVerticesToReplace.push_back(nextVertex);
            continue;
        }

////        int const largestColor(evaluateColors.back());
////        if (R.size() + largestColor <= m_uMaximumCliqueSize) { //// && !selectMinNeighbors) {
////            TesterMISQ::ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
////            P.clear();
////            return;
////        }

////        evaluateColors.pop_back();

        size_t const numLeft = evaluateP.size();
        if (numLeft > 10 && !m_bQuiet) {
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

        X.push_back(nextVertex);

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

    bool const selectMinNeighbors(false); ////depth <= 1);
    while (!P.empty() && !pivot) {
        if (depth == 0) {
            if (!m_bQuiet) {
                cout << "Only " << P.size() << " more vertices to go! " << GetTimeInSeconds(clock() - startTime) << endl;
            }
        }

        int const largestColor(vColors.back());
        if (R.size() + largestColor <= m_uMaximumCliqueSize && !selectMinNeighbors) {
            TesterMISQ::ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
            P.clear();
            return;
        }

        int vertexToChoose(P.back());
        size_t maxNeighborCount(0);

        size_t numLeft = P.size();
        for (; numLeft > 0; --numLeft) {
            if (R.size() + vColors[numLeft-1] <= m_uMaximumCliqueSize) { break; }
        }

        if (selectMinNeighbors) {
            for (size_t index = P.size(); index > numLeft; --index) {
                vMarkedVertices[P[index-1]] = true;
            }

            for (size_t index = P.size(); index > numLeft; --index) {
                int const vertex(P[index-1]);
                size_t neighborCount(0);
                for (int const neighbor : isolates.Neighbors()[vertex]) {
                    if (vMarkedVertices[neighbor]) neighborCount++;
                }
                if (neighborCount > maxNeighborCount) {
                    maxNeighborCount = neighborCount;
                    vertexToChoose = vertex;
                }
            }

            for (size_t index = P.size(); index > numLeft; --index) {
                vMarkedVertices[P[index-1]] = false;
            }
        }

////        cout << "Before subtraction numleft=" << numLeft << endl;
        numLeft = P.size() - numLeft;
////        cout << "After  subtraction numleft=" << numLeft << endl;

        int const nextVertex(vertexToChoose); 
        if (nextVertex != P.back()) {
            P.erase(find(P.begin(), P.end(), nextVertex));
            Color(vVertexOrder, P, vColors);
        } else {
            P.pop_back();
            vColors.pop_back();
        }

        for (int const neighbor : isolates.Neighbors()[nextVertex]) {
            vMarkedVertices[neighbor] = true;
        }

        bool dominates(false);
        for (int const x : X) {
            dominates = true;
            for (int const neighborX : m_AdjacencyArray[x]) {
                if (isolates.GetInGraph().Contains(neighborX) && !vMarkedVertices[neighborX]) {
                    dominates = false;
                    break;
                }
            }

            if (dominates) {
                break;
            }
        }

        for (int const neighbor : isolates.Neighbors()[nextVertex]) {
            vMarkedVertices[neighbor] = false;
        }

        if (dominates) {
            if (!m_bQuiet)
                cout << depth << ": domination check removed " << nextVertex << endl;
            isolates.RemoveVertex(nextVertex);
            vector<int> &vRemovedVerticesToReplace(stackPersistentOther[depth+1]);
            vRemovedVerticesToReplace.push_back(nextVertex);
            continue;
        }

////        if (numLeft > 10) {
////            cout << "depth = " << depth << ", P.size = " << P.size() << ", P.left = " << numLeft << ", neighbors=" << isolates.Neighbors()[nextVertex].Size() << endl;
////        }

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

        X.push_back(nextVertex);

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

