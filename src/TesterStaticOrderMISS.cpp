#include "TesterStaticOrderMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"
#include "Tools.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>

////#define MAX_TWO_NEIGHBORHOOD
////#define MIN_COMPONENT
#define MIN_SUBPROBLEM

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

////    OrderingTools::InitialOrderingReduction(isolates, vVertexOrder);
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
            cout << "Evaluated " << nodeCount << " nodes. " << Tools::GetTimeInSeconds(clock() - startTime) << endl;
            PrintState();
        }
    }

    size_t numLeft = P.size();
    for (; numLeft > 0; --numLeft) {
        if (R.size() + vColors[numLeft-1] <= m_uMaximumCliqueSize) { break; }
    }

    numLeft = P.size() - numLeft;
////    if (depth <=1)
////    cout << depth << ": numleft=" << numLeft << endl;

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

    if (false) { ////!P.empty() && depth <= 1) {
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
                cout << "Only " << evaluateP.size() << " more vertices to go! " << Tools::GetTimeInSeconds(clock() - startTime) << endl;
            }
        }

        int const nextVertex(evaluateP.back()); evaluateP.pop_back();

////        // check if we should return to non-pivot version of the algorithm
////        P.erase(find(P.begin(), P.end(), nextVertex));
////        vVertexOrder.erase(find(vVertexOrder.begin(), vVertexOrder.end(), nextVertex));
        Color(vVertexOrder, P, vColors);
////        cout << __LINE__ << ": Recoloring..." << endl;
        numLeft = P.size();
        for (; numLeft > 0; --numLeft) {
            if (R.size() + vColors[numLeft-1] <= m_uMaximumCliqueSize) { break; }
        }
        numLeft = P.size() - numLeft;
////        if (numLeft < 3.0*evaluateP.size()/4.0) {
////            if (!m_bQuiet) cout << depth << ": Switching out of pivot mode" << endl;
////            pivot = false;
////            break;
////        }

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
////        cout << __LINE__ << endl;
////            if (!m_bQuiet)
////                cout << depth << ": domination check removed " << nextVertex << endl;
            isolates.RemoveVertex(nextVertex);
            vector<int> &vRemovedVerticesToReplace(stackPersistentOther[depth+1]);
            vRemovedVerticesToReplace.push_back(nextVertex);
////        cout << __LINE__ << endl;
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
////
    vector<pair<int,int>> vAddedEdgesUnused;
#ifdef MIN_COMPONENT
    auto ChooseNextVertex = [this, &vAddedEdgesUnused](Isolates3<ArraySet> &theIsolates, vector<int> const &activeSet)
    {
////        cout << "Before choosing vertex: graph contains " << theIsolates.GetInGraph().Size() << " elements" << endl;
        size_t bestComponentSize(0); ///string::npos);
        int    bestVertex(-1);
////        cout << __LINE__ << endl;
////        vector<int> vToConsider(theIsolates.GetInGraph().begin(), theIsolates.GetInGraph().end());
        vector<int> vToConsider(activeSet);
////        vector<int> vToConsider;
////        if (activeSet.size() > 1000) {
////            vToConsider.insert(vToConsider.end(), activeSet.end()-1000, activeSet.end());
////        } else {
////            vToConsider = activeSet;
////        }
            for (int const i : vToConsider) {
                vector<int> v_iNeighbors;
                theIsolates.RemoveVertexAndNeighbors(i, v_iNeighbors);
////                theIsolates.RemoveAllIsolates(0, v_iNeighbors, v_iNeighbors, vAddedEdgesUnused, false);

////                vector<vector<int>> vComponents;
////                GraphTools::ComputeConnectedComponents(theIsolates, vComponents, m_AdjacencyArray.size());
////                size_t uSizeOfLargestComponent(0);
////                for (vector<int> const &vComponent : vComponents) {
////                    uSizeOfLargestComponent = max(vComponent.size(), uSizeOfLargestComponent);
////                }

////                if (i == 58)
////                    cout << "Can remove " << i << " for max component of size " << uSizeOfLargestComponent << endl;

                size_t numRemovedFromActiveSet(0);
                for (int const activeVertex : activeSet) {
                    if (!theIsolates.GetInGraph().Contains(activeVertex)) {
////                        removedFromActiveSet = true;
////                        break;
                        numRemovedFromActiveSet++;
                    }
                }

////                if (uSizeOfLargestComponent < bestComponentSize) { //// && removedFromActiveSet) { // needs to decrease active set by at least one
////                    bestComponentSize = uSizeOfLargestComponent;
                if (numRemovedFromActiveSet > bestComponentSize) { //// && removedFromActiveSet) { // needs to decrease active set by at least one
                    bestComponentSize = numRemovedFromActiveSet;
////                    cout << "best so far: remove " << i << " for " << vComponents.size() << " components, with max component of size " << bestComponentSize << endl;
                    bestVertex = i;
                }

                v_iNeighbors.push_back(i);
                theIsolates.ReplaceAllRemoved(v_iNeighbors);
            }
////        cout << "After  choosing vertex: graph contains " << theIsolates.GetInGraph().Size() << " elements" << endl;
////        cout << "best vertex: " << bestVertex  << " for max component of size " << bestComponentSize << endl;
////        cout << __LINE__ << endl;
        return bestVertex;
    };
#endif // MIN_COMPONENT

#ifdef MIN_SUBPROBLEM
    auto ChooseNextVertex = [this, &vAddedEdgesUnused](Isolates3<ArraySet> &theIsolates, vector<int> const &activeSet)
    {

        size_t best(ULONG_MAX);
        int savedi(-1), savedj(-1), savedk(-1);
        vector<int> vVertices(theIsolates.GetInGraph().begin(), theIsolates.GetInGraph().end());
        vector<int> vNumEdgesInTwoNeighborhood(m_AdjacencyArray.size(), 0);
        for (int const vertex : vVertices) {
            for (int const neighbor : theIsolates.Neighbors()[vertex]) {
                vNumEdgesInTwoNeighborhood[vertex] += theIsolates.Neighbors()[neighbor].Size();
            }
        }

        int const size = vVertices.size();

        ////    sort (vVertices.begin(), vVertices.end(), [this](int const left, int const right) { return m_AdjacencyArray[left].size() > m_AdjacencyArray[right].size(); });
        sort (vVertices.begin(), vVertices.end(), [&vNumEdgesInTwoNeighborhood](int const left, int const right) { return vNumEdgesInTwoNeighborhood[left] > vNumEdgesInTwoNeighborhood[right]; });

        for (size_t i = 0; i < size*0.10; i++) {
            vector<int> v_iNeighbors;
            theIsolates.RemoveVertexAndNeighbors(vVertices[i], v_iNeighbors);
            theIsolates.RemoveAllIsolates(0, v_iNeighbors, v_iNeighbors, vAddedEdgesUnused, false);

            vector<int> vLeftOver1(theIsolates.GetInGraph().begin(), theIsolates.GetInGraph().end());
            vector<int> vTwoNeighborhood1(m_AdjacencyArray.size(), 0);
            for (int const vertex : vLeftOver1) {
                for (int const neighbor : theIsolates.Neighbors()[vertex]) {
                    vTwoNeighborhood1[vertex] += theIsolates.Neighbors()[neighbor].Size();
                }
            }

            sort (vLeftOver1.begin(), vLeftOver1.end(), [&vTwoNeighborhood1](int const left, int const right) { return vTwoNeighborhood1[left] > vTwoNeighborhood1[right]; });
            ////        for (size_t j = i+1; j < i + size*0.01; j++) {
            for (size_t j = 0; j < min(int(vLeftOver1.size()*0.01), 100); j++) {
                vector<int> v_jNeighbors;
                int const secondVertex(vLeftOver1[j]);
                if (theIsolates.GetInGraph().Contains(secondVertex)) {
                    v_jNeighbors.push_back(secondVertex);
                    theIsolates.RemoveVertexAndNeighbors(secondVertex, v_jNeighbors);
                    theIsolates.RemoveAllIsolates(0, v_jNeighbors, v_jNeighbors, vAddedEdgesUnused, false);
                }

                vector<int> vLeftOver2(theIsolates.GetInGraph().begin(), theIsolates.GetInGraph().end());
                vector<int> vTwoNeighborhood2(m_AdjacencyArray.size(), 0);
                for (int const vertex : vLeftOver2) {
                    for (int const neighbor : theIsolates.Neighbors()[vertex]) {
                        vTwoNeighborhood2[vertex] += theIsolates.Neighbors()[neighbor].Size();
                    }
                }

                sort (vLeftOver2.begin(), vLeftOver2.end(), [&vTwoNeighborhood2](int const left, int const right) { return vTwoNeighborhood2[left] > vTwoNeighborhood2[right]; });
                ////            for (size_t k = j+1; k < j+ size*0.05; k++) {
                for (size_t k = 0; k < min(int(vLeftOver2.size()*0.01), 100); k++) {
                    ////                if (loops %10 == 0) { break; }
                    int const thirdVertex(vLeftOver2[k]);
                    vector<int> v_kNeighbors;
                    if (theIsolates.GetInGraph().Contains(thirdVertex)) {
                        v_kNeighbors.push_back(thirdVertex);
                        theIsolates.RemoveVertexAndNeighbors(thirdVertex, v_kNeighbors);
                        theIsolates.RemoveAllIsolates(0, v_kNeighbors, v_kNeighbors, vAddedEdgesUnused, false);
                    }

                    vector<vector<int>> vComponents;
                    GraphTools::ComputeConnectedComponents(theIsolates, vComponents, m_AdjacencyArray.size());
                    size_t biggest(0);
                    for (vector<int> const &vComponent : vComponents) {
                        biggest = max(vComponent.size(), biggest);
                    }


                if (biggest < best) {
                    best = biggest;
                    cout << "best so far: remove " << vVertices[i] << " " << vLeftOver1[j] << " " << vLeftOver2[k] << " for max component of size " << best << endl;
                    cout << "best so far: components=" << vComponents.size() << ", max-component-size=" << best << endl;
                    savedi = vVertices[i];
                    savedj = vLeftOver1[j];
                    savedk = vLeftOver2[k];
                }
                theIsolates.ReplaceAllRemoved(v_kNeighbors);
            }
            theIsolates.ReplaceAllRemoved(v_jNeighbors);
////            if (loops % 100 == 0) break;
        }
        v_iNeighbors.push_back(vVertices[i]);
        theIsolates.ReplaceAllRemoved(v_iNeighbors);

////        int const percentage((int((i*size*0.15 + j*0.10 + k)*100.0/((size*0.01)*(size*0.05)*(size*0.10)))));
////        if (percentage > lastPercentage + 1) {
////            cout << "finished evaluating " << i << "/" << size*0.10 << " first vertices" << endl;
////            lastPercentage = percentage;
////        }
        ////        if (loops % 5000 == 0) { break; }
////            cout << "Removing " << savedi << " " << savedj << " " << savedk << endl;
////            vector<int> vRemovedNeighbors;
////            theIsolates.RemoveVertexAndNeighbors(savedi, vRemovedNeighbors);
////            theIsolates.RemoveVertexAndNeighbors(savedj, vRemovedNeighbors);
////            theIsolates.RemoveVertexAndNeighbors(savedk, vRemovedNeighbors);
////        }
        }

        if (vNumEdgesInTwoNeighborhood[savedi] >= vNumEdgesInTwoNeighborhood[savedj] && vNumEdgesInTwoNeighborhood[savedi] >= vNumEdgesInTwoNeighborhood[savedk]) {
            return savedi;
        }

        if (vNumEdgesInTwoNeighborhood[savedj] >= vNumEdgesInTwoNeighborhood[savedi] && vNumEdgesInTwoNeighborhood[savedj] >= vNumEdgesInTwoNeighborhood[savedk]) {
            return savedj;
        }

        return savedk;
    };

#endif // MIN_SUBPROBLEM

#ifdef MIN_COMPONENT
    auto ChooseNextVertex = [this, &vAddedEdgesUnused](Isolates3<ArraySet> &theIsolates, vector<int> const &activeSet)
    {
////        cout << "Before choosing vertex: graph contains " << theIsolates.GetInGraph().Size() << " elements" << endl;
        size_t bestComponentSize(0); ///string::npos);
        int    bestVertex(-1);
////        cout << __LINE__ << endl;
////        vector<int> vToConsider(theIsolates.GetInGraph().begin(), theIsolates.GetInGraph().end());
        vector<int> vToConsider(activeSet);
////        vector<int> vToConsider;
////        if (activeSet.size() > 1000) {
////            vToConsider.insert(vToConsider.end(), activeSet.end()-1000, activeSet.end());
////        } else {
////            vToConsider = activeSet;
////        }
            for (int const i : vToConsider) {
                vector<int> v_iNeighbors;
                theIsolates.RemoveVertexAndNeighbors(i, v_iNeighbors);
////                theIsolates.RemoveAllIsolates(0, v_iNeighbors, v_iNeighbors, vAddedEdgesUnused, false);

////                vector<vector<int>> vComponents;
////                GraphTools::ComputeConnectedComponents(theIsolates, vComponents, m_AdjacencyArray.size());
////                size_t uSizeOfLargestComponent(0);
////                for (vector<int> const &vComponent : vComponents) {
////                    uSizeOfLargestComponent = max(vComponent.size(), uSizeOfLargestComponent);
////                }

////                if (i == 58)
////                    cout << "Can remove " << i << " for max component of size " << uSizeOfLargestComponent << endl;

                size_t numRemovedFromActiveSet(0);
                for (int const activeVertex : activeSet) {
                    if (!theIsolates.GetInGraph().Contains(activeVertex)) {
////                        removedFromActiveSet = true;
////                        break;
                        numRemovedFromActiveSet++;
                    }
                }

////                if (uSizeOfLargestComponent < bestComponentSize) { //// && removedFromActiveSet) { // needs to decrease active set by at least one
////                    bestComponentSize = uSizeOfLargestComponent;
                if (numRemovedFromActiveSet > bestComponentSize) { //// && removedFromActiveSet) { // needs to decrease active set by at least one
                    bestComponentSize = numRemovedFromActiveSet;
////                    cout << "best so far: remove " << i << " for " << vComponents.size() << " components, with max component of size " << bestComponentSize << endl;
                    bestVertex = i;
                }

                v_iNeighbors.push_back(i);
                theIsolates.ReplaceAllRemoved(v_iNeighbors);
            }
////        cout << "After  choosing vertex: graph contains " << theIsolates.GetInGraph().Size() << " elements" << endl;
////        cout << "best vertex: " << bestVertex  << " for max component of size " << bestComponentSize << endl;
////        cout << __LINE__ << endl;
        return bestVertex;
M_SIZE
#endif // MIN_COMPONENT

////
////    auto NewChooseNextVertex = [this, &vAddedEdgesUnused](
////            Isolates3<ArraySet> &theIsolates,
////            vector<int> &vAdjunctOrdering,
////            size_t const uCurrentCliqueSize,
////            size_t const uMaximumCliqueSize,
////            vector<int> const &activeSet,
////            int const nextVertex)
////    {
////
////////        cout << "Test:    Max clique size=" << uMaximumCliqueSize << endl;
////        int vertexToChoose(-1);
////        size_t minNumLeft(string::npos);
////        vector<int> vNewAdjunctOrdering;
////        vector<int> vNewOrdering;
////        vector<int> vNewColoring;
////        vector<int> vToConsider(activeSet); ////theIsolates.GetInGraph().begin(), theIsolates.GetInGraph().end());
////        for (int const vertex : vToConsider) {
////            vector<int> vCliqueVertices;
////            vector<int> vRemoved;
////////            theIsolates.RemoveVertex(vertex);
////            theIsolates.RemoveVertexAndNeighbors(vertex, vRemoved);
////////            vRemoved.push_back(vertex);
////            vCliqueVertices.push_back(vertex);
////            theIsolates.RemoveAllIsolates(0, vCliqueVertices, vRemoved, vAddedEdgesUnused, false);
////
////            size_t uNewSize(0);
////            vNewAdjunctOrdering.resize(vAdjunctOrdering.size());
////            for (size_t index = 0; index < vAdjunctOrdering.size(); ++index) {
////                int const vertexInOrder(vAdjunctOrdering[index]);
////                if (theIsolates.GetInGraph().Contains(vertexInOrder)) {
////                    vNewAdjunctOrdering[uNewSize++] = vertexInOrder;
////                }
////            }
////
////            vNewAdjunctOrdering.resize(uNewSize);
////            vNewOrdering.resize(vNewAdjunctOrdering.size());
////            vNewColoring.resize(vNewAdjunctOrdering.size());
////
////////            if (vertex == nextVertex) {
////////                cout << "Test   Ordering: ";
////////                for (int const vertexInOrder : vNewAdjunctOrdering) {
////////                    cout << vertexInOrder << " ";
////////                }
////////                cout << endl;
////////            }
////
////
////////            GetNewOrder(vNewAdjunctOrdering, vAdjunctOrdering, vNewOrdering, vertex);
////////            vNewOrdering.resize(vNewAdjunctOrdering.size());
////////            vNewColoring.resize(vNewAdjunctOrdering.size());
////
////            coloringStrategy.Recolor(m_AdjacencyMatrix, vNewAdjunctOrdering, vNewOrdering, vNewColoring, uMaximumCliqueSize, uCurrentCliqueSize + vCliqueVertices.size());
////
////            if (vertex == nextVertex) { ////5 || vertex == 802 || vertex == 45 || vertex == 33) {
////////                cout << "Test:    Vertex " << vertex << " has " << uCurrentCliqueSize << " + " <<  vCliqueVertices.size() << " clique vertices " << endl;
////                cout << "Test   : current clique size is: " << uCurrentCliqueSize << " + " <<  vCliqueVertices.size() << endl;
////                cout << "Test   : max     clique size is: " << uMaximumCliqueSize << endl;
////                cout << "Test   : New P is: ";
////                for (int const pVertex : vNewOrdering) {
////                    cout << pVertex << " ";
////                }
////                cout << endl;
////                cout << "Test   : New Adjunct is: ";
////                for (int const aVertex : vNewAdjunctOrdering) {
////                    cout << aVertex << " ";
////                }
////                cout << endl;
////            }
////
////            size_t numLeft = vNewColoring.size();
////            for (; numLeft > 0; --numLeft) {
////                if (uCurrentCliqueSize + vCliqueVertices.size() + vNewColoring[numLeft-1] <= uMaximumCliqueSize) { break; }
////            }
////            numLeft = vNewColoring.size() - numLeft;
////            if (vertex == nextVertex) {
////                cout << "vertex " << vertex << ": P.left = " << numLeft << endl;
////            }
////            if (numLeft < minNumLeft) {
////                vertexToChoose = vertex;
////                minNumLeft = numLeft;
////////                cout << ", P.left = " << numLeft << endl;
////////                cout << "vertex " << vertex << ": P.left = " << numLeft << endl;
////            }
////
////            theIsolates.ReplaceAllRemoved(vCliqueVertices);
////            theIsolates.ReplaceAllRemoved(vRemoved);
////        }
////
////////        cout << "Choosing next vertex " << vertexToChoose << " will minimize the number of vertices for consideration to " << minNumLeft << endl;
////
////        return vertexToChoose;
////    };
////
////    bool const selectMinNeighbors(false); ////depth <= 1);
    bool firstIteration(true);
    clock_t start(clock());

#if (defined(MIN_COMPONENT) || defined(MAX_TWO_NEIGHBORHOOD) || defined(MIN_SUBPROBLEM))
    bool selectMinComponent(isolates.GetInGraph().Size() > 750); ////depth <= 3 && isolates.GetInGraph().Size() > 950 || firstIteration);
    bool const wasLarge(isolates.GetInGraph().Size() > 900); ////depth <= 3 && isolates.GetInGraph().Size() > 950 || firstIteration);
    bool checkTime(isolates.GetInGraph().Size() > 700);
    int numEvaluated(0);
#else
    bool selectMinComponent(false);
#endif // MIN_COMPONENT
////    bool selectMinComponent(false);
////    bool checkTime(false); ////isolates.GetInGraph().Size() > 700);
    while (!P.empty() && !pivot) {
        // if one criteria is slow, use the other one.
#if (defined(MIN_COMPONENT) || defined(MAX_TWO_NEIGHBORHOOD) || defined(MIN_SUBPROBLEM))
        if (isolates.GetInGraph().Size()<=700) {
            selectMinComponent = false;
            checkTime = false;
        }
        if (wasLarge && isolates.GetInGraph().Size()<=850) {
            selectMinComponent = false;
            checkTime = false;
        }
        if (selectMinComponent) {
            numEvaluated++;
        }
////        if (checkTime) {
////            if ((clock() - start)/CLOCKS_PER_SEC > 10 || firstIteration) {
////            ////            selectMinComponent = !selectMinComponent;
////                selectMinComponent = (isolates.GetInGraph().Size() > 750);
////                selectMinComponent = true;
////            } else {
////                selectMinComponent = false;
////            }
////        }
#endif // 0

        numLeft = P.size();
        for (; numLeft > 0; --numLeft) {
            if (R.size() + vColors[numLeft-1] <= m_uMaximumCliqueSize) { break; }
        }

        numLeft = P.size() - numLeft;

////        if (numEvaluated > numLeft) {
////            selectMinComponent = false;
////            checkTime = false;
////        }

////        start = clock();

        firstIteration = false;
        if (depth == 0) {
            if (!m_bQuiet) {
                cout << "Only " << P.size() << "(real=" << numLeft << ") more vertices to go! " << Tools::GetTimeInSeconds(clock() - startTime) << endl;
            }
        }

        int const largestColor(vColors.back());
        if (R.size() + largestColor <= m_uMaximumCliqueSize /*&& !selectMinNeighbors && !selectMinComponent*/) {
////            cout << __LINE__ << endl;
////            cout << "P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////            cout << "Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
            TesterMISQ::ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
////        cout << __LINE__ << endl;
////        cout << depth << ": P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
            P.clear();
            return;
        }

        int vertexToChoose(P.back());
////        size_t maxNeighborCount(0);
////
////        size_t numLeft = P.size();
////        for (; numLeft > 0; --numLeft) {
////            if (R.size() + vColors[numLeft-1] <= m_uMaximumCliqueSize) { break; }
////        }
////
////        if (selectMinNeighbors) {
////            for (size_t index = P.size(); index > numLeft; --index) {
////                vMarkedVertices[P[index-1]] = true;
////            }
////
////            for (size_t index = P.size(); index > numLeft; --index) {
////                int const vertex(P[index-1]);
////                size_t neighborCount(0);
////                for (int const neighbor : isolates.Neighbors()[vertex]) {
////                    if (vMarkedVertices[neighbor]) neighborCount++;
////                }
////                if (neighborCount > maxNeighborCount) {
////                    maxNeighborCount = neighborCount;
////                    vertexToChoose = vertex;
////                }
////            }
////
////            for (size_t index = P.size(); index > numLeft; --index) {
////                vMarkedVertices[P[index-1]] = false;
////            }
////        } else
        if (selectMinComponent){
////            cout << __LINE__ << endl;
////            cout << "P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////            cout << "Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
            size_t const realNumLeft = P.size() - numLeft;
            vector<int> const activeSet(P.begin() + numLeft, P.end());
            int const proposedVertexToChoose = ChooseNextVertex(isolates, activeSet);
////            int proposedVertexToChoose(-1);
////            if (firstIteration) {
////                firstIteration = false;
////                proposedVertexToChoose = ChooseNextVertex(isolates, activeSet);
////            } else {
////                evenIteration = !evenIteration;
////                if (evenIteration)
////                    proposedVertexToChoose = NewChooseNextVertex(isolates, vVertexOrder, R.size(), m_uMaximumCliqueSize, activeSet, -1);
////            }

            if (proposedVertexToChoose != -1) {
                vertexToChoose = proposedVertexToChoose;
            }
////            cout << __LINE__ << endl;
////            cout << "P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////            cout << "Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
        }

////        cout << "Before subtraction numleft=" << numLeft << endl;
////        numLeft = P.size() - numLeft;
////        if (depth <= 1)
////            cout << depth << ": loop numleft=" << numLeft << endl;

////        cout << __LINE__ << endl;
        int const nextVertex(vertexToChoose); 
////        int const nextVertex(P.back()); 
////        cout << "Removing vertex: " << nextVertex << endl;
        bool const vertexAtEnd(nextVertex == P.back());
////        cout << __LINE__ << endl;
////        cout << depth << ": P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": VO       " << ((find(vVertexOrder.begin(), vVertexOrder.end(), 201) != vVertexOrder.end()) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": P        " << ((find(P.begin(), P.end(), nextVertex) != P.end()) ? "contains " : "does not contain ") << nextVertex << endl;
////        cout << depth << ": VO       " << ((find(vVertexOrder.begin(), vVertexOrder.end(), nextVertex) != vVertexOrder.end()) ? "contains " : "does not contain ") << nextVertex << endl;
////        cout << depth << ": Isolates " << (isolates.GetInGraph().Contains(nextVertex) ? "contains " : "does not contain ") << nextVertex << endl;
        if (!vertexAtEnd) {
            vVertexOrder.erase(find(vVertexOrder.begin(), vVertexOrder.end(), nextVertex));
////        cout << __LINE__ << endl;
////        cout << depth << ": P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": VO       " << ((find(vVertexOrder.begin(), vVertexOrder.end(), 201) != vVertexOrder.end()) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": P        " << ((find(P.begin(), P.end(), nextVertex) != P.end()) ? "contains " : "does not contain ") << nextVertex << endl;
////        cout << depth << ": VO       " << ((find(vVertexOrder.begin(), vVertexOrder.end(), nextVertex) != vVertexOrder.end()) ? "contains " : "does not contain ") << nextVertex << endl;
////        cout << __LINE__ << endl;
            size_t uNewSize(0);
            for (size_t index = 0; index < P.size(); ++index) {
                if (P[index] == nextVertex) continue;
                P[uNewSize] = P[index];
                vColors[uNewSize++] = vColors[index];

            }
            P.resize(uNewSize);
            vColors.resize(uNewSize);

#ifdef RECOLOR
            vector<int> proposedP(vVertexOrder.size());
            vector<int> proposedColors(vVertexOrder.size());
////            cout << __LINE__ << ": Recoloring..." << endl;
            Color(vVertexOrder, proposedP, proposedColors);

            size_t currentNumLeft = P.size();
            for (; currentNumLeft > 0; --currentNumLeft) {
                if (R.size() + vColors[currentNumLeft-1] <= m_uMaximumCliqueSize) { break; }
            }

            currentNumLeft = P.size() - currentNumLeft;

            size_t proposedNumLeft = proposedP.size();
            for (; proposedNumLeft > 0; --proposedNumLeft) {
                if (R.size() + proposedColors[proposedNumLeft-1] <= m_uMaximumCliqueSize) { break; }
            }

            proposedNumLeft = P.size() - proposedNumLeft;

            if (proposedNumLeft < currentNumLeft) {
////                cout << __LINE__ << ": Recoloring..." << endl;
                P = proposedP;
                vColors = proposedColors;
            }
#endif // RECOLOR


////            P.erase(find(P.begin(), P.end(), nextVertex));
////        cout << __LINE__ << endl;
////        cout << depth << ": P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": VO       " << ((find(vVertexOrder.begin(), vVertexOrder.end(), 201) != vVertexOrder.end()) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": P        " << ((find(P.begin(), P.end(), nextVertex) != P.end()) ? "contains " : "does not contain ") << nextVertex << endl;
////        cout << depth << ": VO       " << ((find(vVertexOrder.begin(), vVertexOrder.end(), nextVertex) != vVertexOrder.end()) ? "contains " : "does not contain ") << nextVertex << endl;
////        cout << __LINE__ << endl;
        } else {
////            cout << "Popping off " << P.back() << endl;
////        cout << __LINE__ << endl;
            P.pop_back();
////        cout << __LINE__ << endl;
            vColors.pop_back();
////        cout << __LINE__ << endl;
        }
////        cout << __LINE__ << endl;
////        cout << depth << ": P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": VO       " << ((find(vVertexOrder.begin(), vVertexOrder.end(), 201) != vVertexOrder.end()) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;

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
////        cout << __LINE__ << endl;
////            if (!m_bQuiet)
////                cout << depth << ": domination check removed " << nextVertex << endl;
            isolates.RemoveVertex(nextVertex);
            if (vertexAtEnd) vVertexOrder.erase(find(vVertexOrder.begin(), vVertexOrder.end(), nextVertex));
            vector<int> &vRemovedVerticesToReplace(stackPersistentOther[depth+1]);
            vRemovedVerticesToReplace.push_back(nextVertex);
////        cout << __LINE__ << endl;
////        cout << __LINE__ << endl;
////        cout << depth << ": P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
            continue;
        }

////        cout << __LINE__ << endl;
////        cout << depth << ": P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;

////        if (numLeft > 10) {
////            cout << "depth = " << depth << ", P.size = " << P.size() << ", P.left = " << numLeft << ", neighbors=" << isolates.Neighbors()[nextVertex].Size() << endl;
////        }

        ////        cout << depth << ": Choosing next vertex: " << nextVertex << endl;

////            cout << __LINE__ << endl;
////            cout << "P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////            cout << "Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
        GetNewOrder(vNewVertexOrder, vVertexOrder, P, nextVertex);
////            cout << __LINE__ << endl;
////            cout << "P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////            cout << "Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;

        if (!vNewVertexOrder.empty()) {
            vNewP.resize(vNewVertexOrder.size());
            vNewColors.resize(vNewVertexOrder.size());
            Color(vNewVertexOrder/* evaluation order */, vNewP /* color order */, vNewColors);

////            if (depth <= 0) {
////                cout << "Actual : current clique size is: " << R.size() << endl;
////                cout << "Actual : max     clique size is: " << m_uMaximumCliqueSize << endl;
////                cout << "Actual : New P is: ";
////                for (int const pVertex : vNewP) {
////                    cout << pVertex << " ";
////                }
////                cout << endl;
////                cout << "Actual : New Adjunct is: ";
////                for (int const aVertex : vNewVertexOrder) {
////                    cout << aVertex << " ";
////                }
////                cout << endl;
////            }
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
////                    cout << __LINE__ << endl;
////                    cout << "P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////                    cout << "Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
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
////            cout << __LINE__ << endl;
////            cout << "P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////            cout << "Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
        ProcessOrderAfterRecursion(vVertexOrder, P, vColors, nextVertex);
////        cout << __LINE__ << endl;
////        cout << depth << ": P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////        cout << depth << ": Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;

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

////            cout << __LINE__ << endl;
////            cout << "P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////            cout << "Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
    TesterMISQ::ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
////            cout << __LINE__ << endl;
////            cout << "P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////            cout << "Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
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

