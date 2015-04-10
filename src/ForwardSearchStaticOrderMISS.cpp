#include "ForwardSearchStaticOrderMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"
#include "Tools.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <climits>

////#define MAX_TWO_NEIGHBORHOOD
////#define MIN_COMPONENT
////#define MIN_SUBPROBLEM

using namespace std;

ForwardSearchStaticOrderMISS::ForwardSearchStaticOrderMISS(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
: TesterMISQ(vAdjacencyMatrix, vAdjacencyArray)
, onlyConsider(vAdjacencyArray.size())
, vMarkedVertices(vAdjacencyArray.size(), false)
{
    SetName("forward-search-static-order-miss");
    R.reserve(m_AdjacencyArray.size());

    stackP.resize(m_AdjacencyArray.size() + 1);
    stackColors.resize(m_AdjacencyArray.size() + 1);
    stackOrder.resize(m_AdjacencyArray.size() + 1);
    stackEvaluatedHalfVertices.resize(m_AdjacencyArray.size() + 1);

    // don't reserve for 0-th vectors, they get std::move'd by InitialOrdering
    for (int index = 0; index < stackP.size(); ++index) {
        stackP[index].reserve(m_AdjacencyArray.size());
        stackColors[index].reserve(m_AdjacencyArray.size());
        stackOrder[index].reserve(m_AdjacencyArray.size());
    }
}

void ForwardSearchStaticOrderMISS::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMISR(m_AdjacencyArray, P, vColors, m_uMaximumCliqueSize);
////    OrderingTools::InitialOrderingMISR(m_AdjacencyMatrix, P, vColors, m_uMaximumCliqueSize);

////    OrderingTools::InitialOrderingReduction(isolates, vVertexOrder);
    vVertexOrder = P; ////std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyMatrix, true /* non-decreasing */)); //// = P; //?
}

void ForwardSearchStaticOrderMISS::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
    vector<int> &vCliqueVertices(stackClique[depth+1]); vCliqueVertices.clear(); vCliqueVertices.push_back(chosenVertex);
    vector<int> &vRemoved(stackOther[depth+1]); vRemoved.clear();
    vector<Reduction> &vReductions(stackReductions[depth+1]); vReductions.clear();
    isolates.RemoveVertexAndNeighbors(chosenVertex, vRemoved, vReductions);

#ifdef NO_ISOLATES_P_LEFT_10
    vector<int> const &vColors(stackColors[depth]);
    bool bRemoveIsolates((vColors.size() < 10) || (vColors[vColors.size()-10] + R.size() + isolates.GetFoldedVertexCount() <= m_uMaximumCliqueSize));
#else
    bool bRemoveIsolates(true);
#endif
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
            isolates.RemoveAllIsolates(0/*unused*/, vCliqueVertices, vRemoved, vReductions /* unused */, false /* only consider updated vertices */);
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

void ForwardSearchStaticOrderMISS::RunRecursive(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors)
{
    nodeCount++;
    vector<int> &vNewP(stackP[R.size()+1]);
    vector<int> &vNewColors(stackColors[R.size()+1]);
    vector<int> &vNewVertexOrder(stackOrder[R.size()+1]);

    stackEvaluatedHalfVertices[depth + 1] = true;

////    vector<int> X;

    size_t const uOriginalPSize(P.size());

////    if (nodeCount%10000 == 0) {
    if (nodeCount%10000 == 0) {
        if (!m_bQuiet) {
            cout << "Evaluated " << nodeCount << " nodes. " << Tools::GetTimeInSeconds(clock() - startTime) << endl;
            PrintState();
        }
        if (m_TimeOut > 0 && (clock() - m_StartTime > m_TimeOut)) {
            m_bTimedOut = true;
            return;
        }
    }

    vector<pair<int,int>> vAddedEdgesUnused;

    bool firstIteration(true);
    clock_t start(clock());

    bool selectMinComponent(false);
    bool checkTime(false); ////isolates.GetInGraph().Size() > 700);

    size_t uIndexIntoP(0);

    vector<Reduction> vRemovedVerticesFromP;
    for (int const p : P) {
        isolates.RemoveVertex(p, vRemovedVerticesFromP);
    }

    size_t sizeOfOriginalP(P.size());

    P.clear();

    vector<int> const vOriginalVertexOrder(vVertexOrder);

    while (P.size() != sizeOfOriginalP && !m_bTimedOut) {

        int const nextVertex(vRemovedVerticesFromP[uIndexIntoP].GetVertex());
        P.push_back(nextVertex);
        vector<Reduction> vVertexReduction;
        vVertexReduction.push_back(vRemovedVerticesFromP[uIndexIntoP]);
        isolates.ReplaceAllRemoved(vVertexReduction);

        vVertexOrder.resize(P.size());
        size_t vertexOrderIndex(0);
        for (int const orderedVertex : vOriginalVertexOrder) {
            if (isolates.GetInGraph().Contains(orderedVertex)) {
                vVertexOrder[vertexOrderIndex++] = orderedVertex;
            }
        }

        uIndexIntoP++;

        firstIteration = false;
        if (depth < 3) {
            if (!m_bQuiet) {
                size_t numLeft = P.size();
                for (; numLeft > 0; --numLeft) {
                    if (R.size() + vColors[numLeft-1] + isolates.GetFoldedVertexCount() <= m_uMaximumCliqueSize) { break; }
                }

                numLeft = P.size() - numLeft;

                cout << depth << ": Only " << P.size() << "(real=" << numLeft << ") more vertices to go! " << Tools::GetTimeInSeconds(clock() - startTime) << endl;
            }
        }

            // TODO/DS: implement pruning...
////        int const largestColor(vColors.back());
////        if (R.size() + largestColor + isolates.GetFoldedVertexCount() <= m_uMaximumCliqueSize /*&& !selectMinNeighbors && !selectMinComponent*/) {
////////            cout << __LINE__ << endl;
////////            cout << "P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////////            cout << "Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
////            ForwardSearchStaticOrderMISS::ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
////////        cout << __LINE__ << endl;
////////        cout << depth << ": P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////////        cout << depth << ": Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
////            P.clear();
////            return;
////        }

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

            depth++;
            RunRecursive(vNewP, vNewVertexOrder, cliques, vNewColors);
            depth--;
        } else if (R.size() + isolates.GetFoldedVertexCount() > m_uMaximumCliqueSize) {
            cliques.back().clear();
            cliques.back().insert(cliques.back().end(), R.begin(), R.end());
            ExecuteCallBacks(cliques.back());
            m_uMaximumCliqueSize = R.size() + isolates.GetFoldedVertexCount();
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

////        X.push_back(nextVertex);

////        if (R.size() > m_uMaximumCliqueSize && bPIsEmpty && P.empty()) {
////            cout << "ERROR!" << endl << flush;
////        }

        if (!bPIsEmpty && P.empty()) {
            if (R.size() + isolates.GetFoldedVertexCount() > m_uMaximumCliqueSize) {
                cliques.back().clear();
                cliques.back().insert(cliques.back().end(), R.begin(), R.end());
                ExecuteCallBacks(cliques.back());
                m_uMaximumCliqueSize = R.size() + isolates.GetFoldedVertexCount();
                timeToLargestClique = clock() - startTime;
            }
        }
    }

////            cout << __LINE__ << endl;
////            cout << "P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////            cout << "Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
    ForwardSearchStaticOrderMISS::ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
////            cout << __LINE__ << endl;
////            cout << "P        " << ((find(P.begin(), P.end(), 201) != P.end()) ? "contains " : "does not contain ") << 201 << endl;
////            cout << "Isolates " << (isolates.GetInGraph().Contains(201) ? "contains " : "does not contain ") << 201 << endl;
    P.clear();

    vNewColors.clear();
    vNewP.clear();
}

void ForwardSearchStaticOrderMISS::ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex)
{
////        stackX[depth+1].push_back(chosenVertex);
        std::vector<int> &vCliqueVertices(stackClique[depth+1]);
        std::vector<int> &vRemoved(stackOther[depth+1]);
        std::vector<Reduction> &vReductions(stackReductions[depth+1]);
////        std::vector<Reduction> &vPersistentReductions(stackPersistentReductions[depth+1]);
        // return R back to the state it was at the beginning of the loop.
        // remove reduction clique vertices from R
        for (int const vertex : vCliqueVertices)
            R.pop_back();

////        // one more for first vertex added to clique
////        R.pop_back();

        // put back removed vertices.
////        cout << "Putting back into graph: " << endl;
////        for (int const cliqueVertex : vCliqueVertices) {
////            cout << cliqueVertex << " ";
////        }
////        cout << " | ";
////        for (int const otherVertex : vRemoved) {
////            cout << otherVertex << " ";
////        }
////        cout << endl;
////        isolates.ReplaceAllRemoved(vCliqueVertices);
////        isolates.ReplaceAllRemoved(vRemoved);
        isolates.ReplaceAllRemoved(vReductions);

        vCliqueVertices.clear();
        vRemoved.clear();
        vReductions.clear();

////        if (chosenVertex != -1) isolates.RemoveVertex(chosenVertex, vPersistentReductions);

        // remove vertices
////        vector<int> vTempCliqueVertices;
////        vector<int> vTempRemovedVertices;
////        vector<pair<int,int>> vAddedEdgesUnused;

////        double const density(isolates.GetDensity());
////        size_t const maxDegree(isolates.GetMaxDegree());

}

void ForwardSearchStaticOrderMISS::ProcessOrderBeforeReturn(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors)
{
////    vector<int> &vCliqueVerticesToReplace(stackPersistentClique[depth+1]);
////    vector<int> &vRemovedVerticesToReplace(stackPersistentOther[depth+1]);
////    std::vector<Reduction> &vPersistentReductions(stackPersistentReductions[depth+1]);
////
////////    isolates.ReplaceAllRemoved(vCliqueVerticesToReplace);
////////    isolates.ReplaceAllRemoved(vRemovedVerticesToReplace);
////    isolates.ReplaceAllRemoved(vPersistentReductions);
////
////    for (int const cliqueVertex : vCliqueVerticesToReplace) {
////        R.pop_back();
////    }
////
////    vCliqueVerticesToReplace.clear();
////    vRemovedVerticesToReplace.clear();
////    vPersistentReductions.clear();
////    vFoldedVertexCounts[depth+1] = 0;
}
