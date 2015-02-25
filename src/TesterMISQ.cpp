#include "TesterMISQ.h"
#include "OrderingTools.h"
#include "GraphTools.h"
#include "Tools.h"

#include <cmath>
#include <iostream>

using namespace std;

TesterMISQ::TesterMISQ(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
: LightWeightMISQ(vAdjacencyMatrix)
////, m_AdjacencyMatrix(vAdjacencyMatrix)
, m_AdjacencyArray(vAdjacencyArray)
////, coloringStrategy(m_AdjacencyMatrix)
////, m_uMaximumCliqueSize(0)
////, stackP(vAdjacencyMatrix.size())
////, stackColors(vAdjacencyMatrix.size())
////, stackOrder(vAdjacencyMatrix.size())
////, stackX(vAdjacencyMatrix.size() + 1)
, stackClique(vAdjacencyMatrix.size() + 1)
, stackOther(vAdjacencyMatrix.size() + 1)
, stackPersistentClique(vAdjacencyMatrix.size() + 1)
, stackPersistentOther(vAdjacencyMatrix.size() + 1)
////, nodeCount(0)
////, depth(-1)
, isolates(vAdjacencyArray)
////, startTime(clock())
////, m_bInvert(0)
{
    SetName("reduction-misq");
    stackEvaluatedHalfVertices.resize(vAdjacencyArray.size(), false);
}

////void TesterMISQ::SetInvert(bool const invert)
////{
////    m_bInvert = invert;
////}

void TesterMISQ::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
////    P = std::move(GraphTools::OrderVerticesByDegree(isolates.GetInGraph(), isolates.Neighbors(), true /* non-decreasing*/));
    P = std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyArray, true /* non-decreasing*/));

    size_t maxDegree(0);
#ifdef SPARSE
    for (SparseArraySet const &neighborSet : isolates.Neighbors()) {
#else
    for (ArraySet const &neighborSet : isolates.Neighbors()) {
#endif //SPARSE
        maxDegree = max(maxDegree, neighborSet.Size());
    }

    vColors.reserve(P.size());
    vColors.clear();
    for (int degree = 1; degree <= maxDegree; degree++ ) {
        vColors.push_back(degree);
    }

    vColors.resize(P.size(), maxDegree + 1);

    vVertexOrder = P;
}

////void Contains(vector<int> const &vVertices, int const vertex, int const lineNo)
////{
////    if (find(vVertices.begin(), vVertices.end(), vertex) != vVertices.end()) {
////        cout << lineNo << ": vector contains " << vertex << endl << flush;
////    }
////}

size_t TesterMISQ::ComputeConnectedComponents(vector<vector<int>> &vComponents)
{
    ArraySet remaining = isolates.GetInGraph();

    Set currentSearch;
    Set evaluated;

    size_t componentCount(0);
    vComponents.clear();

    if (!remaining.Empty()) {
        int const startVertex = *remaining.begin();
        currentSearch.Insert(startVertex);
        remaining.Remove(startVertex);
        componentCount++;
        vComponents.resize(componentCount);
    }

    while (!remaining.Empty() && !currentSearch.Empty()) {
        int const nextVertex(*currentSearch.begin());
        evaluated.Insert(nextVertex);
        vComponents[componentCount - 1].push_back(nextVertex);
        currentSearch.Remove(nextVertex);
        remaining.Remove(nextVertex);
        for (int const neighbor : isolates.Neighbors()[nextVertex]) {
            if (!evaluated.Contains(neighbor)) {
                currentSearch.Insert(neighbor);
            }
        }

        if (currentSearch.Empty() && !remaining.Empty()) {
            int const startVertex = *remaining.begin();
            currentSearch.Insert(startVertex);
            remaining.Remove(startVertex);
            componentCount++;
            vComponents.resize(componentCount);
        }
    }

    return componentCount;
}


void TesterMISQ::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
    vector<int> &vCliqueVertices(stackClique[depth+1]); vCliqueVertices.clear(); vCliqueVertices.push_back(chosenVertex);
    vector<int> &vRemoved(stackOther[depth+1]); vRemoved.clear();
    vector<pair<int,int>> vAddedEdgesUnused;

////    vector<int> const &vColors(stackColors[depth]);
////

////    bool bRemoveIsolates(false);
////    size_t index = vColors.size();
////    for (; index > 0; --index) {
////        if (R.size() + vColors[index] <= m_uMaximumCliqueSize) { bRemoveIsolates = (P.size() - index < 10); break; }
////    }

////    bool const &bRemoveIsolates(stackEvaluatedHalfVertices[depth+1]);
    isolates.RemoveVertexAndNeighbors(chosenVertex, vRemoved);
////    cout << "Size of clique before reduction: " << R.size() << endl;

////    double const density(isolates.GetDensity());
////    size_t const maxDegree(isolates.GetMaxDegree());

////    if (bRemoveIsolates)
////        isolates.RemoveAllIsolates(0/*unused*/, vCliqueVertices, vRemoved, vAddedEdgesUnused /* unused */, false /* only consider updated vertices */, onlyConsider);
        isolates.RemoveAllIsolates(0/*unused*/, vCliqueVertices, vRemoved, vAddedEdgesUnused /* unused */, false /* only consider updated vertices */);

////    cout << __LINE__ << ": density=" << density << ", max-degree=" << maxDegree << ", clique-vertices=" << vCliqueVertices.size() << ", other-removed=" << vRemoved.size()  << ", percent-total-removed=" << (vCliqueVertices.size() + vRemoved.size())/static_cast<double>(P.size())*100 << "%" << endl;

////    cout << __LINE__ << ": Removed " << vCliqueVertices.size() + vRemoved.size() << " vertices " << endl;
    vNewVertexOrder.resize(P.size());
    size_t uNewIndex(0);
    for (int const candidate : P) {
////    if (!m_bInvert && m_AdjacencyMatrix[chosenVertex][candidate]) vNewVertexOrder[uNewIndex++] = candidate;
        if (isolates.GetInGraph().Contains(candidate)) vNewVertexOrder[uNewIndex++] = candidate;
    }
    vNewVertexOrder.resize(uNewIndex);

    R.insert(R.end(), vCliqueVertices.begin(), vCliqueVertices.end());

////    vector<vector<int>> vComponents;
////    cerr << "size of subgraph             : " << isolates.GetInGraph().Size() << endl << flush;
////    cerr << "# connected components       : " << ComputeConnectedComponents(vComponents) << endl << flush;
////    cerr << "size of connected components : ";
////    cout << "[ ";
////    for (vector<int> const& vComponent : vComponents) {
////        cout << vComponent.size() << " ";
////    }
////    cout << "]" << endl << flush;
}

void TesterMISQ::ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex)
{
////        stackX[depth+1].push_back(chosenVertex);
        std::vector<int> &vCliqueVertices(stackClique[depth+1]);
        std::vector<int> &vRemoved(stackOther[depth+1]);
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
        isolates.ReplaceAllRemoved(vCliqueVertices);
        isolates.ReplaceAllRemoved(vRemoved);

        vCliqueVertices.clear();
        vRemoved.clear();

        if (chosenVertex != -1) isolates.RemoveVertex(chosenVertex);

        // remove vertices
        vector<int> vTempCliqueVertices;
        vector<int> vTempRemovedVertices;
        vector<pair<int,int>> vAddedEdgesUnused;

////        double const density(isolates.GetDensity());
////        size_t const maxDegree(isolates.GetMaxDegree());
#ifndef REMOVE_ISOLATES_BEFORE_ONLY
#ifdef ALWAYS_REMOVE_ISOLATES_AFTER
        bool const bRemoveIsolates(true);
#else
        bool const bRemoveIsolates(depth <= 2); ////stackEvaluatedHalfVertices[depth+1]);
#endif //ALWAYS_REMOVE_ISOLATES_AFTER
////////        cout << "Size of clique before reduction: " << R.size() << endl;
        if (bRemoveIsolates)
            isolates.RemoveAllIsolates(0 /*unused*/,vTempCliqueVertices, vTempRemovedVertices, vAddedEdgesUnused, chosenVertex == -1 /* either consider all (true) or consider only changed vertices (false) */);
#endif //REMOVE_ISOLATES_BEFORE_ONLY
////        cout << __LINE__ << ": Removed " << vTempCliqueVertices.size() + vTempRemovedVertices.size() << " vertices " << endl;
////        cout << __LINE__ << ": density=" << density << ", max-degree=" << maxDegree << ", clique-vertices=" << vCliqueVertices.size() << ", other-removed=" << vRemoved.size()  << ", percent-total-removed=" << (vCliqueVertices.size() + vRemoved.size())/static_cast<double>(P.size())*100 << "%" << endl;
////
        R.insert(R.end(), vTempCliqueVertices.begin(), vTempCliqueVertices.end());

////        Contains(vTempCliqueVertices, 5975, __LINE__);
////        Contains(vTempCliqueVertices, 4202, __LINE__);

        vector<int> &vCliqueVerticesToReplace(stackPersistentClique[depth+1]);
        vector<int> &vRemovedVerticesToReplace(stackPersistentOther[depth+1]);

        vCliqueVerticesToReplace.insert(vCliqueVerticesToReplace.end(), vTempCliqueVertices.begin(), vTempCliqueVertices.end());
        if (chosenVertex != -1) vRemovedVerticesToReplace.push_back(chosenVertex);
        vRemovedVerticesToReplace.insert(vRemovedVerticesToReplace.end(), vTempRemovedVertices.begin(), vTempRemovedVertices.end());

////        if (P.size() != isolates.GetInGraph().Size()) {
////            cout << "Node " << recursionNode << ": Consistancy Error" << endl;
////            cout << "P: " << endl;
////            for (int const vertexInP : P) {
////                cout << vertexInP << "(" << isolates.GetInGraph().Contains(vertexInP) << ") ";
////            }
////            cout << endl;
////        }

        if (true) { //P.size() != isolates.GetInGraph().Size()) { // if false?
////            cout << __LINE__ << ": This should not be triggered!" << endl;
            size_t uNewIndex(0);
            // pull vertices out of P and vColors
            for (size_t index = 0; index < P.size(); ++index) {
                if (isolates.GetInGraph().Contains(P[index])) {
                    P[uNewIndex] = P[index];
                    vColors[uNewIndex] = vColors[index];
                    uNewIndex++;
                }
            }
            P.resize(uNewIndex);
            vColors.resize(uNewIndex);

            // pull vertices out of vVertexOrder
            uNewIndex = 0;
            for (size_t index = 0; index < vVertexOrder.size(); ++index) {
                if (isolates.GetInGraph().Contains(vVertexOrder[index])) {
                    vVertexOrder[uNewIndex++] = vVertexOrder[index];
                }
            }
            vVertexOrder.resize(uNewIndex);
////            P.resize(uNewIndex);
////            vColors.resize(uNewIndex);
////            Color(vVertexOrder/* evaluation order */, P /* color order */, vColors);

            // TODO/DS verify that we don't need to recolor,
            // removing unneeded colors should suffice
#ifndef REMOVE_ISOLATES_BEFORE_ONLY
            if (bRemoveIsolates) {
#ifdef RECOLOR
////                cout << __LINE__ << ": Recoloring..." << endl;
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
////                    cout << __LINE__ << ": Recoloring..." << endl;
                    P = proposedP;
                    vColors = proposedColors;
                }
#else
                Color(vVertexOrder, P, vColors);
#endif //RECOLOR
            }
#endif //REMOVE_ISOLATES_BEFORE_ONLY
        }

}

void TesterMISQ::ProcessOrderBeforeReturn(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors)
{
    vector<int> &vCliqueVerticesToReplace(stackPersistentClique[depth+1]);
    vector<int> &vRemovedVerticesToReplace(stackPersistentOther[depth+1]);

    isolates.ReplaceAllRemoved(vCliqueVerticesToReplace);
    isolates.ReplaceAllRemoved(vRemovedVerticesToReplace);

    for (int const cliqueVertex : vCliqueVerticesToReplace) {
        R.pop_back();
    }

    vCliqueVerticesToReplace.clear();
    vRemovedVerticesToReplace.clear();
}

////void TesterMISQ::PrintState() const
////{
////    cout << "(";
////    for (size_t index = 0; index < stackP.size(); ++index) {
////        cout << stackP[index];
////        if (index+1 != stackP.size()) cout << ", ";
////    }
////    cout << ")" << endl << flush;
////}
