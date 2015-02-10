#include "LightWeightReductionDominationMISQ.h"
#include "OrderingTools.h"
#include "GraphTools.h"
#include "Tools.h"
#include "Reducer.h"

#include <cmath>
#include <iostream>

using namespace std;

#define REMOVE_INITIAL_ISOLATES

LightWeightReductionDominationMISQ::LightWeightReductionDominationMISQ(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
: LightWeightMISQ(vAdjacencyMatrix)
////, m_AdjacencyMatrix(vAdjacencyMatrix)
, m_AdjacencyArray(vAdjacencyArray)
////, coloringStrategy(m_AdjacencyMatrix)
////, m_uMaximumCliqueSize(0)
////, stackP(vAdjacencyMatrix.size())
////, stackColors(vAdjacencyMatrix.size())
////, stackOrder(vAdjacencyMatrix.size())
, stackX(vAdjacencyMatrix.size() + 1)
, stackClique(vAdjacencyMatrix.size() + 1)
, stackOther(vAdjacencyMatrix.size() + 1)
, stackPersistentClique(vAdjacencyMatrix.size() + 1)
, stackPersistentOther(vAdjacencyMatrix.size() + 1)
////, nodeCount(0)
////, depth(-1)
////, isolates(vAdjacencyArray)
, reducer(vAdjacencyArray)
////, startTime(clock())
////, m_bInvert(0)
{
    SetName("reduction-domination-misq");
}

////void LightWeightReductionDominationMISQ::SetInvert(bool const invert)
////{
////    m_bInvert = invert;
////}

void LightWeightReductionDominationMISQ::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
////    P = std::move(GraphTools::OrderVerticesByDegree(isolates.GetInGraph(), isolates.Neighbors(), true /* non-decreasing*/));
    P = std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyArray, true /* non-decreasing*/));

    size_t maxDegree(0);
////#ifdef SPARSE
////    for (SparseArraySet const &neighborSet : isolates.Neighbors()) {
////#else
////    for (ArraySet const &neighborSet : isolates.Neighbors()) {
////#endif //SPARSE
    for (vector<int> const &vNeighbors : m_AdjacencyArray) {
        maxDegree = max(maxDegree, vNeighbors.size());
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

void LightWeightReductionDominationMISQ::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
    vector<int> &vCliqueVertices(stackClique[depth+1]); vCliqueVertices.clear(); vCliqueVertices.push_back(chosenVertex);
    vector<int> &vRemoved(stackOther[depth+1]); vRemoved.clear();
////    isolates.RemoveVertexAndNeighbors(chosenVertex, vRemoved);
////    isolates.RemoveAllIsolates(0/*unused*/, vCliqueVertices, vRemoved, vAddedEdgesUnused /* unused */, false /* only consider updated vertices */);
    reducer.RemoveVertexAndNeighbors(chosenVertex, vRemoved);

    vector<int> const &oldX(stackX[depth]);
    vector<int> &newX(stackX[depth+1]); newX.clear();

////    if (depth == 0) {
////    cout << depth << ": X: Old       :";
////    for (int const x : oldX) {
////        cout << x << " ";
////    }
////    cout << endl;
////    }

#if 1
    vector<bool> vMarkedVertices(m_AdjacencyArray.size(), false);
    for (int const neighbor : m_AdjacencyArray[chosenVertex]) {
        vMarkedVertices[neighbor] = true;
    }

    for (int const vertexInX : oldX) {
        if (!vMarkedVertices[vertexInX]) newX.push_back(vertexInX);
    }

    for (int const neighbor : m_AdjacencyArray[chosenVertex]) {
        vMarkedVertices[neighbor] = false;
    }
#endif // 0

////    if (depth == 0) {
////    cout << depth << ": X: New       :";
////    for (int const x : newX) {
////        cout << x << " ";
////    }
////    cout << endl;
////    }

////    cout << "Size of clique before reduction: " << R.size() << endl;

////    cout << depth << ": ";
    reducer.Reduce(newX, vCliqueVertices, vRemoved);
////    cout << __LINE__ << ": Removed " << vCliqueVertices.size() + vRemoved.size() << " vertices " << endl;

    vNewVertexOrder.resize(P.size());
    size_t uNewIndex(0);
    for (int const candidate : P) {
////    if (!m_bInvert && m_AdjacencyMatrix[chosenVertex][candidate]) vNewVertexOrder[uNewIndex++] = candidate;
////        if (isolates.GetInGraph().Contains(candidate)) vNewVertexOrder[uNewIndex++] = candidate;
        if (reducer.InRemainingGraph(candidate)) vNewVertexOrder[uNewIndex++] = candidate;
    }
    vNewVertexOrder.resize(uNewIndex);

    R.insert(R.end(), vCliqueVertices.begin(), vCliqueVertices.end());
}

void LightWeightReductionDominationMISQ::ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex)
{
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
    ////        isolates.ReplaceAllRemoved(vCliqueVertices);
    ////        isolates.ReplaceAllRemoved(vRemoved);
    reducer.ReplaceVertices(vCliqueVertices);
    reducer.ReplaceVertices(vRemoved);

    vCliqueVertices.clear();
    vRemoved.clear();

    ////        if (chosenVertex != -1) isolates.RemoveVertex(chosenVertex);

    // remove vertices
    vector<int> vTempCliqueVertices;
    vector<int> vTempRemovedVertices;

    ////        isolates.RemoveAllIsolates(0 /*unused*/,vTempCliqueVertices, vTempRemovedVertices, vAddedEdgesUnused, chosenVertex == -1 /* consider all or only changed vertices */);
    if (chosenVertex == -1) {
        reducer.InitialReduce(vTempCliqueVertices);
////        cout << __LINE__ << ": Removed " << vTempCliqueVertices.size() << " + other vertices " << endl;
    } else {
        vector<int> &X(stackX[depth]);
        X.push_back(chosenVertex);
        reducer.RemoveVertex(chosenVertex);
////        cout << "Size of clique before reduction: " << R.size() << endl;
////        cout << depth << " (persistent) : ";
        reducer.Reduce(X, vTempCliqueVertices, vTempRemovedVertices);
////        cout << __LINE__ << ": Removed " << vTempCliqueVertices.size() + vTempRemovedVertices.size() << " vertices " << endl;
////        if (depth == 0) {
////        cout << depth << ": P: AfterRemoval of " << chosenVertex << ":";
////        for (int const p : P) {
////            if (p != chosenVertex) cout << p << " ";
////        }
////        cout << endl;
////
////        cout << depth << ": X: AfterRemoval of " << chosenVertex << ":";
////        for (int const x : X) {
////            cout << x << " ";
////        }
////        cout << endl;
////        }
    }

    R.insert(R.end(), vTempCliqueVertices.begin(), vTempCliqueVertices.end());

////    Contains(vTempCliqueVertices, 5975, __LINE__);
////    Contains(vTempCliqueVertices, 4202, __LINE__);

    vector<int> &vCliqueVerticesToReplace(stackPersistentClique[depth+1]);
    vector<int> &vRemovedVerticesToReplace(stackPersistentOther[depth+1]);

    vCliqueVerticesToReplace.insert(vCliqueVerticesToReplace.end(), vTempCliqueVertices.begin(), vTempCliqueVertices.end());
    if (chosenVertex != -1) vRemovedVerticesToReplace.push_back(chosenVertex);
    vRemovedVerticesToReplace.insert(vRemovedVerticesToReplace.end(), vTempRemovedVertices.begin(), vTempRemovedVertices.end());

////    if (P.size() != isolates.GetInGraph().Size()) {
////        cout << "Node " << recursionNode << ": Consistancy Error" << endl;
////        cout << "P: " << endl;
////        for (int const vertexInP : P) {
////            cout << vertexInP << "(" << isolates.GetInGraph().Contains(vertexInP) << ") ";
////        }
////        cout << endl;
////    }

////    if (P.size() != isolates.GetInGraph().Size()) { // if false?
    if (P.size() != reducer.RemainingGraphSize()) { // if false?
        ////            cout << __LINE__ << ": This should not be triggered!" << endl;
        size_t uNewIndex(0);
        // pull vertices out of P and vColors
        for (size_t index = 0; index < P.size(); ++index) {
            ////                if (isolates.GetInGraph().Contains(P[index])) {
            if (reducer.InRemainingGraph(P[index])) {
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
            ////                if (isolates.GetInGraph().Contains(vVertexOrder[index])) {
            if (reducer.InRemainingGraph(vVertexOrder[index])) {
                vVertexOrder[uNewIndex++] = vVertexOrder[index];
            }
        }
        vVertexOrder.resize(uNewIndex);
////        P.resize(uNewIndex);
////        vColors.resize(uNewIndex);
////        Color(vVertexOrder/* evaluation order */, P /* color order */, vColors);

    }

    Color(vVertexOrder, P, vColors);
}

void LightWeightReductionDominationMISQ::ProcessOrderBeforeReturn(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors)
{
    vector<int> &vCliqueVerticesToReplace(stackPersistentClique[depth+1]);
    vector<int> &vRemovedVerticesToReplace(stackPersistentOther[depth+1]);

////    isolates.ReplaceAllRemoved(vCliqueVerticesToReplace);
////    isolates.ReplaceAllRemoved(vRemovedVerticesToReplace);
    reducer.ReplaceVertices(vCliqueVerticesToReplace);
    reducer.ReplaceVertices(vRemovedVerticesToReplace);

    for (int const cliqueVertex : vCliqueVerticesToReplace) {
        R.pop_back();
    }

    vCliqueVerticesToReplace.clear();
    vRemovedVerticesToReplace.clear();
}

////void LightWeightReductionDominationMISQ::PrintState() const
////{
////    cout << "(";
////    for (size_t index = 0; index < stackP.size(); ++index) {
////        cout << stackP[index];
////        if (index+1 != stackP.size()) cout << ", ";
////    }
////    cout << ")" << endl << flush;
////}
