#include "LightWeightReductionSparseMISQ.h"
#include "OrderingTools.h"
#include "GraphTools.h"
#include "Tools.h"

#include <cmath>
#include <iostream>

using namespace std;

#define REMOVE_INITIAL_ISOLATES

LightWeightReductionSparseMISQ::LightWeightReductionSparseMISQ(vector<vector<int>> const &vAdjacencyArray)
: MaxSubgraphAlgorithm("reduction-sparse-misq")
, m_AdjacencyArray(vAdjacencyArray)
, coloringStrategy(vAdjacencyArray)
, stackClique(vAdjacencyArray.size() + 1)
, stackOther(vAdjacencyArray.size() + 1)
, stackPersistentClique(vAdjacencyArray.size() + 1)
, stackPersistentOther(vAdjacencyArray.size() + 1)
, reducer(vAdjacencyArray)
////, m_bInvert(0)
{
    R.reserve(m_AdjacencyArray.size());

    stackP.resize(m_AdjacencyArray.size() + 1);
    stackColors.resize(m_AdjacencyArray.size() + 1);
    stackOrder.resize(m_AdjacencyArray.size() + 1);

    stackEvaluatedHalfVertices.resize(vAdjacencyArray.size(), false);
}

////void LightWeightReductionSparseMISQ::SetInvert(bool const invert)
////{
////    m_bInvert = invert;
////}

void LightWeightReductionSparseMISQ::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
////    P = std::move(GraphTools::OrderVerticesByDegree(isolates.GetInGraph(), isolates.Neighbors(), true /* non-decreasing*/));
    P = std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyArray, true /* non-decreasing*/));

    size_t maxDegree(0);
    for (vector<int> const &vNeighbors: m_AdjacencyArray) {
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

void LightWeightReductionSparseMISQ::Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors)
{
    coloringStrategy.Color(m_AdjacencyArray, vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors);
}

void LightWeightReductionSparseMISQ::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
    vector<int> &vCliqueVertices(stackClique[depth+1]); vCliqueVertices.clear(); vCliqueVertices.push_back(chosenVertex);
    vector<int> &vRemoved(stackOther[depth+1]); vRemoved.clear();
    reducer.RemoveVertexAndNeighbors(chosenVertex, vRemoved);
    vector<int> vUnused;

////    bool const &bRemoveIsolates(stackEvaluatedHalfVertices[depth+1]);
////    if (bRemoveIsolates)
    reducer.Reduce(vUnused, vCliqueVertices, vRemoved);
    vNewVertexOrder.resize(P.size());
    size_t uNewIndex(0);
    for (int const candidate : P) {
////    if (!m_bInvert && m_AdjacencyMatrix[chosenVertex][candidate]) vNewVertexOrder[uNewIndex++] = candidate;
        if (reducer.InRemainingGraph(candidate)) vNewVertexOrder[uNewIndex++] = candidate;
    }
    vNewVertexOrder.resize(uNewIndex);

    R.insert(R.end(), vCliqueVertices.begin(), vCliqueVertices.end());
}

void LightWeightReductionSparseMISQ::ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex)
{
        std::vector<int> &vCliqueVertices(stackClique[depth+1]);
        std::vector<int> &vRemoved(stackOther[depth+1]);
        // return R back to the state it was at the beginning of the loop.
        // remove reduction clique vertices from R
        for (int const vertex : vCliqueVertices)
            R.pop_back();

        reducer.ReplaceVertices(vCliqueVertices);
        reducer.ReplaceVertices(vRemoved);

        vCliqueVertices.clear();
        vRemoved.clear();

        if (chosenVertex != -1) reducer.RemoveVertex(chosenVertex);

        // remove vertices
        vector<int> vTempCliqueVertices;
        vector<int> vTempRemovedVertices;
        vector<pair<int,int>> vAddedEdgesUnused;
        vector<int> vUnused;

        bool const &bRemoveIsolates(stackEvaluatedHalfVertices[depth+1]);
        if (chosenVertex == -1) {
            reducer.InitialReduce(vTempCliqueVertices);
        } else {
            if (bRemoveIsolates)
                reducer.Reduce(vUnused, vTempCliqueVertices, vTempRemovedVertices);
        }

        R.insert(R.end(), vTempCliqueVertices.begin(), vTempCliqueVertices.end());

        vector<int> &vCliqueVerticesToReplace(stackPersistentClique[depth+1]);
        vector<int> &vRemovedVerticesToReplace(stackPersistentOther[depth+1]);

        vCliqueVerticesToReplace.insert(vCliqueVerticesToReplace.end(), vTempCliqueVertices.begin(), vTempCliqueVertices.end());
        if (chosenVertex != -1) vRemovedVerticesToReplace.push_back(chosenVertex);
        vRemovedVerticesToReplace.insert(vRemovedVerticesToReplace.end(), vTempRemovedVertices.begin(), vTempRemovedVertices.end());

        if (true) { ////P.size() != reducer.RemainingGraphSize()) { // if false?
////            cout << __LINE__ << ": This should not be triggered!" << endl;
            size_t uNewIndex(0);
            // pull vertices out of P and vColors
            for (size_t index = 0; index < P.size(); ++index) {
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
                if (reducer.InRemainingGraph(vVertexOrder[index])) {
                    vVertexOrder[uNewIndex++] = vVertexOrder[index];
                }
            }
            vVertexOrder.resize(uNewIndex);
////            P.resize(uNewIndex);
            vColors.resize(uNewIndex);
            if (bRemoveIsolates)
                Color(vVertexOrder/* evaluation order */, P /* color order */, vColors);
        }

////        if (bRemoveIsolates)
////        Color(vVertexOrder, P, vColors);
}

void LightWeightReductionSparseMISQ::ProcessOrderBeforeReturn(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors)
{
    vector<int> &vCliqueVerticesToReplace(stackPersistentClique[depth+1]);
    vector<int> &vRemovedVerticesToReplace(stackPersistentOther[depth+1]);

    reducer.ReplaceVertices(vCliqueVerticesToReplace);
    reducer.ReplaceVertices(vRemovedVerticesToReplace);

    for (int const cliqueVertex : vCliqueVerticesToReplace) {
        R.pop_back();
    }

    vCliqueVerticesToReplace.clear();
    vRemovedVerticesToReplace.clear();
}
