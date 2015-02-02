#include "LightWeightReductionStaticOrderMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

LightWeightReductionStaticOrderMISS::LightWeightReductionStaticOrderMISS(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
: LightWeightReductionMISQ(vAdjacencyMatrix, vAdjacencyArray)
{
    SetName("reduction-static-order-miss");
}

void LightWeightReductionStaticOrderMISS::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMISR(m_AdjacencyMatrix, P, vColors, m_uMaximumCliqueSize);
    vVertexOrder = P; ////std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyMatrix, true /* non-decreasing */)); //// = P; //?
}

void LightWeightReductionStaticOrderMISS::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
    vector<int> &vCliqueVertices(stackClique[depth+1]); vCliqueVertices.clear(); vCliqueVertices.push_back(chosenVertex);
    vector<int> &vRemoved(stackOther[depth+1]); vRemoved.clear();
    vector<pair<int,int>> vAddedEdgesUnused;
    isolates.RemoveVertexAndNeighbors(chosenVertex, vRemoved);
    isolates.RemoveAllIsolates(0/*unused*/, vCliqueVertices, vRemoved, vAddedEdgesUnused /* unused */, false /* only consider updated vertices */);
    vNewVertexOrder.resize(vVertexOrder.size());
    size_t uNewIndex(0);
    for (int const candidate : vVertexOrder) {
////    if (!m_bInvert && m_AdjacencyMatrix[chosenVertex][candidate]) vNewVertexOrder[uNewIndex++] = candidate;
        if (isolates.GetInGraph().Contains(candidate)) vNewVertexOrder[uNewIndex++] = candidate;
    }
    vNewVertexOrder.resize(uNewIndex);

    R.insert(R.end(), vCliqueVertices.begin(), vCliqueVertices.end());
}

////void LightWeightReductionStaticOrderMISS::PostProcessOrder(std::vector<int> &vVertexOrder, int const chosenVertex)
////{
////    if (chosenVertex == -1) return;
////    ////        vVertexOrder.erase(find(vVertexOrder.begin(), vVertexOrder.end(), chosenVertex));
////    // try searching from end, might be faster in general...
////    size_t indexAfterVertex(0);
////    for (indexAfterVertex = vVertexOrder.size(); indexAfterVertex >= 1; indexAfterVertex--) {
////        if (vVertexOrder[indexAfterVertex-1] == chosenVertex) {
////            break;
////        }
////    }
////
////    for (; indexAfterVertex < vVertexOrder.size(); indexAfterVertex++) {
////        vVertexOrder[indexAfterVertex-1] = vVertexOrder[indexAfterVertex];
////    }
////    vVertexOrder.pop_back();
////    R.pop_back();
////}
