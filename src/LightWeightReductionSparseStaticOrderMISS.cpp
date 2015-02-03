#include "LightWeightReductionSparseStaticOrderMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

LightWeightReductionSparseStaticOrderMISS::LightWeightReductionSparseStaticOrderMISS(vector<vector<int>> const &vAdjacencyArray)
: LightWeightReductionSparseMISQ(vAdjacencyArray)
{
    SetName("reduction-sparse-static-order-miss");
}

void LightWeightReductionSparseStaticOrderMISS::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMISR(m_AdjacencyArray, P, vColors, m_uMaximumCliqueSize);
    vVertexOrder = std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyArray, true /* non-decreasing */)); //// = P; //?
}

void LightWeightReductionSparseStaticOrderMISS::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
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

////void LightWeightReductionSparseStaticOrderMISS::ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex)
////{
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
////}
