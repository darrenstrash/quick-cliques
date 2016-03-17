#include "LightWeightSparseStaticOrderMCS.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

LightWeightSparseStaticOrderMCS::LightWeightSparseStaticOrderMCS(vector<vector<int>> const &vAdjacencyArray)
: LightWeightSparseMCQ(vAdjacencyArray)
{
    SetName("static-order-mcs");
}

void LightWeightSparseStaticOrderMCS::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMCR(m_AdjacencyArray, P, vColors, m_uMaximumCliqueSize);
    vVertexOrder = P; //// = std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyMatrix, false /* non-increasing */)); //// = P; //?
}

void LightWeightSparseStaticOrderMCS::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
    vNewVertexOrder.resize(P.size());

    for (int const neighbor : m_AdjacencyArray[chosenVertex]) {
        vMarkedVertices[neighbor] = true;
    }

    {
        size_t uNewIndex(0);
        for (int const candidate : vVertexOrder) {
            if (vMarkedVertices[candidate]) vNewVertexOrder[uNewIndex++] = candidate;
        }
        vNewVertexOrder.resize(uNewIndex);
    }

    for (int const neighbor : m_AdjacencyArray[chosenVertex]) {
        vMarkedVertices[neighbor] = false;
    }

    R.push_back(chosenVertex);
}

void LightWeightSparseStaticOrderMCS::ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex)
{
    if (chosenVertex == -1) return;
#if 1
    vVertexOrder.erase(find(vVertexOrder.begin(), vVertexOrder.end(), chosenVertex));
#else
    // try searching from end, might be faster in general...
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
#endif // 0
    R.pop_back();
}
