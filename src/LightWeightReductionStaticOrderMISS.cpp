#include "LightWeightReductionStaticOrderMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

LightWeightReductionStaticOrderMISS::LightWeightReductionStaticOrderMISS(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
: LightWeightReductionMISQ(vAdjacencyMatrix, vAdjacencyArray)
, onlyConsider(vAdjacencyMatrix.size())
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
