#include "ComparisonStaticOrderMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

ComparisonStaticOrderMISS::ComparisonStaticOrderMISS(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
: ComparisonMISQ(vAdjacencyMatrix, vAdjacencyArray)
, onlyConsider(vAdjacencyMatrix.size())
{
    SetName("reduction-static-order-miss");
}

void ComparisonStaticOrderMISS::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMISR(m_AdjacencyMatrix, P, vColors, m_uMaximumCliqueSize);
    vVertexOrder = P; ////std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyMatrix, true /* non-decreasing */)); //// = P; //?
}

void ComparisonStaticOrderMISS::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
    vector<int> &vCliqueVertices(stackClique[depth+1]); vCliqueVertices.clear(); vCliqueVertices.push_back(chosenVertex);
    vector<int> &vRemoved(stackOther[depth+1]); vRemoved.clear();
    vector<pair<int,int>> vAddedEdgesUnused;
    isolates.RemoveVertexAndNeighbors(chosenVertex, vRemoved);

////    vector<int> const &vColors(stackColors[depth]);
////
////    for (size_t index = vColors.size() + 1; index > 0; --index) {
////        if (R.size() + vColors[index-1] <= m_uMaximumCliqueSize) break;
////        onlyConsider.Insert(P[index-1]);
////    }

////    double const density(isolates.GetDensity());
////    size_t const maxDegree(isolates.GetMaxDegree());
////    bool const &bRemoveIsolates(stackEvaluatedHalfVertices[depth + 1]);
////    if (bRemoveIsolates) {
////        if (onlyConsider.Size() > 50) {
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

    // test the two algorithms
    if (!vNewVertexOrder.empty() && depth < 3) {
        vector<int> vNewP;
        vector<int> vNewColors;

        vNewP.resize(vNewVertexOrder.size());
        vNewColors.resize(vNewVertexOrder.size());
        list<list<int>> testCliques(1);
        Color(vNewVertexOrder/* evaluation order */, vNewP /* color order */, vNewColors);
        vector<int> testColors = vNewColors;
        vector<int> testOrder  = vNewVertexOrder;
        vector<int> testP      = vNewP;

        algorithm1.SetR(R); algorithm1.SetMaximumCliqueSize(m_uMaximumCliqueSize);
        size_t index = vNewColors.empty() ? 0 : vNewColors.size();
        for (; index > 0; --index) {
            if (R.size() + vNewColors[index-1] <= m_uMaximumCliqueSize) break;
        }

        if (testColors.empty() || (testColors.back() + R.size() <= m_uMaximumCliqueSize)) return;

        if (index > 0) index--;


        cout << "Depth     = " << depth << endl;
        cout << "P.size    = " << vNewP.size() << endl;
        cout << "P.left    = " << vNewP.size() - index << endl;
        cout << "MaxDegree = " << isolates.GetMaxDegree() << endl;
        cout << "Density   = " << isolates.GetDensity() << endl;
        clock_t start(clock());
        algorithm1.SetIsolates(isolates);
        algorithm1.RunRecursive(testP, testOrder, testCliques, testColors);
        cout << "    CC   Algorithm finished in " << (double)(clock()-start)/(double)(CLOCKS_PER_SEC) << " seconds" << endl << flush;

        testColors = vNewColors;
        testOrder  = vNewVertexOrder;
        testP      = vNewP;

        algorithm2.SetIsolates(isolates);
        clock_t start2 = clock();
        algorithm2.SetR(R); algorithm2.SetMaximumCliqueSize(m_uMaximumCliqueSize);
        algorithm2.RunRecursive(testP, testOrder, testCliques, testColors);
        cout << "    Test Algorithm finished in " << (double)(clock()-start2)/(double)(CLOCKS_PER_SEC) << " seconds" << endl << flush;
    }
}
