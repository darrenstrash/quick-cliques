#include "ComparisonMISQ.h"
#include "OrderingTools.h"
#include "GraphTools.h"
#include "Tools.h"

#include <cmath>
#include <iostream>

using namespace std;

#define REMOVE_INITIAL_ISOLATES

ComparisonMISQ::ComparisonMISQ(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
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
, algorithm1(vAdjacencyMatrix, vAdjacencyArray)
, algorithm2(vAdjacencyMatrix, vAdjacencyArray)
{
    SetName("reduction-misq");
    stackEvaluatedHalfVertices.resize(vAdjacencyArray.size(), false);

    algorithm1.SetQuiet(true);
    algorithm2.SetQuiet(true);
}

////void ComparisonMISQ::SetInvert(bool const invert)
////{
////    m_bInvert = invert;
////}

void ComparisonMISQ::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
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

void ComparisonMISQ::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
    vector<int> &vCliqueVertices(stackClique[depth+1]); vCliqueVertices.clear(); vCliqueVertices.push_back(chosenVertex);
    vector<int> &vRemoved(stackOther[depth+1]); vRemoved.clear();
    vector<pair<int,int>> vAddedEdgesUnused;

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

    // test the two algorithms
        if (!vNewVertexOrder.empty()) {
            vector<int> vNewP;
            vector<int> vNewColors;

            vNewP.resize(vNewVertexOrder.size());
            vNewColors.resize(vNewVertexOrder.size());
            vector<int> testColors = vNewColors;
            vector<int> testOrder  = vNewVertexOrder;
            vector<int> testP      = vNewP;
            list<list<int>> testCliques(1);
            Color(vNewVertexOrder/* evaluation order */, vNewP /* color order */, vNewColors);
            algorithm1.SetR(R); algorithm1.SetMaximumCliqueSize(m_uMaximumCliqueSize);

            clock_t start(clock());
            algorithm1.SetIsolates(isolates);
            algorithm1.RunRecursive(testP, testOrder, testCliques, testColors);
            cout << "Algorithm1 finished in " << (double)(clock()-start)/(double)(CLOCKS_PER_SEC) << " seconds" << endl;

            testColors = vNewColors;
            testOrder  = vNewVertexOrder;
            testP      = vNewP;

            start = clock();
            algorithm2.SetIsolates(isolates);
            algorithm2.SetR(R); algorithm2.SetMaximumCliqueSize(m_uMaximumCliqueSize);
            algorithm2.RunRecursive(testP, testOrder, testCliques, testColors);
            cout << "Algorithm2 finished in " << (double)(clock()-start)/(double)(CLOCKS_PER_SEC) << " seconds" << endl;
        }

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

void ComparisonMISQ::ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex)
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
        bool const bRemoveIsolates(depth <= 2); ////stackEvaluatedHalfVertices[depth+1]);
////////        cout << "Size of clique before reduction: " << R.size() << endl;
        if (bRemoveIsolates)
            isolates.RemoveAllIsolates(0 /*unused*/,vTempCliqueVertices, vTempRemovedVertices, vAddedEdgesUnused, chosenVertex == -1 /* either consider all (true) or consider only changed vertices (false) */);
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

            if (bRemoveIsolates)
                Color(vVertexOrder, P, vColors);
        }

}

void ComparisonMISQ::ProcessOrderBeforeReturn(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors)
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

////void ComparisonMISQ::PrintState() const
////{
////    cout << "(";
////    for (size_t index = 0; index < stackP.size(); ++index) {
////        cout << stackP[index];
////        if (index+1 != stackP.size()) cout << ", ";
////    }
////    cout << ")" << endl << flush;
////}
