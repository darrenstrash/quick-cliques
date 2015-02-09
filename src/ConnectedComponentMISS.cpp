#include "ConnectedComponentMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"
#include "Tools.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

ConnectedComponentMISS::ConnectedComponentMISS(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
: LightWeightReductionMISQ(vAdjacencyMatrix, vAdjacencyArray)
, m_vSubgraphClique()
{
    SetName("connected-component-miss");
}

ConnectedComponentMISS::~ConnectedComponentMISS()
{
////    cout << "SubgraphClique has size " << m_vSubgraphClique.size() << endl << flush;
}

void ConnectedComponentMISS::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMISR(m_AdjacencyMatrix, P, vColors, m_uMaximumCliqueSize);
    vVertexOrder = P; ////std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyMatrix, true /* non-decreasing */)); //// = P; //?
}

void ConnectedComponentMISS::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
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

void ConnectedComponentMISS::RunRecursive(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors)
{
    nodeCount++;
    vector<int> &vNewP(stackP[depth+1]);
    vector<int> &vNewColors(stackColors[depth+1]);
    vector<int> &vNewVertexOrder(stackOrder[depth+1]);

////    cout << depth << ": P: ";
////    for (int const p : P) {
////        cout << p << " ";
////    }
////    cout << endl;

#if 1
    if (P.size() > 200) {
        vector<vector<int>> vComponents;
        GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());
        size_t uRemainingVertices(P.size());
        size_t addedToR(0);
        if (vComponents.size() > 1) {

            list<int> bestClique(R.begin(), R.end());

            vector<int> const savedR = std::move(R); //// m_vSubgraphClique);
            size_t            savedCliqueSize(m_uMaximumCliqueSize);
            list<int>         savedClique(cliques.back().begin(), cliques.back().end());

////            cerr << "# connected components       : " << vComponents.size() << endl << flush;
////            cerr << "size of connected components : ";
////            cout << "[ ";
////            for (vector<int> const& vComponent : vComponents) {
////                cout << vComponent.size() << " ";
////            }
////            cout << "]" << endl << flush;
            for (vector<int> const &vComponent : vComponents) {
                m_uMaximumCliqueSize = 0;
                R.clear();
                cliques.back().clear();
                isolates.SetConnectedComponent(vComponent);

                size_t uNewSize(0);
                vNewVertexOrder.resize(vVertexOrder.size());
                for (size_t index = 0; index < vVertexOrder.size(); ++index) {
                    if (isolates.GetInGraph().Contains(vVertexOrder[index])) vNewVertexOrder[uNewSize++] = vVertexOrder[index];
                }
                vNewVertexOrder.resize(uNewSize);
                vNewP.resize(uNewSize);
                vNewColors.resize(uNewSize);
////                cout << depth << ": vertexOrder size = " << vNewVertexOrder.size() << endl;

                // TODO/DS: add delta value to the colors so that we get better pruning.
                // delta = P.Size() - vComponent.size()
                Color(vNewVertexOrder, vNewP, vNewColors); //, delta) //// like this

////                for (size_t index = 0; index < P.size(); ++index) {
////                    if (isolates.GetInGraph().Contains(vVertexOrder[index])) {
////                        vNewP[uNewSize] = P[index];
////                        vNewColors[uNewSize] = vColors[index];
////                        uNewSize++;
////                    }
////                }
////                vNewP.resize(uNewSize);
////                vNewColors.resize(uNewSize);

                // TODO/DS: need a flag to keep track of largest clique found found in subcall, so we can add it to R.
                // otherwise, we will miss out on larger cliques entirely!

                //save current best clique and startingDepth; // like this
////                m_uSubgraphDepth = depth;
////                m_vSubgraphClique = R;

                //bestSubClique = empty vector; cliqueStartingdepth = depth; // update as recursion updates
                depth++;
                RunRecursive(vNewP, vNewVertexOrder, cliques, vColors);
                depth--;

////                cout << depth << ": Found clique of size " << cliques.back().size() << endl;
////                cout << depth << ": or              size " << m_uMaximumCliqueSize  << endl;

                bestClique.insert(bestClique.end(), cliques.back().begin(), cliques.back().end());

                if (bestClique.size() > savedCliqueSize) {
                    savedCliqueSize = bestClique.size();
                    savedClique.clear();
                    savedClique.insert(savedClique.end(), bestClique.begin(), bestClique.end());
                    ExecuteCallBacks(cliques.back());
                }

////                // TODO/DS: add new clique to R.
////                for (int const vertex : m_vSubgraphClique) {
////                    R.push_back(vertex);
////                } // like this

////                restore bestsubclique and startingdepth // restore, like this

                // TODO/DS: if we can prune, break;
            }

            // TODO/DS: remove newly found vertices (independent sets) from R.

            // put R, clique size, and graph back the way they were.
            R = savedR;
            cliques.back() = bestClique;
            m_uMaximumCliqueSize = savedCliqueSize;
            isolates.SetConnectedComponent(P); 

////            m_vSubgraphClique = savedClique;

            return;
        }
    }
#endif // 0

    if (nodeCount%10000 == 0) {
        cout << "Evaluated " << nodeCount << " nodes. " << GetTimeInSeconds(clock() - startTime) << endl;
        ////        PrintState();
    }

    while (!P.empty()) {

        if (depth == 0 && !m_bQuiet) {
            cout << "Only " << P.size() << " more vertices to go! " << GetTimeInSeconds(clock() - startTime) << endl;
        }

        int const largestColor(vColors.back());
        if (R.size() + largestColor <= m_uMaximumCliqueSize) {
            ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
            return;
        }

        vColors.pop_back();
        int const nextVertex(P.back()); P.pop_back();
////        cout << depth << ": adding vertex " << nextVertex << endl;
        GetNewOrder(vNewVertexOrder, vVertexOrder, P, nextVertex);

        if (!vNewVertexOrder.empty()) {
            vNewP.resize(vNewVertexOrder.size());
            vNewColors.resize(vNewVertexOrder.size());
            Color(vNewVertexOrder/* evaluation order */, vNewP /* color order */, vNewColors);
            depth++;
            RunRecursive(vNewP, vNewVertexOrder, cliques, vNewColors);
            depth--;
        } else if (R.size() > m_uMaximumCliqueSize) {
            cliques.back().clear();
            cliques.back().insert(cliques.back().end(), R.begin(), R.end());
            ExecuteCallBacks(cliques.back());
            m_uMaximumCliqueSize = R.size();
        }

        bool bPIsEmpty(P.empty());
        ProcessOrderAfterRecursion(vVertexOrder, P, vColors, nextVertex);

////        if (R.size() > m_uMaximumCliqueSize && bPIsEmpty && P.empty()) {
////            cout << "ERROR!" << endl << flush;
////        }
////        if (R.size() > m_vSubgraphClique.size()) {
////            m_vSubgraphClique = R;
////        }

        if (!bPIsEmpty && P.empty()) {
            if (R.size() > m_uMaximumCliqueSize) {
                cliques.back().clear();
                cliques.back().insert(cliques.back().end(), R.begin(), R.end());
                ExecuteCallBacks(cliques.back());
                m_uMaximumCliqueSize = R.size();
            }
        }
    }

    ProcessOrderBeforeReturn(vVertexOrder, P, vColors);

    vNewColors.clear();
    vNewP.clear();
}
