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
, m_vStackDelta()
{
    SetName("connected-component-miss");

    m_vStackDelta.resize(vAdjacencyMatrix.size() + 1, 0);
}

ConnectedComponentMISS::~ConnectedComponentMISS()
{
////    cout << "SubgraphClique has size " << m_vSubgraphClique.size() << endl << flush;
}

void ConnectedComponentMISS::Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors)
{
////    vVerticesToReorder = vVertexOrder;
////    coloringStrategy.ColorWithoutReorder(m_AdjacencyMatrix, vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors);
    coloringStrategy.Color(m_AdjacencyMatrix, vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors);
}

void ConnectedComponentMISS::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    OrderingTools::InitialOrderingMISR(m_AdjacencyMatrix, P, vColors, m_uMaximumCliqueSize);
////    OrderingTools::InitialOrderingConnectedComponent(m_AdjacencyMatrix, vVertexOrder, vColors);
    vVertexOrder = P; ////std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyMatrix, true /* non-decreasing */)); //// = P; //?

////    vColors.resize(vVertexOrder.size(), 0);
////    P.resize(vVertexOrder.size(), 0);
////    Color(vVertexOrder, P, vColors); //, delta) //// like this
}

void ConnectedComponentMISS::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex)
{
////    cout << depth << ": GetNewOrder goes from size " << vVertexOrder.size();
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

////    cout << " to size " << vNewVertexOrder.size() << endl;
}

void ConnectedComponentMISS::RunRecursive(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors)
{
    RunRecursive(P, vVertexOrder, cliques, vColors, true /* separate connected components */);
}

void ConnectedComponentMISS::RunRecursive(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors, bool const checkForConnectedComponents)
{
    nodeCount++;
    vector<int> &vNewP(stackP[depth+1]);
    vector<int> &vNewColors(stackColors[depth+1]);
    vector<int> &vNewVertexOrder(stackOrder[depth+1]);
    m_vStackDelta[depth+1] = m_vStackDelta[depth];

////    cout << depth << ": P: ";
////    for (int const p : P) {
////        cout << p << " ";
////    }
////    cout << endl;

#if 1
    if (checkForConnectedComponents && P.size() > 1000) {
        vector<vector<int>> vComponents;
        GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());
        size_t uRemainingVertices(P.size());
        size_t addedToR(0);
        if (vComponents.size() > 1) {

            cout << depth << ": Found connected components before recursion" << endl;

            sort(vComponents.begin(), vComponents.end(), [](vector<int> const &left, vector<int> const &right) { return right < left;});

            list<int> bestClique(R.begin(), R.end());

            vector<int> const savedR(R);

           //save current best clique and startingDepth; // like this
           vector<int> vSavedClique(std::move(m_vSubgraphClique));
           m_vSubgraphClique = R;
           list<int>         savedClique(cliques.back().begin(), cliques.back().end());

////            cerr << "# connected components       : " << vComponents.size() << endl << flush;
////            cerr << "size of connected components : ";
////            cout << "[ ";
////            for (vector<int> const& vComponent : vComponents) {
////                cout << vComponent.size() << " ";
////            }
////            cout << "]" << endl << flush;
            size_t uNumberOfVerticesEvaluated(0);

            vector<int> vDelta(vComponents.size(), 0);

            // since each connected component is isolated, each vertex contributes
            // to an increase in the color of vertices in other components.
            vDelta[vComponents.size()-1] = 0; // the last component only has colors within the component.
            for (size_t index = vComponents.size()-1; index >= 1; --index) {
                vDelta[index-1] = vDelta[index] + vComponents[index].size();
            }

            for (size_t componentIndex = 0; componentIndex < vComponents.size(); ++componentIndex) {
                vector<int> const &vComponent(vComponents[componentIndex]);

                //// computing maximum cliques without knowing the global max clique makes us
                //// traverse many more nodes than are necessary
////                m_uMaximumCliqueSize = 0;
////                R.clear();
////                cliques.back().clear();
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

                // add delta value to the colors so that we get better pruning.
                // delta = P.Size() - vComponent.size()
                Color(vNewVertexOrder, vNewP, vNewColors); //, delta) //// like this

////                for (int &color : vNewColors) {
////                    color += vDelta[componentIndex];
////                }

                m_vStackDelta[depth+1] += vDelta[componentIndex];

////                cout << "Colors (" << vNewColors.size() << "): ";
////                for (int index = max(0,static_cast<int>(vNewColors.size() - 11)); index < static_cast<int>(vNewColors.size()); ++index) {
////                    cout << vNewColors[index] << " ";
////                }
////                cout << endl;

                // TODO/DS: need to recognize when returning because of color pruning.
                // then we can avoid evaluating the rest of the connected components entirely.

                // TODO/DS: need a flag to keep track of largest clique found found in subcall, so we can add it to R.
                // after completing recursion on a component. Otherwise, we will miss out on larger cliques entirely!

                depth++;
                RunRecursive(vNewP, vNewVertexOrder, cliques, vNewColors, false);
                depth--;

                m_vStackDelta[depth+1] -= vDelta[componentIndex];

////                cout << depth << ": Found clique of size " << cliques.back().size() << endl;
////                cout << depth << ": or              size " << m_uMaximumCliqueSize  << endl;

////                bestClique.insert(bestClique.end(), cliques.back().begin(), cliques.back().end());

////                if (bestClique.size() > savedCliqueSize) {
////                    savedCliqueSize = bestClique.size();
////                    savedClique.clear();
////                    savedClique.insert(savedClique.end(), bestClique.begin(), bestClique.end());
////                    ExecuteCallBacks(cliques.back());
////                }

                R = m_vSubgraphClique;

                uNumberOfVerticesEvaluated += vComponent.size();

        if (depth == 0 && !m_bQuiet) {
            cout << "Only " << P.size() - uNumberOfVerticesEvaluated << " more vertices to go! " << GetTimeInSeconds(clock() - startTime) << endl;
            cout << "Current best clique: " << m_uMaximumCliqueSize << endl;
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
            R = std::move(savedR); // could also pop_back() // might be faster.

            // only put back the clique if the newfound clique is smaller than
            // the already-known one.
            if (m_vSubgraphClique.size() < vSavedClique.size())
                m_vSubgraphClique = vSavedClique;
////            cliques.back() = bestClique;
////            m_uMaximumCliqueSize = savedCliqueSize;
            isolates.SetConnectedComponent(P); 

////            m_vSubgraphClique = savedClique;

            return;
        }
    }
#endif // 0

    if (nodeCount%10000 == 0) {
        cout << "Evaluated " << nodeCount << " nodes. " << GetTimeInSeconds(clock() - startTime) << endl;
        PrintState();
    }

    while (!P.empty()) {

////        cout << depth << ": P has " << P.size() << " elements" << endl;

        if (depth == 0 && !m_bQuiet) {
            cout << "Only " << P.size() << " more vertices to go! " << GetTimeInSeconds(clock() - startTime) << endl;
            cout << "Current best clique: " << m_uMaximumCliqueSize << endl;
        }

////        cout << "P.size=" << P.size() << ", vColors.size=" << vColors.size() << endl << flush;
        int largestColorInSubgraph(0);
        for (int const color : vColors) {
            largestColorInSubgraph = max(color, largestColorInSubgraph);
        }
////        int const largestColor(vColors.back() + m_vStackDelta[depth+1]);
        int const largestColor(largestColorInSubgraph + m_vStackDelta[depth+1]);
        if (R.size() + largestColor <= m_uMaximumCliqueSize) {
////            cout << depth << ": Pruned: "; 
////            cout << "Colors (" << vColors.size() << "): ";
////            for (int index = max(0,static_cast<int>(vColors.size() - 11)); index < static_cast<int>(vColors.size()); ++index) {
////                cout << vColors[index] << " ";
////            }
////            cout << endl;
////            cout << R.size() << " + " << largestColor << " <= " << m_uMaximumCliqueSize << endl << flush;
////            cout << "delta=" << m_vStackDelta[depth] << endl << flush;
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

////            cout << "Colors (" << vNewColors.size() << ": ";
////            for (int index = max(0,static_cast<int>(vNewColors.size() - 11)); index < static_cast<int>(vNewColors.size()); ++index) {
////                cout << vNewColors[index] << " ";
////            }
////            cout << endl;
            RunRecursive(vNewP, vNewVertexOrder, cliques, vNewColors, true);
            depth--;
        } else if (R.size() > m_uMaximumCliqueSize) {
            cliques.back().clear();
            cliques.back().insert(cliques.back().end(), R.begin(), R.end());
            ExecuteCallBacks(cliques.back());
            m_uMaximumCliqueSize = R.size();
        }

        bool bPIsEmpty(P.empty());
////        cout << depth << ": P has " << P.size() << " elements before postprocessing" << endl;
        ProcessOrderAfterRecursion(vVertexOrder, P, vColors, nextVertex);
////        cout << depth << ": P has " << P.size() << " elements after  postprocessing" << endl;

////        if (R.size() > m_uMaximumCliqueSize && bPIsEmpty && P.empty()) {
////            cout << "ERROR!" << endl << flush;
////        }
////        if (R.size() > m_vSubgraphClique.size()) {
////            m_vSubgraphClique = R;
////        }

    // if the graph has become disconnected, we can evaluate the connected components separately
#if 1
    if (P.size() > 1000) {
        vector<vector<int>> vComponents;
        GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());
        size_t uRemainingVertices(P.size());
        size_t addedToR(0);
        if (vComponents.size() > 1) {

////            cout << depth << ": Found connected components after  recursion" << endl;

            sort(vComponents.begin(), vComponents.end(), [](vector<int> const &left, vector<int> const &right) { return right < left;});

            list<int> bestClique(R.begin(), R.end());

            vector<int> const savedR(R);

           //save current best clique and startingDepth; // like this
           vector<int> vSavedClique(std::move(m_vSubgraphClique));
           m_vSubgraphClique = R;
           list<int>         savedClique(cliques.back().begin(), cliques.back().end());

////            cerr << "# connected components       : " << vComponents.size() << endl << flush;
////            cerr << "size of connected components : ";
////            cout << "[ ";
////            for (vector<int> const& vComponent : vComponents) {
////                cout << vComponent.size() << " ";
////            }
////            cout << "]" << endl << flush;
            size_t uNumberOfVerticesEvaluated(0);

            vector<int> vDelta(vComponents.size(), 0);

            // since each connected component is isolated, each vertex contributes
            // to an increase in the color of vertices in other components.
            vDelta[vComponents.size()-1] = 0; // the last component only has colors within the component.
            for (size_t index = vComponents.size()-1; index >= 1; --index) {
                vDelta[index-1] = vDelta[index] + vComponents[index].size();
            }

            for (size_t componentIndex = 0; componentIndex < vComponents.size(); ++componentIndex) {
                vector<int> const &vComponent(vComponents[componentIndex]);

                //// computing maximum cliques without knowing the global max clique makes us
                //// traverse many more nodes than are necessary
////                m_uMaximumCliqueSize = 0;
////                R.clear();
////                cliques.back().clear();
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

                // add delta value to the colors so that we get better pruning.
                // delta = P.Size() - vComponent.size()
                Color(vNewVertexOrder, vNewP, vNewColors); //, delta) //// like this

////                for (int &color : vNewColors) {
////                    color += vDelta[componentIndex];
////                }

                m_vStackDelta[depth+1] += vDelta[componentIndex];

////                cout << "Colors (" << vNewColors.size() << "): ";
////                for (int index = max(0,static_cast<int>(vNewColors.size() - 11)); index < static_cast<int>(vNewColors.size()); ++index) {
////                    cout << vNewColors[index] << " ";
////                }
////                cout << endl;

                // TODO/DS: need to recognize when returning because of color pruning.
                // then we can avoid evaluating the rest of the connected components entirely.

                // TODO/DS: need a flag to keep track of largest clique found found in subcall, so we can add it to R.
                // after completing recursion on a component. Otherwise, we will miss out on larger cliques entirely!

                depth++;
                RunRecursive(vNewP, vNewVertexOrder, cliques, vNewColors, false);
                depth--;

                m_vStackDelta[depth+1] -= vDelta[componentIndex];

////                cout << depth << ": Found clique of size " << cliques.back().size() << endl;
////                cout << depth << ": or              size " << m_uMaximumCliqueSize  << endl;

////                bestClique.insert(bestClique.end(), cliques.back().begin(), cliques.back().end());

////                if (bestClique.size() > savedCliqueSize) {
////                    savedCliqueSize = bestClique.size();
////                    savedClique.clear();
////                    savedClique.insert(savedClique.end(), bestClique.begin(), bestClique.end());
////                    ExecuteCallBacks(cliques.back());
////                }

                R = m_vSubgraphClique;

                uNumberOfVerticesEvaluated += vComponent.size();

        if (depth == 0 && !m_bQuiet) {
            cout << "Only " << P.size() - uNumberOfVerticesEvaluated << " more vertices to go! " << GetTimeInSeconds(clock() - startTime) << endl;
            cout << "Current best clique: " << m_uMaximumCliqueSize << endl;
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
            R = std::move(savedR); // could also pop_back() // might be faster.

            // only put back the clique if the newfound clique is smaller than
            // the already-known one.
            if (m_vSubgraphClique.size() < vSavedClique.size())
                m_vSubgraphClique = vSavedClique;
////            cliques.back() = bestClique;
////            m_uMaximumCliqueSize = savedCliqueSize;
            isolates.SetConnectedComponent(P); 

////            m_vSubgraphClique = savedClique;

            ProcessOrderBeforeReturn(vVertexOrder, P, vColors);
            return;
        }
    }
#endif // 0


        if (!bPIsEmpty && P.empty()) {
            if (R.size() > m_vSubgraphClique.size()) {
                m_vSubgraphClique = R;
            }

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
