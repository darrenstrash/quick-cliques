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
{
    SetName("connected-component-miss");
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

#if 1
void ConnectedComponentMISS::RunRecursive(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors)
{
    nodeCount++;
    vector<int> &vNewP(stackP[R.size()+1]);
    vector<int> &vNewColors(stackColors[R.size()+1]);
    vector<int> &vNewVertexOrder(stackOrder[R.size()+1]);

////    cout << depth << ": P: ";
////    for (int const p : P) {
////        cout << p << " ";
////    }
////    cout << endl;

#if 1
    if (P.size() > 100) {
        vector<vector<int>> vComponents;
        GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());
        size_t uRemainingVertices(P.size());
        size_t addedToR(0);
        if (vComponents.size() > 1) {
////            cerr << "# connected components       : " << vComponents.size() << endl << flush;
////            cerr << "size of connected components : ";
////            cout << "[ ";
////            for (vector<int> const& vComponent : vComponents) {
////                cout << vComponent.size() << " ";
////            }
////            cout << "]" << endl << flush;
            for (vector<int> const &vComponent : vComponents) {
                map<int,int> vertexRemap;
                map<int,int> reverseMap;
                int newVertex(0);
                vector<vector<char>> adjacencyMatrix(vComponent.size(), vector<char>());
                vector<vector<int>> adjacencyArray(vComponent.size(), vector<int>());
                for (int const vertex : vComponent) {
                    vertexRemap[vertex] = newVertex++;
                    reverseMap[newVertex-1] = vertex;
                    adjacencyMatrix[newVertex-1].resize(vComponent.size());
                }

////        cout << "Matrix resized to " << vComponent.size() << "x" << vComponent.size() << endl << flush;

                for (int const vertex : vComponent) {
                    int const newVertex(vertexRemap[vertex]);
                    for (int const neighbor : m_AdjacencyArray[vertex]) {
                        if (vertexRemap.find(neighbor) != vertexRemap.end()) {
                            int const newNeighbor(vertexRemap[neighbor]);
                            ////                    cout << "Adding edge " << newVertex << "," << newNeighbor << endl << flush;
                            adjacencyMatrix[newVertex][newNeighbor] = 1;
                            adjacencyArray[newVertex].push_back(newNeighbor);
                        }
                    }
                }

                list<list<int>> result;

                ConnectedComponentMISS algorithm(adjacencyMatrix, adjacencyArray);
                algorithm.SetQuiet(true);
                algorithm.SetNodeCount(nodeCount);
                algorithm.Run(result);
                if (!result.empty()) {
                    for (int const vertex : result.back()) {
                        R.push_back(reverseMap[vertex]);
                    }
                }

                nodeCount = algorithm.GetNodeCount();

                addedToR += result.back().size();

                uRemainingVertices -= vComponents.size();
                // prune if remaining graph is too small to yield a new result.
                if (R.size() + uRemainingVertices < m_uMaximumCliqueSize) {
                    for (size_t index = 0; index < addedToR; ++index) {
                        R.pop_back();
                    }
                    return;
                }

                if (R.size() > m_uMaximumCliqueSize) {
                    cliques.back().clear();
                    cliques.back().insert(cliques.back().end(), R.begin(), R.end());
                    ExecuteCallBacks(cliques.back());
                    m_uMaximumCliqueSize = R.size();
                }
            }


            for (size_t index = 0; index < addedToR; ++index) {
                R.pop_back();
            }
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
#endif // 0
