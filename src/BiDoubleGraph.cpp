#include "BiDoubleGraph.h"
#include "GraphTools.h"

#include <set>
#include <vector>
#include <list>
#include <iostream>

#define UNMATCHED_VERTEX -1

using namespace std;

BiDoubleGraph::BiDoubleGraph(vector<vector<int>> adjacencyList)
: m_AdjacencyList()
, m_Stack()
////, m_InStack(adjacencyList.size()*2, false)
, m_Evaluated(adjacencyList.size()*2, -1)
{
    m_Stack.reserve(adjacencyList.size()*2);
    m_AdjacencyList = std::move(GraphTools::ComputeBiDoubleGraph(adjacencyList));
}

BiDoubleGraph::~BiDoubleGraph()
{
}

vector<int> const &BiDoubleGraph::Neighbors(int const vertex) const
{
    return m_AdjacencyList[vertex];
}

bool BiDoubleGraph::InLeftSide(int const vertex) const
{
    return vertex < m_AdjacencyList.size()/2;
}

void BiDoubleGraph::PushOnStack(int const vertex)
{
    m_Stack.push_back(vertex);
    m_Evaluated[vertex] = m_iCurrentRound;
}

int BiDoubleGraph::PopOffStack()
{
    int const top(m_Stack.back());
    m_Stack.pop_back();
    return top;
}

bool BiDoubleGraph::IsEvaluated(int const vertex) const
{
    return m_Evaluated[vertex] == m_iCurrentRound;
}

bool BiDoubleGraph::ComputeResidualPath(vector<int> const &vMatching, vector<int> &vPath)
{
////    m_InStack.clear(); m_InStack.resize(m_AdjacencyList.size(), false);
////    m_Evaluated.clear(); m_Evaluated.resize(m_AdjacencyList.size(), false);

    vPath.clear();
    m_Stack.clear();
    m_iCurrentRound++;
    if (m_iCurrentRound == 0) {
        m_Evaluated.clear(); m_Evaluated.resize(m_AdjacencyList.size(), -1);
    }

    size_t numVerticesInLeftSide(vMatching.size()/2);

    // depth first search, starting from imaginary start vertex
    // i.e., search starts from left-side.
    for (size_t index = 0; index < numVerticesInLeftSide; ++index) {
        ////        for (size_t index = numVerticesInLeftSide; index > 0; --index) {
        // only insert vertices without edges in matching, otherwise
        // imaginary first vertex has residual capacity 0 to that vertex.
        if (vMatching[index] != UNMATCHED_VERTEX) continue;
        PushOnStack(index);
    }

    vector<int> vPreviousVertexOnPath(m_AdjacencyList.size(), UNMATCHED_VERTEX);
    int endVertex(UNMATCHED_VERTEX);

    bool foundPath(false);
    while (!m_Stack.empty() && !foundPath) {
        int const vertex = PopOffStack();
        for (int const neighbor : m_AdjacencyList[vertex]) {
            // evaluate neighbor if the edge to that neighbor has residual capacity.
            if (IsEvaluated(neighbor)) continue;

            // forward edge with residual capacity
            if (InLeftSide(vertex) && vMatching[vertex] != neighbor) {
                vPreviousVertexOnPath[neighbor] = vertex;
                PushOnStack(neighbor);

                if (!InLeftSide(neighbor) && vMatching[neighbor] == UNMATCHED_VERTEX) { //found path
                    foundPath = true;
                    endVertex = neighbor;
                    break;
                }
            }

            // backward edge that we can "undo" by pushing flow back...
            else if (InLeftSide(neighbor) && vMatching[neighbor] == vertex) {
                vPreviousVertexOnPath[neighbor] = vertex;
                PushOnStack(neighbor);
            }
        }
    }

    if (endVertex == UNMATCHED_VERTEX) return false;
    vPath.push_back(endVertex);
    while (vPreviousVertexOnPath[endVertex] != UNMATCHED_VERTEX) {
        vPath.push_back(vPreviousVertexOnPath[endVertex]);
        endVertex = vPreviousVertexOnPath[endVertex];
    }

    std::reverse(vPath.begin(), vPath.end());

    ////        cout << "Path through residual graph: ";
    ////        for (int const vertex : vPath) {
    ////            cout << vertex << " ";
    ////        }
    ////        cout << endl;
    return true;
}

void BiDoubleGraph::ComputeMaximumMatching(vector<int> &vMatching)
{
    vector<int> path;
    path.reserve(vMatching.size());
    ComputeResidualPath(vMatching, path);
    size_t iterations(0);
    while (!path.empty()) {
        iterations++;
////        cout << "Found path of length " << path.size() << ", iteration: " << iterations << endl << flush;
        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const backwardEdge(vertex1 > vertex2);
            if (backwardEdge) {
                vMatching[vertex1] = UNMATCHED_VERTEX;
                vMatching[vertex2] = UNMATCHED_VERTEX;
////                cout << "Remove backward edge " << vertex1 << "->" << vertex2 << " from vMatching" << endl << flush;
            }
////            if (vMatching[vertex1] != UNMATCHED_VERTEX) {
////                vMatching[vMatching[vertex1]] = UNMATCHED_VERTEX;
////            }
////            if (vMatching[vertex2] != UNMATCHED_VERTEX) {
////                vMatching[vMatching[vertex2]] = UNMATCHED_VERTEX;
////            }
////            vMatching[vertex1] = UNMATCHED_VERTEX;
////            vMatching[vertex2] = UNMATCHED_VERTEX;
        }

        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const forwardEdge(vertex1 < vertex2);
            if (forwardEdge) {
                vMatching[vertex1] = vertex2;
                vMatching[vertex2] = vertex1;

////                cout << "Add    forward  edge " << vertex1 << "->" << vertex2 << " to   vMatching" << endl << flush;
            }
        }

        ComputeResidualPath(vMatching, path);
    }

#ifdef VERIFY
    bool bVerified(true);
    for (int vertex = 0; vertex < vMatching.size(); ++vertex) {
        if (vMatching[vertex] != UNMATCHED_VERTEX) {
            if (vMatching[vMatching[vertex]] != vertex) {
                cout << "ERROR! mismatch: " << vertex << " -> " << vMatching[vertex] << " -> " << vMatching[vMatching[vertex]] << endl << flush;
                bVerified = false;
                break;
            }
        }
    }


    cout << "Verification of matching " << (bVerified ? "passed!" : "failed!") << endl << flush;
#endif // VERIFY
}


//
//  Versions that track the vertices in the graph.
//  TODO/DS: Refactor to remove vInGraph and setInGraph.
//

bool BiDoubleGraph::ComputeResidualPath(vector<int> const &vMatching, vector<int> &vPath, vector<bool> const &vInGraph, set<int> const &setInGraph)
{
////    m_InStack.clear(); m_InStack.resize(m_AdjacencyList.size(), false);
////    m_Evaluated.clear(); m_Evaluated.resize(m_AdjacencyList.size(), false);

    vPath.clear();
    m_Stack.clear();
    m_iCurrentRound++;
    if (m_iCurrentRound == 0) {
        m_Evaluated.clear(); m_Evaluated.resize(m_AdjacencyList.size(), -1);
    }

    auto inGraph = [&vInGraph] (int const vertex) {
        return vInGraph[vertex];
    };

////    size_t numVerticesInLeftSide(vMatching.size()/2);

    // depth first search, starting from imaginary start vertex
    // i.e., search starts from left-side.
////    for (size_t index = 0; index < numVerticesInLeftSide; ++index) {
    for (int const vertex : setInGraph) {
        ////        for (size_t index = numVerticesInLeftSide; index > 0; --index) {
        // only insert vertices without edges in matching, otherwise
        // imaginary first vertex has residual capacity 0 to that vertex.
        if (vMatching[vertex] != UNMATCHED_VERTEX || !InLeftSide(vertex)) continue;
        PushOnStack(vertex);
    }

    vector<int> vPreviousVertexOnPath(m_AdjacencyList.size(), UNMATCHED_VERTEX);
    int endVertex(UNMATCHED_VERTEX);

    bool foundPath(false);
    while (!m_Stack.empty() && !foundPath) {
        int const vertex = PopOffStack();
        for (int const neighbor : m_AdjacencyList[vertex]) {
            // evaluate neighbor if the edge to that neighbor has residual capacity.
            if (!inGraph(neighbor) || IsEvaluated(neighbor)) continue;

            // forward edge with residual capacity
            if (InLeftSide(vertex) && vMatching[vertex] != neighbor) {
                vPreviousVertexOnPath[neighbor] = vertex;
                PushOnStack(neighbor);

                if (!InLeftSide(neighbor) && vMatching[neighbor] == UNMATCHED_VERTEX) { //found path
                    foundPath = true;
                    endVertex = neighbor;
                    break;
                }
            }

            // backward edge that we can "undo" by pushing flow back...
            else if (InLeftSide(neighbor) && vMatching[neighbor] == vertex) {
                vPreviousVertexOnPath[neighbor] = vertex;
                PushOnStack(neighbor);
            }
        }
    }

    if (endVertex == UNMATCHED_VERTEX) return false;
    vPath.push_back(endVertex);
    while (vPreviousVertexOnPath[endVertex] != UNMATCHED_VERTEX) {
        vPath.push_back(vPreviousVertexOnPath[endVertex]);
        endVertex = vPreviousVertexOnPath[endVertex];
    }

    std::reverse(vPath.begin(), vPath.end());

    ////        cout << "Path through residual graph: ";
    ////        for (int const vertex : vPath) {
    ////            cout << vertex << " ";
    ////        }
    ////        cout << endl;
    return true;
}

void BiDoubleGraph::ComputeMaximumMatching(vector<int> &vMatching, vector<bool> const &vInGraph, set<int> const &setInGraph)
{
    vector<int> path;
    path.reserve(vMatching.size());
    ComputeResidualPath(vMatching, path, vInGraph, setInGraph);
    size_t iterations(0);
    while (!path.empty()) {
        iterations++;
////        cout << "Found path of length " << path.size() << ", iteration: " << iterations << endl << flush;
        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const backwardEdge(vertex1 > vertex2);
            if (backwardEdge) {
                vMatching[vertex1] = UNMATCHED_VERTEX;
                vMatching[vertex2] = UNMATCHED_VERTEX;
////                cout << "Remove backward edge " << vertex1 << "->" << vertex2 << " from vMatching" << endl << flush;
            }
////            if (vMatching[vertex1] != UNMATCHED_VERTEX) {
////                vMatching[vMatching[vertex1]] = UNMATCHED_VERTEX;
////            }
////            if (vMatching[vertex2] != UNMATCHED_VERTEX) {
////                vMatching[vMatching[vertex2]] = UNMATCHED_VERTEX;
////            }
////            vMatching[vertex1] = UNMATCHED_VERTEX;
////            vMatching[vertex2] = UNMATCHED_VERTEX;
        }

        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const forwardEdge(vertex1 < vertex2);
            if (forwardEdge) {
                vMatching[vertex1] = vertex2;
                vMatching[vertex2] = vertex1;

////                cout << "Add    forward  edge " << vertex1 << "->" << vertex2 << " to   vMatching" << endl << flush;
            }
        }

        ComputeResidualPath(vMatching, path, vInGraph, setInGraph);
    }

#ifdef VERIFY
    bool bVerified(true);
    for (int vertex = 0; vertex < vMatching.size(); ++vertex) {
        if (vMatching[vertex] != UNMATCHED_VERTEX) {
            if (vMatching[vMatching[vertex]] != vertex) {
                cout << "ERROR! mismatch: " << vertex << " -> " << vMatching[vertex] << " -> " << vMatching[vMatching[vertex]] << endl << flush;
                bVerified = false;
                break;
            }
        }
    }


    cout << "Verification of matching " << (bVerified ? "passed!" : "failed!") << endl << flush;
#endif // VERIFY
}
