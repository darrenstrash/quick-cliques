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
{
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

bool BiDoubleGraph::ComputeResidualPath(vector<int> const &vMatching, vector<int> &vPath) const
{
    vector<bool> inStack(m_AdjacencyList.size(), false);
    vector<bool> evaluated(m_AdjacencyList.size(), false);
    list<int> stack;

    vPath.clear();

    size_t numVerticesInLeftSide(vMatching.size()/2);

    // depth first search, starting from imaginary start vertex
    // i.e., search starts from left-side.
    for (size_t index = 0; index < numVerticesInLeftSide; ++index) {
        ////        for (size_t index = numVerticesInLeftSide; index > 0; --index) {
        // only insert vertices without edges in matching, otherwise
        // imaginary first vertex has residual capacity 0 to that vertex.
        if (vMatching[index] != UNMATCHED_VERTEX) continue;
        stack.push_back(index);
        inStack[index] = true;
    }

    vector<int> vPreviousVertexOnPath(m_AdjacencyList.size(), UNMATCHED_VERTEX);
    int endVertex(UNMATCHED_VERTEX);

    bool foundPath(false);
    while (!stack.empty() && !foundPath) {
        int const vertex = stack.back(); stack.pop_back();
        evaluated[vertex] = true;
        inStack[vertex] = false;
        for (int const neighbor : m_AdjacencyList[vertex]) {
            // evaluate neighbor if the edge to that neighbor has residual capacity.
            if (evaluated[neighbor] || inStack[neighbor]) continue;

            // forward edge with residual capacity
            if (InLeftSide(vertex) && vMatching[vertex] != neighbor) {
                vPreviousVertexOnPath[neighbor] = vertex;
                stack.push_back(neighbor);
                inStack[neighbor] = true;

                if (!InLeftSide(neighbor) && vMatching[neighbor] == UNMATCHED_VERTEX) { //found path
                    foundPath = true;
                    endVertex = neighbor;
                    break;
                }
            }
            // backward edge that we can "undo" by pushing flow back...
            else if (InLeftSide(neighbor) && vMatching[neighbor] == vertex) {
                ////                else if (!InLeftSide(vertex) && vMatching[neighbor] == vertex) {
                vPreviousVertexOnPath[neighbor] = vertex;
                stack.push_back(neighbor);
                inStack[neighbor] = true;
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

void BiDoubleGraph::ComputeMaximumMatching(vector<int> &vMatching) const
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
