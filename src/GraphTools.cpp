#include "GraphTools.h"
#include "SparseArraySet.h"
#include "ArraySet.h"
#include "Isolates2.h"
#include "Isolates3.h"

#include <set>
#include <vector>
#include <list>
#include <map>
#include <iostream>

using namespace std;

void GraphTools::ComputeInducedSubgraph(vector<vector<int>> &graph, set<int> const &vertices, vector<vector<int>> &subgraph, map<int,int> &remapping)
{
    subgraph.clear();
    remapping.clear();

    cout << "Forming induced subgraph on " << vertices.size() << " vertices." << endl;

    map<int,int> forwardMapping;

    int vertexIndex(0);
    auto mappedVertex = [&vertexIndex, &remapping, &forwardMapping](int const vertex)
    {
        if (forwardMapping.find(vertex) == forwardMapping.end()) {
            forwardMapping[vertex] = vertexIndex;
            remapping[vertexIndex] = vertex;
            vertexIndex++;
        }
        return forwardMapping[vertex];
    };

    for (int const vertex : vertices) {
        mappedVertex(vertex);
    }

    subgraph.resize(vertices.size());

    for (int vertex = 0; vertex < graph.size(); ++vertex) {
        if (vertices.find(vertex) == vertices.end()) continue;

        vector<int> const &neighbors(graph[vertex]);
        int const newVertex = mappedVertex(vertex);
        cout << newVertex << " : ";
        for (int const neighbor : neighbors) {
            if (vertices.find(neighbor) == vertices.end()) continue;
            int const newNeighbor = mappedVertex(neighbor);
            subgraph[newVertex].push_back(newNeighbor);
////            subgraph[newNeighbor].push_back(newVertex);
            cout << newNeighbor << " ";
        }
        cout << endl;
    }
}

vector<int> GraphTools::OrderVerticesByDegree(vector<vector<int>> const &adjacencyList, bool const ascending)
{
    vector<int> vOrderedVertices(adjacencyList.size(), -1);

    size_t maxDegree(0);
    for (vector<int> const &neighbors : adjacencyList) {
        maxDegree = max(maxDegree, neighbors.size());
    }

    vector<list<int>> vlVerticesByDegree(maxDegree + 1, list<int>());

    for (size_t vertex = 0; vertex < adjacencyList.size(); ++vertex) {
////        std::cout << "maxDegree=" << maxDegree << ", degree=" << adjacencyList[vertex].size() << endl << flush;
        vlVerticesByDegree[adjacencyList[vertex].size()].push_back(vertex);
    }

    if (ascending) {
        size_t index(0);
        for (size_t degree = 0; degree <= maxDegree; ++degree) {
            for (int const vertex : vlVerticesByDegree[degree]) {
                vOrderedVertices[index++] = vertex;
            }
        }
    } else {
        size_t index(0);
        for (size_t degree = 0; degree <= maxDegree; ++degree) {
            for (int const vertex : vlVerticesByDegree[maxDegree - degree]) {
                vOrderedVertices[index++] = vertex;
            }
        }
    }

    return vOrderedVertices;
}

vector<int> GraphTools::OrderVerticesByDegree(vector<vector<char>> const &adjacencyMatrix, vector<int> const &vDegree, bool const ascending)
{
    vector<int> vOrderedVertices(adjacencyMatrix.size(), -1);

    size_t maxDegree(0);
    for (int const degree : vDegree) {
        maxDegree = max(maxDegree, static_cast<size_t>(degree));
    }

    vector<list<int>> vlVerticesByDegree(maxDegree + 1, list<int>());

    for (size_t vertex = 0; vertex < adjacencyMatrix.size(); ++vertex) {
////        std::cout << "maxDegree=" << maxDegree << ", degree=" << adjacencyList[vertex].size() << endl << flush;
        vlVerticesByDegree[vDegree[vertex]].push_back(vertex);
    }

    if (ascending) {
        size_t index(0);
        for (size_t degree = 0; degree <= maxDegree; ++degree) {
            for (int const vertex : vlVerticesByDegree[degree]) {
                vOrderedVertices[index++] = vertex;
            }
        }
    } else {
        size_t index(0);
        for (size_t degree = 0; degree <= maxDegree; ++degree) {
            for (int const vertex : vlVerticesByDegree[maxDegree - degree]) {
                vOrderedVertices[index++] = vertex;
            }
        }
    }

    return vOrderedVertices;
}

vector<int> GraphTools::OrderVerticesByDegree(vector<vector<char>> const &adjacencyMatrix, bool const ascending)
{
    vector<int> vDegree(adjacencyMatrix.size(), 0);
    for (size_t u = 0; u < adjacencyMatrix.size(); ++u) {
        for (size_t v = 0; v < adjacencyMatrix.size(); ++v) {
            if (adjacencyMatrix[u][v]) vDegree[u]++;
        }
    }
    return std::move(OrderVerticesByDegree(adjacencyMatrix, vDegree, ascending));
}

// NOT CURRENTLY USED.

////void GraphTools::RemoveVertices(vector<vector<int>> &adjacencyList, vector<int> const &vVertices)
////{
////    for (int const vertex : vVertices) {
////        if (vertex >= adjacencyList.size() || vertex < 0) {
////            continue;
////        }
////
////        vector<int> &neighbors(adjacencyList[vertex]);
////        for (int const neighbor : neighbors) {
////
////        }
////    }
////}

vector<int> GraphTools::OrderVerticesByDegree(ArraySet const &inGraph, vector<ArraySet> const &neighborSets, bool const ascending)
{
    vector<int> vOrderedVertices(inGraph.Size(), -1);

    size_t maxDegree(0);
    for (ArraySet const &neighborSet : neighborSets) {
        maxDegree = max(maxDegree, neighborSet.Size());
    }

    vector<list<int>> vlVerticesByDegree(maxDegree + 1, list<int>());

    for (size_t vertex = 0; vertex < neighborSets.size(); ++vertex) {
        if (!inGraph.Contains(vertex)) continue;
////        std::cout << "maxDegree=" << maxDegree << ", degree=" << adjacencyList[vertex].size() << endl << flush;
        vlVerticesByDegree[neighborSets[vertex].Size()].push_back(vertex);
    }

    if (ascending) {
        size_t index(0);
        for (size_t degree = 0; degree <= maxDegree; ++degree) {
            for (int const vertex : vlVerticesByDegree[degree]) {
                vOrderedVertices[index++] = vertex;
            }
        }
    } else {
        size_t index(0);
        for (size_t degree = 0; degree <= maxDegree; ++degree) {
            for (int const vertex : vlVerticesByDegree[maxDegree - degree]) {
                vOrderedVertices[index++] = vertex;
            }
        }
    }

    return vOrderedVertices;
}

vector<int> GraphTools::OrderVerticesByDegree(ArraySet const &inGraph, vector<SparseArraySet> const &neighborSets, bool const ascending)
{
    vector<int> vOrderedVertices(inGraph.Size(), -1);

    size_t maxDegree(0);
    for (SparseArraySet const &neighborSet : neighborSets) {
        maxDegree = max(maxDegree, neighborSet.Size());
    }

    vector<list<int>> vlVerticesByDegree(maxDegree + 1, list<int>());

    for (size_t vertex = 0; vertex < neighborSets.size(); ++vertex) {
        if (!inGraph.Contains(vertex)) continue;
////        std::cout << "maxDegree=" << maxDegree << ", degree=" << adjacencyList[vertex].size() << endl << flush;
        vlVerticesByDegree[neighborSets[vertex].Size()].push_back(vertex);
    }

    if (ascending) {
        size_t index(0);
        for (size_t degree = 0; degree <= maxDegree; ++degree) {
            for (int const vertex : vlVerticesByDegree[degree]) {
                vOrderedVertices[index++] = vertex;
            }
        }
    } else {
        size_t index(0);
        for (size_t degree = 0; degree <= maxDegree; ++degree) {
            for (int const vertex : vlVerticesByDegree[maxDegree - degree]) {
                vOrderedVertices[index++] = vertex;
            }
        }
    }

    return vOrderedVertices;
}

template<typename IsolatesType>
void GraphTools::ComputeConnectedComponents(IsolatesType const &isolates, vector<vector<int>> &vComponents, size_t const uNumVertices) {
    ArraySet remaining = isolates.GetInGraph();

    ArraySet currentSearch(uNumVertices);
    vector<bool> evaluated(uNumVertices, 0);

    size_t componentCount(0);
    vComponents.clear();

    if (!remaining.Empty()) {
        int const startVertex = *remaining.begin();
        currentSearch.Insert(startVertex);
        remaining.Remove(startVertex);
        componentCount++;
        vComponents.resize(componentCount);
    }

    while (!remaining.Empty() && !currentSearch.Empty()) {
        int const nextVertex(*currentSearch.begin());
        evaluated[nextVertex] = true;
        vComponents[componentCount - 1].push_back(nextVertex);
        currentSearch.Remove(nextVertex);
        remaining.Remove(nextVertex);
        for (int const neighbor : isolates.Neighbors()[nextVertex]) {
            if (!evaluated[neighbor]) {
                currentSearch.Insert(neighbor);
            }
        }

        if (currentSearch.Empty() && !remaining.Empty()) {
            int const startVertex = *remaining.begin();
            currentSearch.Insert(startVertex);
            remaining.Remove(startVertex);
            componentCount++;
            vComponents.resize(componentCount);
        }
    }
}

vector<vector<int>> GraphTools::ComputeBiDoubleGraph(vector<vector<int>> const &adjacencyArray)
{
    vector<vector<int>> biDoubleGraph(adjacencyArray.size()*2);

    int const size(adjacencyArray.size());
    for (int vertex = 0; vertex < adjacencyArray.size(); ++vertex) {
        vector<int> const &neighbors(adjacencyArray[vertex]);
        for (int const neighbor : neighbors) {
            biDoubleGraph[vertex].push_back(neighbor + size);
        }
    }

    return biDoubleGraph;
}


// compute size of unweighted maximum matchings in biDouble graphs (which are bipartite).
// first half of biDOuble graph is the first set of the bipartite graph.
// simplified implementation of Ford-Fulkerson algorithm.
int GraphTools::ComputeMaximumMatchingSize(vector<vector<int>> const &biDoubleGraph)
{
    // each edge in matching is indexed by first vertex, and contains the second vertex
    // initially there are no edges in matching (thus invalid second vertex -1).
    vector<int> matching(biDoubleGraph.size(), -1);

    // residual graph is represented the same as edge in matching.
    // we can iterate over neighbors (excluding neighbor in matching
    // to find residual paths.
////    vector<int> residualGraph(matching);

    auto inRightSide = [&biDoubleGraph] (int const vertex) {
        return vertex < biDoubleGraph.size()/2;
    };

    // finds a path from start to finish with excess capacity.
    auto findPath = [&inRightSide, &biDoubleGraph] (vector<int> const &matching) {
        vector<bool> inStack(matching.size());
        vector<bool> evaluated(matching.size());
        list<int> stack;

        // depth first search, starting from imaginary start vertex
        for (size_t index = 0; index < matching.size()/2; ++index) {
            // only insert vertices without edges in matching, otherwise
            // imaginary first vertex has residual capacity 0 to that vertex.
            if (matching[index] != -1) continue;
            stack.push_back(index);
            inStack[index] = true;
        }

        vector<int> vPreviousVertexOnPath(biDoubleGraph.size(), -1);
        int endVertex(-1);

        bool foundPath(false);
        while (!stack.empty() && !foundPath) {
            int const vertex = stack.back(); stack.pop_back();
            evaluated[vertex] = true;
            inStack[vertex] = false;
            for (int const neighbor : biDoubleGraph[vertex]) {
                // evaluate neighbor if the edge to that neighbor has residual capacity.
                if (evaluated[neighbor] || inStack[neighbor]) continue;

                stack.push_back(neighbor);
                inStack[neighbor] = true;

                // forward edge with residual capacity
                if (inRightSide(vertex) && matching[vertex] != neighbor) {
                    vPreviousVertexOnPath[neighbor] = vertex;
                    if (!inRightSide(neighbor) && matching[neighbor] == -1) { //found path
                        foundPath = true;
                        endVertex = neighbor;
                        break;
                    }
                }
                // backward edge that we can "undo" by pushing flow back...
                else if (!inRightSide(vertex) && matching[neighbor] == vertex) {
                    vPreviousVertexOnPath[neighbor] = vertex;
                    stack.push_back(neighbor);
                }
            }
        }

        vector<int> vPath;
        if (endVertex == -1) return vPath;
        vPath.push_back(endVertex);
        while (vPreviousVertexOnPath[endVertex] != -1) {
            vPath.push_back(vPreviousVertexOnPath[endVertex]);
            endVertex = vPreviousVertexOnPath[endVertex];
        }

        std::reverse(vPath.begin(), vPath.end());

        return vPath;
    };

    vector<int> path;
    path = findPath(matching);
    while (!path.empty()) {
        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const forwardEdge(vertex1 < vertex2);
            if (forwardEdge) {
                matching[vertex1] = vertex2;
                matching[vertex2] = vertex1;
            }
        }

        path = findPath(matching);
    }

    int count(0);
    for (int const vertex : matching) {
        if (vertex == -1) continue;
        count++;
    }

    return count/2;
}

bool GraphTools::TestMatchingCount()
{
    cout << "Maximum Matching: ";
    vector<vector<int>> vAdjacencyList(6);
    vAdjacencyList[0].push_back(5);
    vAdjacencyList[5].push_back(0);
    if (ComputeMaximumMatchingSize(vAdjacencyList) != 1) {
        cout << "FAILED: MaximumMatching one-edge test" << endl << flush;
        return false;
    }

    vAdjacencyList[0].push_back(4);
    vAdjacencyList[4].push_back(0);
    if (ComputeMaximumMatchingSize(vAdjacencyList) != 1) {
        cout << "FAILED: MaximumMatching two-edge-one-path test" << endl << flush;
        return false;
    }

    vAdjacencyList[1].push_back(4);
    vAdjacencyList[4].push_back(1);
    if (ComputeMaximumMatchingSize(vAdjacencyList) != 2) {
        cout << "FAILED: MaximumMatching three-edge-two-path test" << endl << flush;
        return false;
    }

    vAdjacencyList[2].push_back(4);
    vAdjacencyList[4].push_back(2);
    if (ComputeMaximumMatchingSize(vAdjacencyList) != 2) {
        cout << "FAILED: MaximumMatching four-edge-two-path test" << endl << flush;
        return false;
    }

    vAdjacencyList[0] = {3,4};
    vAdjacencyList[1] = {};
    vAdjacencyList[2] = {4,5};
    vAdjacencyList[3] = {0};
    vAdjacencyList[4] = {0,2};
    vAdjacencyList[5] = {2};

    if (ComputeMaximumMatchingSize(vAdjacencyList) != 2) {
        cout << "FAILED: MaximumMatching \"sigma\" case" << endl << flush;
        return false;
    }

    vAdjacencyList[0] = {3,4,5};
    vAdjacencyList[1] = {3,4,5};
    vAdjacencyList[2] = {3,4,5};
    vAdjacencyList[3] = {0,1,2};
    vAdjacencyList[4] = {0,1,2};
    vAdjacencyList[5] = {0,1,2};

    if (ComputeMaximumMatchingSize(vAdjacencyList) != 3) {
        cout << "FAILED: MaximumMatching full bipartite case" << endl << flush;
        return false;
    }

    cout << "PASSED!" << endl;
    return true;
}

template
void GraphTools::ComputeConnectedComponents<Isolates2<ArraySet>>(Isolates2<ArraySet> const &isolates, vector<vector<int>> &vComponents, size_t const uNumVertices);

template
void GraphTools::ComputeConnectedComponents<Isolates3<ArraySet>>(Isolates3<ArraySet> const &isolates, vector<vector<int>> &vComponents, size_t const uNumVertices);

