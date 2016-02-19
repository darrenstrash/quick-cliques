#include "GraphTools.h"
#include "SparseArraySet.h"
#include "ArraySet.h"
#include "Isolates2.h"
#include "Isolates3.h"
#include "Isolates4.h"
#include "IsolatesWithMatrix.h"
#include "FastIsolates.h"

#include <set>
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include <algorithm>

using namespace std;

#define UNMATCHED_VERTEX -1

void GraphTools::ComputeInducedSubgraph(vector<vector<int>> const &graph, set<int> const &vertices, vector<vector<int>> &subgraph, map<int,int> &remapping)
{
    subgraph.clear();
    remapping.clear();

////    cout << "Forming induced subgraph on " << vertices.size() << " vertices." << endl;

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
////        cout << newVertex << " : ";
        for (int const neighbor : neighbors) {
            if (vertices.find(neighbor) == vertices.end()) continue;
            int const newNeighbor = mappedVertex(neighbor);
            subgraph[newVertex].push_back(newNeighbor);
////            subgraph[newNeighbor].push_back(newVertex);
////            cout << newNeighbor << " ";
        }
////        cout << endl;
    }
}

template <typename IsolatesType>
void GraphTools::ComputeInducedSubgraphIsolates(IsolatesType const &isolates, set<int> const &vertices, vector<vector<int>> &subgraph, map<int,int> &remapping)
{
    subgraph.clear();
    remapping.clear();

////    cout << "Forming induced subgraph on " << vertices.size() << " vertices." << endl;

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

    for (int const vertex : isolates.GetInGraph()) {
        if (vertices.find(vertex) == vertices.end()) continue;

        int const newVertex = mappedVertex(vertex);
////        cout << newVertex << " : ";
        for (int const neighbor : isolates.Neighbors()[vertex]) {
            if (vertices.find(neighbor) == vertices.end()) continue;
            int const newNeighbor = mappedVertex(neighbor);
            subgraph[newVertex].push_back(newNeighbor);
////            subgraph[newNeighbor].push_back(newVertex);
////            cout << newNeighbor << " ";
        }
////        cout << endl;
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
            biDoubleGraph[vertex + size].push_back(neighbor);
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

    auto inLeftSide = [&biDoubleGraph] (int const vertex) {
        return vertex < biDoubleGraph.size()/2;
    };

    // finds a path from start to finish with excess capacity.
    auto findPath = [&inLeftSide, &biDoubleGraph] (vector<int> const &matching) {
        vector<bool> inStack(matching.size(), false);
        vector<bool> evaluated(matching.size(), false);
        list<int> stack;

        // depth first search, starting from imaginary start vertex
        for (size_t index = 0; index < matching.size()/2; ++index) {
////        for (size_t index = matching.size()/2; index > 0; --index) {
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
////            for (int index = biDoubleGraph[vertex].size(); index > 0; --index) {
////                int const neighbor(biDoubleGraph[vertex][index-1]);
                // evaluate neighbor if the edge to that neighbor has residual capacity.
                if (evaluated[neighbor] || inStack[neighbor]) continue;

                // forward edge with residual capacity
                if (inLeftSide(vertex) && matching[vertex] != neighbor) {
                    vPreviousVertexOnPath[neighbor] = vertex;
                    stack.push_back(neighbor);
                    inStack[neighbor] = true;

                    if (!inLeftSide(neighbor) && matching[neighbor] == -1) { //found path
                        foundPath = true;
                        endVertex = neighbor;
                        break;
                    }
                }
                // backward edge that we can "undo" by pushing flow back...
                else if (inLeftSide(neighbor) && matching[neighbor] == vertex) {
////                else if (!inLeftSide(vertex) && matching[neighbor] == vertex) {
                    vPreviousVertexOnPath[neighbor] = vertex;
                    stack.push_back(neighbor);
                    inStack[neighbor] = true;
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

////        cout << "Path through residual graph: ";
////        for (int const vertex : vPath) {
////            cout << vertex << " ";
////        }
////        cout << endl;
        return vPath;
    };

    vector<int> path;
    path = findPath(matching);
    while (!path.empty()) {
        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            if (matching[vertex1] != -1) {
                matching[matching[vertex1]] = -1;
            }
            if (matching[vertex2] != -1) {
                matching[matching[vertex2]] = -1;
            }
            matching[vertex1] = -1;
            matching[vertex2] = -1;
        }

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

set<int> GraphTools::ComputeLeftMIS(vector<vector<int>> const &biDoubleGraph)
{
    // each edge in matching is indexed by first vertex, and contains the second vertex
    // initially there are no edges in matching (thus invalid second vertex in matching).
    vector<int> matching(biDoubleGraph.size(), UNMATCHED_VERTEX);

    // residual graph is represented the same as edge in matching.
    // we can iterate over neighbors (excluding neighbor in matching
    // to find residual paths.
////    vector<int> residualGraph(matching);

    cout << "Computing BiDoubleMIS..." << endl << flush;

    auto inLeftSide = [&biDoubleGraph] (int const vertex) {
        return (vertex < biDoubleGraph.size()/2);
    };

    // finds a path from start to finish with excess capacity.
    auto findPath = [&inLeftSide, &biDoubleGraph] (vector<int> const &matching) {
        vector<bool> inStack(matching.size(), false);
        vector<bool> evaluated(matching.size(), false);
        list<int> stack;

        // depth first search, starting from imaginary start vertex
        // i.e., search starts from left-side.
        for (size_t index = 0; index < matching.size()/2; ++index) {
////        for (size_t index = matching.size()/2; index > 0; --index) {
            // only insert vertices without edges in matching, otherwise
            // imaginary first vertex has residual capacity 0 to that vertex.
            if (matching[index] != UNMATCHED_VERTEX) continue;
            stack.push_back(index);
            inStack[index] = true;
        }

        vector<int> vPreviousVertexOnPath(biDoubleGraph.size(), UNMATCHED_VERTEX);
        int endVertex(UNMATCHED_VERTEX);

        bool foundPath(false);
        while (!stack.empty() && !foundPath) {
            int const vertex = stack.back(); stack.pop_back();
            evaluated[vertex] = true;
            inStack[vertex] = false;
            for (int const neighbor : biDoubleGraph[vertex]) {
                // evaluate neighbor if the edge to that neighbor has residual capacity.
                if (evaluated[neighbor] || inStack[neighbor]) continue;

                // forward edge with residual capacity
                if (inLeftSide(vertex) && matching[vertex] != neighbor) {
                    vPreviousVertexOnPath[neighbor] = vertex;
                    stack.push_back(neighbor);
                    inStack[neighbor] = true;

                    if (!inLeftSide(neighbor) && matching[neighbor] == UNMATCHED_VERTEX) { //found path
                        foundPath = true;
                        endVertex = neighbor;
                        break;
                    }
                }
                // backward edge that we can "undo" by pushing flow back...
                else if (inLeftSide(neighbor) && matching[neighbor] == vertex) {
////                else if (!inLeftSide(vertex) && matching[neighbor] == vertex) {
                    vPreviousVertexOnPath[neighbor] = vertex;
                    stack.push_back(neighbor);
                    inStack[neighbor] = true;
                }
            }
        }

        vector<int> vPath;
        if (endVertex == UNMATCHED_VERTEX) return vPath;
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
        return vPath;
    };

    vector<int> path;
    path = findPath(matching);
    size_t iterations(0);
    while (!path.empty()) {
        iterations++;
////        cout << "Found path of length " << path.size() << ", iteration: " << iterations << endl << flush;
        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const backwardEdge(vertex1 > vertex2);
            if (backwardEdge) {
                matching[vertex1] = UNMATCHED_VERTEX;
                matching[vertex2] = UNMATCHED_VERTEX;
////                cout << "Remove backward edge " << vertex1 << "->" << vertex2 << " from matching" << endl << flush;
            }
////            if (matching[vertex1] != UNMATCHED_VERTEX) {
////                matching[matching[vertex1]] = UNMATCHED_VERTEX;
////            }
////            if (matching[vertex2] != UNMATCHED_VERTEX) {
////                matching[matching[vertex2]] = UNMATCHED_VERTEX;
////            }
////            matching[vertex1] = UNMATCHED_VERTEX;
////            matching[vertex2] = UNMATCHED_VERTEX;
        }

        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const forwardEdge(vertex1 < vertex2);
            if (forwardEdge) {
                matching[vertex1] = vertex2;
                matching[vertex2] = vertex1;

////                cout << "Add    forward  edge " << vertex1 << "->" << vertex2 << " to   matching" << endl << flush;
            }
        }

        path = findPath(matching);
    }

    // From each unmatched vertex, calculate the alternating paths.
    // -- Paths that alternate matched edges and non-matched edges.
    // Mark all vertices that are connected by an alternating path.

#ifdef VERIFY
    bool bVerified(true);
    for (int vertex = 0; vertex < matching.size(); ++vertex) {
        if (matching[vertex] != UNMATCHED_VERTEX) {
            if (matching[matching[vertex]] != vertex) {
                cout << "ERROR! mismatch: " << vertex << " -> " << matching[vertex] << " -> " << matching[matching[vertex]] << endl << flush;
                bVerified = false;
                break;
            }
        }
    }

    cout << "Verification of matching " << (bVerified ? "passed!" : "failed!") << endl << flush;
#endif // VERIFY

    enum LastEdge {NO_LAST_EDGE, MATCHED_LAST_EDGE, UNMATCHED_LAST_EDGE, BOTH_LAST_EDGE, NULL_LAST_EDGE};

    // TODO/DS: Make more efficient. Stop search when after marking a vertex in each matchings.
    // i.e., number of vertices marked = number of matchings.
    vector<LastEdge> marked(matching.size(), NO_LAST_EDGE);

    // first, mark unmatched vertices, and their neighbors.
    vector<int> matchedVertices;

    vector<int> previousVertex(matching.size(), -1);


    // TODO/DS: Check that lambda works 
    auto isMatchedEdge = [&matching] (int const vertex, int const neighbor) {
        return matching[neighbor] == vertex; 
    };

    int coverVertexCount(0);
    //Mark all unmatched vertices on left-hand side
    //mark their matched neighbors.
    for (int vertex = 0; vertex < matching.size()/2; ++vertex) {
        if (matching[vertex] == UNMATCHED_VERTEX) {
////            cout << "Mark " << vertex << " as unmatched!" << endl << flush;
            marked[vertex] = LastEdge::NULL_LAST_EDGE;
            coverVertexCount++;
            for (int const neighbor : biDoubleGraph[vertex]) {
                if (marked[neighbor] == LastEdge::NO_LAST_EDGE)
                    coverVertexCount++;
                marked[neighbor] = LastEdge::UNMATCHED_LAST_EDGE;
                previousVertex[neighbor] = vertex;
////                if (isMatchedEdge(vertex, neighbor)) {
////                    cout << "ERROR! (" << __LINE__ << ")" << endl << flush;
////                }
            }
        } else {
           matchedVertices.push_back(vertex); 
        }
    }

    cout << "Matching       has " << matchedVertices.size() << " vertices" << endl << flush;
    cout << "Cover count     is " << coverVertexCount                      << endl << flush;


    int iteration(0);
    while (iteration < matching.size()) {
        int const savedCoverCount(coverVertexCount);
        //cout << "Iteration: " << iteration << "/" << matching.size() << endl << flush;
        iteration++;
////        for (int const matchedVertex : matchedVertices) {
        for (int matchedVertex = 0; matchedVertex < matching.size(); matchedVertex++) {
            if (marked[matchedVertex] == NO_LAST_EDGE) continue;

            LastEdge const edgeIntoVertex(marked[matchedVertex]);

            // If the last edge was matched, follow unmatched edges
            if (edgeIntoVertex == LastEdge::MATCHED_LAST_EDGE) {
////                if (!inLeftSide(matchedVertex))
////                    cout << "ERROR! (" << __LINE__ << ") Vertex " << matchedVertex << " has   matched last edge " << endl << flush;
                for (int const neighbor : biDoubleGraph[matchedVertex]) {
                    if (!isMatchedEdge(matchedVertex, neighbor)) {
                        if (marked[neighbor] == LastEdge::NO_LAST_EDGE) {
                            marked[neighbor] = LastEdge::UNMATCHED_LAST_EDGE;
                            coverVertexCount++;
                        }
                    }
////                    if (marked[neighbor] == LastEdge::MATCHED_LAST_EDGE)
////                        cout << "ERROR! Vertex " << neighbor << " on " << (inLeftSide(neighbor)? "left": "right") << " has   matched last edge." << endl << flush;
                }

            // If last edge was unmatched, follow matched edges.
            } else if (edgeIntoVertex == LastEdge::UNMATCHED_LAST_EDGE) {
////                if (inLeftSide(matchedVertex))
////                    cout << "ERROR! (" << __LINE__ << ") Vertex " << matchedVertex << " has unmatched last edge " << endl << flush;
                for (int const neighbor : biDoubleGraph[matchedVertex]) {
                    if (isMatchedEdge(matchedVertex, neighbor)) {
                        if (marked[neighbor] == NO_LAST_EDGE) {
                            marked[neighbor] = LastEdge::MATCHED_LAST_EDGE;
                            coverVertexCount++;
                        }
////                        else if (marked[neighbor] == LastEdge::UNMATCHED_LAST_EDGE)
////                            cout << "ERROR! Vertex " << neighbor << " on " << (inLeftSide(neighbor)? "left": "right") << " has unmatched last edge." << endl << flush;
                    }
                }
            }
        }

        if (savedCoverCount != coverVertexCount) {
            cout << "Cover increased to " << coverVertexCount << endl << flush;
        }
        // TODO/DS: doesn't work, I'm not sure why...
////        else {
////            break;
////        }
    }

    int matchVertices(0);
    for (int const vertex : matching) {
        if (vertex != UNMATCHED_VERTEX) {
            matchVertices++;
        }
    }

    matchVertices >>= 1;

    vector<bool> mvc(matching.size(), false);
    size_t numMVCVertices(0);
    set<int> misToReturn;
    size_t numTotalMISVertices(0);

    // put vertices in mvc
    for (int vertex = 0; vertex < matching.size(); ++vertex) {
        if ((inLeftSide(vertex) && (marked[vertex] == LastEdge::NO_LAST_EDGE)) ||
            (!inLeftSide(vertex) && (marked[vertex] != LastEdge::NO_LAST_EDGE))) {
#ifdef VERIFY
            if (matching[vertex] == UNMATCHED_VERTEX) {
                cout << "ERROR! vertex " << vertex << " is not in matching!" << endl << flush;
            }
#endif //VERIFY
            mvc[vertex] = true;
            numMVCVertices++;

            if (inLeftSide(vertex) && marked[vertex] == LastEdge::UNMATCHED_LAST_EDGE) {
                cout << "ERROR!" << endl << flush;
            }
            if (!inLeftSide(vertex) && marked[vertex] == LastEdge::MATCHED_LAST_EDGE) {
                cout << "ERROR!" << endl << flush;
            }
        }
        else {
            numTotalMISVertices++;
            if (inLeftSide(vertex)) {
                misToReturn.insert(vertex);
            }

            // TODO/DS: I don't know if this belongs or not...
////            else {
////                misToReturn.insert(vertex - matching.size()/2);
////            }
        }
    }

#ifdef VERIFY
    for (int const vertex : misToReturn) {
        for (int const neighbor : biDoubleGraph[vertex]) {
            if (misToReturn.find(neighbor) != misToReturn.end()) {
                cout << "ERROR! Failed independent set validation." << endl << flush;
            }
        }
    }
#endif // VERIFY


    cout << "MVC             Vertices: " << numMVCVertices << endl << flush;
    cout << "MIS             Vertices: " << misToReturn.size() << endl << flush;
    cout << "Cover and Match Vertices: " << coverVertexCount << "/" << matchVertices << endl << flush;

    if (numMVCVertices != matchVertices) {
        cout << "ERROR! MVC != size of matching" << endl << flush;
    }

    // first mark unmatched vertices
    // mark them as having the last edge matched, because the only option
    // for the next edge is an unmatched edge.
////    for (int vertex = 0; vertex < matching.size(); ++vertex) {
////        if (matching[vertex] == UNMATCHED_VERTEX) {
////            marked[vertex] = LastEdge::MATCHED_LAST_EDGE;
////        }
////    }
////    for (int vertex = 0; vertex < matching.size(); ++vertex) {

////    for (int const vertex : matching) {
////        if (inLeftSide(vertex) && vertex == UNMATCHED_VERTEX) {
////            misToReturn.push_back(vertex);
////        }
////    }

    return misToReturn;
}

set<int> GraphTools::ComputeBiDoubleMIS(vector<vector<int>> const &biDoubleGraph, vector<bool> const &vInGraph, set<int> const &setInGraph)
{
    // each edge in matching is indexed by first vertex, and contains the second vertex
    // initially there are no edges in matching (thus invalid second vertex in matching).
    vector<int> matching(biDoubleGraph.size(), UNMATCHED_VERTEX);

    // residual graph is represented the same as edge in matching.
    // we can iterate over neighbors (excluding neighbor in matching
    // to find residual paths.
////    vector<int> residualGraph(matching);

////    cout << "Computing BiDoubleMIS..." << endl << flush;

    auto inLeftSide = [&biDoubleGraph] (int const vertex) {
        return (vertex < biDoubleGraph.size()/2);
    };

    auto inGraph = [&vInGraph] (int const vertex) {
        return vInGraph[vertex];
    };

    // finds a path from start to finish with excess capacity.
    auto findPath = [&inLeftSide, &biDoubleGraph, &inGraph, &setInGraph] (vector<int> const &matching) {
        vector<bool> inStack(matching.size(), false);
        vector<bool> evaluated(matching.size(), false);
        list<int> stack;

        // depth first search, starting from imaginary start vertex
        // i.e., search starts from left-side.
        for (int const vertex : setInGraph) {
////        for (size_t index = matching.size()/2; index > 0; --index) {
            // only insert vertices without edges in matching, otherwise
            // imaginary first vertex has residual capacity 0 to that vertex.
            if (matching[vertex] != UNMATCHED_VERTEX || !inLeftSide(vertex)) continue;
            stack.push_back(vertex);
            inStack[vertex] = true;
        }

        vector<int> vPreviousVertexOnPath(biDoubleGraph.size(), UNMATCHED_VERTEX);
        int endVertex(UNMATCHED_VERTEX);

        bool foundPath(false);
        while (!stack.empty() && !foundPath) {
            int const vertex = stack.back(); stack.pop_back();
            evaluated[vertex] = true;
            inStack[vertex] = false;
            for (int const neighbor : biDoubleGraph[vertex]) {
                // evaluate neighbor if the edge to that neighbor has residual capacity.
                if (!inGraph(neighbor) || evaluated[neighbor] || inStack[neighbor]) continue;

                // forward edge with residual capacity
                if (inLeftSide(vertex) && matching[vertex] != neighbor) {
                    vPreviousVertexOnPath[neighbor] = vertex;
                    stack.push_back(neighbor);
                    inStack[neighbor] = true;

                    if (!inLeftSide(neighbor) && matching[neighbor] == UNMATCHED_VERTEX) { //found path
                        foundPath = true;
                        endVertex = neighbor;
                        break;
                    }
                }
                // backward edge that we can "undo" by pushing flow back...
                else if (inLeftSide(neighbor) && matching[neighbor] == vertex) {
////                else if (!inLeftSide(vertex) && matching[neighbor] == vertex) {
                    vPreviousVertexOnPath[neighbor] = vertex;
                    stack.push_back(neighbor);
                    inStack[neighbor] = true;
                }
            }
        }

        vector<int> vPath;
        if (endVertex == UNMATCHED_VERTEX) return vPath;
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
        return vPath;
    };

    vector<int> path;
    path = findPath(matching);
    size_t iterations(0);
    while (!path.empty()) {
        iterations++;
////        cout << "Found path of length " << path.size() << ", iteration: " << iterations << endl << flush;
        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const backwardEdge(vertex1 > vertex2);
            if (backwardEdge) {
                matching[vertex1] = UNMATCHED_VERTEX;
                matching[vertex2] = UNMATCHED_VERTEX;
////                cout << "Remove backward edge " << vertex1 << "->" << vertex2 << " from matching" << endl << flush;
            }
////            if (matching[vertex1] != UNMATCHED_VERTEX) {
////                matching[matching[vertex1]] = UNMATCHED_VERTEX;
////            }
////            if (matching[vertex2] != UNMATCHED_VERTEX) {
////                matching[matching[vertex2]] = UNMATCHED_VERTEX;
////            }
////            matching[vertex1] = UNMATCHED_VERTEX;
////            matching[vertex2] = UNMATCHED_VERTEX;
        }

        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const forwardEdge(vertex1 < vertex2);
            if (forwardEdge) {
                matching[vertex1] = vertex2;
                matching[vertex2] = vertex1;

////                cout << "Add    forward  edge " << vertex1 << "->" << vertex2 << " to   matching" << endl << flush;
            }
        }

        path = findPath(matching);
    }

    // From each unmatched vertex, calculate the alternating paths.
    // -- Paths that alternate matched edges and non-matched edges.
    // Mark all vertices that are connected by an alternating path.

#ifdef VERIFY
    bool bVerified(true);
    for (int const vertex : setInGraph) {
        if (matching[vertex] != UNMATCHED_VERTEX) {
            if (matching[matching[vertex]] != vertex) {
                cout << "ERROR! mismatch: " << vertex << " -> " << matching[vertex] << " -> " << matching[matching[vertex]] << endl << flush;
                bVerified = false;
                break;
            }
        }
    }

    cout << "Verification of matching " << (bVerified ? "passed!" : "failed!") << endl << flush;
#endif // VERIFY

    enum LastEdge {NO_LAST_EDGE, MATCHED_LAST_EDGE, UNMATCHED_LAST_EDGE, BOTH_LAST_EDGE, NULL_LAST_EDGE};

    // TODO/DS: Make more efficient. Stop search when after marking a vertex in each matchings.
    // i.e., number of vertices marked = number of matchings.
    vector<LastEdge> marked(matching.size(), NO_LAST_EDGE);

    // first, mark unmatched vertices, and their neighbors.
    vector<int> matchedVertices;

    vector<int> previousVertex(matching.size(), -1);

    auto isMatchedEdge = [&matching] (int const vertex, int const neighbor) {
        return matching[neighbor] == vertex; 
    };

    int coverVertexCount(0);
    //Mark all unmatched vertices on left-hand side
    //mark their matched neighbors.
    for (int const vertex : setInGraph) {
        if (!inLeftSide(vertex)) continue;
        if (matching[vertex] == UNMATCHED_VERTEX) {
////            cout << "Mark " << vertex << " as unmatched!" << endl << flush;
            marked[vertex] = LastEdge::NULL_LAST_EDGE;
            coverVertexCount++;
            for (int const neighbor : biDoubleGraph[vertex]) {
                if (!inGraph(neighbor)) continue;
                if (marked[neighbor] == LastEdge::NO_LAST_EDGE)
                    coverVertexCount++;
                marked[neighbor] = LastEdge::UNMATCHED_LAST_EDGE;
                previousVertex[neighbor] = vertex;
////                if (isMatchedEdge(vertex, neighbor)) {
////                    cout << "ERROR! (" << __LINE__ << ")" << endl << flush;
////                }
            }
        } else {
           matchedVertices.push_back(vertex); 
        }
    }

////    cout << "Matching       has " << matchedVertices.size() << " vertices" << endl << flush;
////    cout << "Cover count     is " << coverVertexCount                      << endl << flush;


    int iteration(0);
    while (iteration < setInGraph.size()) {
        int const savedCoverCount(coverVertexCount);
        //cout << "Iteration: " << iteration << "/" << matching.size() << endl << flush;
        iteration++;
////        for (int const matchedVertex : matchedVertices) {
        for (int const matchedVertex : setInGraph) {
            if (marked[matchedVertex] == NO_LAST_EDGE) continue;

            LastEdge const edgeIntoVertex(marked[matchedVertex]);

            // If the last edge was matched, follow unmatched edges
            if (edgeIntoVertex == LastEdge::MATCHED_LAST_EDGE) {
////                if (!inLeftSide(matchedVertex))
////                    cout << "ERROR! (" << __LINE__ << ") Vertex " << matchedVertex << " has   matched last edge " << endl << flush;
                for (int const neighbor : biDoubleGraph[matchedVertex]) {
                    if (!inGraph(neighbor)) continue;
                    if (!isMatchedEdge(matchedVertex, neighbor)) {
                        if (marked[neighbor] == LastEdge::NO_LAST_EDGE) {
                            marked[neighbor] = LastEdge::UNMATCHED_LAST_EDGE;
                            coverVertexCount++;
                        }
                    }
////                    if (marked[neighbor] == LastEdge::MATCHED_LAST_EDGE)
////                        cout << "ERROR! Vertex " << neighbor << " on " << (inLeftSide(neighbor)? "left": "right") << " has   matched last edge." << endl << flush;
                }

            // If last edge was unmatched, follow matched edges.
            } else if (edgeIntoVertex == LastEdge::UNMATCHED_LAST_EDGE) {
////                if (inLeftSide(matchedVertex))
////                    cout << "ERROR! (" << __LINE__ << ") Vertex " << matchedVertex << " has unmatched last edge " << endl << flush;
                for (int const neighbor : biDoubleGraph[matchedVertex]) {
                    if (!inGraph(neighbor)) continue;
                    if (isMatchedEdge(matchedVertex, neighbor)) {
                        if (marked[neighbor] == NO_LAST_EDGE) {
                            marked[neighbor] = LastEdge::MATCHED_LAST_EDGE;
                            coverVertexCount++;
                        }
////                        else if (marked[neighbor] == LastEdge::UNMATCHED_LAST_EDGE)
////                            cout << "ERROR! Vertex " << neighbor << " on " << (inLeftSide(neighbor)? "left": "right") << " has unmatched last edge." << endl << flush;
                    }
                }
            }
        }

////        if (savedCoverCount != coverVertexCount) {
////            cout << "Cover increased to " << coverVertexCount << endl << flush;
////        }
        // TODO/DS: doesn't work, I'm not sure why...
////        else {
////            break;
////        }
    }

    int matchVertices(0);
    for (int const vertex : setInGraph) {
        if (matching[vertex] != UNMATCHED_VERTEX) {
            matchVertices++;
        }
    }

    matchVertices >>= 1;

    vector<bool> mvc(matching.size(), false);
    set<int> misToReturn;

    size_t numMVCVertices(0);
    size_t numTotalMISVertices(0);

    // put vertices in mvc
    for (int const vertex : setInGraph) {
        if ((inLeftSide(vertex) && (marked[vertex] == LastEdge::NO_LAST_EDGE)) ||
                (!inLeftSide(vertex) && (marked[vertex] != LastEdge::NO_LAST_EDGE))) {
#ifdef VERIFY
            if (matching[vertex] == UNMATCHED_VERTEX) {
                cout << "ERROR! vertex " << vertex << " is not in matching!" << endl << flush;
            }
#endif //VERIFY
            mvc[vertex] = true;
            numMVCVertices++;

            if (inLeftSide(vertex) && marked[vertex] == LastEdge::UNMATCHED_LAST_EDGE) {
                cout << "ERROR!" << endl << flush;
            }
            if (!inLeftSide(vertex) && marked[vertex] == LastEdge::MATCHED_LAST_EDGE) {
                cout << "ERROR!" << endl << flush;
            }
        }
        else {
            numTotalMISVertices++;
            if (inLeftSide(vertex)) {
                misToReturn.insert(vertex);
            }

            // TODO/DS: I don't know if this belongs or not...
            else {
                misToReturn.insert(vertex);
            }
        }
    }

#ifdef VERIFY
    for (int vertex=0; vertex < mvc.size(); ++vertex) {
        if (!inGraph(vertex)) continue;
        for (int const neighbor : biDoubleGraph[vertex]) {
            if (!inGraph(neighbor)) continue;
            if (!mvc[neighbor] && !mvc[vertex]) {
                cout << "ERROR! MVC is not a vertex cover!, edge (" << vertex << "," << neighbor << ") is not covered." << endl << flush;
            }
        }
    }

    for (int const vertex : misToReturn) {
        if (!inGraph(vertex)) {
            cout << "ERROR! vertex " << vertex << " in MIS is not in graph!" << endl;
        }
        for (int const neighbor : biDoubleGraph[vertex]) {
            if (misToReturn.find(neighbor) != misToReturn.end()) {
                cout << "ERROR! Failed independent set validation." << endl << flush;
            }
        }
    }
#endif // VERIFY


////    cout << "Graph Vertices: " << setInGraph.size() << endl << flush;          
////    cout << "Match Vertices: " << matchVertices << endl << flush;          
////    cout << "MVC   Vertices: " << numMVCVertices << endl << flush;
////    cout << "MIS   Vertices: " << misToReturn.size() << endl << flush;

    if (numMVCVertices + numTotalMISVertices != setInGraph.size()) {
        cout << "ERROR! MVC + MIS != Graph" << endl << flush;
    }

    if (numMVCVertices != matchVertices) {
        cout << "ERROR! MVC != size of matching" << endl << flush;
    }

    // first mark unmatched vertices
    // mark them as having the last edge matched, because the only option
    // for the next edge is an unmatched edge.
////    for (int vertex = 0; vertex < matching.size(); ++vertex) {
////        if (matching[vertex] == UNMATCHED_VERTEX) {
////            marked[vertex] = LastEdge::MATCHED_LAST_EDGE;
////        }
////    }
////    for (int vertex = 0; vertex < matching.size(); ++vertex) {

////    for (int const vertex : matching) {
////        if (inLeftSide(vertex) && vertex == UNMATCHED_VERTEX) {
////            misToReturn.push_back(vertex);
////        }
////    }

    return misToReturn;
}

set<int> GraphTools::ComputeLeftMIS(vector<vector<int>> const &biDoubleGraph, vector<bool> const &vInGraph, set<int> const &setInGraph)
{
    set<int> mis(GraphTools::ComputeBiDoubleMIS(biDoubleGraph, vInGraph, setInGraph));
    for (set<int>::iterator it = mis.begin(); it != mis.end(); ++it) {
        if (*it >= biDoubleGraph.size()/2) {
            mis.erase(it, mis.end());
            break;
        }
    }

    return mis;
}


void GraphTools::PrintGraphInEdgesFormat(vector<vector<int>> const &adjacencyArray)
{
    cout << adjacencyArray.size() << endl;
    size_t edges(0);
    for (vector<int> const &neighborList : adjacencyArray) {
        edges+= neighborList.size();
    }
    cout << edges << endl;

    for (size_t index = 0; index < adjacencyArray.size(); ++index) {
        for (int const neighbor : adjacencyArray[index]) {
            cout << index << "," << neighbor << endl;
        }
    }
}

void GraphTools::PrintGraphInSNAPFormat(vector<vector<int>> const &adjacencyArray)
{
    for (size_t index = 0; index < adjacencyArray.size(); ++index) {
        for (int const neighbor : adjacencyArray[index]) {
            cout << (index+1) << " " << (neighbor+1) << endl << flush;
        }
    }
}

template
void GraphTools::ComputeInducedSubgraphIsolates<Isolates4<SparseArraySet>>(Isolates4<SparseArraySet> const &isolates, set<int> const &vertices, vector<vector<int>> &subgraph, map<int,int> &remapping);

template
void GraphTools::ComputeConnectedComponents<Isolates2<ArraySet>>(Isolates2<ArraySet> const &isolates, vector<vector<int>> &vComponents, size_t const uNumVertices);

template
void GraphTools::ComputeConnectedComponents<Isolates3<ArraySet>>(Isolates3<ArraySet> const &isolates, vector<vector<int>> &vComponents, size_t const uNumVertices);

template
void GraphTools::ComputeConnectedComponents<Isolates4<SparseArraySet>>(Isolates4<SparseArraySet> const &isolates, vector<vector<int>> &vComponents, size_t const uNumVertices);
template
void GraphTools::ComputeConnectedComponents<Isolates4<ArraySet>>(Isolates4<ArraySet> const &isolates, vector<vector<int>> &vComponents, size_t const uNumVertices);

template
void GraphTools::ComputeConnectedComponents<IsolatesWithMatrix<ArraySet>>(IsolatesWithMatrix<ArraySet> const &isolates, vector<vector<int>> &vComponents, size_t const uNumVertices);

#if 0
template
void GraphTools::ComputeConnectedComponents<FastIsolates<SparseArraySet>>(FastIsolates<SparseArraySet> const &isolates, vector<vector<int>> &vComponents, size_t const uNumVertices);
#endif // 0
