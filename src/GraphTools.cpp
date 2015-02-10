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

template
void GraphTools::ComputeConnectedComponents<Isolates2<ArraySet>>(Isolates2<ArraySet> const &isolates, vector<vector<int>> &vComponents, size_t const uNumVertices);

template
void GraphTools::ComputeConnectedComponents<Isolates3<ArraySet>>(Isolates3<ArraySet> const &isolates, vector<vector<int>> &vComponents, size_t const uNumVertices);

