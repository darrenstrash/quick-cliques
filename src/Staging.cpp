// local includes
#include "Staging.h"
#include "CliqueTools.h"
#include "Isolates4.h"
#include "LightWeightReductionMISR.h"
#include "LightWeightMISQ.h"
#include "LightWeightStaticOrderMISS.h"
#include "LightWeightFullMISS.h"
#include "LightWeightReductionStaticOrderMISS.h"
#include "LightWeightReductionFullMISS.h"
#include "LightWeightReductionSparseFullMISS.h"
#include "LightWeightReductionSparseStaticOrderMISS.h"
#include "TesterMISS.h"
#include "ConnectedComponentMISS.h"
#include "OrderingTools.h"
#include "GraphTools.h"
#include "Tools.h"

#include "IsolatesIndependentSetColoringStrategy.h"

// system includes
#include <vector>
#include <list>
#include <set>
#include <iostream>
#include <algorithm>
#include <climits>
#include <random>

////#define TWO_LEVEL

using namespace std;

Staging::Staging(std::vector<std::vector<int>> &adjacencyList)
 : Algorithm("staging")
 , m_AdjacencyList(adjacencyList)
{
}

Staging::~Staging()
{
}

size_t ComputeConnectedComponents(Isolates4<SparseArraySet> const &isolates, vector<vector<int>> &vComponents)
{
    ArraySet remaining = isolates.GetInGraph();

    Set currentSearch;
    Set evaluated;

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
        evaluated.Insert(nextVertex);
        vComponents[componentCount - 1].push_back(nextVertex);
        currentSearch.Remove(nextVertex);
        remaining.Remove(nextVertex);
        for (int const neighbor : isolates.Neighbors()[nextVertex]) {
            if (!evaluated.Contains(neighbor)) {
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

    return componentCount;
}

vector<vector<int>> ComputeSubgraphOfSize(Isolates4<SparseArraySet> const &isolates, size_t const graphSize, size_t const subgraphSize)
{
    ArraySet remaining = isolates.GetInGraph();

    map<int,int> vertexRemap;

    Set currentSearch;
    vector<bool> evaluated(graphSize, false);

    size_t componentCount(0);

    if (!remaining.Empty()) {
        int const startVertex = *remaining.begin();
        currentSearch.Insert(startVertex);
        remaining.Remove(startVertex);
        componentCount++;
    }

    size_t numEvaluated(0);

    while (!remaining.Empty() && numEvaluated < subgraphSize && !currentSearch.Empty()) {
        int const nextVertex(*currentSearch.begin());
        evaluated[nextVertex] = true;
        vertexRemap[nextVertex] = numEvaluated++;
        currentSearch.Remove(nextVertex);
        remaining.Remove(nextVertex);
        for (int const neighbor : isolates.Neighbors()[nextVertex]) {
            if (!evaluated[neighbor]) {
                currentSearch.Insert(neighbor);
            }
        }

////        if (currentSearch.Empty() && !remaining.Empty()) {
////            int const startVertex = *remaining.begin();
////            currentSearch.Insert(startVertex);
////            remaining.Remove(startVertex);
////            componentCount++;
////            vComponents.resize(componentCount);
////        }
    }

    vector<vector<int>> vAdjacencyArray(vertexRemap.size());

    for (pair<int,int> const &mapPair : vertexRemap) {
        int const oldVertex(mapPair.first);
        int const newVertex(mapPair.second);
        for (int const neighbor : isolates.Neighbors()[oldVertex]) {
            if (vertexRemap.find(neighbor) != vertexRemap.end()) {
////                cout << "vAdjacencyArray.size=" << vAdjacencyArray.size() << endl;
////                cout << "newVertex           =" << newVertex << endl; 
                vAdjacencyArray[newVertex].push_back(vertexRemap[neighbor]);
            }
        }
    }

    Isolates4<SparseArraySet> subgraphIsolates(vAdjacencyArray);
    vector<int> vRemoved;
    vector<int> vIsolates;
    set<int>    setRemoved;
////    vector<pair<int,int>> vAddedEdges;
////    subgraphIsolates.RemoveAllIsolates(0, vIsolates, vRemoved, vAddedEdges, true /* consider all vertices for reduction */);
    vector<Reduction> vReductions;
    subgraphIsolates.RemoveAllIsolates(0, vIsolates, vRemoved, vReductions, true /* consider all vertices for reduction */);

    if (subgraphIsolates.GetInGraph().Size() == vAdjacencyArray.size()) {
        return vAdjacencyArray;
    }

    vAdjacencyArray.clear();
    vAdjacencyArray.resize(subgraphIsolates.GetInGraph().Size());
    vertexRemap.clear();
    size_t uNewIndex = 0;

    for (int const vertex : subgraphIsolates.GetInGraph()) {
        vertexRemap[vertex] = uNewIndex++;
    }

    for (pair<int,int> const &mapPair : vertexRemap) {
        int const oldVertex(mapPair.first);
        int const newVertex(mapPair.second);
        for (int const neighbor : subgraphIsolates.Neighbors()[oldVertex]) {
            if (vertexRemap.find(neighbor) != vertexRemap.end()) {
////                cout << "vAdjacencyArray.size=" << vAdjacencyArray.size() << endl;
////                cout << "newVertex           =" << newVertex << endl; 
                vAdjacencyArray[newVertex].push_back(vertexRemap[neighbor]);
            }
        }
    }

    return vAdjacencyArray;
}

void PrintGraphInEdgesFormat(vector<vector<int>> const &adjacencyArray)
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

void PrintSubgraphInEdgesFormat(vector<vector<int>> const &adjacencyArray, map<int,int> const &vertexMap)
{
    cout << vertexMap.size() << endl;
    size_t edges(0);
    for (int index = 0; index < adjacencyArray.size(); ++index) {
        if (vertexMap.find(index) == vertexMap.end()) continue;
        for (int const neighbor : adjacencyArray[index]) {
            if (vertexMap.find(neighbor) != vertexMap.end()) edges++;
        }
    }
    cout << edges << endl;

    for (int index = 0; index < adjacencyArray.size(); ++index) {
        if (vertexMap.find(index) == vertexMap.end()) continue;
        for (int const neighbor : adjacencyArray[index]) {
            if (vertexMap.find(neighbor) == vertexMap.end()) continue;
            cout << vertexMap.at(index) << "," << vertexMap.at(neighbor) << endl;
        }
    }
}

void PrintSubgraphInEdgesFormat(Isolates4<SparseArraySet> const &isolates)
{
    map<int,int> vertexMap;
    cout << isolates.GetInGraph().Size() << endl;
    size_t edges(0);
    int newVertex(0);
    for (int const vertex : isolates.GetInGraph()) {
        vertexMap[vertex] = newVertex++;
        for (int const neighbor : isolates.Neighbors()[vertex]) {
            edges++;
        }
    }
    cout << edges << endl;

    for (int const vertex : isolates.GetInGraph()) {
        if (vertexMap.find(vertex) == vertexMap.end()) continue;
        for (int const neighbor : isolates.Neighbors()[vertex]) {
            if (vertexMap.find(neighbor) == vertexMap.end()) continue;
            cout << vertexMap.at(vertex) << "," << vertexMap.at(neighbor) << endl;
        }
    }
}


void GetVertexRemap(vector<int> const &vVertices, Isolates4<SparseArraySet> const &isolates, map<int,int> &vertexRemap)
{
////    map<int,int> vertexRemap;
////    map<int,int> reverseMap;
    int newVertex(0);
////    vector<vector<char>> componentMatrix(vVertices.size(), vector<char>());
////    vector<vector<int>> componentArray(vVertices.size(), vector<int>());
    for (int const vertex : vVertices) {
        vertexRemap[vertex] = newVertex++;
////        reverseMap[newVertex-1] = vertex;
////        componentMatrix[newVertex-1].resize(vVertices.size());
    }
}

void ComputeOnConnectedComponent(vector<int> const &vVertices, Isolates4<SparseArraySet> const &isolates, vector<vector<int>> const &adjacencyArray, list<int> &realClique, list<list<int>> cliques, bool const bSetCliqueSize)
{
    map<int,int> vertexRemap;
    map<int,int> reverseMap;
    int newVertex(0);
    vector<vector<char>> componentMatrix(vVertices.size(), vector<char>());
    vector<vector<int>> componentArray(vVertices.size(), vector<int>());
    for (int const vertex : vVertices) {
        vertexRemap[vertex] = newVertex++;
        reverseMap[newVertex-1] = vertex;
        componentMatrix[newVertex-1].resize(vVertices.size());
    }

    ////        cout << "Matrix resized to " << vVertices.size() << "x" << vVertices.size() << endl << flush;

    cout << "Graph, size=" << vVertices.size() << ":" << endl << flush;
    for (int const vertex : vVertices) {
        int const newVertex(vertexRemap[vertex]);
        componentMatrix[newVertex][newVertex] = 1;
        ////                    cout << newVertex << ":";
        for (int const neighbor : isolates.Neighbors()[vertex]) {
            if (vertexRemap.find(neighbor) != vertexRemap.end()) {
                int const newNeighbor(vertexRemap[neighbor]);
                ////                    cout << "Adding edge " << newVertex << "," << newNeighbor << endl << flush;
                componentMatrix[newVertex][newNeighbor] = 1;
                componentArray[newVertex].push_back(newNeighbor);
                ////                            cout << " " << newNeighbor;
            }
        }
        ////                    cout << endl << flush;
    }
    ////                cout << endl << flush;

////    PrintGraphInEdgesFormat(componentArray);
////    PrintSubgraphInEdgesFormat(adjacencyArray, vertexRemap);

    list<list<int>> componentResult;

    ////                cout << "Start algorithm: " << endl << flush;
    TesterMISS algorithm(componentMatrix, componentArray);
////    LightWeightFullMISS algorithm(componentMatrix);
////    algorithm.SetQuiet(true);
////    algorithm.SetQuiet(false);

    if (bSetCliqueSize) {
        if (cliques.back().size() >= realClique.size())
            algorithm.SetMaximumCliqueSize(cliques.back().size() - realClique.size());
////        algorithm.SetMaximumCliqueSize(0);
    }

    algorithm.Run(componentResult);

    if (!componentResult.empty() && !componentResult.back().empty()) {
        for (int const vertex : componentResult.back()) {
            realClique.push_back(reverseMap[vertex]);
        }
    }
    ////                cout << "    Component has independent set of size " << componentResult.back().size() << endl;
}

void ComputeOnConnectedComponent(vector<int> const &vVertices, Isolates4<SparseArraySet> const &isolates, vector<vector<int>> const &adjacencyArray, list<int> &realClique, list<list<int>> cliques, map<int,int> const &vertexRemap, bool const bSetCliqueSize)
{
////    map<int,int> vertexRemap;
    map<int,int> reverseMap;
    int newVertex(0);
    vector<vector<char>> componentMatrix(vertexRemap.size(), vector<char>());
    vector<vector<int>> componentArray(vertexRemap.size(), vector<int>());
    for (pair<int,int> const &mapPair : vertexRemap) {
        int const vertex(mapPair.first);
        int const newVertex(mapPair.second);
////        vertexRemap[vertex] = newVertex++;
        reverseMap[newVertex] = vertex;
        componentMatrix[newVertex].resize(vertexRemap.size());
    }

    ////        cout << "Matrix resized to " << vVertices.size() << "x" << vVertices.size() << endl << flush;

    cout << "Graph, size=" << vVertices.size() << ":" << endl << flush;
    for (int const vertex : vVertices) {
        int const newVertex(vertexRemap.find(vertex)->second);
        componentMatrix[newVertex][newVertex] = 1;
        ////                    cout << newVertex << ":";
        for (int const neighbor : isolates.Neighbors()[vertex]) {
            if (vertexRemap.find(neighbor) != vertexRemap.end()) {
                int const newNeighbor(vertexRemap.find(neighbor)->second);
                ////                    cout << "Adding edge " << newVertex << "," << newNeighbor << endl << flush;
                componentMatrix[newVertex][newNeighbor] = 1;
                componentArray[newVertex].push_back(newNeighbor);
                ////                            cout << " " << newNeighbor;
            }
        }
        ////                    cout << endl << flush;
    }
    ////                cout << endl << flush;

////    PrintGraphInEdgesFormat(componentArray);
////    PrintSubgraphInEdgesFormat(adjacencyArray, vertexRemap);

    list<list<int>> componentResult;

    ////                cout << "Start algorithm: " << endl << flush;
    TesterMISS algorithm(componentMatrix, componentArray);
////    algorithm.SetQuiet(true);
    algorithm.SetQuiet(false);

    if (bSetCliqueSize) {
        if (cliques.back().size() >= realClique.size())
            algorithm.SetMaximumCliqueSize(cliques.back().size() - realClique.size());
////        algorithm.SetMaximumCliqueSize(0);
    }

    algorithm.Run(componentResult);

    if (!componentResult.empty() && !componentResult.back().empty()) {
        for (int const vertex : componentResult.back()) {
            realClique.push_back(reverseMap[vertex]);
        }
    }
    ////                cout << "    Component has independent set of size " << componentResult.back().size() << endl;
}

void Staging::BuildSingleEdgeIsolateGraph(ArraySet &removed, ArraySet &isolatedSet)
{
    vector<int> vGraphVertices(m_AdjacencyList.size());
    for (int i = 0; i < vGraphVertices.size(); ++i) {
        vGraphVertices[i] = i;
    }

    sort (vGraphVertices.begin(), vGraphVertices.end(), [this](int const &left, int const &right) { return m_AdjacencyList[left].size() < m_AdjacencyList[right].size(); });

    cout << "min-degree=" << m_AdjacencyList[vGraphVertices.front()].size() << ", max-degree=" << m_AdjacencyList[vGraphVertices.back()].size() << endl;

    ArraySet adjacentToIsolated(m_AdjacencyList.size());

    ArraySet remaining(m_AdjacencyList.size());
    for (int const vertex : vGraphVertices) {
        remaining.Insert(vertex);
    }

    // select first edge to be removed = min degree with min degree neighbor.

    int const firstVertex(vGraphVertices[0]);
    int secondVertex(-1);
    size_t smallestDegree(ULONG_MAX);
    for (int const neighbor : m_AdjacencyList[firstVertex]) {
        if (m_AdjacencyList[neighbor].size() < smallestDegree) {
            secondVertex = neighbor;
            smallestDegree = m_AdjacencyList[neighbor].size();
        }
    }

    removed.Insert(firstVertex);
    removed.Insert(secondVertex);
    remaining.Remove(firstVertex);
    remaining.Remove(secondVertex);
    isolatedSet.Insert(firstVertex);
    for (int const neighbor : m_AdjacencyList[firstVertex]) {
        if (!removed.Contains(neighbor)) {
            adjacentToIsolated.Insert(neighbor);
        }
    }

    while (!remaining.Empty()) {
        cout << "Must remove " << adjacentToIsolated.Size() << " vertices to reduce " << removed.Size() << " vertices" << endl;

        int neighborToRemove(-1);
        int isolateToRemove(-1);
        size_t neighbors(ULONG_MAX);
////        for (int const vertex : adjacentToIsolated) {
        for (int const vertex : remaining) {
            bool goodVertex(true);
            int tempNeighborToRemove(-1);
            for (int const neighbor : m_AdjacencyList[vertex]) {
                if (removed.Contains(neighbor)) { goodVertex = false; break; } // can't remove as an isolate
                tempNeighborToRemove = neighbor;
            }

            if (goodVertex && neighbors > m_AdjacencyList[vertex].size()) {
                isolateToRemove = vertex;
                neighborToRemove = tempNeighborToRemove;
                neighbors = m_AdjacencyList[vertex].size();
            }
        }

        if (isolateToRemove != -1) { // found a neighbor we could remove...
            adjacentToIsolated.Remove(isolateToRemove);
            if (adjacentToIsolated.Contains(neighborToRemove)) {
                adjacentToIsolated.Remove(neighborToRemove);
            }

            removed.Insert(isolateToRemove);
            removed.Insert(neighborToRemove);

            remaining.Remove(neighborToRemove);
            remaining.Remove(isolateToRemove);
            isolatedSet.Insert(isolateToRemove);

            for (int const neighbor : m_AdjacencyList[isolateToRemove]) {
                if (!removed.Contains(neighbor) && !adjacentToIsolated.Contains(neighbor)) {
                    adjacentToIsolated.Insert(neighbor);
                }
            }
        } else { // no more to remove...
            break;
        }
////        else { // no neighbor we could remove, so remove a different one
////            bool goodVertex(true)
////            for (int const vertex : remaining) {
////                for (int const neighbor : m_AdjacencyList[vertex]) {
////                    if (remove.Contains(neighbor)) { goodVertex = false; break; }
////                }
////
////                if (goodVertex == true) {
////                }
////            }
////        }
    }
    cout << "Must remove " << adjacentToIsolated.Size() << " vertices to reduce " << removed.Size() << " vertices" << endl;
}


void Staging::Run()
{
    cout << "Running Staging..." << endl << flush;

    //CliqueTools::FindMaximalIndependentSetInCliqueGraph(m_AdjacencyList);

#if 0
    cout << "Applying Original Reductions..." << endl << flush;

    set<int> isolates;
    set<int> removed;
    set<int> remaining;
    set<int> inGraph;

    vector<set<int>> neighbors(m_AdjacencyList.size());
    for (size_t u = 0; u < m_AdjacencyList.size(); ++u) {
        neighbors[u].insert(m_AdjacencyList[u].begin(), m_AdjacencyList[u].end());
        remaining.insert(u);
        inGraph.insert(u);
    }

    vector<bool> vMarkedVertices(m_AdjacencyList.size(), false);

    cout << "Removing isolates..." << endl;
    RemoveAllIsolates(neighbors, remaining, isolates, removed, vMarkedVertices);
    for (int const vertex : removed) {
        inGraph.erase(vertex);
    }
    cout << "# vertices remaining in graph: " << inGraph.size() << "/" << m_AdjacencyList.size() << endl << flush;

    // as of now, this loop helps find a really small Independent set, this is to help limit the recursion depth.
    // need to expand this into a branch-and-bound, pick the vertices that constrain the search the most, to keep
    // the search at a reasonable depth. And try to "peel off" as many vertices as possible through reductions.
    while (removed.size() != m_AdjacencyList.size()) {
        // find vertex whose removal maximizes the number of vertices moved to IS by reduction.
        int vertexWithMaxReductions(-1);

        remaining = inGraph;

        int maxIsolates(-1);
        cout << "Testing remaining " << inGraph.size() << " vertices, to maximize isolate removal" << endl << flush;
        for (int const vertex : inGraph) {
            set<int> tempRemaining(remaining);
            set<int> tempRemoved(removed);
            set<int> tempIsolates(isolates);

            RemoveVertexAndNeighbors(vertex, neighbors, tempRemaining, tempIsolates, tempRemoved);
            RemoveAllIsolates(neighbors, tempRemaining, tempIsolates, tempRemoved, vMarkedVertices);
            if (static_cast<int>(tempIsolates.size() - isolates.size()) > maxIsolates) {
                maxIsolates = static_cast<int>(tempIsolates.size() - isolates.size());
                vertexWithMaxReductions = vertex;
            }

            vector<int> newlyRemovedVertices(tempRemoved.size() - removed.size(), -1);

            vector<int>::iterator it = set_difference(tempRemoved.begin(), tempRemoved.end(), removed.begin(), removed.end(), newlyRemovedVertices.begin());

            newlyRemovedVertices.resize(it - newlyRemovedVertices.begin());

            ReplaceRemovedVertices(m_AdjacencyList, neighbors, inGraph, newlyRemovedVertices);
        }

        cout << "Removing vertex " << vertexWithMaxReductions << " maximizes isolate removal (" << maxIsolates << " isolates)." << endl << flush;
        //cout << "Removed " << removed.size() << "/" << m_AdjacencyList.size() << " vertices" << endl << flush;

        cout << "Removing next vertex..." << endl << flush;
        // don't actually remove anything, just run like before. // TODO/DS actually remove vertex!
        set<int> tempRemoved(removed);
        if (vertexWithMaxReductions == -1) {
            //RemoveAllIsolates    (neighbors, remaining, isolates, tempRemoved, vMarkedVertices); // TODO: remove, does nothing
            RemoveMaxDegreeVertex(neighbors, remaining, isolates, tempRemoved, vMarkedVertices); // TODO: remove max degree instead?
        } else {
            RemoveVertexAndNeighbors(vertexWithMaxReductions, neighbors, remaining, isolates, tempRemoved);
            RemoveAllIsolates(neighbors, remaining, isolates, tempRemoved, vMarkedVertices);
        }

        vector<int> newlyRemovedVertices(tempRemoved.size() - removed.size(), -1);
        vector<int>::iterator it = set_difference(tempRemoved.begin(), tempRemoved.end(), removed.begin(), removed.end(), newlyRemovedVertices.begin()); // TODO: This is an expensive diff, just grab the added vertices before we get to this point.
        newlyRemovedVertices.resize(it - newlyRemovedVertices.begin());

        for (int const removedVertex : newlyRemovedVertices) {
            //cout << "Removing " << removedVertex << " from graph " << endl;
            inGraph.erase(removedVertex);
        }
        //        inGraph.erase(tempRemoved.begin(), tempRemoved.end());
        removed.insert(newlyRemovedVertices.begin(), newlyRemovedVertices.end());

        cout << "# vertices remaining in graph: " << inGraph.size() << "/" << m_AdjacencyList.size() << endl << flush;
    }

    cout << "Removed " << removed.size() << " vertices." << endl << flush;
    cout << "Found independent set of size: " << isolates.size() << endl << flush;
#else
    cerr << "Applying New      Reductions..." << endl << flush;

    Isolates4<SparseArraySet> isolates(m_AdjacencyList);

    cerr << "Removing isolates..." << endl;
    vector<int> vRemoved;
    vector<int> vIsolates;
    set<int>    setRemoved;
////    vector<pair<int,int>> vAddedEdges;
    vector<Reduction> vReductions;
    isolates.RemoveAllIsolates(0, vIsolates, vRemoved, vReductions, true /* consider all vertices for reduction */);
    cout << "Removed " << vIsolates.size() << " isolates, graph has " << isolates.GetInGraph().Size() << "/" << m_AdjacencyList.size() << " vertices remaining" << endl;

    cout << "Number to add to final MIS: " << vIsolates.size() + isolates.GetFoldedVertexCount() << endl;

    PrintSubgraphInEdgesFormat(isolates);

#if 0
    vRemoved.insert(vRemoved.end(), vIsolates.begin(), vIsolates.end());
    setRemoved.insert(vRemoved.begin(), vRemoved.end());

    vector<vector<int>> vComponents;

    cerr << "# vertices remaining in graph: " << m_AdjacencyList.size() - vRemoved.size() << "/" << m_AdjacencyList.size() << endl << flush;
    cerr << "# connected components       : " << ComputeConnectedComponents(isolates, vComponents) << endl << flush;
    cerr << "size of connected components : ";
    cout << "[ ";

    size_t maxComponentSize(0);
    size_t maxComponentIndex(0);

    for (size_t index = 0; index < vComponents.size(); ++index) {
        vector<int> const& vComponent(vComponents[index]);
        cout << vComponent.size() << " ";
        if (vComponent.size() > maxComponentSize) {
            maxComponentSize = vComponent.size();
            maxComponentIndex = index;
        }
    }

    cout << "]" << endl << flush;

    cout << "Max component size=" << maxComponentSize << endl;

    isolates.SetConnectedComponent(vComponents[maxComponentIndex]);

    vector<vector<int>> const subgraphAdjacencyList = ComputeSubgraphOfSize(isolates, m_AdjacencyList.size(), 3000);

    cout << "Subgraph size=" << subgraphAdjacencyList.size() << endl << flush;

    vector<vector<char>> subgraphAdjacencyMatrix(subgraphAdjacencyList.size());
    for (size_t index = 0; index < subgraphAdjacencyList.size(); ++index) {
        subgraphAdjacencyMatrix[index].resize(subgraphAdjacencyList.size(), 0);
        for (int const neighbor : subgraphAdjacencyList[index]) {
            subgraphAdjacencyMatrix[index][neighbor] = 1;
        }
    }

    list<list<int>> cliques;
////    LightWeightReductionStaticOrderMISS algorithm(subgraphAdjacencyMatrix, subgraphAdjacencyList);
////    LightWeightReductionFullMISS algorithm(subgraphAdjacencyMatrix, subgraphAdjacencyList);
////    LightWeightFullMISS algorithm(subgraphAdjacencyMatrix);
////    LightWeightReductionSparseFullMISS algorithm(subgraphAdjacencyList);
////    LightWeightReductionSparseStaticOrderMISS algorithm(subgraphAdjacencyList);
////    LightWeightStaticOrderMISS algorithm(subgraphAdjacencyMatrix); ////, subgraphAdjacencyList);

    cout << subgraphAdjacencyList.size() << endl;
    size_t edges(0);
    for (vector<int> const &neighborList : subgraphAdjacencyList) {
        edges+= neighborList.size();
    }
    cout << edges << endl;

    for (size_t index = 0; index < subgraphAdjacencyList.size(); ++index) {
        for (int const neighbor : subgraphAdjacencyList[index]) {
            cout << index << "," << neighbor << endl;
        }
    }

    auto printCliqueSize = [](list<int> const &clique) {
        cout << "Found clique of size " << clique.size() << endl << flush;
    };
////    algorithm.AddCallBack(printCliqueSize);
////    algorithm.Run(cliques);

#if 0
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
            for (int const neighbor : m_AdjacencyList[vertex]) {
                if (vertexRemap.find(neighbor) != vertexRemap.end()) {
                    int const newNeighbor(vertexRemap[neighbor]);
////                    cout << "Adding edge " << newVertex << "," << newNeighbor << endl << flush;
                    adjacencyMatrix[newVertex][newNeighbor] = 1;
                    adjacencyArray[newVertex].push_back(newNeighbor);
                }
            }
        }

        list<list<int>> result;

////        LightWeightReductionMISR algorithm(adjacencyMatrix, adjacencyArray);
////        LightWeightStaticOrderMISS algorithm(adjacencyMatrix);
////        LightWeightMISQ algorithm(adjacencyMatrix);
        LightWeightReductionStaticOrderMISS algorithm(adjacencyMatrix, adjacencyArray);
        algorithm.Run(result);
        if (!result.empty())
            for (int const vertex : result.back()) {
                realClique.push_back(reverseMap[vertex]);
            }
    }

    cout << "Total clique size: " << realClique.size() << endl << flush;
#endif // 0
#else
#if 0
    size_t numVertices(0);
    size_t numEdges(0);

    map<int,int> vertexRemap;
    for (int const vertex : isolates.GetInGraph()) {
        if (!isolates.Neighbors()[vertex].Empty()) {
            numEdges += isolates.Neighbors()[vertex].Size();
            vertexRemap[vertex] = numVertices;
            numVertices++;
        }
    }

    cout << numVertices << endl;
    cout << numEdges << endl;

    for (int const vertex : isolates.GetInGraph()) {
        for (int const neighbor : isolates.Neighbors()[vertex]) {
            cout << vertexRemap[vertex] << "," << vertexRemap[neighbor] << endl;
        }
    }
#endif // 0
#endif // 0

#if 0
    vector<vector<char>> adjacencyMatrix(m_AdjacencyList.size());
    for (size_t index = 0; index < m_AdjacencyList.size(); ++index) {
        adjacencyMatrix[index].resize(m_AdjacencyList.size(), 0);
        for (int const neighbor : m_AdjacencyList[index]) {
            adjacencyMatrix[index][neighbor] = 1;
        }
    }

    // is there some small set of vertices, such that if I remove them and their neighbors, the connected components all have
    // small size (don't really care about balanced, could be 500 200 1 1 1 1 1 1 or 500 400 1 1 1 1.

    size_t best(string::npos);
    size_t const size(adjacencyMatrix.size());
    size_t loops(0);

    vector<int> vVertices;
    for (int i = 0; i < m_AdjacencyList.size(); ++i) {
        vVertices.push_back(i);
    }

    vector<int> vNumEdgesInTwoNeighborhood(m_AdjacencyList.size(), 0);
    for (int vertex = 0; vertex < m_AdjacencyList.size(); ++vertex) {
        for (int const neighbor : m_AdjacencyList[vertex]) {
            vNumEdgesInTwoNeighborhood[vertex] += m_AdjacencyList[neighbor].size();
        }
    }


////    sort (vVertices.begin(), vVertices.end(), [this](int const left, int const right) { return m_AdjacencyList[left].size() > m_AdjacencyList[right].size(); });
    sort (vVertices.begin(), vVertices.end(), [&vNumEdgesInTwoNeighborhood](int const left, int const right) { return vNumEdgesInTwoNeighborhood[left] > vNumEdgesInTwoNeighborhood[right]; });

    std::random_device generator;
    std::uniform_int_distribution<int> distribution1(0,(int)(size * 0.01));
    std::uniform_int_distribution<int> distribution5(0,(int)(size * 0.05));
    std::uniform_int_distribution<int> distribution10(0,(int)(size * 0.10));
    std::uniform_int_distribution<int> distribution20(0,(int)(size * 0.20));
    std::uniform_int_distribution<int> distribution30(0,(int)(size * 0.30));
    std::uniform_int_distribution<int> distribution50(0,(int)(size * 0.50));
////    int dice_roll = distribution(generator);  // generates number in the range 1..6

#endif // 0

    vector<pair<int,int>> vAddedEdgesUnused;

#if 0
    for (int initialIndex = 0; initialIndex < vVertices.size(); ++initialIndex) {
        vector<int> vInitialRemoved;
        int const initialVertex(vVertices[initialIndex]);
        vInitialRemoved.push_back(initialVertex);
        isolates.RemoveVertexAndNeighbors(initialVertex, vInitialRemoved);
        int innerLoops(0);
        while (innerLoops < 200) {
            vector<int> vRemoved;
            for (int i=0; i < 2; ++i) {
                int vertexToRemove(-1);
                int const generated(distribution5(generator));
////                int const firstGenerated(initialIndex >= generated5 ? (initialIndex - generated5) : generated5 );
                int const firstGenerated(generated);
                switch (i) {
                    case 0: vertexToRemove = vVertices[firstGenerated]; break;
                    case 1: vertexToRemove = vVertices[firstGenerated + distribution1(generator)]; break;
                    case 2: vertexToRemove = vVertices[distribution10(generator)]; break;
                    case 3: vertexToRemove = vVertices[distribution20(generator)]; break;
                };
                if (!isolates.GetInGraph().Contains(vertexToRemove)) {
                    i--;
                    continue;
                }
                vRemoved.push_back(vertexToRemove);
                isolates.RemoveVertexAndNeighbors(vertexToRemove, vRemoved);
                isolates.RemoveVertex(vertexToRemove);
                isolates.RemoveAllIsolates(0, vRemoved, vRemoved, vAddedEdgesUnused, false);
            }

            vector<vector<int>> vComponents;
            ComputeConnectedComponents(isolates, vComponents);
            size_t biggest(0);
            for (vector<int> const &vComponent : vComponents) {
                biggest = max(vComponent.size(), biggest);
            }
            if (biggest < best) {
                best = biggest;
                cout << "best so far: components=" << vComponents.size() << ", max-component-size=" << best << endl;
                innerLoops = 0; // this might be a good vertex to keep considering...
            }
            innerLoops++;
            isolates.ReplaceAllRemoved(vRemoved);
        }

        isolates.ReplaceAllRemoved(vInitialRemoved);
    }
#endif //0
#if 0
    int savedi(0), savedj(0), savedk(0);
#if 1
#if 0
    for (size_t i = 0; i < size; i++) {
        vector<int> v_iNeighbors;
        isolates.RemoveVertexAndNeighbors(vVertices[i], v_iNeighbors);
        isolates.RemoveAllIsolates(0, v_iNeighbors, v_iNeighbors, vAddedEdgesUnused, false);

        vector<vector<int>> vComponents;
        ComputeConnectedComponents(isolates, vComponents);
        size_t biggest(0);
        for (vector<int> const &vComponent : vComponents) {
            biggest = max(vComponent.size(), biggest);
        }

        if (i % int(size*0.1) == 0) {
            cout << "finished evaluating " << (i*100.0/size) << "%" << endl;
        }

        if (biggest < best) {
            best = biggest;
            cout << "best so far: remove " << i << " for max component of size " << best << endl;
            cout << "best so far: components=" << vComponents.size() << ", max-component-size=" << best << endl;
            savedi = i;
        }

        v_iNeighbors.push_back(i);
        isolates.ReplaceAllRemoved(v_iNeighbors);
////        if (loops % 5000 == 0) { break; }
////            cout << "Removing " << savedi << " " << savedj << " " << savedk << endl;
////            vector<int> vRemovedNeighbors;
////            isolates.RemoveVertexAndNeighbors(savedi, vRemovedNeighbors);
////            isolates.RemoveVertexAndNeighbors(savedj, vRemovedNeighbors);
////            isolates.RemoveVertexAndNeighbors(savedk, vRemovedNeighbors);
////        }
    }
#endif // 0

    int lastPercentage = 0;
#if 0
////    for (size_t i = 0; i < size*0.01; i++) {
    for (size_t i = 0; i < size*0.05; i++) {
        vector<int> v_iNeighbors;
        isolates.RemoveVertexAndNeighbors(vVertices[i], v_iNeighbors);
////        for (size_t j = i+1; j < i + size*0.10; j++) {
        for (size_t j = i+1; j < i + size*0.5; j++) {
            vector<int> v_jNeighbors;
            isolates.RemoveVertexAndNeighbors(vVertices[j], v_jNeighbors);
            isolates.RemoveAllIsolates(0, v_jNeighbors, v_jNeighbors, vAddedEdgesUnused, false);

                vector<vector<int>> vComponents;
                ComputeConnectedComponents(isolates, vComponents);
                size_t biggest(0);
                for (vector<int> const &vComponent : vComponents) {
                    biggest = max(vComponent.size(), biggest);
                }

////                int const percentage((int((i*size*0.10 + j)*100.0/((size*0.01)*(size*0.10)))));
////                int const percentage((int((i*size + j)*100.0/((size)*(size)))));
                int const percentage((int((i*size + j)*100.0/((size*0.05)*(size*0.5)))));
                if (percentage > lastPercentage) {
                    cout << "finished evaluating " << percentage << "%" << endl;
                    lastPercentage = percentage;
                }

                if (biggest < best) {
                    best = biggest;
                    cout << "best so far: remove " << i << " " << j << " for max component of size " << best << endl;
                    cout << "best so far: components=" << vComponents.size() << ", max-component-size=" << best << endl;
                    savedi = i;
                    savedj = j;
                }
                v_jNeighbors.push_back(j);
                isolates.ReplaceAllRemoved(v_jNeighbors);
                ////            if (loops % 100 == 0) break;
        }
        v_iNeighbors.push_back(i);
        isolates.ReplaceAllRemoved(v_iNeighbors);
    }
#endif // 0
#endif // 0
////        if (loops % 5000 == 0) { break; }
////            cout << "Removing " << savedi << " " << savedj << " " << savedk << endl;
////            vector<int> vRemovedNeighbors;
////            isolates.RemoveVertexAndNeighbors(savedi, vRemovedNeighbors);
////            isolates.RemoveVertexAndNeighbors(savedj, vRemovedNeighbors);
////            isolates.RemoveVertexAndNeighbors(savedk, vRemovedNeighbors);
////        }
////    }

    lastPercentage = 0;
////    for (size_t i = 0; i < size; i++) {
    for (size_t i = 0; i < size*0.10; i++) {
        vector<int> v_iNeighbors;
        isolates.RemoveVertexAndNeighbors(vVertices[i], v_iNeighbors);
        isolates.RemoveAllIsolates(0, v_iNeighbors, v_iNeighbors, vAddedEdgesUnused, false);

        vector<int> vLeftOver1(isolates.GetInGraph().begin(), isolates.GetInGraph().end());
        vector<int> vTwoNeighborhood1(m_AdjacencyList.size(), 0);
        for (int const vertex : vLeftOver1) {
            for (int const neighbor : isolates.Neighbors()[vertex]) {
                vTwoNeighborhood1[vertex] += isolates.Neighbors()[neighbor].Size();
            }
        }

        sort (vLeftOver1.begin(), vLeftOver1.end(), [&vTwoNeighborhood1](int const left, int const right) { return vTwoNeighborhood1[left] > vTwoNeighborhood1[right]; });
////        for (size_t j = i+1; j < i + size*0.01; j++) {
        for (size_t j = 0; j < min(int(vLeftOver1.size()*0.01), 100); j++) {
            vector<int> v_jNeighbors;
            int const secondVertex(vLeftOver1[j]);
            if (isolates.GetInGraph().Contains(secondVertex)) {
                v_jNeighbors.push_back(secondVertex);
                isolates.RemoveVertexAndNeighbors(secondVertex, v_jNeighbors);
                isolates.RemoveAllIsolates(0, v_jNeighbors, v_jNeighbors, vAddedEdgesUnused, false);
            }

            vector<int> vLeftOver2(isolates.GetInGraph().begin(), isolates.GetInGraph().end());
            vector<int> vTwoNeighborhood2(m_AdjacencyList.size(), 0);
            for (int const vertex : vLeftOver2) {
                for (int const neighbor : isolates.Neighbors()[vertex]) {
                    vTwoNeighborhood2[vertex] += isolates.Neighbors()[neighbor].Size();
                }
            }

            sort (vLeftOver2.begin(), vLeftOver2.end(), [&vTwoNeighborhood2](int const left, int const right) { return vTwoNeighborhood2[left] > vTwoNeighborhood2[right]; });
////            for (size_t k = j+1; k < j+ size*0.05; k++) {
            for (size_t k = 0; k < min(int(vLeftOver2.size()*0.01), 100); k++) {
                loops++;
////                if (loops %10 == 0) { break; }
                int const thirdVertex(vLeftOver2[k]);
                vector<int> v_kNeighbors;
                if (isolates.GetInGraph().Contains(thirdVertex)) {
                    v_kNeighbors.push_back(thirdVertex);
                    isolates.RemoveVertexAndNeighbors(thirdVertex, v_kNeighbors);
                    isolates.RemoveAllIsolates(0, v_kNeighbors, v_kNeighbors, vAddedEdgesUnused, false);
                }

                vector<vector<int>> vComponents;
                ComputeConnectedComponents(isolates, vComponents);
                size_t biggest(0);
                for (vector<int> const &vComponent : vComponents) {
                    biggest = max(vComponent.size(), biggest);
                }


                if (biggest < best) {
                    best = biggest;
                    cout << "best so far: remove " << vVertices[i] << " " << vLeftOver1[j] << " " << vLeftOver2[k] << " for max component of size " << best << endl;
                    cout << "best so far: components=" << vComponents.size() << ", max-component-size=" << best << endl;
                    savedi = i;
                    savedj = j;
                    savedk = k;

                    if (best == 1470) {
                        PrintGraphInEdgesFormat(ComputeSubgraphOfSize(isolates, m_AdjacencyList.size(), 3000));
                    }
                }
                isolates.ReplaceAllRemoved(v_kNeighbors);
            }
            isolates.ReplaceAllRemoved(v_jNeighbors);
////            if (loops % 100 == 0) break;
        }
        v_iNeighbors.push_back(vVertices[i]);
        isolates.ReplaceAllRemoved(v_iNeighbors);

////        int const percentage((int((i*size*0.15 + j*0.10 + k)*100.0/((size*0.01)*(size*0.05)*(size*0.10)))));
////        if (percentage > lastPercentage + 1) {
            cout << "finished evaluating " << i << "/" << size*0.10 << " first vertices" << endl;
////            lastPercentage = percentage;
////        }
        ////        if (loops % 5000 == 0) { break; }
////            cout << "Removing " << savedi << " " << savedj << " " << savedk << endl;
////            vector<int> vRemovedNeighbors;
////            isolates.RemoveVertexAndNeighbors(savedi, vRemovedNeighbors);
////            isolates.RemoveVertexAndNeighbors(savedj, vRemovedNeighbors);
////            isolates.RemoveVertexAndNeighbors(savedk, vRemovedNeighbors);
////        }
    }
#endif // 0

#if 0
    vector<int> const criticalSet = CliqueTools::ComputeMaximumCriticalIndependentSet(m_AdjacencyList);
    cout << "Critical set (" << criticalSet.size() << " elements):" << endl;
    vector<int> vRemovedUnused;
    for (int const vertex : criticalSet) {
        cout << vertex << " ";
        isolates.RemoveVertexAndNeighbors(vertex, vRemovedUnused);
    }
    cout << endl;
    cout << "After Removing Critical Set, graph has " << isolates.GetInGraph().Size() << " vertices remaining." << endl;

    set<int> const setVertices(isolates.GetInGraph().begin(), isolates.GetInGraph().end());

    map<int,int> mapUnused;

    vector<vector<int>> subgraphAdjacencyList;
    GraphTools::ComputeInducedSubgraph(m_AdjacencyList, setVertices, subgraphAdjacencyList, mapUnused);

    cout << "Subgraph size=" << subgraphAdjacencyList.size() << endl << flush;

    vector<vector<char>> subgraphAdjacencyMatrix(subgraphAdjacencyList.size());
    for (size_t index = 0; index < subgraphAdjacencyList.size(); ++index) {
        subgraphAdjacencyMatrix[index].resize(subgraphAdjacencyList.size(), 0);
        for (int const neighbor : subgraphAdjacencyList[index]) {
            subgraphAdjacencyMatrix[index][neighbor] = 1;
        }
    }

    TesterMISS algorithm(subgraphAdjacencyMatrix, subgraphAdjacencyList);
////    LightWeightFullMISS algorithm(subgraphAdjacencyMatrix);
    list<list<int>> cliques;
    algorithm.Run(cliques);
    cout << "Actual largest independentSet: " << (cliques.back().size() + vIsolates.size() + criticalSet.size()) << endl << flush;
#endif // 0

#if 0
////    vector<vector<int>> const biDoubleArray(GraphTools::ComputeBiDoubleGraph(m_AdjacencyList));
////    PrintGraphInEdgesFormat(biDoubleArray);

    vector<vector<int>> const biDoubleArray(m_AdjacencyList);
    PrintGraphInEdgesFormat(biDoubleArray);

    vector<vector<char>> biDoubleMatrix(biDoubleArray.size());
    for (size_t index = 0; index < biDoubleArray.size(); ++index) {
        biDoubleMatrix[index].resize(biDoubleArray.size(), 0);
        for (int const neighbor : biDoubleArray[index]) {
            biDoubleMatrix[index][neighbor] = 1;
        }
    }

    TesterMISS algorithm(biDoubleMatrix, biDoubleArray);
    list<list<int>> cliques;
    algorithm.Run(cliques);

    int const minVertexCoverSize(GraphTools::ComputeMaximumMatchingSize(biDoubleArray));

    cout << "Flow-computed min vertex cover size : " << minVertexCoverSize << endl;
    cout << "Flow-computed independent set size  : " << (biDoubleArray.size() - minVertexCoverSize) << endl;
#endif //0

#if 0
    // try removing high-impact reduction vertices at first, then run regular algorithm on resulting connected
    // components.
////    vector<pair<int,int>> vAddedEdgesUnused;
////    vector<int> vRemoved1;
    vector<int> vCliqueVerticesPersistent1;
    vector<int> vRemovedPersistent1;
    list<list<int>> cliques;
    cliques.push_back(list<int>());


    auto ChooseNextVertex = [&vAddedEdgesUnused](Isolates4<SparseArraySet> &theIsolates, bool const firstLevel)
    {
////        cout << "Before choosing vertex: graph contains " << theIsolates.GetInGraph().Size() << " elements" << endl;
        size_t bestComponentSize(string::npos);
        int    bestVertex(-1);
        vector<int> vToConsider(theIsolates.GetInGraph().begin(), theIsolates.GetInGraph().end());
        if (!firstLevel) {
            for (int const i : vToConsider) {
                vector<int> v_iNeighbors;
                theIsolates.RemoveVertexAndNeighbors(i, v_iNeighbors);
                theIsolates.RemoveAllIsolates(0, v_iNeighbors, v_iNeighbors, vAddedEdgesUnused, false);

                vector<vector<int>> vComponents;
                ComputeConnectedComponents(theIsolates, vComponents);
                size_t uSizeOfLargestComponent(0);
                for (vector<int> const &vComponent : vComponents) {
                    uSizeOfLargestComponent = max(vComponent.size(), uSizeOfLargestComponent);
                }

////                if (i == 58)
////                    cout << "Can remove " << i << " for max component of size " << uSizeOfLargestComponent << endl;

                if (uSizeOfLargestComponent < bestComponentSize) {
                    bestComponentSize = uSizeOfLargestComponent;
////                    cout << "best so far: remove " << i << " for " << vComponents.size() << " components, with max component of size " << bestComponentSize << endl;
                    bestVertex = i;
                }

                v_iNeighbors.push_back(i);
                theIsolates.ReplaceAllRemoved(v_iNeighbors);
            }
        } else {
            sort (vToConsider.begin(), vToConsider.end(), [&theIsolates](int const left, int const right) { return theIsolates.Neighbors()[left].Size() > theIsolates.Neighbors()[right].Size(); });

            for (size_t i = 0; i < vToConsider.size() /**0.20*/; i++) {
                vector<int> v_iNeighbors;
                theIsolates.RemoveVertexAndNeighbors(vToConsider[i], v_iNeighbors);
                theIsolates.RemoveAllIsolates(0, v_iNeighbors, v_iNeighbors, vAddedEdgesUnused, false);
                v_iNeighbors.push_back(vToConsider[i]);
////                if (i == 46)
                size_t savedSize(theIsolates.GetInGraph().Size());
////                    cout << "After removing " << i << " there are " << theIsolates.GetInGraph().Size() << " elements remaining" << endl;
                for (size_t j = i+1; j < vToConsider.size()*0.10; j++) {
                    if (!theIsolates.GetInGraph().Contains(vToConsider[j])) continue;
                    if (savedSize != theIsolates.GetInGraph().Size()) {
                        cout << "Consistency Error: Before removing " << vToConsider[i] << ", " << vToConsider[j] << " there are " << theIsolates.GetInGraph().Size() << "!=" << savedSize << " elements remaining" << endl;
                    }
////                    cout << "Before removing " << i << ", " << j << " there are " << theIsolates.GetInGraph().Size() << " elements remaining" << endl;
////                    if (i == 0 && j == 72) {
////                        cout << "Before removing" << vToConsider[i] << ", " << vToConsider[j] << ":    " << endl;
////                        set<int> inGraph(theIsolates.GetInGraph().begin(), theIsolates.GetInGraph().end());
////                        for (int const neighbor : inGraph) {
////                            cout << neighbor << " ";
////                        }
////                        cout << endl;
////                    }

////                    if (i == 46 && j == 58)
                    vector<int> v_jNeighbors;
                    theIsolates.RemoveVertexAndNeighbors(vToConsider[j], v_jNeighbors);
                    theIsolates.RemoveAllIsolates(0, v_jNeighbors, v_jNeighbors, vAddedEdgesUnused, false);
                    v_jNeighbors.push_back(vToConsider[j]);
////                    if (i == 0 && j == 72) {
////                        cout << "119" << ((find(v_jNeighbors.begin(), v_jNeighbors.end(), 119) != v_jNeighbors.end())?" is":" is not") << " removed during reduction " << endl;
////                    }
////                    cout << "Removed: " << endl;
////                    for (int const neighbor : v_jNeighbors) {
////                        cout << neighbor << " ";
////                    }
////                    cout << endl;
////                    if (i == 46 && j == 58)
////                        cout << "After removing " << i << ", " << j << " there are " << theIsolates.GetInGraph().Size() << " elements remaining" << endl;

                    vector<vector<int>> vComponents;
                    ComputeConnectedComponents(theIsolates, vComponents);
                    size_t uSizeOfLargestComponent(0);
                    for (vector<int> const &vComponent : vComponents) {
                        uSizeOfLargestComponent = max(vComponent.size(), uSizeOfLargestComponent);
                    }

                    if (uSizeOfLargestComponent < bestComponentSize) {
                        bestComponentSize = uSizeOfLargestComponent;
                        cout << "best so far: remove " << vToConsider[i] << " " << vToConsider[j] << " for " << vComponents.size() << " components, with max component of size " << bestComponentSize << endl;
                        bestVertex = vToConsider[i];
                    }
                    theIsolates.ReplaceAllRemoved(v_jNeighbors);

                    if (savedSize != theIsolates.GetInGraph().Size()) {
                        cout << "After  removing" << i << ", " << j << ":    " << endl;
                        set<int> inGraphAfter(theIsolates.GetInGraph().begin(), theIsolates.GetInGraph().end());
                        for (int const neighbor : inGraphAfter) {
                            cout << neighbor << " ";
                        }
                        cout << endl;
                        cout << "Consistency Error: After  removing " << i << ", " << j << " there are " << theIsolates.GetInGraph().Size() << "!=" << savedSize << " elements remaining" << endl;
                    }
                    ////            if (loops % 100 == 0) break;
                }
                theIsolates.ReplaceAllRemoved(v_iNeighbors);
            }

////            for (size_t i = 0; i < vToConsider.size()*0.01; i++) {
////                vector<int> v_iNeighbors;
////                theIsolates.RemoveVertexAndNeighbors(vToConsider[i], v_iNeighbors);
////                for (size_t j = i+1; j < vToConsider.size()*0.01; j++) {
////                    vector<int> v_jNeighbors;
////                    theIsolates.RemoveVertexAndNeighbors(vToConsider[j], v_jNeighbors);
////                    for (size_t k = j+1; k < vToConsider.size()*0.01; k++) {
////                        ////                if (loops %10 == 0) { break; }
////                        vector<int> v_kNeighbors;
////                        theIsolates.RemoveVertexAndNeighbors(vToConsider[k], v_kNeighbors);
////                        theIsolates.RemoveAllIsolates(0, v_kNeighbors, v_kNeighbors, vAddedEdgesUnused, false);
////
////                        vector<vector<int>> vComponents;
////                        ComputeConnectedComponents(theIsolates, vComponents);
////                        size_t uSizeOfLargestComponent(0);
////                        for (vector<int> const &vComponent : vComponents) {
////                            uSizeOfLargestComponent = max(vComponent.size(), uSizeOfLargestComponent);
////                        }
////
////                        if (uSizeOfLargestComponent < bestComponentSize) {
////                            bestComponentSize = uSizeOfLargestComponent;
////                            cout << "best so far: remove " << vToConsider[i] << " " << vToConsider[j] << " " << vToConsider[k] << " for " << vComponents.size() << " components, max component of size " << bestComponentSize << endl;
////                            bestVertex = vToConsider[i];
////                        }
////                        theIsolates.ReplaceAllRemoved(v_kNeighbors);
////                    }
////                    v_jNeighbors.push_back(vToConsider[j]);
////                    theIsolates.ReplaceAllRemoved(v_jNeighbors);
////                    ////            if (loops % 100 == 0) break;
////                }
////                v_iNeighbors.push_back(vToConsider[i]);
////                theIsolates.ReplaceAllRemoved(v_iNeighbors);
////                ////        if (loops % 5000 == 0) { break; }
////                ////            cout << "Removing " << savedi << " " << savedj << " " << savedk << endl;
////                ////            vector<int> vRemovedNeighbors;
////                ////            isolates.RemoveVertexAndNeighbors(savedi, vRemovedNeighbors);
////                ////            isolates.RemoveVertexAndNeighbors(savedj, vRemovedNeighbors);
////                ////            isolates.RemoveVertexAndNeighbors(savedk, vRemovedNeighbors);
////                ////        }
////            }
            }

////        cout << "After  choosing vertex: graph contains " << theIsolates.GetInGraph().Size() << " elements" << endl;
        ////        cout << "best vertex: " << bestVertex  << " for max component of size " << bestComponentSize << endl;
        return bestVertex;
    };

    IndependentSetColoringStrategy coloringStrategy(adjacencyMatrix);

    auto NewChooseNextVertex = [&adjacencyMatrix, &vAddedEdgesUnused, &coloringStrategy](
            Isolates4<SparseArraySet> &theIsolates,
            vector<int> const &vAdjunctOrdering,
            size_t const uCurrentCliqueSize,
            size_t const uMaximumCliqueSize,
            int const nextVertex)
    {

        cout << "Test:    Max clique size=" << uMaximumCliqueSize << endl;
        int vertexToChoose(-1);
        size_t minNumLeft(string::npos);
        vector<int> vNewAdjunctOrdering;
        vector<int> vNewOrdering;
        vector<int> vNewColoring;
        vector<int> vToConsider(theIsolates.GetInGraph().begin(), theIsolates.GetInGraph().end());
        for (int const vertex : vToConsider) {
            vector<int> vCliqueVertices;
            vector<int> vRemoved;
            theIsolates.RemoveVertex(vertex);
////            theIsolates.RemoveVertexAndNeighbors(vertex, vRemoved);
            vRemoved.push_back(vertex);
            theIsolates.RemoveAllIsolates(0, vCliqueVertices, vRemoved, vAddedEdgesUnused, false);

            size_t uNewSize(0);
            vNewAdjunctOrdering.resize(vAdjunctOrdering.size());
            for (size_t index = 0; index < vAdjunctOrdering.size(); ++index) {
                int const vertexInOrder(vAdjunctOrdering[index]);
                if (theIsolates.GetInGraph().Contains(vertexInOrder)) {
                    vNewAdjunctOrdering[uNewSize++] = vertexInOrder;
                }
            }

            vNewAdjunctOrdering.resize(uNewSize);
            vNewOrdering.resize(vNewAdjunctOrdering.size());
            vNewColoring.resize(vNewAdjunctOrdering.size());

////            if (vertex == nextVertex) {
////                cout << "Test   Ordering: ";
////                for (int const vertexInOrder : vNewAdjunctOrdering) {
////                    cout << vertexInOrder << " ";
////                }
////                cout << endl;
////            }

            if (vertex == nextVertex) { ////5 || vertex == 802 || vertex == 45 || vertex == 33) {
                cout << "Test:    Vertex " << vertex << " has " << uCurrentCliqueSize << " + " <<  vCliqueVertices.size() << " clique vertices " << endl;
            }
            coloringStrategy.Recolor(adjacencyMatrix, vNewAdjunctOrdering, vNewOrdering, vNewColoring, uCurrentCliqueSize + vCliqueVertices.size(), uMaximumCliqueSize);

            size_t numLeft = vNewColoring.size();
            for (; numLeft > 0; --numLeft) {
                if (uCurrentCliqueSize + vCliqueVertices.size() + vNewColoring[numLeft-1] <= uMaximumCliqueSize) { break; }
            }
            numLeft = vNewColoring.size() - numLeft;
            if (vertex == nextVertex) {
                cout << "vertex " << vertex << ": P.left = " << numLeft << endl;
            }
            if (numLeft < minNumLeft) {
                vertexToChoose = vertex;
                minNumLeft = numLeft;
////                cout << ", P.left = " << numLeft << endl;
            }

            theIsolates.ReplaceAllRemoved(vCliqueVertices);
            theIsolates.ReplaceAllRemoved(vRemoved);
        }

        cout << "Choosing next vertex " << vertexToChoose << " will minimize the number of vertices for consideration to " << minNumLeft << endl;

        return vertexToChoose;
    };


    vector<int> vOrdering;
    vector<int> vColoring;

    size_t cliqueSize(0);
////    OrderingTools::InitialOrderingReduction(isolates, vOrdering, vColoring);
    OrderingTools::InitialOrderingMISR(m_AdjacencyList, isolates, vOrdering, vColoring, cliqueSize);

    vector<int> vAdjunctOrdering = vOrdering;

    if (cliqueSize > 0) {
        cliques.back() = list<int>(vOrdering.begin(), vOrdering.begin() + cliqueSize);
    }

////    coloringStrategy.Recolor(adjacencyMatrix, vAdjunctOrdering, vOrdering, vColoring, cliqueSize, cliqueSize);

////    bool firstIteration(true);
////    size_t count1(0);

    map<int,int> remap;

    while (isolates.GetInGraph().Size() > 913) {
////    while (count1++ < 5 && !isolates.GetInGraph().Empty()) {
#ifdef TWO_LEVEL
        int nextVertex1 = ChooseNextVertex(isolates, true /* first level selection*/);
////        int nextVertex1 = ChooseNextVertex(isolates, false /* only select single vertex best */);
////        int const nextVertex1 = NewChooseNextVertex(isolates, vAdjunctOrdering, vCliqueVerticesPersistent1.size(), cliques.back().size(), -1);

#else
        int nextVertex1 = ChooseNextVertex(isolates, false /* only select single vertex best */);
#endif // TWO_LEVEL
////        int const unusedVertex = NewChooseNextVertex(isolates, vAdjunctOrdering, vCliqueVerticesPersistent1.size(), cliques.back().size(), nextVertex1);
////        if (firstIteration) {
////            firstIteration = false;
////            nextVertex1 = ChooseNextVertex(isolates, false /* only select single vertex best */);
////        } else {
////            nextVertex1 = NewChooseNextVertex(isolates, vAdjunctOrdering, vCliqueVerticesPersistent1.size(), cliques.back().size(), nextVertex1);
////        }

        if (isolates.GetInGraph().Size() == 924) {
            list<int> realClique;
            realClique.insert(realClique.end(), vCliqueVerticesPersistent1.begin(), vCliqueVerticesPersistent1.end());

            vector<int> const vRemainingVertices(isolates.GetInGraph().begin(), isolates.GetInGraph().end());
            GetVertexRemap(vRemainingVertices, isolates, remap);
        }

        if (nextVertex1 == -1) break;
        cout << "Next vertex: " << nextVertex1 << endl;
        if (vColoring.back() + vCliqueVerticesPersistent1.size() <= cliques.back().size()) break;

        cout << "Vertices remaining: " << isolates.GetInGraph().Size() << endl << flush;

////        {
////        size_t numLeft = vColoring.size();
////        for (; numLeft > 0; --numLeft) {
////            if (vCliqueVerticesPersistent1.size() + vColoring[numLeft-1] <= cliques.back().size()) { break; }
////        }
////        numLeft = vColoring.size() - numLeft;
////        cout << ", P.before.left = " << numLeft << endl;
////        }

        vector<int> vRemoved1;
        vector<int> vCliqueVertices1; vCliqueVertices1.push_back(nextVertex1);
        isolates.RemoveVertexAndNeighbors(nextVertex1, vRemoved1);
        ////        vRemoved1.push_back(nextVertex1);
        isolates.RemoveAllIsolates(0, vCliqueVertices1, vRemoved1, vAddedEdgesUnused, false);
        ////            cout << "After removing isolates there are " << isolates.GetInGraph().Size() << " elements remaining " << endl;

        list<int> currentClique(vCliqueVertices1.begin(), vCliqueVertices1.end());
        currentClique.insert(currentClique.end(), vCliqueVerticesPersistent1.begin(), vCliqueVerticesPersistent1.end());

        size_t numLeft = vColoring.size();
        for (; numLeft > 0; --numLeft) {
            if (currentClique.size() + vColoring[numLeft-1] <= cliques.back().size()) { break; }
        }
        numLeft = vColoring.size() - numLeft;
        cout << ", P.after.left = " << numLeft << endl;

#ifdef TWO_LEVEL

        vector<int> vCliqueVerticesPersistent2;
        vector<int> vRemovedPersistent2;

        vector<int> const vToReplace(isolates.GetInGraph().begin(), isolates.GetInGraph().end());

        if (currentClique.size() > cliques.back().size()) {
            cliques.back().clear();
            cliques.back() = currentClique;
            cout << "Found independent set of size: " << currentClique.size() << endl;
        }

        vector<int> vNewAdjunctOrdering = vAdjunctOrdering;
        vector<int> vNewOrdering = vOrdering;
        vector<int> vNewColoring = vColoring;
        ////        while (!isolates.GetInGraph().Empty()) {
////        size_t count2(0);
////        while (!isolates.GetInGraph().Empty() && count2++ < 5) {
        while (isolates.GetInGraph().Size() > 700) {
////        if (!isolates.GetInGraph().Empty()) {

            vector<vector<int>> vComponents;
            ComputeConnectedComponents(isolates, vComponents);
            if (vComponents.size() > 1) {
                cout << "Could break graph into " << vComponents.size() << "components:" << endl;
                cout << "[ ";

                for (size_t index = 0; index < vComponents.size(); ++index) {
                    vector<int> const& vComponent(vComponents[index]);
                    cout << vComponent.size() << " ";
                }
                cout << "]" << endl;
            }

            list<int> newRealClique(vCliqueVertices1.begin(), vCliqueVertices1.end());
            newRealClique.insert(newRealClique.end(), vCliqueVerticesPersistent1.begin(), vCliqueVerticesPersistent1.end());
            ////            newRealClique.insert(newRealClique.end(), vCliqueVertices2.begin(), vCliqueVertices2.end());
            newRealClique.insert(newRealClique.end(), vCliqueVerticesPersistent2.begin(), vCliqueVerticesPersistent2.end());

            if (vNewColoring.empty()) {
                if (newRealClique.size() > cliques.back().size()) {
                    cliques.back().clear();
                    cliques.back() = newRealClique;
                    cout << "Found independent set of size: " << newRealClique.size() << endl;
                }
                break;
            }

            if (vNewColoring.back() + newRealClique.size() <= cliques.back().size()) break;

            int const nextVertex2 = ChooseNextVertex(isolates, false /* second level selection */);
////            int const nextVertex2 = NewChooseNextVertex(isolates, vNewAdjunctOrdering, newRealClique.size(), cliques.back().size(), -1);
            if (nextVertex2 == -1) {
                break;
            }

            ////            cout << "    Next vertex: " << nextVertex2 << endl;
            cout << "Vertices remaining (inner loop): " << isolates.GetInGraph().Size() << endl << flush;

            vector<int> vRemoved2;
            vector<int> vCliqueVertices2; vCliqueVertices2.push_back(nextVertex2);
            isolates.RemoveVertexAndNeighbors(nextVertex2, vRemoved2);
            isolates.RemoveAllIsolates(0, vCliqueVertices2, vRemoved2, vAddedEdgesUnused, false);
            ////            cout << "After removing isolates there are " << isolates.GetInGraph().Size() << " elements remaining " << endl;
#endif // TWO_LEVEL

            list<int> realClique;
            realClique.insert(realClique.end(), vCliqueVerticesPersistent1.begin(), vCliqueVerticesPersistent1.end());
            realClique.insert(realClique.end(), vCliqueVertices1.begin(), vCliqueVertices1.end());
#ifdef TWO_LEVEL
            realClique.insert(realClique.end(), vCliqueVerticesPersistent2.begin(), vCliqueVerticesPersistent2.end());
            realClique.insert(realClique.end(), vCliqueVertices2.begin(), vCliqueVertices2.end());
#endif // TWO_LEVEL

            vector<vector<int>> vNewComponents;
            ComputeConnectedComponents(isolates, vNewComponents);

            sort (vNewComponents.begin(), vNewComponents.end(), [this](vector<int> const &left, vector<int> const &right) { return left.size() < right.size(); });

////            cout << "Components:" << endl;
////            for (size_t index = 0; index < vNewComponents.size(); ++index) {
////                if (vNewComponents[index].size() == 0) continue;
////                ComputeOnConnectedComponent(vNewComponents[index], isolates, m_AdjacencyList, realClique, cliques, index == vNewComponents.size()-1);
////            }

            if (realClique.size() > cliques.back().size()) {
                cliques.back().clear();
                cliques.back() = realClique;
                cout << "Found independent set of size: " << realClique.size() << endl;
            }

#ifdef TWO_LEVEL
            isolates.ReplaceAllRemoved(vCliqueVertices2);
            isolates.ReplaceAllRemoved(vRemoved2);
            isolates.RemoveVertex(nextVertex2);
            isolates.RemoveAllIsolates(0, vCliqueVerticesPersistent2, vRemovedPersistent2, vAddedEdgesUnused, false);

            size_t uNewSize(0);
            vNewAdjunctOrdering.resize(vAdjunctOrdering.size());
            for (size_t index = 0; index < vAdjunctOrdering.size(); ++index) {
                int const vertex(vAdjunctOrdering[index]);
                if (isolates.GetInGraph().Contains(vertex)) {
                    vNewAdjunctOrdering[uNewSize++] = vertex;
                }
            }
            vNewAdjunctOrdering.resize(uNewSize);
            vNewOrdering.resize(vNewAdjunctOrdering.size());
            vNewColoring.resize(vNewAdjunctOrdering.size());

            list<int> newClique(vCliqueVertices1.begin(), vCliqueVertices1.end());
            newClique.insert(newClique.end(), vCliqueVerticesPersistent1.begin(), vCliqueVerticesPersistent1.end());
////            newClique.insert(newClique.end(), vCliqueVertices2.begin(), vCliqueVertices2.end());
            newClique.insert(newClique.end(), vCliqueVerticesPersistent2.begin(), vCliqueVerticesPersistent2.end());

            ////        OrderingTools::InitialOrderingMISR(m_AdjacencyList, isolates, vOrdering, vColoring, cliqueSize);
            coloringStrategy.Recolor(adjacencyMatrix, vNewAdjunctOrdering, vNewOrdering, vNewColoring, newClique.size(), cliques.back().size());
        }

        list<int> realClique;
        realClique.insert(realClique.end(), vCliqueVertices1.begin(), vCliqueVertices1.end());
        realClique.insert(realClique.end(), vCliqueVerticesPersistent1.begin(), vCliqueVerticesPersistent1.end());
        realClique.insert(realClique.end(), vCliqueVerticesPersistent2.begin(), vCliqueVerticesPersistent2.end());
        if (!isolates.GetInGraph().Empty() && vNewColoring.back() + realClique.size() > cliques.back().size()) {

            vector<int> const vRemainingVertices(isolates.GetInGraph().begin(), isolates.GetInGraph().end());
            ComputeOnConnectedComponent(vRemainingVertices, isolates, m_AdjacencyList, realClique, cliques, true);

            if (realClique.size() > cliques.back().size()) {
                cliques.back().clear();
                cliques.back() = realClique;
                cout << "Found independent set of size: " << realClique.size() << endl;
            }
        }

////        isolates.ReplaceAllRemoved(vCliqueVerticesPersistent2);
////        isolates.ReplaceAllRemoved(vRemovedPersistent2);
        isolates.ReplaceAllRemoved(vToReplace);
#endif // TWO_LEVEL

        isolates.ReplaceAllRemoved(vCliqueVertices1);
        isolates.ReplaceAllRemoved(vRemoved1);
        isolates.RemoveVertex(nextVertex1);
        isolates.RemoveAllIsolates(0, vCliqueVerticesPersistent1, vRemovedPersistent1, vAddedEdgesUnused, false);

        size_t uNewSize(0);
        for (size_t index = 0; index < vAdjunctOrdering.size(); ++index) {
            int const vertex(vAdjunctOrdering[index]);
            if (isolates.GetInGraph().Contains(vertex)) {
                vAdjunctOrdering[uNewSize++] = vertex;
            }
        }
        vAdjunctOrdering.resize(uNewSize);
        vOrdering.resize(vAdjunctOrdering.size());
        vColoring.resize(vAdjunctOrdering.size());

////        cout << "Actual Ordering: ";
////        for (int const vertexInOrder : vAdjunctOrdering) {
////            cout << vertexInOrder << " ";
////        }
////        cout << endl;

////        OrderingTools::InitialOrderingMISR(m_AdjacencyList, isolates, vOrdering, vColoring, cliqueSize);
        currentClique = list<int>(vCliqueVerticesPersistent1.begin(), vCliqueVerticesPersistent1.end());
        coloringStrategy.Recolor(adjacencyMatrix, vAdjunctOrdering, vOrdering, vColoring, currentClique.size(), cliques.back().size());

        cout << "Actual:  Max clique size=" << cliques.back().size() << endl;
        cout << "Actual:  Vertex " << nextVertex1 << " has " << currentClique.size() << " clique vertices " << endl;
    }

    if (!isolates.GetInGraph().Empty() && vColoring.back() + vCliqueVerticesPersistent1.size() > cliques.back().size()) {
        list<int> realClique;
        realClique.insert(realClique.end(), vCliqueVerticesPersistent1.begin(), vCliqueVerticesPersistent1.end());

        vector<int> const vRemainingVertices(isolates.GetInGraph().begin(), isolates.GetInGraph().end());
        ComputeOnConnectedComponent(vRemainingVertices, isolates, m_AdjacencyList, realClique, cliques, remap, true);

        if (realClique.size() > cliques.back().size()) {
            cliques.back().clear();
            cliques.back() = realClique;
            cout << "Found independent set of size: " << realClique.size() << endl;
        }
    }

    cout << "Found independent set of size: " << cliques.back().size() << endl << flush;
#endif // 0

#if 0
int savedi(0), savedj(0), savedk(0);
for (size_t i = distribution(generator); i < size-1; i = distribution(generator)) {
    vector<int> v_iNeighbors;
    isolates.RemoveVertexAndNeighbors(i, v_iNeighbors);
    for (size_t j = distribution(generator); ; j = distribution(generator)) {
        vector<int> v_jNeighbors;
        isolates.RemoveVertexAndNeighbors(j, v_jNeighbors);
        for (size_t k = distribution(generator); ; k = distribution(generator)) {
            loops++;
            if (loops %10 == 0) {
                break;
            }
            vector<int> v_kNeighbors;
            isolates.RemoveVertexAndNeighbors(k, v_kNeighbors);
            vector<vector<int>> vComponents;
            ComputeConnectedComponents(isolates, vComponents);
            size_t biggest(0);
            for (vector<int> const &vComponent : vComponents) {
                biggest = max(vComponent.size(), biggest);
            }

            if (biggest < best) {
                best = biggest;
                cout << "best so far: remove " << i << " " << j << " " << k << " for max component of size " << best << endl;
                savedi = i;
                savedj = j;
                savedk = k;
            }
                isolates.ReplaceAllRemoved(v_kNeighbors);
            }
            v_jNeighbors.push_back(j);
            isolates.ReplaceAllRemoved(v_jNeighbors);
            if (loops % 100 == 0) break;
        }
        v_iNeighbors.push_back(i);
        isolates.ReplaceAllRemoved(v_iNeighbors);
        if (loops % 5000 == 0) {
            cout << "Removing " << savedi << " " << savedj << " " << savedk << endl;
            vector<int> vRemovedNeighbors;
            isolates.RemoveVertexAndNeighbors(savedi, vRemovedNeighbors);
            isolates.RemoveVertexAndNeighbors(savedj, vRemovedNeighbors);
            isolates.RemoveVertexAndNeighbors(savedk, vRemovedNeighbors);
        }
    }
#endif // 0

#if 0
    // read iperm file
    vector<int> vOrdering;
    vector<int> vColoring;

////    OrderingTools::InitialOrderingConnectedComponent("/Users/strash/graphs/kernels/biogrid-yeast.955.kernel.graph.iperm", adjacencyMatrix, vOrdering, vColoring);
////    OrderingTools::InitialOrderingConnectedComponent("/Users/strash/graphs/kernels/facebook.kernel.graph.iperm", adjacencyMatrix, vOrdering, vColoring);
    OrderingTools::InitialOrderingConnectedComponent("/Users/strash/graphs/kernels/facebook.2754.kernel.graph.iperm", adjacencyMatrix, vOrdering, vColoring);

    // analyze the size of the connected components and the degrees of the removed vertices.
    // need the degrees to be high, and the connected components to be roughly the same size.

    vector<int> vNeighborsRemoved;
////    vector<pair<int,int>> vAddedEdgesUnused;
    for (size_t index = vOrdering.size(); index > 0; --index) {
        if (isolates.GetInGraph().Contains(vOrdering[index-1]))
            isolates.RemoveVertexAndNeighbors(vOrdering[index-1], vNeighborsRemoved);
////            isolates.RemoveVertex(vOrdering[index-1]);
////            isolates.RemoveVertex(vOrdering[index-1]);
////    for (size_t index = 0; index < vOrdering.size();  ++index) {
////        if (isolates.GetInGraph().Contains(vOrdering[index]))
////            isolates.RemoveVertexAndNeighbors(vOrdering[index], vNeighborsRemoved);
////        cout << "After removing " << vOrdering.size() - (index - 1) << " vertices: " << endl;
        else continue;
        vector<vector<int>> vComponents;
        ComputeConnectedComponents(isolates, vComponents);
        cerr << "# vertices remaining in graph: " << isolates.GetInGraph().Size() << "/" << m_AdjacencyList.size() << endl << flush;
        cerr << "# connected components       : " << vComponents.size() << endl << flush;
        cerr << "size of connected components : ";
        cout << "[ ";

        for (size_t index = 0; index < vComponents.size(); ++index) {
            vector<int> const& vComponent(vComponents[index]);
            cout << vComponent.size() << " ";
        }
        cout << "]" << endl;

        isolates.RemoveAllIsolates(0, vNeighborsRemoved, vNeighborsRemoved, vAddedEdgesUnused, false);
    }

    // can we use the separation to find an approximate solution fast? basically we would
    // be within - (total size of separators). And really, how many levels could there be until we get something
    // that is reasonably solvable? Certainly no more than 5 or 6 for 4000-vertex graphs, right?

    // what could throw a wrench in this is if the connected components have many vertices of large degree. Hopefully
    // that's not the case, because they should be in the separator.

    // as of now, this loop helps find a really small Independent set, this is to help limit the recursion depth.
    // need to expand this into a branch-and-bound, pick the vertices that constrain the search the most, to keep
    // the search at a reasonable depth. And try to "peel off" as many vertices as possible through reductions.

////    list<list<int>> result;
////    TesterMISS algorithm(adjacencyMatrix, m_AdjacencyList);
////    algorithm.Run(result);
#endif // 0

#if 0
    while (vRemoved.size() != m_AdjacencyList.size()) {
        // find vertex whose removal maximizes the number of vertices moved to IS by reduction.
        int vertexWithMaxReductions(-1);

        int maxIsolates(-1);
        cout << "Testing remaining " << m_AdjacencyList.size() - vRemoved.size() << " vertices, to maximize isolate removal" << endl << flush;

        int const vertexToRemove(isolates.NextVertexToRemove());

        if (vertexToRemove == -1) {
            cout << "Something is WRONG: NextVertexToRemove returned -1" << endl << flush;
            exit(1);
        }

        cout << "Removing next vertex..." << endl << flush;
        vector<int> vNextRemoved;
        vector<int> vNextIsolates;
        vector<pair<int,int>> vNextAddedEdges;
        isolates.RemoveVertexAndNeighbors(vertexToRemove, vNextRemoved);
        isolates.RemoveAllIsolates(vIsolates.size(), vNextIsolates, vNextRemoved, vNextAddedEdges);
        cout << "Newly added edges:" << endl;
        for (pair<int,int> const &edge : vNextAddedEdges) {
            cout << "    " << edge.first << "," << edge.second << endl;
        }

        independentSet.push_back(vertexToRemove);
        for (int const vertex : vNextIsolates) {
            independentSet.push_back(vertex);
        }

        vNextRemoved.insert(vNextRemoved.end(), vNextIsolates.begin(), vNextIsolates.end());

        for (int const vertex : vNextRemoved) {
            if (setRemoved.find(vertex) != setRemoved.end()) {
                cout << "Vertex " << vertex << " is being removed twice!" << endl << flush;
            }
            setRemoved.insert(vertex);
        }

        vRemoved.insert(vRemoved.end(), vNextRemoved.begin(), vNextRemoved.end());
        vAddedEdges.insert(vAddedEdges.end(), vNextAddedEdges.begin(), vNextAddedEdges.end());

        cout << "# vertices remaining in graph: " << m_AdjacencyList.size() - vRemoved.size() << "/" << m_AdjacencyList.size() << endl << flush;
    }

    ////ExecuteCallBacks(independentSet);

    set<int> testIndependentSet;
    list<int> trueIndependentSet;
////    testIndependentSet.insert(independentSet.begin(), independentSet.end());
    // build up the real independent set using alternatives from path isolate removal.
    for (list<int>::const_reverse_iterator cit = independentSet.rbegin(); cit != independentSet.rend(); ++cit) {
        int const vertex(*cit);
        bool useAlternative(false);
        for (int const neighbor : m_AdjacencyList[vertex]) {
            if (testIndependentSet.find(neighbor) != testIndependentSet.end()) {
                cout << "Need to repair independent set: vertex " << vertex << " has neighbor " << neighbor << " in the independent set." << endl;
                useAlternative = true;
                break;
            }
        }

        if (!useAlternative) {
            testIndependentSet.insert(vertex);
            trueIndependentSet.push_front(vertex);
            continue;
        }

        int const alternative(isolates.GetAlternativeVertex(vertex));
        if (alternative == -1) {
            cout << "Unable to repair independent set: vertex " << vertex << " has no alternative" << endl;
            continue;
        }

        for (int const neighbor2 : m_AdjacencyList[alternative]) {
            if (testIndependentSet.find(neighbor2) != testIndependentSet.end()) {
                cout << "Unable to repair independent set: " << vertex << "'s alternative " << alternative << " does not fix independent set. " << alternative << " has neighbor " << neighbor2 << " in independent set." << endl;
                break;
            }
        }

        trueIndependentSet.push_front(alternative);
        testIndependentSet.insert(alternative);
    }

    ExecuteCallBacks(trueIndependentSet);

    cout << "Removed " << vRemoved.size() << " vertices." << endl << flush;
    cout << "Found independent set of size: " << isolates.size() << endl << flush;
#endif // 0

#endif // 0

    //// print graph size and kernel size.
#if 0
    size_t edges(0);
    for (vector<int> const &neighbors : m_AdjacencyList) {
        edges += neighbors.size();
    }
    edges >>= 1;

    cout << "Begin removal of rest of graph..." << endl << flush;

    // remove remaining vertices
    while (!isolates.GetInGraph().Empty()) {
        int const nextVertex = (*isolates.GetInGraph().begin());
        isolates.RemoveVertexAndNeighbors(nextVertex, vRemoved, vReductions);
////        isolates.RemoveVertex(nextVertex, vReductions);
        isolates.RemoveAllIsolates(0, vIsolates, vRemoved, vReductions, true /* consider all vertices for reduction */);
    }

    size_t setSize(0);
    for (Reduction const &reduction : vReductions) {
        if (reduction.GetType() != DOMINATED_VERTEX) {
            setSize++;
        }
    }
    cout << "Found set of size: " << setSize << endl << flush;

    cout << "Checking undo reductions..." << endl << flush;
    if (true)
    {
        isolates.ReplaceAllRemoved(vReductions);
        for (size_t vertex = 0; vertex < m_AdjacencyList.size(); ++vertex) {
            if (!isolates.GetInGraph().Contains(vertex)) {
                cout << "ERROR during sanity check." << endl << flush;
                cout << "Vertex " << vertex << " is excluded from graph after undoing all reductions." << endl;
                exit(1);
            }

            if (isolates.Neighbors()[vertex].Size() != m_AdjacencyList[vertex].size()) {
                cout << "ERROR during sanity check." << endl << flush;
                cout << "Vertex " << vertex << " does not have same number of neighbors after undoing all reductions." << endl;
                exit(1);
            }

            for (int const neighbor : m_AdjacencyList[vertex]) {
                if (!isolates.Neighbors()[vertex].Contains(neighbor)) {
                    cout << "ERROR during sanity check." << endl << flush;
                    cout << "Vertex " << vertex << " is missing neighbor " << neighbor << " after undoing all reductions." << endl;
                    exit(1);
                }
            }
        }
    }

    vIsolates.clear();
    vRemoved.clear();
    vReductions.clear();
    isolates.RemoveAllIsolates(0, vIsolates, vRemoved, vReductions, true /* consider all vertices for reduction */);

    cout << "Removed " << vIsolates.size() << " isolates, graph has " << isolates.GetInGraph().Size() << "/" << m_AdjacencyList.size() << " vertices remaining" << endl;

    vector<vector<int>> vComponents;
    GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyList.size());

    size_t maxComponentSize(0);
    size_t maxComponentIndex(0);
    for (size_t index = 0; index < vComponents.size(); ++index) {
        vector<int> const &vComponent(vComponents[index]);
        if (vComponent.size() > maxComponentSize) {
            maxComponentSize = max(maxComponentSize, vComponent.size());
            maxComponentIndex = index;
        }
    }

    cout << "latex:" << m_AdjacencyList.size() << " & " << edges << " & " << isolates.GetInGraph().Size() << " & " << vComponents.size() << " & " << maxComponentSize << endl;

    size_t independentSetSize(isolates.size());
    vector<bool> vIndependentVertices(m_AdjacencyList.size(), false);

    int const testVertex(2745);
    int const testVertex2(168585);

    for (int const vertex : vIsolates) {
        vIndependentVertices[vertex] = true;
        if (vertex == testVertex) {
            cout << "Test vertex " << testVertex << " is removed in initial isolates" << endl << flush;
        }
        if (vertex == testVertex2) {
            cout << "Test vertex " << testVertex2 << " is removed in initial isolates" << endl << flush;
        }
    }

////    cout << "Neighbors of test vertex " << testVertex << ":" << endl << flush;
////    for (int const neighbor : isolates.Neighbors()[testVertex]) {
////        cout << neighbor << ": ";
////        for (int const nNeighbor : isolates.Neighbors()[neighbor]) {
////            cout << nNeighbor << " ";
////        }
////        cout << endl << flush;
////    }
////    cout << endl << flush;

////    if (maxComponentSize < 1000) {
    for (vector<int> const vComponent : vComponents) {
        list<int> clique;
        list<list<int>> cliques;
        cout << "Running connect component algorithm..." << endl;
////        ComputeOnConnectedComponent(vComponents[maxComponentIndex], isolates, m_AdjacencyList, clique, cliques, false /* don't set clique size */);
        ComputeOnConnectedComponent(vComponent, isolates, m_AdjacencyList, clique, cliques, false /* don't set clique size */);
        independentSetSize += clique.size();
        for (int const vertex : clique) {
            vIndependentVertices[vertex] = true;
            if (vertex == testVertex) {
                cout << "Test vertex " << testVertex << " is added in clique finding." << endl << flush;
            }
            if (vertex == testVertex2) {
                cout << "Test vertex " << testVertex2 << " is added in clique finding." << endl << flush;
            }
        }

#if 1
        list<int> vertexSet;
        for (int vertex = 0; vertex < vIndependentVertices.size(); ++vertex) {
            if (vIndependentVertices[vertex])
                vertexSet.push_back(vertex);
        }

////        if (vIndependentVertices[testVertex]) {
////            cout << "Test vertex " << testVertex << " is in independent set" << endl << flush;
////        }

        if (!CliqueTools::IsIndependentSet(m_AdjacencyList, vertexSet, false /* be quiet */)) {
            cout << "Is not an independent set!" << endl << flush;
            if (vComponent.size() == 111) {
                map<int,int> remap;
                GetVertexRemap(vComponent, isolates, remap);
                PrintSubgraphInEdgesFormat(m_AdjacencyList, remap);
            }
        }
#endif // 0
    }
    cout << "Original graph size           : " << m_AdjacencyList.size() << endl << flush;
    cout << "New      graph size           : " << isolates.GetInGraph().Size() << endl << flush;
    cout << "Removed  vertices             : " << vIsolates.size() + vRemoved.size() << endl << flush;
    cout << "Recomputed new graph size     : " << m_AdjacencyList.size() - (vIsolates.size() + vRemoved.size()) << endl << flush;
    cout << "Computed independent set size : " << independentSetSize  << endl << flush;
    // include folded vertices into independent set:

////    cout << "Neighbors of test vertex " << testVertex << ":" << endl << flush;
////    for (int const neighbor : isolates.Neighbors()[testVertex]) {
////        cout << neighbor << ": ";
////        for (int const nNeighbor : isolates.Neighbors()[neighbor]) {
////            cout << nNeighbor << " ";
////        }
////        cout << endl << flush;
////    }
////    cout << endl << flush;
////    if (vIndependentVertices[testVertex] && vIndependentVertices[testVertex2]) {
////        cout << "Test vertex " << testVertex << " is in independent set with test vertex " << testVertex2 << endl << flush;
////    }
    int actualAdditions(0);
    for (size_t index = vReductions.size(); index > 0; index--) {
        Reduction const &reduction(vReductions[index-1]);
        if (reduction.GetType() != FOLDED_VERTEX) continue;
        int const newVertex(reduction.GetVertex());
        int const removedVertex1(reduction.GetNeighbors()[0]);
        int const removedVertex2(reduction.GetNeighbors()[1]);
 #if 0
        list<int> vertexSet;
        cout << "Adding vertex fold #" << actualAdditions << ":" << vertexFold.newVertex << " <- " << vertexFold.removedVertex1 << " , " << vertexFold.removedVertex2  << endl << flush;
        for (int vertex = 0; vertex < vIndependentVertices.size(); ++vertex) {
            if (vIndependentVertices[vertex])
                vertexSet.push_back(vertex);
        }

////        if (vIndependentVertices[testVertex]) {
////            cout << "Test vertex " << testVertex << " is in independent set" << endl << flush;
////        }

        if (!CliqueTools::IsIndependentSet(m_AdjacencyList, vertexSet, false /* be quiet */)) {
            cout << "Is not an independent set!" << endl << flush;
        }
#endif // 0

        if (vIndependentVertices[newVertex]) {
            actualAdditions++;
            if (vIndependentVertices[removedVertex1] || vIndependentVertices[removedVertex2]) {
                cout << "ERROR: Something is wrong..." << endl << flush;
            }
            vIndependentVertices[removedVertex1] = true;
            vIndependentVertices[removedVertex2] = true;
            vIndependentVertices[newVertex] = false;
////            cout << "Removing " << vertexFold.newVertex << ", replacing with " << vertexFold.removedVertex1 << " , " << vertexFold.removedVertex2 << endl << flush;
        } else {
            actualAdditions++;
            if (vIndependentVertices[removedVertex1] || vIndependentVertices[removedVertex2]) {
                cout << "ERROR: Something is wrong..." << endl << flush;
            }
            vIndependentVertices[newVertex] = true;
////            cout << "Adding " << vertexFold.newVertex << " to independent set" << endl << flush;
        }
    }

    list<int> vertexSet;
    for (int vertex = 0; vertex < vIndependentVertices.size(); ++vertex) {
        if (vIndependentVertices[vertex])
            vertexSet.push_back(vertex);
    }

    cout << "Final Check: " << endl << flush;
    if (!CliqueTools::IsIndependentSet(m_AdjacencyList, vertexSet, false /* be quiet */)) {
        cout << "Is not an independent set!" << endl << flush;
    }

    cout << "Actual    additions           : " << actualAdditions << endl << flush;
    cout << "New      independent set size : " << independentSetSize + actualAdditions << endl << flush;

#if 1
    cout << endl << "Processing dominated vertices..." << endl << flush;
    for (size_t index = vRemoved.size(); index > 0; index--) {
        int const vertex(vRemoved[index-1]);
        if (vIndependentVertices[vertex]) continue;
        bool noNeighbor(true);
        for (int const neighbor : m_AdjacencyList[vertex]) {
            if (vIndependentVertices[neighbor]) {
                noNeighbor = false;
                break;
            }
        }

        if (noNeighbor) {
            actualAdditions++;
            vIndependentVertices[vertex] = true;
        }
    }

    cout << "Final Check(2): " << endl << flush;
    if (!CliqueTools::IsIndependentSet(m_AdjacencyList, vertexSet, false /* be quiet */)) {
        cout << "Is not an independent set!" << endl << flush;
    }

    cout << "New       additions           : " << actualAdditions << endl << flush;
    cout << "New      independent set size : " << independentSetSize + actualAdditions << endl << flush;
#endif // 0
#endif //1

// try to build up a graph from isolated vertices. remaining graph needs to be removed to empty the graph
// assume there aren't any isolates vertices
#if 0

    // Attempt #1, build up graph with only degree 1 vertices.
    cout << "Computing upper bound: " << endl << flush;

    vector<int> vOriginalAdjunctOrdering;
    vector<int> vOriginalOrdering;
    vector<int> vOriginalColoring;

    size_t cliqueSize(0);

////    OrderingTools::InitialOrderingMISR(m_AdjacencyList, vOriginalAdjunctOrdering, vOriginalColoring, cliqueSize);
    OrderingTools::InitialOrderingMISR(m_AdjacencyList, vOriginalAdjunctOrdering, vOriginalColoring, cliqueSize);

    cout << "Original upper bound from ordering: " << vOriginalColoring.back() << endl << flush;
    cout << "Original lower bound from ordering: " << cliqueSize << endl << flush;

    ArraySet removed(m_AdjacencyList.size());
    ArraySet isolatedSet(m_AdjacencyList.size());
#if 0
    BuildSingleEdgeIsolateGraph(removed, isolatedSet);
#else
    vector<int> independentSet(vOriginalAdjunctOrdering.begin(), vOriginalAdjunctOrdering.begin() + cliqueSize);
    ArraySet forbidden(m_AdjacencyList.size());
    for (int const independentVertex : independentSet) {
        isolatedSet.Insert(independentVertex);
        removed.Insert(independentVertex);
        for (int const neighbor : m_AdjacencyList[independentVertex]) {
            if (!removed.Contains(neighbor) && !forbidden.Contains(neighbor)) {
                removed.Insert(neighbor);
                break;
            }
        }

        for (int const neighbor : m_AdjacencyList[independentVertex]) {
            forbidden.Insert(neighbor);
        }
    }
#endif // 0

    for (int const vertex : removed) {
        isolates.RemoveVertex(vertex);
    }

    IsolatesIndependentSetColoringStrategy<Isolates4<SparseArraySet>> coloringStrategy(isolates, m_AdjacencyList.size());

    coloringStrategy.Recolor(isolates, vOriginalAdjunctOrdering, vOriginalOrdering, vOriginalColoring, cliqueSize, 0);
    cout << "Aftering recoloring: " << vOriginalColoring.back() << endl << flush;

    vector<int> vAdjunctOrdering;
    vector<int> vColoring;
    vector<int> vOrderedVertices;

    vector<int> vNewIsolates;
    vector<int> vNewRemoved;
////    isolates.RemoveAllIsolates(0, vNewIsolates, vNewRemoved, vAddedEdgesUnused, true /* consider all vertices for removal */);

    for (int const vertex : vOriginalAdjunctOrdering) {
        if (isolates.GetInGraph().Contains(vertex)) {
            vAdjunctOrdering.push_back(vertex);
        }
    }

    vColoring.resize(vAdjunctOrdering.size());
    vOrderedVertices.resize(vAdjunctOrdering.size());

    coloringStrategy.Recolor(isolates, vAdjunctOrdering, vOrderedVertices, vColoring, cliqueSize, 0);

    cout << "Estimate: " << isolatedSet.Size() << "+" << vNewIsolates.size() << "+" <<  vColoring.back() << "=" << (isolatedSet.Size() + vNewIsolates.size() + vColoring.back()) << endl;

    cout << "Remaining vertices: " << isolates.GetInGraph().Size() << "/" << m_AdjacencyList.size() << endl << flush;

    vector<vector<int>> const subgraphList = ComputeSubgraphOfSize(isolates, m_AdjacencyList.size(), 3000);

    vector<vector<char>> subgraphMatrix(subgraphList.size());
    for (size_t index = 0; index < subgraphList.size(); ++index) {
        subgraphMatrix[index].resize(subgraphList.size(), 0);
        for (int const neighbor : subgraphList[index]) {
            subgraphMatrix[index][neighbor] = 1;
        }
    }

    list<list<int>> cliques;
    TesterMISS algorithm(subgraphMatrix, subgraphList);
////    algorithm.SetQuiet(true);
////    algorithm.SetOnlyVertex(vOrdering[splitPoint]);
    algorithm.Run(cliques);
#else
#if 0
    vector<int> vOrdering;
    vector<int> vColoring;
    size_t cliqueSize(0);
    OrderingTools::InitialOrderingMISR(m_AdjacencyList, vOrdering, vColoring, cliqueSize);

//    size_t pivot(0);
//    for (size_t afterIndex = vColoring.size(); afterIndex > 0; afterIndex--) {
//        if (vColoring[afterIndex-1] <= cliqueSize) { pivot = afterIndex-1; break; }
//    }

    map<int,int> mapUnused;
////    set<int>     setVertices;
    size_t const startVertex(2000);
////    size_t const startVertex(0);
////    set<int>     setVertices(vOrdering.begin(), vOrdering.begin() + startVertex);
    set<int>     setVertices(vOrdering.begin() + startVertex, vOrdering.end());
    clock_t startTime(clock());
////    for (size_t splitPoint = startVertex/*vOrdering.size()/2*/; splitPoint < vOrdering.size(); splitPoint++) {
    for (size_t splitPoint = startVertex/*vOrdering.size()/2*/; splitPoint > 0; splitPoint--) {
        setVertices.insert(vOrdering[splitPoint]);
        vector<vector<int>> subgraphAdjacencyList;
        GraphTools::ComputeInducedSubgraph(m_AdjacencyList, setVertices, subgraphAdjacencyList, mapUnused);

        cout << "Subgraph size=" << subgraphAdjacencyList.size() << " " << Tools::GetTimeInSeconds(clock() - startTime) << endl << flush;

        vector<vector<char>> subgraphAdjacencyMatrix(subgraphAdjacencyList.size());
        for (size_t index = 0; index < subgraphAdjacencyList.size(); ++index) {
            subgraphAdjacencyMatrix[index].resize(subgraphAdjacencyList.size(), 0);
            for (int const neighbor : subgraphAdjacencyList[index]) {
                subgraphAdjacencyMatrix[index][neighbor] = 1;
            }
        }

        list<list<int>> cliques;

        Isolates3<ArraySet> newIsolates(subgraphAdjacencyList);
        vector<int> vRemoved;
        newIsolates.RemoveAllIsolates(0, vRemoved, vRemoved, vAddedEdgesUnused, true);
        vector<vector<int>> vComponents;
        GraphTools::ComputeConnectedComponents(newIsolates, vComponents, subgraphAdjacencyList.size());
        cerr << "# vertices remaining in graph: " << newIsolates.GetInGraph().Size() << "/" << subgraphAdjacencyList.size() << endl << flush;
        cerr << "# connected components       : " << vComponents.size() << endl << flush;
        cerr << "size of connected components : ";
        cout << "[ ";

        for (size_t index = 0; index < vComponents.size(); ++index) {
            vector<int> const& vComponent(vComponents[index]);
            cout << vComponent.size() << " ";
        }
        cout << "]" << endl;

        TesterMISS algorithm(subgraphAdjacencyMatrix, subgraphAdjacencyList);
////        algorithm.SetQuiet(true);
////        algorithm.SetOnlyVertex(vOrdering[splitPoint]);
        algorithm.Run(cliques);


        if (cliques.back().size() > cliqueSize) {
            cliqueSize = cliques.back().size();
            cout << "At " << splitPoint << "/" << vOrdering.size() << ", found a better independent set: size=" << cliqueSize << endl;
        }
        break;
    }

#endif // 0
#endif // 1
}
