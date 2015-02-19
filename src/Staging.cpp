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
#include "ConnectedComponentMISS.h"

// system includes
#include <vector>
#include <list>
#include <set>
#include <iostream>
#include <algorithm>
#include <climits>

using namespace std;

Staging::Staging(std::vector<std::vector<int>> &adjacencyList)
 : Algorithm("staging")
 , m_AdjacencyList(adjacencyList)
{
}

Staging::~Staging()
{
}

bool RemoveIsolatedClique(int const vertex, vector<set<int>> &neighbors, set<int> &remaining, set<int> &isolates, set<int> &removed, vector<bool> &vMarkedVertices)
{
    size_t neighborCount(0);
    bool superSet(true);
    // TODO/DS: Put back!
    for (int const neighbor : neighbors[vertex]) {
        if (neighbors[neighbor].size() < neighbors[vertex].size()) {
            superSet = false; break;
        }
    }

    ////    if (vertex == 21952) {
    ////        cout << "vertex" << vertex << " has " << neighbors[vertex].size() << " neighbors." << endl << flush;
    ////    }

    if (!superSet) return false;

////    bool const debug(true);//vertex==194);
////    if (debug) {
////        cout << "vertex " << vertex << " has neighbors: ";
////    }
////    if (debug) {
////        for (int const neighbor : neighbors[vertex]) {
////            cout << neighbor << " " << flush;
////        }
////        cout << endl << flush;
////    }

    for (int const neighbor : neighbors[vertex]) {
////        if (debug) {
////            cout << "Considering neighbor " <<  neighbor  << flush;
////        }
        neighborCount++;

        for (int const nNeighbor : neighbors[neighbor]) {
            vMarkedVertices[nNeighbor] = true;
        }
        vMarkedVertices[neighbor] = true;

        for (int const neighbor2 : neighbors[vertex]) {
////            if (debug && superSet && !vMarkedVertices[neighbor2]) {
////                cout << "(missing neighbor " << neighbor2 << ") ";
////            }
            superSet = superSet && vMarkedVertices[neighbor2];
        }

        for (int const nNeighbor : neighbors[neighbor]) {
            vMarkedVertices[nNeighbor] = false;
        }
        vMarkedVertices[neighbor] = false;

        if (!superSet) {
////            if (debug) {
////                cout << "[x] " << endl;
////            }
            return false;
        }

////        if (debug)
////            cout << endl << flush;

        ////            cout << "Iteration(" << neighborCount << ": Intersection=";
        ////            for (int const vertex : intersection) {
        ////                cout << vertex << " ";
        ////            }
        ////            cout << endl << flush;
    }

    ////        cout << "Iteration(" << neighborCount << ": Intersection=";
    ////        for (int const vertex : intersection) {
    ////            cout << vertex << " ";
    ////        }
    ////        cout << endl << flush;

    if (superSet) {
////        cout << "Removing Clique: [" << vertex;
////        for (int const neighbor : neighbors[vertex]) {
////            cout << " " << neighbor;
////        }
////        cout << "]" << endl;
        removed.insert(neighbors[vertex].begin(), neighbors[vertex].end());
        removed.insert(vertex);
        isolates.insert(vertex);
        for (int const neighbor : neighbors[vertex]) {
            remaining.erase(neighbor);
            ////                cout << "   Expunging neighbor " << neighbor << " with " << neighbors[neighbor].size() << " neighbors" << endl << flush;
            for (int const nNeighbor : neighbors[neighbor]) {
                if (removed.find(nNeighbor) == removed.end()) {
                    remaining.insert(nNeighbor);
                }

                if (nNeighbor != vertex) {
                    neighbors[nNeighbor].erase(neighbor);
                }
            }
            neighbors[neighbor].clear();
        }
        neighbors[vertex].clear();

        ////    if (vertex == 21952) {
        ////        cout << "vertex" << vertex << " is being removed." << endl << flush;
        ////    }

        return true;
    }
    return false;
}

bool RemoveIsolatedPath(int const vertex, vector<set<int>> &neighbors, set<int> &remaining, set<int> &isolates, set<int> &removed, vector<bool> &vMarkedVertices)
{
    if (neighbors[vertex].size() != 2) return false;

    for (int const neighbor : neighbors[vertex]) {
        if (neighbors[neighbor].size() == 2) {
            removed.insert(vertex);
            removed.insert(neighbor);
            isolates.insert(vertex);

            int endVertex(-1);

            // remove from other neighbor's list...
            for (int const nNeighbor : neighbors[neighbor]) {
                if (nNeighbor != vertex) {
                    endVertex = nNeighbor;
                    neighbors[nNeighbor].erase(neighbor);
                    remaining.insert(nNeighbor);
                    break;
                }
            }

            for (int const otherNeighbor : neighbors[vertex]) {
                if (otherNeighbor != neighbor) {
                    neighbors[otherNeighbor].erase(vertex);
                    remaining.insert(otherNeighbor);

                    neighbors[otherNeighbor].insert(endVertex);
                    neighbors[endVertex].insert(otherNeighbor);
                    break;
                }
            }

            neighbors[vertex].clear();
            neighbors[neighbor].clear();

            remaining.erase(vertex);
            remaining.erase(neighbor);
////            cout << "Removing Path: [" << vertex << " " << neighbor << "]" << endl << flush;

            return true;
        }
    }
    return false;
}

bool RemoveVertexAndNeighbors(int const vertex, vector<set<int>> &neighbors, set<int> &remaining, set<int> &isolates, set<int> &removedVertices)
{
    //    cout << __LINE__ << ": Removing vertex " << vertex << endl << flush;
    removedVertices.insert(vertex);
    remaining.erase(vertex);
    isolates.insert(vertex);

    for (int const neighbor : neighbors[vertex]) {
        //        cout << __LINE__ << ":     Removing neighbor " << neighbor << endl << flush;
        removedVertices.insert(neighbor);
        remaining.erase(neighbor);
        for (int const nNeighbor : neighbors[neighbor]) {
            //            cout << __LINE__ << ":         Removing from neighbor's neighbor "<< nNeighbor << endl << flush;
            if (nNeighbor == vertex) continue;
            neighbors[nNeighbor].erase(neighbor);
            //            cout << __LINE__ << ":         Done removing" << endl << flush;
            if (removedVertices.find(nNeighbor) == removedVertices.end())
                remaining.insert(nNeighbor);
        }
        //        cout << __LINE__ << ":     Done Removing neighbor" << neighbor << endl << flush;
        neighbors[neighbor].clear();
        //        cout << __LINE__ << ": Cleared neighbors: " << neighbor << endl << flush;
    }

    //    cout << __LINE__ << ": Done Removing vertex " << vertex << endl << flush;
    neighbors[vertex].clear();
    //    cout << __LINE__ << ": Cleared neighbors: " << vertex << endl << flush;
    return true;
}


bool RemoveMinDegreeVertex(vector<set<int>> &neighbors, set<int> &remaining, set<int> &isolates, set<int> &removed, vector<bool> &vMarkedVertices)
{
    if (removed.size() == neighbors.size()) return false;

    int vertex(-1);
    int minDegree(INT_MAX);
    for (int i = 0; i < neighbors.size(); ++i) {
        if (neighbors[i].size() < minDegree && removed.find(i) == removed.end()) {
            minDegree = neighbors[i].size();
            vertex = i;
        }
    }

    //cout << "Removing vertex " << vertex << " with degree " << minDegree << endl << flush;

    RemoveVertexAndNeighbors(vertex, neighbors, remaining, isolates, removed);

    return true;
}

bool RemoveMaxDegreeVertex(vector<set<int>> &neighbors, set<int> &remaining, set<int> &isolates, set<int> &removed, vector<bool> &vMarkedVertices)
{
    if (removed.size() == neighbors.size()) return false;

    int vertex(-1);
    int maxDegree(INT_MIN);
    for (int i = 0; i < neighbors.size(); ++i) {
        if (static_cast<int>(neighbors[i].size()) > maxDegree && removed.find(i) == removed.end()) {
            maxDegree = static_cast<int>(neighbors[i].size());
            vertex = i;
        }
    }

    ////    cout << "Removing vertex " << vertex << " with degree " << maxDegree << endl << flush;

    RemoveVertexAndNeighbors(vertex, neighbors, remaining, isolates, removed);

    return true;
}

bool RemoveAllIsolates(vector<set<int>> &neighbors, set<int> &remaining, set<int> &isolates, set<int> &removed, vector<bool> &vMarkedVertices)
{
    //cout << "Removing all isolates." << endl << flush;
    int isolateSize(isolates.size());
    while (!remaining.empty()) {
        ////            if (remaining.size() %10000 == 0)
        //cout << "Remaining: " << remaining.size() << ", Removed: " << removed.size() << endl << flush;
        set<int>::iterator sit = remaining.begin();
        int const vertex = *sit;
        remaining.erase(sit);

////        cout << "Attempting to remove vertex " << vertex << endl << flush;

        bool reduction = RemoveIsolatedClique(vertex, neighbors, remaining, isolates, removed, vMarkedVertices);
        if (!reduction) {
            reduction = RemoveIsolatedPath(vertex, neighbors, remaining, isolates, removed, vMarkedVertices);
        }
    }

////    cout << "Removed " << isolates.size() - isolateSize << " isolates." << endl << flush;
    return true;
}


template <typename C> bool ReplaceRemovedVertices(vector<vector<int>> const &adjacencyArray, vector<set<int>> &neighbors, set<int> const &inGraph, C const &removedVertices)
{
    //cout << "Replacing all removed vertices." << endl << flush;
    for (typename C::value_type const removedVertex : removedVertices) {
        ////            if (remaining.size() %10000 == 0)
        ////        cout << "Remaining: " << remaining.size() << ", Removed: " << removed.size() << endl << flush;
        for (int const neighbor : adjacencyArray[removedVertex]) {
            if (inGraph.find(neighbor) == inGraph.end()) continue;
            neighbors[removedVertex].insert(neighbor);
            neighbors[neighbor].insert(removedVertex);
        }
    }
    //cout << "Replaced vertices..." << endl << flush;
    return true;
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
    vector<pair<int,int>> vAddedEdges;
    subgraphIsolates.RemoveAllIsolates(0, vIsolates, vRemoved, vAddedEdges, true /* consider all vertices for reduction */);

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
    vector<pair<int,int>> vAddedEdges;
    isolates.RemoveAllIsolates(0, vIsolates, vRemoved, vAddedEdges, true /* consider all vertices for reduction */);

#if 1
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

    vector<vector<int>> const subgraphAdjacencyList = ComputeSubgraphOfSize(isolates, m_AdjacencyList.size(), 1024);

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
    LightWeightReductionFullMISS algorithm(subgraphAdjacencyMatrix, subgraphAdjacencyList);
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
    algorithm.AddCallBack(printCliqueSize);
    algorithm.Run(cliques);

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

    // as of now, this loop helps find a really small Independent set, this is to help limit the recursion depth.
    // need to expand this into a branch-and-bound, pick the vertices that constrain the search the most, to keep
    // the search at a reasonable depth. And try to "peel off" as many vertices as possible through reductions.

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

#endif // 1
}
