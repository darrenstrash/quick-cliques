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
        for (int const neighbor : adjacencyArray[vertex]) {
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

    list<list<int>> componentResult;

    ////                cout << "Start algorithm: " << endl << flush;
    TesterMISS algorithm(componentMatrix, componentArray);
    algorithm.SetQuiet(true);

////    if (bSetCliqueSize) {
////        algorithm.SetMaximumCliqueSize(cliques.back().size() - realClique.size());
////    }

    algorithm.Run(componentResult);

    if (!componentResult.empty() && !componentResult.back().empty()) {
        for (int const vertex : componentResult.back()) {
            realClique.push_back(reverseMap[vertex]);
        }
    }
    ////                cout << "    Component has independent set of size " << componentResult.back().size() << endl;
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

    sort (vVertices.begin(), vVertices.end(), [this](int const left, int const right) { return m_AdjacencyList[left].size() > m_AdjacencyList[right].size(); });

    std::random_device generator;
    std::uniform_int_distribution<int> distribution(0,(int)(size *0.30));
////    int dice_roll = distribution(generator);  // generates number in the range 1..6

    vector<pair<int,int>> vAddedEdgesUnused;

#if 0
    while (loops < 10000) {
        vector<int> vRemoved;
        for (int i=0; i < 3; ++i) {
            int vertexToRemove = vVertices[distribution(generator)];
            if (!isolates.GetInGraph().Contains(vertexToRemove)) {
                i--;
                continue;
            }
            vRemoved.push_back(vertexToRemove);
            isolates.RemoveVertexAndNeighbors(vertexToRemove, vRemoved);
////            isolates.RemoveVertex(vertexToRemove);
        }

        isolates.RemoveAllIsolates(0, vRemoved, vRemoved, vAddedEdgesUnused, false);
        vector<vector<int>> vComponents;
        ComputeConnectedComponents(isolates, vComponents);
        size_t biggest(0);
        for (vector<int> const &vComponent : vComponents) {
            biggest = max(vComponent.size(), biggest);
        }
        if (biggest < best) {
            best = biggest;
            cout << "best so far: max component of size " << best << endl;
        }
        loops++;
        isolates.ReplaceAllRemoved(vRemoved);
    }
#endif //0
#if 0
    int savedi(0), savedj(0), savedk(0);
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

        if (biggest < best) {
            best = biggest;
            cout << "best so far: remove " << i << " for max component of size " << best << endl;
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

    for (size_t i = 0; i < size*0.10; i++) {
        vector<int> v_iNeighbors;
        isolates.RemoveVertexAndNeighbors(vVertices[i], v_iNeighbors);
        for (size_t j = i+1; j < size*0.10; j++) {
            vector<int> v_jNeighbors;
            isolates.RemoveVertexAndNeighbors(vVertices[j], v_jNeighbors);
            isolates.RemoveAllIsolates(0, v_jNeighbors, v_jNeighbors, vAddedEdgesUnused, false);

                vector<vector<int>> vComponents;
                ComputeConnectedComponents(isolates, vComponents);
                size_t biggest(0);
                for (vector<int> const &vComponent : vComponents) {
                    biggest = max(vComponent.size(), biggest);
                }

                if (biggest < best) {
                    best = biggest;
                    cout << "best so far: remove " << i << " " << j << " for max component of size " << best << endl;
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
////        if (loops % 5000 == 0) { break; }
////            cout << "Removing " << savedi << " " << savedj << " " << savedk << endl;
////            vector<int> vRemovedNeighbors;
////            isolates.RemoveVertexAndNeighbors(savedi, vRemovedNeighbors);
////            isolates.RemoveVertexAndNeighbors(savedj, vRemovedNeighbors);
////            isolates.RemoveVertexAndNeighbors(savedk, vRemovedNeighbors);
////        }
////    }

    for (size_t i = 0; i < size*0.05; i++) {
        vector<int> v_iNeighbors;
        isolates.RemoveVertexAndNeighbors(vVertices[i], v_iNeighbors);
        for (size_t j = i+1; j < size*0.05; j++) {
            vector<int> v_jNeighbors;
            isolates.RemoveVertexAndNeighbors(vVertices[j], v_jNeighbors);
            for (size_t k = j+1; k < size*0.05; k++) {
                loops++;
////                if (loops %10 == 0) { break; }
                vector<int> v_kNeighbors;
                isolates.RemoveVertexAndNeighbors(vVertices[k], v_kNeighbors);
                isolates.RemoveAllIsolates(0, v_kNeighbors, v_kNeighbors, vAddedEdgesUnused, false);

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
////            if (loops % 100 == 0) break;
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


    // try removing high-impact reduction vertices at first, then run regular algorithm on resulcting connected
    // components.
////    vector<pair<int,int>> vAddedEdgesUnused;
////    vector<int> vRemoved1;
    vector<int> vCliqueVerticesPersistent1;
    vector<int> vRemovedPersistent1;
    list<list<int>> cliques;
    cliques.push_back(list<int>());


    auto ChooseNextVertex = [&vAddedEdgesUnused](Isolates4<SparseArraySet> &theIsolates)
    {
        size_t bestComponentSize(string::npos);
        int    bestVertex(-1);
        vector<int> vToConsider(theIsolates.GetInGraph().begin(), theIsolates.GetInGraph().end());
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

            if (uSizeOfLargestComponent < bestComponentSize) {
                bestComponentSize = uSizeOfLargestComponent;
////                cout << "best so far: remove " << i << " for max component of size " << bestComponentSize << endl;
                bestVertex = i;
            }

            v_iNeighbors.push_back(i);
            theIsolates.ReplaceAllRemoved(v_iNeighbors);
        }

////        cout << "best vertex: " << bestVertex  << " for max component of size " << bestComponentSize << endl;
        return bestVertex;
    };

    IndependentSetColoringStrategy coloringStrategy(adjacencyMatrix);

    vector<int> vOrdering;
    vector<int> vColoring;

    size_t cliqueSize(0);
    OrderingTools::InitialOrderingMISR(m_AdjacencyList, isolates, vOrdering, vColoring, cliqueSize);

    vector<int> vAdjunctOrdering = vOrdering;

    if (cliqueSize > 0) {
        cliques.back() = list<int>(vOrdering.begin(), vOrdering.begin() + cliqueSize);
    }

////    coloringStrategy.Recolor(adjacencyMatrix, vAdjunctOrdering, vOrdering, vColoring, cliqueSize, cliqueSize);

    while (isolates.GetInGraph().Size() > 913) {
////    while (!isolates.GetInGraph().Empty()) {
        int const nextVertex1 = ChooseNextVertex(isolates);
        if (nextVertex1 == -1) break;
////        if (vColoring.back() <= cliques.back().size()) break;

        cout << "Vertices remaining: " << isolates.GetInGraph().Size() << endl << flush;

        vector<int> vRemoved1;
        vector<int> vCliqueVertices1; vCliqueVertices1.push_back(nextVertex1);
        isolates.RemoveVertexAndNeighbors(nextVertex1, vRemoved1);
////        vRemoved1.push_back(nextVertex1);
        isolates.RemoveAllIsolates(0, vCliqueVertices1, vRemoved1, vAddedEdgesUnused, false);

        vector<int> currentClique(vCliqueVertices1.begin(), vCliqueVertices1.end());
        currentClique.insert(currentClique.end(), vCliqueVerticesPersistent1.begin(), vCliqueVerticesPersistent1.end());

        size_t numLeft = vColoring.size();
        for (; numLeft > 0; --numLeft) {
            if (currentClique.size() + vColoring[numLeft-1] <= cliques.back().size()) { break; }
        }
        numLeft = vColoring.size() - numLeft;
        cout << ", P.left = " << numLeft << endl;

#ifdef TWO_LEVEL
        vector<int> vCliqueVerticesPersistent2;
        vector<int> vRemovedPersistent2;

        vector<int> const vToReplace(isolates.GetInGraph().begin(), isolates.GetInGraph().end());

        while (!isolates.GetInGraph().Empty()) {
            int const nextVertex2 = ChooseNextVertex(isolates);
            if (nextVertex2 == -1) {
                break;
            }
////            cout << "Vertices remaining (inner loop): " << isolates.GetInGraph().Size() << endl << flush;

            vector<int> vRemoved2;
            vector<int> vCliqueVertices2; vCliqueVertices2.push_back(nextVertex2);
            isolates.RemoveVertexAndNeighbors(nextVertex2, vRemoved2);
            isolates.RemoveAllIsolates(0, vCliqueVertices2, vRemoved2, vAddedEdgesUnused, false);
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
            for (size_t index = 0; index < vNewComponents.size(); ++index) {
                ComputeOnConnectedComponent(vNewComponents[index], isolates, m_AdjacencyList, realClique, cliques, index == vNewComponents.size()-1);
            }

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
        }

////        isolates.ReplaceAllRemoved(vCliqueVerticesPersistent2);
////        isolates.ReplaceAllRemoved(vRemovedPersistent2);
        isolates.ReplaceAllRemoved(vToReplace);
#endif // TWO_LEVEL

        isolates.ReplaceAllRemoved(vCliqueVertices1);
        isolates.ReplaceAllRemoved(vRemoved1);
        isolates.RemoveVertex(nextVertex1);
        isolates.RemoveAllIsolates(0, vCliqueVerticesPersistent1, vRemovedPersistent1, vAddedEdgesUnused, false);

        vAdjunctOrdering.erase(find(vAdjunctOrdering.begin(), vAdjunctOrdering.end(), nextVertex1));
        vOrdering.resize(vAdjunctOrdering.size());
        vColoring.resize(vAdjunctOrdering.size());

        coloringStrategy.Recolor(adjacencyMatrix, vAdjunctOrdering, vOrdering, vColoring, currentClique.size(), cliques.back().size());
    }

    list<int> realClique;
    realClique.insert(realClique.end(), vCliqueVerticesPersistent1.begin(), vCliqueVerticesPersistent1.end());

    vector<int> const vRemainingVertices(isolates.GetInGraph().begin(), isolates.GetInGraph().end());
    ComputeOnConnectedComponent(vRemainingVertices, isolates, m_AdjacencyList, realClique, cliques, true);

    if (realClique.size() > cliques.back().size()) {
        cliques.back().clear();
        cliques.back() = realClique;
        cout << "Found independent set of size: " << realClique.size() << endl;
    }

    cout << "Found independent set of size: " << cliques.back().size() << endl << flush;

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
    vector<pair<int,int>> vAddedEdgesUnused;
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

#endif // 1
}
