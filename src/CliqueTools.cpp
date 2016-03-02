
// local includes
#include "CliqueTools.h"
#include "GraphTools.h"
#include "MatchingTools.h"
#include "CliqueGraphAlgorithm.h"
#include "DegeneracyVertexSets.h"
#include "DegeneracyIndependentSets.h"
#include "MaximumCliqueAlgorithm.h"
#include "MaximalCliqueAlgorithm.h"

#include "LightWeightReductionSparseStaticOrderMISS.h"
#include "LightWeightReductionSparseFullMISS.h"

// system includes
#include <vector>
#include <list>
#include <set>
#include <iostream>

////#define TESTING

using namespace std;

void CliqueTools::ComputeCliqueGraph(vector<vector<int>> &adjacencyList, vector<vector<int>> &cliqueGraphAdjacencyList, vector<vector<int>> &vertexToClique)
{
    CliqueGraphAlgorithm cliqueGraphAlgorithm(new DegeneracyVertexSets(adjacencyList), adjacencyList);

    list<list<int>> cliques;
    cliqueGraphAlgorithm.Run(cliques);

    cliqueGraphAlgorithm.MoveCliqueGraph(cliqueGraphAdjacencyList, vertexToClique);
}

void CliqueTools::FindMaximalIndependentSetInCliqueGraph(vector<vector<int>> &adjacencyList)
{
    vector<vector<int>> cliqueGraph;
    vector<vector<int>> vertexToClique;
    ComputeCliqueGraph(adjacencyList, cliqueGraph, vertexToClique);

    list<list<int>> cliques;
    {
        MaximalCliqueAlgorithm algorithm(new DegeneracyIndependentSets(cliqueGraph));
        algorithm.Run(cliques);

        if (cliques.empty()) {
            cout << __LINE__ << ": ERROR: clique algorithm should always return a clique in a non-empty graph." << endl << flush;
            return;
        }
    }

    set<int> subgraphVertices;

////#ifdef DEBUG
    cout << "# clique vertices: " << vertexToClique.size() << endl;
////#endif

    list<int> const &clique(cliques.front());
    for (int const vertex : clique) {
////#ifdef DEBUG
        cout << vertex << ":";
////#endif
        subgraphVertices.insert(vertexToClique[vertex].begin(), vertexToClique[vertex].end());
        for (int const realVertex : vertexToClique[vertex]) {
            cout << realVertex << " ";
        }
        cout << endl;
    }
    cliques.clear();

    vector<vector<int>> inducedSubgraph;
    map<int,int> vertexMap;
    GraphTools::ComputeInducedSubgraph(adjacencyList, subgraphVertices, inducedSubgraph, vertexMap);

    MaximumCliqueAlgorithm algorithm(new DegeneracyIndependentSets(inducedSubgraph));
    algorithm.Run(cliques);
    if (!cliques.empty()) {
        cout << "Found independent set of size: " << cliques.front().size() << endl;
    }

}

void CliqueTools::DecomposeIntoDisjointCliques(vector<vector<int>> &adjacencyList, vector<vector<int>> &cliquesToReturn)
{
    set<int> verticesToExclude;
    list<list<int>> cliques;

    int iterationNumber(0);
    while (iterationNumber++ < adjacencyList.size() && verticesToExclude.size() != adjacencyList.size()) {
        MaximumCliqueAlgorithm algorithm(new DegeneracyVertexSets(adjacencyList, verticesToExclude));
        algorithm.SetQuiet(true);
        algorithm.Run(cliques);
        if (cliques.empty()) {
            cout << __LINE__ << ": ERROR: clique algorithm should always return a clique in a non-empty graph." << endl << flush;
            break;
        }

#ifdef DEBUG 
        cout << "#Maximum cliques found: " << cliques.size() << endl;
        cout << "Size of first   clique: " << cliques.front().size() << endl;
        cout << "    {";
        for (int const vertex : cliques.front()) {
            cout << vertex << " ";
        }

        cout << "}" << endl;
#endif

        vector<int> clique(cliques.front().begin(), cliques.front().end());
        verticesToExclude.insert(clique.begin(), clique.end());
        cliquesToReturn.emplace_back(std::move(clique));
#ifdef DEBUG
        cout << "Exclude : ";
        for (int const vertex : verticesToExclude) {
            cout << vertex << " ";
        }
        cout << endl;
        cout << "At most " << (adjacencyList.size() - verticesToExclude.size()) << " iterations to go!" << endl;
#endif
    }

    if (verticesToExclude.size() != adjacencyList.size()) {
        cout << __LINE__ << ": ERROR: Did not finish decomposing graph into cliques. Missed " << adjacencyList.size() - verticesToExclude.size() << " vertices." << endl << flush;
    }
 
}

bool CliqueTools::IsMaximalClique(vector<vector<int>> &adjacencyArray, list<int> const&clique, bool const verbose)
{
    size_t const cliqueSize(clique.size());
    vector<bool> vMarkedVertices(adjacencyArray.size(), false);

    for (int const vertex : clique) {
        vMarkedVertices[vertex] = true;
    }

    // first check that it is a clique
    #if 0
    for (int const vertex : clique) {
        size_t neighborsInClique(0)
        for (int const neighbor : adjacencyArray[vertex]) {
            if (vMarkedVertices[neighbor]) neighborsInClique++;
        }
        if (neighborsInClique != cliqueSize - 1) {
            if (verbose) {
                cout << "Maximal clique test failed: " << vertex << " is not in the clique!" << endl;
            }
            return false;
    }
    #else
    for (size_t vertex = 0; vertex < adjacencyArray.size(); ++vertex) {
        bool const inClique(vMarkedVertices[vertex]);
        size_t neighborsInClique(0);
        for (int const neighbor : adjacencyArray[vertex]) {
            if (vMarkedVertices[neighbor]) neighborsInClique++;
        }

        if (inClique && neighborsInClique != cliqueSize - 1) {
            if (verbose) {
                cout << "Maximal clique test failed: " << vertex << " should not be in clique!" << endl;
            }
            return false;
        }

        if (!inClique && neighborsInClique == cliqueSize) {
            if (verbose) {
                cout << "Maximal clique test failed: " << vertex << " can be added to clique!" << endl;
            }
            return false;
        }
    }
    #endif
    return true;
}

bool CliqueTools::IsClique(vector<vector<char>> &adjacencyMatrix, list<int> const&clique, bool const verbose)
{
    for (int const vertex : clique) {
        for (int const otherVertex : clique) {
            if (vertex == otherVertex) continue;
            if (!adjacencyMatrix[vertex][otherVertex]) {
                cout << "clique test failed: " << vertex << " should not be in the same set with " << otherVertex << endl << flush;
                return false;
            }
        }
    }

    return true;
}

bool CliqueTools::IsIndependentSet(vector<vector<char>> &adjacencyMatrix, list<int> const&clique, bool const verbose)
{
    for (int const vertex : clique) {
        for (int const otherVertex : clique) {
            if (vertex == otherVertex) continue;
            if (adjacencyMatrix[vertex][otherVertex]) {
                cout << "Independent set test failed: " << vertex << " should not be in the same set with " << otherVertex << endl << flush;
                return false;
            }
        }
    }

    return true;
}


bool CliqueTools::IsMaximalIndependentSet(vector<vector<int>> &adjacencyArray, list<int> const &vertexSet, bool const verbose)
{
    size_t const setSize(vertexSet.size());
    vector<bool> vMarkedVertices(adjacencyArray.size(), false);

    for (size_t const vertex : vertexSet) {
        vMarkedVertices[vertex] = true;
    }

    // first check that it is a clique
    #if 0
    for (int const vertex : clique) {
        size_t neighborsInClique(0)
        for (int const neighbor : adjacencyArray[vertex]) {
            if (vMarkedVertices[neighbor]) neighborsInClique++;
        }
        if (neighborsInClique != cliqueSize - 1) {
            if (verbose) {
                cout << "Maximal clique test failed: " << vertex << " is not in the clique!" << endl;
            }
            return false;
    }
    #else
    for (size_t vertex = 0; vertex < adjacencyArray.size(); ++vertex) {
        bool const inSet(vMarkedVertices[vertex]);
        size_t neighborsInSet(0);
        for (int const neighbor : adjacencyArray[vertex]) {
            if (vMarkedVertices[neighbor]) {
                neighborsInSet++;
                if (inSet) {
                    cout << "Maximal independent set test failed: " << vertex << " should not be in set with neighbor " << neighbor << "!" << endl;
                }
            }
        }

        if (inSet && neighborsInSet != 0) {
////            if (verbose) {
////                cout << "Maximal independent set test failed: " << vertex << " should not be in set!" << endl;
////            }
            return false;
        }

        if (!inSet && neighborsInSet == 0) {
            if (verbose) {
                cout << "Maximal independent set test failed: " << vertex << " can be added to set!" << endl;
            }
            return false;
        }
    }
    #endif
    return true;
}

bool CliqueTools::IsIndependentSet(vector<vector<int>> &adjacencyArray, list<int> const &vertexSet, bool const verbose)
{
    size_t const setSize(vertexSet.size());
    vector<bool> vMarkedVertices(adjacencyArray.size(), false);

    for (size_t const vertex : vertexSet) {
        vMarkedVertices[vertex] = true;
    }

    // first check that it is a clique
    #if 0
    for (int const vertex : clique) {
        size_t neighborsInClique(0)
        for (int const neighbor : adjacencyArray[vertex]) {
            if (vMarkedVertices[neighbor]) neighborsInClique++;
        }
        if (neighborsInClique != cliqueSize - 1) {
            if (verbose) {
                cout << "Maximal clique test failed: " << vertex << " is not in the clique!" << endl;
            }
            return false;
    }
    #else
    for (size_t vertex = 0; vertex < adjacencyArray.size(); ++vertex) {
        bool const inSet(vMarkedVertices[vertex]);
        size_t neighborsInSet(0);
        for (int const neighbor : adjacencyArray[vertex]) {
            if (vMarkedVertices[neighbor]) {
                neighborsInSet++;
                if (inSet) {
                    cout << "Maximal independent set test failed: " << vertex << " should not be in set with neighbor " << neighbor << "!" << endl;
                }
            }
        }

        if (inSet && neighborsInSet != 0) {
////            if (verbose) {
////                cout << "Maximal independent set test failed: " << vertex << " should not be in set!" << endl;
////            }
            return false;
        }

////        if (!inSet && neighborsInSet == 0) {
////            if (verbose) {
////                cout << "Maximal independent set test failed: " << vertex << " can be added to set!" << endl;
////            }
////            return false;
////        }
    }
    #endif
    return true;
}

// Ageev's algorithm (See Butenko and Trukhanov)
set<int> CliqueTools::ComputeCriticalIndependentSet(vector<vector<int>> const &adjacencyList)
{
    clock_t start_time(clock());
    vector<vector<int>> const biDoubleGraph(std::move(GraphTools::ComputeBiDoubleGraph(adjacencyList)));

    int const graphSize(biDoubleGraph.size());
    cout << "bi-double graph size=" << biDoubleGraph.size() << endl;
    set<int> setRemovedVertices;
    set<int> setRemainingVertices;
    set<int> setVerticesToEvaluate;
    for (int vertex = 0; vertex < biDoubleGraph.size(); ++vertex) {
        setRemainingVertices.insert(vertex);
        setVerticesToEvaluate.insert(vertex);
    }

////    GraphTools::PrintGraphInSNAPFormat(biDoubleGraph);

#ifndef TESTING
    set<int> criticalSet(std::move(MatchingTools::ComputeLeftMIS(biDoubleGraph)));
    cout << "Critical             set found: " << criticalSet.size() << endl << flush;

#else
    clock_t start_time1(clock());
    set<int> criticalSet(std::move(MatchingTools::ComputeLeftMIS(biDoubleGraph)));
    clock_t end_time1(clock());
    cout << "Time to compute critical set  : " << Tools::GetTimeInSeconds(end_time1 - start_time1) << endl << flush;
    cout << "Critical             set found: " << criticalSet.size() << endl << flush;
    clock_t start_time2(clock());
    set<int> const criticalSet2(std::move(MatchingTools::ComputeCriticalSet(adjacencyList)));
    clock_t end_time2(clock());
    cout << "Critical           set 2 found: " << criticalSet2.size() << endl << flush;
    cout << "Time to compute critical set 2: " << Tools::GetTimeInSeconds(end_time2 - start_time2) << endl << flush;

    if (criticalSet.size() != criticalSet2.size()) {
        cout << "ERROR! Critical sets are different!" << endl << flush;
    }
    for (int const vertex : criticalSet) {
        if (criticalSet2.find(vertex) == criticalSet2.end()) {
            cout << "ERROR! Critical set 2 does not contain " << vertex << "!" << endl << flush;
        }
    }
#endif // TESTING

    set<int> toRemove;
    for (int const vertex : criticalSet) {
        for (int const neighbor : adjacencyList[vertex]) {
            set<int>::iterator const it(criticalSet.find(neighbor));
            if (it != criticalSet.end()) {
                toRemove.insert(neighbor);
            }
        }
    }

    set<int> criticalIndependentSet(std::move(criticalSet));
    for (int const vertexToRemove: toRemove) {
        criticalIndependentSet.erase(criticalIndependentSet.find(vertexToRemove));
    }

    clock_t end_time(clock());

    cout << "Critical independent set found: " << criticalIndependentSet.size() << endl << flush;
////    size_t const independenceNumber(graphSize - maximumMatching.size());
    cout << "Time to compute  with matching: " << Tools::GetTimeInSeconds(end_time - start_time) << endl << flush;

    return criticalIndependentSet;
}

// Iterate critical independent set computation (this is a critical set reduction)
// returns vertices in kernel.
set<int> CliqueTools::IterativelyRemoveCriticalIndependentSets(vector<vector<int>> const &adjacencyList, set<int> &independentVertices)
{
    independentVertices.clear();
    clock_t start_time(clock());
    vector<vector<int>> biDoubleGraph(std::move(GraphTools::ComputeBiDoubleGraph(adjacencyList)));

    BiDoubleGraph biDouble(adjacencyList);
    vector<int> matching(biDouble.Size(), -1);

    int const biDoubleGraphSize(biDoubleGraph.size());
////    cout << "bi-double graph size=" << biDoubleGraph.size() << endl;
    set<int> setRemainingBiDoubleVertices;
    vector<bool> vInBiDoubleGraph(biDoubleGraph.size(), true);
#ifdef DEBUG
    size_t numberOfSingletons(0);
    cout << "Singletons: ";
    for (int vertex=0; vertex < adjacencyList.size(); ++vertex) {
        vector<int> const & neighbors(adjacencyList[vertex]);
        if (neighbors.empty()) {
            cout << vertex << " ";
            numberOfSingletons++;
        }
    }
    cout << endl << flush;
    cout << "Number of Singletons: " << numberOfSingletons << endl << flush;
#endif //DEBUG
    for (int vertex = 0; vertex < biDoubleGraph.size(); ++vertex) {
        setRemainingBiDoubleVertices.insert(vertex);
    }

////    GraphTools::PrintGraphInSNAPFormat(biDoubleGraph);

    // remove a vertex and its twin, if they exist
    auto removeFromBiDouble = [&biDoubleGraph, &matching, &setRemainingBiDoubleVertices, &vInBiDoubleGraph, &biDoubleGraphSize] (int const vertex) {
        set<int>::iterator it = setRemainingBiDoubleVertices.find(vertex);
        if (it != setRemainingBiDoubleVertices.end()) {
            biDoubleGraph[vertex].clear();
            biDoubleGraph[vertex + biDoubleGraphSize/2].clear();

            vInBiDoubleGraph[vertex] = false;
            vInBiDoubleGraph[vertex + biDoubleGraphSize/2] = false;

            setRemainingBiDoubleVertices.erase(it);
            it = setRemainingBiDoubleVertices.find(vertex + biDoubleGraphSize/2);
            setRemainingBiDoubleVertices.erase(it);
            if (matching[vertex] != -1) {
                matching[matching[vertex]] = -1;
                matching[vertex] = -1;
            }
            if (matching[vertex + biDoubleGraph.size()/2] != -1) {
                matching[matching[vertex + biDoubleGraph.size()/2]] = -1;
                matching[vertex + biDoubleGraph.size()/2] = -1;
            }
            return true;
        }
        return false;
    };

    int previousGraphSize(adjacencyList.size());
    int newGraphSize(-1);

    while (previousGraphSize != newGraphSize) {
#ifdef DEBUG
        vector<vector<int>> biDoubleSubgraph;
        map<int,int> unusedRemapping;
        GraphTools::ComputeInducedSubgraph(biDoubleGraph, setRemainingBiDoubleVertices, biDoubleSubgraph, unusedRemapping);
        cout << "---Start BiDouble Subgraph---" << endl << flush;
        GraphTools::PrintGraphInSNAPFormat(biDoubleSubgraph);
        cout << "---End   BiDouble Subgraph---" << endl << flush;
#endif // DEBUG
        
        previousGraphSize = newGraphSize;
        set<int> criticalSet(std::move(MatchingTools::ComputeLeftMISOptimizedWithMatching(biDouble, matching, vInBiDoubleGraph, setRemainingBiDoubleVertices)));
////        cout << "Critical          set found: " << criticalSet.size() << endl << flush;

        set<int> toRemove;
        for (int const vertex : criticalSet) {
            for (int const neighbor : adjacencyList[vertex]) {
                set<int>::iterator const it(criticalSet.find(neighbor));
                if (it != criticalSet.end()) {
                    toRemove.insert(neighbor);
                }
            }
        }

        set<int> criticalIndependentSet(std::move(criticalSet));
        for (int const vertexToRemove: toRemove) {
            criticalIndependentSet.erase(criticalIndependentSet.find(vertexToRemove));
        }

////        cout << "Num critical indepset nodes: " << criticalIndependentSet.size() << endl << flush;

        independentVertices.insert(criticalIndependentSet.begin(), criticalIndependentSet.end());
#if DEBUG
        cout << "indepset: ";
        for (int const vertex : criticalIndependentSet) {
            cout << vertex << " ";
        }
        cout << endl;
#endif

        size_t removedFromBiDoubleCount(0);
        size_t sumOfNeighbors(0);

        for (int const vertex : criticalIndependentSet) {
            if (removeFromBiDouble(vertex)) {
                removedFromBiDoubleCount++;
            }
            for (int const neighbor : adjacencyList[vertex]) {
                sumOfNeighbors++;
                if (removeFromBiDouble(neighbor)) {
                    removedFromBiDoubleCount++;
                }
            }
        }

////        cout << "Num vertices removed        : " << removedFromBiDoubleCount << endl << flush;
////        cout << "Sum of neighbors            : " << sumOfNeighbors << endl << flush;

        newGraphSize = setRemainingBiDoubleVertices.size();
        if (newGraphSize == 0) {
////            cout << "Kernel is empty, breaking..." << endl << flush;
            break;
        }
    }

    clock_t end_time(clock());

////    cout << "Critical set kernel size: " << setRemainingBiDoubleVertices.size()/2 << endl << flush;
////    cout << "Independent set    nodes: " << independentSetNodes.size() << endl << flush;
////    size_t const independenceNumber(graphSize - maximumMatching.size());
////    cout << "Time to compute         : " << Tools::GetTimeInSeconds(end_time - start_time) << endl << flush;

    return setRemainingBiDoubleVertices;
}

// Larson algorithm
set<int> CliqueTools::IterativelyRemoveMaximumCriticalIndependentSets(vector<vector<int>> const &adjacencyList, set<int> &independentVertices)
{
    clock_t start_time(clock());
    vector<vector<int>> biDoubleGraph(std::move(GraphTools::ComputeBiDoubleGraph(adjacencyList)));

    BiDoubleGraph biDouble(adjacencyList);
    vector<int> matching(biDouble.Size(), -1);

    int const biDoubleGraphSize(biDoubleGraph.size());
////    cout << "bi-double graph size=" << biDoubleGraph.size() << endl;
    set<int> setRemainingBiDoubleVertices;
    set<int> setRemainingVertices;
    vector<bool> vInBiDoubleGraph(biDoubleGraph.size(), true);
#ifdef DEBUG
    size_t numberOfSingletons(0);
    cout << "Singletons: ";
    for (int vertex=0; vertex < adjacencyList.size(); ++vertex) {
        vector<int> const & neighbors(adjacencyList[vertex]);
        if (neighbors.empty()) {
            cout << vertex << " ";
            numberOfSingletons++;
        }
    }
    cout << endl << flush;
    cout << "Number of Singletons: " << numberOfSingletons << endl << flush;
#endif //DEBUG
    for (int vertex = 0; vertex < biDoubleGraph.size(); ++vertex) {
        setRemainingBiDoubleVertices.insert(vertex);
        if (vertex < adjacencyList.size()) setRemainingVertices.insert(vertex);
    }

////    GraphTools::PrintGraphInSNAPFormat(biDoubleGraph);

    // remove a vertex and its twin, if they exist
    auto removeFromBiDouble = [&biDoubleGraph, &matching, &setRemainingBiDoubleVertices, &setRemainingVertices, &vInBiDoubleGraph, &biDoubleGraphSize] (int const vertex) {
        set<int>::iterator it = setRemainingBiDoubleVertices.find(vertex);
        if (it != setRemainingBiDoubleVertices.end()) {
////            biDoubleGraph[vertex].clear();
////            biDoubleGraph[vertex + biDoubleGraphSize/2].clear();

            vInBiDoubleGraph[vertex] = false;
            vInBiDoubleGraph[vertex + biDoubleGraphSize/2] = false;

            setRemainingVertices.erase(setRemainingVertices.find(vertex));
            setRemainingBiDoubleVertices.erase(it);
            it = setRemainingBiDoubleVertices.find(vertex + biDoubleGraphSize/2);
            setRemainingBiDoubleVertices.erase(it);
            if (matching[vertex] != -1) {
                matching[matching[vertex]] = -1;
                matching[vertex] = -1;
            }
            if (matching[vertex + biDoubleGraph.size()/2] != -1) {
                matching[matching[vertex + biDoubleGraph.size()/2]] = -1;
                matching[vertex + biDoubleGraph.size()/2] = -1;
            }
            return true;
        }
        return false;
    };

    int previousGraphSize(-1);
    int newGraphSize(adjacencyList.size());

    set<int> independentSetNodes;

    size_t misNodeCount(0);

    while (previousGraphSize != newGraphSize) {
#ifdef DEBUG
        vector<vector<int>> biDoubleSubgraph;
        map<int,int> unusedRemapping;
        GraphTools::ComputeInducedSubgraph(biDoubleGraph, setRemainingBiDoubleVertices, biDoubleSubgraph, unusedRemapping);
        cout << "---Start BiDouble Subgraph---" << endl << flush;
        GraphTools::PrintGraphInSNAPFormat(biDoubleSubgraph);
        cout << "---End   BiDouble Subgraph---" << endl << flush;
#endif // DEBUG
        
        previousGraphSize = newGraphSize;

        set<int> maximumCriticalIndependentSet;

        set<int> setSavedRemainingVertices(setRemainingVertices);

        vector<int> removed;
        // compute maximum Critical independent set
        for (int const vertex : setSavedRemainingVertices) {
////            cout << "Add vertex " << vertex << " to MCIS? ";

            // not in biDoubleGraph
            if (setRemainingBiDoubleVertices.find(vertex) == setRemainingBiDoubleVertices.end()) {
////                cout << "NO, SKIP" << endl << flush;
                continue;
            }

            // Compute actual independent set, not just independent set on left side
#ifndef TESTING
            set<int> const criticalSet(std::move(MatchingTools::ComputeBiDoubleMISOptimizedWithMatching(biDouble, matching, vInBiDoubleGraph, setRemainingBiDoubleVertices)));
#else
            set<int> const criticalSet(std::move(MatchingTools::ComputeBiDoubleMIS(biDoubleGraph, vInBiDoubleGraph, setRemainingBiDoubleVertices)));
#endif // TESTING

            vector<bool> vNewInBiDoubleGraph(vInBiDoubleGraph);
            set<int> setNewRemainingBiDoubleVertices(setRemainingBiDoubleVertices);

            auto removeFromNewBiDouble = [&biDoubleGraph, &matching, &vNewInBiDoubleGraph, &setNewRemainingBiDoubleVertices] (int const vertexToRemove) {
                set<int>::iterator it = setNewRemainingBiDoubleVertices.find(vertexToRemove);
                if (it != setNewRemainingBiDoubleVertices.end()) {
                    setNewRemainingBiDoubleVertices.erase(it);
                    setNewRemainingBiDoubleVertices.erase(setNewRemainingBiDoubleVertices.find(vertexToRemove + biDoubleGraph.size()/2));

                    vNewInBiDoubleGraph[vertexToRemove] = false;
                    vNewInBiDoubleGraph[vertexToRemove + biDoubleGraph.size()/2] = false;
                    if (matching[vertexToRemove] != -1) {
                        matching[matching[vertexToRemove]] = -1;
                        matching[vertexToRemove] = -1;
                    }
                    if (matching[vertexToRemove + biDoubleGraph.size()/2] != -1) {
                        matching[matching[vertexToRemove + biDoubleGraph.size()/2]] = -1;
                        matching[vertexToRemove + biDoubleGraph.size()/2] = -1;
                    }
                    return true;
                }
                return false;
            };

            removeFromNewBiDouble(vertex);

            // remove neighbors from bidouble graph
            for (int const neighbor : adjacencyList[vertex]) {
                removeFromNewBiDouble(neighbor);
            }

            // compute new critical set size

#ifndef TESTING
            set<int> const newCriticalSet(std::move(MatchingTools::ComputeBiDoubleMISOptimizedWithMatching(biDouble, matching, vNewInBiDoubleGraph, setNewRemainingBiDoubleVertices)));
#else
            set<int> const newCriticalSet(std::move(MatchingTools::ComputeBiDoubleMIS(biDoubleGraph, vNewInBiDoubleGraph, setNewRemainingBiDoubleVertices)));
            set<int> const newCriticalSet2(std::move(MatchingTools::ComputeBiDoubleMISOptimizedWithMatching(biDouble, matching, vNewInBiDoubleGraph, setNewRemainingBiDoubleVertices)));

            if (newCriticalSet.size() != newCriticalSet2.size()) {
                cout << "ERROR! Critical sets are different!" << endl << flush;
            }
            for (int const vertex : newCriticalSet) {
                if (newCriticalSet2.find(vertex) == newCriticalSet2.end()) {
                    cout << "ERROR! Critical set 2 does not contain " << vertex << "!" << endl << flush;
                }
            }
#endif // TESTING

            bool const addToIndependentSet(criticalSet.size() == (newCriticalSet.size() + 2));
////            cout << "set size comparison = " << criticalSet.size() << ", " << newCriticalSet.size() << endl << flush;

#ifdef VERIFY
            map<int,int> unusedMapping;
            vector<vector<int>> subgraph1;
            GraphTools::ComputeInducedSubgraph(biDoubleGraph, setRemainingBiDoubleVertices, subgraph1, unusedMapping);
            LightWeightReductionSparseFullMISS algorithm1(subgraph1);
            list<list<int>> cliques1;
            algorithm1.Run(cliques1);
            cout << "First graph, MIS=" << cliques1.back().size() << endl << flush;

            vector<vector<int>> subgraph2;
            GraphTools::ComputeInducedSubgraph(biDoubleGraph, setNewRemainingBiDoubleVertices, subgraph2, unusedMapping);
            LightWeightReductionSparseFullMISS algorithm2(subgraph2);
            list<list<int>> cliques2;
            algorithm2.Run(cliques2);
            cout << "Secnd graph, MIS=" << cliques2.back().size() << endl << flush;
#endif // VERIFY

            if (addToIndependentSet) {
////                cout << "YES" << endl << flush;
                maximumCriticalIndependentSet.insert(vertex);
                misNodeCount++;
                // remove vertex and neighbors from bidoublegraph
                //setRemainingBiDoubleVertices = std::move(setNewRemainingBiDoubleVertices);
                removeFromBiDouble(vertex);
                removed.push_back(vertex);
                for (int const neighbor : adjacencyList[vertex]) {
                    removeFromBiDouble(neighbor);
                    removed.push_back(neighbor);
                }
            } else {
////                cout << "NO" << endl << flush;
                removeFromBiDouble(vertex);
                removed.push_back(vertex);
            }
////        cout << "Critical          set found: " << criticalSet.size() << endl << flush;
////            cout << "Vertices Remaining     : " << setRemainingVertices.size() << endl << flush;
////            cout << "BiDouble Remaining     : " << setRemainingBiDoubleVertices.size() << endl << flush;
        }

////        cout << "Num critical indepset nodes: " << maximumCriticalIndependentSet.size() << endl << flush;
////        cout << "Critical  #  indepset nodes: " << misNodeCount << endl << flush;

        size_t removedFromBiDoubleCount(0);
        size_t sumOfNeighbors(0);

        for (int const vertex : removed) {
            setRemainingVertices.insert(vertex);
            setRemainingBiDoubleVertices.insert(vertex);
            setRemainingBiDoubleVertices.insert(vertex + adjacencyList.size());
            vInBiDoubleGraph[vertex] = true;
            vInBiDoubleGraph[vertex + adjacencyList.size()] = true;
        }

        for (int const vertex : maximumCriticalIndependentSet) {
            if (removeFromBiDouble(vertex)) {
                removedFromBiDoubleCount++;
            }
            for (int const neighbor : adjacencyList[vertex]) {
                sumOfNeighbors++;
                if (removeFromBiDouble(neighbor)) {
                    removedFromBiDoubleCount++;
                }
            }
        }

////        cout << "Num vertices removed        : " << removedFromBiDoubleCount << endl << flush;
////        cout << "Sum of neighbors            : " << sumOfNeighbors << endl << flush;

        independentVertices.insert(maximumCriticalIndependentSet.begin(), maximumCriticalIndependentSet.end());

        newGraphSize = setRemainingBiDoubleVertices.size();
        if (newGraphSize == 0) {
////            cout << "Kernel is empty, breaking..." << endl << flush;
            break;
        }

    }

    clock_t end_time(clock());

////    cout << "Critical set kernel size: " << setRemainingBiDoubleVertices.size()/2 << endl << flush;
////    cout << "Independent set    nodes: " << independentSetNodes.size() << endl << flush;
////    size_t const independenceNumber(graphSize - maximumMatching.size());
////    cout << "Time to compute         : " << Tools::GetTimeInSeconds(end_time - start_time) << endl << flush;

    return setRemainingBiDoubleVertices;
}

