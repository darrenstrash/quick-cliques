
// local includes
#include "CliqueTools.h"
#include "GraphTools.h"
#include "CliqueGraphAlgorithm.h"
#include "DegeneracyVertexSets.h"
#include "DegeneracyIndependentSets.h"
#include "MaximumCliqueAlgorithm.h"
#include "MaximalCliqueAlgorithm.h"

// system includes
#include <vector>
#include <list>
#include <set>
#include <iostream>

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

// the algorithm by Larson (A NOTE ON CRITICAL INDEPENDENCE REDUCTIONS)
vector<int> CliqueTools::ComputeMaximumCriticalIndependentSet(vector<vector<int>> adjacencyList)
{
    vector<int> criticalSet;
    vector<vector<int>> biDoubleGraph(std::move(GraphTools::ComputeBiDoubleGraph(adjacencyList)));

    size_t graphSize(biDoubleGraph.size());
    set<int> setRemovedVertices;
    for (int vertex = 0; vertex < adjacencyList.size(); ++vertex) {
        if (setRemovedVertices.find(vertex)!= setRemovedVertices.end()) continue;
        int const independenceNumber(graphSize - GraphTools::ComputeMaximumMatchingSize(biDoubleGraph));
        cout << "a(   graph)=" << independenceNumber << endl;

        vector<vector<int>> subgraph(biDoubleGraph);
        int const dualVertex(vertex + adjacencyList.size());

        size_t subgraphSize(graphSize);

        // remove vertex and all its neighbors
        for (int const neighbor : subgraph[vertex]) {
            subgraphSize--;
            for (int const nNeighbor : subgraph[neighbor]) {
                subgraph[nNeighbor].erase(find(subgraph[nNeighbor].begin(), subgraph[nNeighbor].end(), neighbor));
            }
            subgraph[neighbor].clear();
        }
        subgraph[vertex].clear();
        subgraphSize--;

        // remove dual vertex and all its neighbors
        for (int const dualNeighbor : subgraph[dualVertex]) {
            for (int const dualnNeighbor : subgraph[dualNeighbor]) {
                subgraph[dualnNeighbor].erase(find(subgraph[dualnNeighbor].begin(), subgraph[dualnNeighbor].end(), dualNeighbor));
            }
            subgraph[dualNeighbor].clear();
            subgraphSize--;
        }
        subgraph[dualVertex].clear();
        subgraphSize--;

        int const subgraphIndependenceNumber(subgraphSize - GraphTools::ComputeMaximumMatchingSize(subgraph));
        cout << "a(subgraph)=" << subgraphIndependenceNumber << endl;

        if (independenceNumber == subgraphIndependenceNumber) {
            criticalSet.push_back(vertex);
            biDoubleGraph = subgraph;
            graphSize = subgraphSize;

            setRemovedVertices.insert(vertex);
            setRemovedVertices.insert(biDoubleGraph[dualVertex].begin(), biDoubleGraph[dualVertex].end());
            cout << vertex << " is in" << endl;
        } else {
            // remove vertex from graph
            for (int const neighbor : subgraph[vertex]) {
                subgraph[neighbor].erase(find(subgraph[neighbor].begin(), subgraph[neighbor].end(), vertex));
            }
            biDoubleGraph[vertex].clear();
            graphSize--;

            // remove dual vertex from graph
            for (int const dualNeighbor : subgraph[dualVertex]) {
                subgraph[dualNeighbor].erase(find(subgraph[dualNeighbor].begin(), subgraph[dualNeighbor].end(), dualVertex));
            }
            biDoubleGraph[vertex].clear();
            graphSize--;

            setRemovedVertices.insert(vertex);
            cout << vertex << " is out" << endl;
        }
    }

    return criticalSet;
}
vector<int> CliqueTools::ComputeBipartiteMaximumIndependentSet(vector<vector<int>> const &biDoubleGraph)
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

    vector<int> independentSet;
    int count(0);
    for (int vertex = 0; vertex < biDoubleGraph.size()/2; ++vertex) {
        if (vertex == -1) independentSet.push_back(vertex);
    }

    return independentSet;
}

