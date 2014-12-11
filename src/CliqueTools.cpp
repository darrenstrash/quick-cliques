
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
