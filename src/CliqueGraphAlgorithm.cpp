/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

// local includes
#include "CliqueGraphAlgorithm.h"
#include "MinimumCliqueAlgorithm.h"
#include "MaximumCliqueAlgorithm.h"
#include "MaximalCliqueAlgorithm.h"
#include "DegeneracyIndependentSets.h"
#include "Tools.h"
#include "Algorithm.h"

// system includes
#include <limits.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <list>
#include <vector>
#include <map>

using namespace std;

/*! \file CliqueGraphAlgorithm.cpp

    \brief This file contains the main algorithm for listing all cliques
           according to the algorithm of Tomita et al. (TCS 2006) with one
           crucial change: the input graph is represented as an adjacency list.

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly

    See the main algorithm's description in 
    http://dx.doi.org/10.1016/j.tcs.2006.06.015, and a description of this
    variant of the algorithm in http://dx.doi.org/10.1007/978-3-642-20662-7_31

    This is a recursive backtracking algorithm that maintains three 
    sets of vertices, R, a partial clique, P, the common neighbors
    of vertices in R that are candidates to add to the partial clique,
    and X, the set of common neighbors of R that have been listed 
    in a maximal clique with R already.

    The algorithm recursively adds vertices to R from P, then 
    updates the sets P and X to be the new common neighbors of R
    and recurses. When P and X are empty, R is a maximal clique,
    and is reported.

    Updating the sets P and X is done by iterating over the
    neighbors of the new vertex v added to R and testing
    if they are in P or X. Neighbors of v remain in their
    respective sets.

*/

static unsigned long largestDifference(0);
static unsigned long numLargeJumps;
static unsigned long stepsSinceLastReportedClique(0);

CliqueGraphAlgorithm::CliqueGraphAlgorithm(VertexSets *pSets, vector<vector<int>> &adjacencyArray)
 : Algorithm("compute-clique-graph-" + pSets->GetName())
 , m_pSets(pSets)
 , m_NeighborSets()
 , m_CliqueGraph()
 , m_CliquesWithVertex(pSets->GetGraphSize())
 , m_AdjacencyArray(adjacencyArray)
{
}

CliqueGraphAlgorithm::~CliqueGraphAlgorithm()
{
    delete m_pSets; m_pSets = nullptr;
}

size_t maxMin2Neighborhood(0);

/*! \brief List all maximal cliques in a given graph using the algorithm
           by Tomita et al. (TCS 2006), modified to use an adjacency list
           representation of the graph instead of an adjacency matrix. 
 
    \param cliques A linked list of cliques to return. <b>(only available when compiled 
                   with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

    \return The number of maximal cliques of the input graph.
*/

long CliqueGraphAlgorithm::Run(list<list<int>> &cliques)
{
    long cliqueCount = 0;
    list<int> partialClique;

    m_pSets->Initialize();

    m_CliquesWithVertex.resize(m_pSets->GetGraphSize());

    while (m_pSets->GetNextTopLevelPartition()) {
        m_pSets->GetTopLevelPartialClique(partialClique);
        RunRecursive(cliqueCount, cliques, partialClique);
        partialClique.clear();
    }

////    cerr << "Largest Difference : " << largestDifference << endl;
////    cerr << "Num     Differences: " << numLargeJumps << endl;

    m_CliquesWithVertex.clear();

    cerr << "Num     Vertices   : " << m_NeighborSets.size() << endl;

    size_t numEdges(0);

    for (set<int> const &neighbors : m_NeighborSets) {
        numEdges += neighbors.size();
    }

    numEdges = numEdges/2;
    cerr << "Num     Edges      : " << numEdges << endl;

    cout << "Maximum min-two-neighborhood: " << maxMin2Neighborhood << endl;

    m_CliqueGraph.resize(m_NeighborSets.size());
    for (int i = 0; i < m_NeighborSets.size(); ++i) {
        m_CliqueGraph[i].insert(m_CliqueGraph[i].end(), m_NeighborSets[i].begin(), m_NeighborSets[i].end());
        m_NeighborSets[i].clear();
    }

#if 0
    cout << "Graph:" << endl;
    cout << m_NeighborSets.size() << endl;
    cout << numEdges*2 << endl;
    for (int vertex = 0; vertex < m_AdjacencyList.size(); ++vertex) {
        for (int const neighbor : m_AdjacencyList[vertex]) {
            cout << vertex << "," << neighbor;
            cout << "{ ";
            for (int const cliqueVertex : m_VerticesInClique[vertex]) {
                cout << cliqueVertex << " ";
            }
            cout << "} { ";
            for (int const cliqueVertex : m_VerticesInClique[neighbor]) {
                cout << cliqueVertex << " ";
            }
            cout << "}" << endl;
        }
    }
#endif

    m_NeighborSets.clear();

#if 1
    int const cliqueGraphDegeneracy(computeDegeneracy(m_CliqueGraph, m_CliqueGraph.size()));
    cout << "Degeneracy of clique graph: " << cliqueGraphDegeneracy << endl;

    Algorithm *pAlgorithm = new MaximalCliqueAlgorithm(new DegeneracyIndependentSets(m_CliqueGraph));

    pAlgorithm->Run(cliques);

    delete pAlgorithm; pAlgorithm = nullptr;

    vector<int> vReducedVertices;
    map<int,int> vMapReducedToOriginal;

    vector<vector<int>> vSuperVertices;

    if (!cliques.empty()) {
        cout << "Choosing one vertex from each line produces an independent set : " << endl;
        std::list<int> const &clique = cliques.front();
        for (int const cliqueVertex : clique) {
            cout << cliqueVertex << " : ";
            vector<int> verticesInClique;
            for (int const realVertex : m_VerticesInClique[cliqueVertex]) {
                vMapReducedToOriginal[realVertex] = vReducedVertices.size();
                vReducedVertices.push_back(realVertex);
                cout << realVertex << " ";
                verticesInClique.push_back(realVertex);
            }
            vSuperVertices.emplace_back(std::move(verticesInClique));
            cout << endl;
        }
    }

////    int largest(0);

#if 0
    vector<int>  vWhichPicked; vWhichPicked.reserve(vSuperVertices.size());
    vector<bool> vWhichEliminated(vSuperVertices.size(), false);
    vector<bool> vWhichVertexEliminated(m_pSets->GetGraphSize(), false);

    cout << "End  :";
    for (vector<int> const &vValue : vSuperVertices) {
        cout << vValue.size() - 1 << " ";
    }
    cout << endl;

    // consider all possible sets.
    // TODO: support skipping eliminated neighbors between sets.
    int bottom(-1);
    int current(0);
    vWhichPicked.push_back(-1);
    while (bottom != vSuperVertices.size() - 1) {

        if (current == vSuperVertices.size()) {
            current--;
            vWhichPicked.pop_back();
        }

////        cout << "State: ";
////        for (int const value : vWhichPicked) {
////            cout << value << " ";
////        }
////
////        cout << endl;

        if (vWhichPicked[current] == vSuperVertices[current].size()-1) { // if just evaluated last element of this set
            if (bottom == current - 1) { // move the bottom up.
                bottom = current;
                current++;
                vWhichPicked.push_back(-1); // pick the first element in the next set
            } else { // go lower
                current--;
                vWhichPicked.pop_back();
            }
        } else {
            vWhichPicked[current++]++;  // pick the next guy and increase current
            vWhichPicked.push_back(-1);
        }
    }

    cout << "Done with States" << endl;
#endif

    vector<vector<int>> reducedGraph(vReducedVertices.size());

    for (pair<int,int> mappedPair : vMapReducedToOriginal) {
        int const originalVertex(mappedPair.first);
        int const newVertex(mappedPair.second);

        for (int const originalNeighbor : m_AdjacencyArray[originalVertex]) {
            map<int,int>::const_iterator cit = vMapReducedToOriginal.find(originalNeighbor);
            if (cit != vMapReducedToOriginal.end()) {
                reducedGraph[newVertex].push_back(cit->second); // insert mapped neighbor into vertex's list
            }
        }
    }

    cout << "Begin running on reduced graph..." << endl;

    cliques.clear();

    pAlgorithm = new MaximumCliqueAlgorithm(new DegeneracyIndependentSets(reducedGraph));
    pAlgorithm->Run(cliques);

    if (!cliques.empty()) {
        for (int const vertex : cliques.front()) {
            cout << vReducedVertices[vertex] << " ";
        }
    }
    cout << endl;

    cout << "Found clique of size " << endl;

#endif

    return cliqueCount;
}

/*! \brief Recursively list all maximal cliques containing all of
           all vertices in R, some vertices in P and no vertices in X.

    \param cliqueCount A pointer to the number of maximal cliques computed 
                       thus far.

    \param cliques A linked list of cliques to return. <b>(only available when compiled 
                   with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

    \param partialClique A linked list storing R, the partial clique for this
                         recursive call. 
*/

static unsigned long recursionNode(0);

void CliqueGraphAlgorithm::RunRecursive(long &cliqueCount, list<list<int>> &cliques, list<int> &partialClique)
{
    int const currentRecursionNode(recursionNode++);

////    cout << currentRecursionNode << endl;
////    if (partialClique.empty()) {
////        cout << "Another vertex down...only " << m_pSets->SizeOfP() << " more to go" << endl;
////    }

////    m_pSets->PrintSummary(__LINE__);

    stepsSinceLastReportedClique++;

    vector<int> dominatedVertices;
    m_pSets->RemoveDominatedVertices(dominatedVertices);

    if (m_pSets->PIsEmpty() && !partialClique.empty()) {
        size_t min2Neighborhood(ULONG_MAX);
        for (int const vertex : partialClique) {
            size_t twoNeighborhood(m_AdjacencyArray[vertex].size());
            if (twoNeighborhood < min2Neighborhood)
                min2Neighborhood = twoNeighborhood;
        }

        if (min2Neighborhood > maxMin2Neighborhood)
            maxMin2Neighborhood = min2Neighborhood;
    }

    // if X is empty and P is empty, return partial clique as maximal
    if (m_pSets->XAndPAreEmpty()) {
        cliqueCount++;
////        if (cliqueCount %50 == 0)
////            cout << "Found clique #" << cliqueCount << endl;

        if (stepsSinceLastReportedClique > partialClique.size()) {
            numLargeJumps++;
            //cout << "steps: " << stepsSinceLastReportedClique << ">" << partialClique.size() << endl;
            if (largestDifference < (stepsSinceLastReportedClique - partialClique.size())) {
                largestDifference = stepsSinceLastReportedClique - partialClique.size();
            }
        }

        stepsSinceLastReportedClique = 0;

        processClique( 
                       #ifdef RETURN_CLIQUES_ONE_BY_ONE
                       cliques,
                       #endif
                       partialClique );

        int const cliqueIndex = cliqueCount-1;
        set<int> neighbors;

////        bool const debug(find(partialClique.begin(), partialClique.end(), 26) != partialClique.end());
////
////        if (debug) {
////            cout << cliqueIndex << ": {";
////        }
////
        // find all my neighbors.
        for (int const vertex : partialClique) {
////            if (debug) cout << vertex << " ";
            neighbors.insert(m_CliquesWithVertex[vertex].begin(), m_CliquesWithVertex[vertex].end());
////            if (cliqueIndex == 11) {
////                cout << "11 has neighbors: ";
////                for (int neighbor : m_CliquesWithVertex[vertex]) {
////                    cout << neighbor << " ";
////                }
////                cout << endl;
////            }

            m_CliquesWithVertex[vertex].push_back(cliqueIndex);

            //m_CliquesWithVertex[vertex].push_back(cliqueIndex);
        }

////        if (debug) cout << "}" << endl;
////
////        if (debug) {
////            cout << "I have neighbors: ";
////        }

        // put me in my neighbors' lists.
        for (int const neighbor : neighbors) {
////            if (debug) cout << neighbor << " ";
            m_NeighborSets[neighbor].insert(cliqueIndex);
        }

////        if (debug) cout << endl;

        // put my neighbors in my list
        m_NeighborSets.emplace_back(std::move(neighbors));
        vector<int> vClique(partialClique.begin(), partialClique.end());
        m_VerticesInClique.push_back(vClique);

////        if (debug) {
////            cout << "0 has neighbors " << endl;
////            for (int const neighbor : m_NeighborSets[0]) {
////                cout << neighbor << " ";
////            }
////            cout << endl;
////        }

        m_pSets->ReturnDominatedVertices(dominatedVertices);
        return;
    }

    // avoid work if P is empty.
    if (m_pSets->PIsEmpty()) {
        m_pSets->ReturnDominatedVertices(dominatedVertices);
        return;
    }

    vector<int> vNonNeighborsOfPivot = std::move(m_pSets->ChoosePivot());

    // add candiate vertices to the partial clique one at a time and 
    // search for maximal cliques
    if (!vNonNeighborsOfPivot.empty()) {
        for (int const vertex : vNonNeighborsOfPivot) {
            m_pSets->MoveFromPToR(vertex);

            // add vertex into partialClique, representing R.
            partialClique.push_back(vertex);
            list<int>::iterator vertexLink(partialClique.end());
            --vertexLink;

#ifdef PRINT_CLIQUES_TOMITA_STYLE
            printf("%d ", vertex);
#endif

            // recursively compute maximal cliques with new sets R, P and X
            RunRecursive(cliqueCount, cliques, partialClique);

#ifdef PRINT_CLIQUES_TOMITA_STYLE
            printf("b ");
#endif

            // remove vertex from partialCliques
            partialClique.erase(vertexLink);
            m_pSets->MoveFromRToX(vertex);
        }

        // swap vertices that were moved to X back into P, for higher recursive calls.
        m_pSets->ReturnVerticesToP(vNonNeighborsOfPivot);
    }

    m_pSets->ReturnDominatedVertices(dominatedVertices);

    stepsSinceLastReportedClique++;
}
