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

#include <limits.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <iostream>

#include "Tools.h"
#include <list>
#include <vector>
#include "MemoryManager.h"
#include "MaximumCliqueAlgorithm.h"
#include "Algorithm.h"

using namespace std;

/*! \file MaximumCliqueAlgorithm.cpp

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

MaximumCliqueAlgorithm::MaximumCliqueAlgorithm(VertexSets *pSets)
 : Algorithm(pSets->GetName())
 , m_pSets(pSets)
 , m_uMaximumCliqueSize(0)
{
}

MaximumCliqueAlgorithm::~MaximumCliqueAlgorithm()
{
    delete m_pSets; m_pSets = nullptr;
}

/*! \brief List all maximal cliques in a given graph using the algorithm
           by Tomita et al. (TCS 2006), modified to use an adjacency list
           representation of the graph instead of an adjacency matrix. 
 
    \param cliques A linked list of cliques to return. <b>(only available when compiled 
                   with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

    \return The number of maximal cliques of the input graph.
*/

long MaximumCliqueAlgorithm::Run(list<list<int>> &cliques)
{
    long cliqueCount = 0;
    list<int> partialClique;

    m_pSets->Initialize();

    while (m_pSets->GetNextTopLevelPartition()) {
        m_pSets->GetTopLevelPartialClique(partialClique);
        RunRecursive(cliqueCount, cliques, partialClique);
        partialClique.clear();
    }

    if (!GetQuiet()) {
        cerr << "Largest Difference : " << largestDifference << endl;
        cerr << "Num     Differences: " << numLargeJumps << endl;
        cerr << "Largest Clique     : " << m_uMaximumCliqueSize << endl;
    }

    return  cliqueCount;
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

void MaximumCliqueAlgorithm::RunRecursive(long &cliqueCount, list<list<int>> &cliques, list<int> &partialClique)
{
    int const currentRecursionNode(recursionNode++);
////
////    cout << currentRecursionNode << endl;
////    if (partialClique.empty()) {
////        cout << "Another vertex down...only " << m_pSets->SizeOfP() << " more to go" << endl;
////    }

////    cout << "Partial Clique: ";
////
////    for (int const vertex : partialClique) {
////        cout << vertex << " ";
////    }
////
////    cout << endl << flush;

////    vector<int> dominatedVertices;
////    m_pSets->RemoveDominatedVertices(dominatedVertices);

    if (partialClique.size() + m_pSets->RemainingSizeEstimate() <= m_uMaximumCliqueSize) {
////        m_pSets->ReturnDominatedVertices(dominatedVertices);
        return;
    }

////    if (partialClique.size() > m_uMaximumCliqueSize) {
////        if (!GetQuiet())
////            cout << "Found clique with size " << partialClique.size() << endl;
////        m_uMaximumCliqueSize = partialClique.size();
////        ExecuteCallBacks(partialClique);
////    }

////    m_Sets.PrintSummary(__LINE__);

    stepsSinceLastReportedClique++;

////    ExecuteCallBacks(partialClique);

    // if X is empty and P is empty, return partial clique as maximal
    if (m_pSets->PIsEmpty()) { //XAndPAreEmpty()) {
        cliqueCount++;
        if (!GetQuiet())
            cout << "Found clique with size " << partialClique.size() << endl;

        if (stepsSinceLastReportedClique > partialClique.size()) {
            numLargeJumps++;
            //cout << "steps: " << stepsSinceLastReportedClique << ">" << partialClique.size() << endl;
            if (largestDifference < (stepsSinceLastReportedClique - partialClique.size())) {
                largestDifference = stepsSinceLastReportedClique - partialClique.size();
            }
        }

        m_uMaximumCliqueSize = partialClique.size();

        stepsSinceLastReportedClique = 0;

        cliques.clear();
        cliques.push_back(partialClique);

        ExecuteCallBacks(partialClique);

        processClique( 
                       #ifdef RETURN_CLIQUES_ONE_BY_ONE
                       cliques,
                       #endif
                       partialClique );

////        m_pSets->ReturnDominatedVertices(dominatedVertices);
        return;
    }

    // avoid work if P is empty.
    if (m_pSets->PIsEmpty()) {
////        m_pSets->ReturnDominatedVertices(dominatedVertices);
        return;
    }

    vector<int> vVerticesToEvaluate = std::move(m_pSets->ChoosePivotNonConst());
    vector<int> vEvaluatedVertices;

    // add candiate vertices to the partial clique one at a time and 
    // search for maximal cliques
    int vertex = m_pSets->GetNextVertexToEvaluate(vVerticesToEvaluate);
    while (vertex != -1) {
        vEvaluatedVertices.push_back(vertex);
        if (currentRecursionNode == 0) {
            cout << "Only " << vVerticesToEvaluate.size() << " more vertices to go!" << endl << flush;
        }
////        if (currentRecursionNode == 15) {
////            cout << currentRecursionNode << ":" << endl << flush;
////            cout << "vertex " << vertex << " is " << (m_pSets->InP(vertex) ? " " : "not ") << " in P" << endl << flush;
////            cout << "Adding vertex " << vertex << " to clique " << endl << flush;
////        }

        m_pSets->MoveFromPToR(partialClique, vertex);

#ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("%d ", vertex);
#endif

////        m_pSets->ApplyReduction(partialClique);

        // recursively compute maximal cliques with new sets R, P and X
        RunRecursive(cliqueCount, cliques, partialClique);

#ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
#endif

////        m_pSets->UndoReduction(partialClique);

        m_pSets->MoveFromRToX(partialClique, vertex);

        if (partialClique.size() + m_pSets->RemainingSizeEstimate() <= m_uMaximumCliqueSize) {
            break;
        }
        vertex = m_pSets->GetNextVertexToEvaluate(vVerticesToEvaluate);
    }

    // swap vertices that were moved to X back into P, for higher recursive calls.
    m_pSets->ReturnVerticesToP(partialClique, vEvaluatedVertices);

    stepsSinceLastReportedClique++;

////    m_pSets->ReturnDominatedVertices(dominatedVertices);
    return;
}
