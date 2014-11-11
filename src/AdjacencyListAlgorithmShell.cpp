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
#include "AdjacencyListAlgorithmShell.h"
#include "MaximalCliquesAlgorithm.h"

using namespace std;

/*! \file AdjacencyListAlgorithmShell.cpp

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

AdjacencyListAlgorithmShell::AdjacencyListAlgorithmShell(AdjacencyListVertexSets &sets)
 : MaximalCliquesAlgorithm("shell-adjlist")
 , m_Sets(sets)
{
}

AdjacencyListAlgorithmShell::~AdjacencyListAlgorithmShell()
{
}

/*! \brief List all maximal cliques in a given graph using the algorithm
           by Tomita et al. (TCS 2006), modified to use an adjacency list
           representation of the graph instead of an adjacency matrix. 
 
    \param cliques A linked list of cliques to return. <b>(only available when compiled 
                   with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

    \return The number of maximal cliques of the input graph.
*/

long AdjacencyListAlgorithmShell::Run(list<list<int>> &cliques)
{
    long cliqueCount = 0;
    list<int> partialClique;

    RunRecursive(cliqueCount, cliques, partialClique);

    cerr << "Largest Difference : " << largestDifference << endl;
    cerr << "Num     Differences: " << numLargeJumps << endl;

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

void AdjacencyListAlgorithmShell::RunRecursive(long &cliqueCount, list<list<int>> &cliques, list<int> &partialClique)
{
////    int const currentRecursionNode(recursionNode++);
////
////    cout << currentRecursionNode << endl;

////    m_Sets.PrintSummary(__LINE__);

    stepsSinceLastReportedClique++;

    // if X is empty and P is empty, return partial clique as maximal
    if (m_Sets.XAndPAreEmpty()) {
        cliqueCount++;

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

        return;
    }

    // avoid work if P is empty.
    if (m_Sets.PIsEmpty())
        return;

    vector<int> vNonNeighborsOfPivot = std::move(m_Sets.ChoosePivot());

    // add candiate vertices to the partial clique one at a time and 
    // search for maximal cliques
    if (!vNonNeighborsOfPivot.empty()) {
        for (int const vertex : vNonNeighborsOfPivot) {
            m_Sets.MoveFromPToR(vertex);

            // add vertex into partialClique, representing R.
            partialClique.push_back(vertex);
            list<int>::iterator vertexLink(partialClique.end());
            --vertexLink;

#ifdef PRINT_CLIQUES_TOMITA_STYLE
            printf("%d ", vertex);
#endif

            // recursively compute maximal cliques with new sets R, P and X
            RunRecursive( cliqueCount, cliques, partialClique);

#ifdef PRINT_CLIQUES_TOMITA_STYLE
            printf("b ");
#endif

            // remove vertex from partialCliques
            partialClique.erase(vertexLink);
            m_Sets.MoveFromRToX(vertex);
        }

        // swap vertices that were moved to X back into P, for higher recursive calls.
        m_Sets.ReturnVerticesToP(vNonNeighborsOfPivot);
    }

    stepsSinceLastReportedClique++;
}
