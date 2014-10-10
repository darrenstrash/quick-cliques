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

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "Tools.h"
#include <list>
#include <vector>
#include "MemoryManager.h"
#include "HybridAlgorithm.h"

using namespace std;

/*! \file hybrid.cpp

    \brief Execute the algorithm in HybridAlgorithm.cpp
           and print the number of cliques found and wall clock
           execution time.

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

int main()
{

#ifdef MEMORY_DEBUG
    fprintf(stderr, "WARNING: MEMORY_DEBUG is defined, timing will be off.\n");
#endif

    int n; // number of vertices
    int m; // 2x number of edges

    //TODO: use std::move, to eliminate copy
    vector<list<int>> const adjacencyList = readInGraphAdjList(&n,&m);

    int** adjList = (int**)Calloc(n, sizeof(int*));
    int* degree = (int*)Calloc(n, sizeof(int));

    for(int i = 0; i < n; i++)
    {
        degree[i] = adjacencyList[i].size();
        adjList[i] = (int*)Calloc(degree[i], sizeof(int));
        int j = 0;
        for(int const neighbor : adjacencyList[i])
        {
            adjList[i][j++] = neighbor;
        }
    }

    #ifdef RETURN_CLIQUES_ONE_BY_ONE
    list<list<int>> cliques;
    #endif

    runAndPrintStatsListList( &listAllMaximalCliquesHybrid,
                              "hybrid",
                              adjacencyList, adjList, 
                              #ifdef RETURN_CLIQUES_ONE_BY_ONE
                              cliques,
                              #endif
                              degree, n );

    // Free up memory from adjacency list.

    #ifdef RETURN_CLIQUES_ONE_BY_ONE
    destroyCliqueResults(cliques);
    #endif

    int i = 0;
    while(i < n)
    {
        Free(adjList[i]);
        i++;
    }

    Free(degree);
    Free(adjList); 

    return 0;
}
