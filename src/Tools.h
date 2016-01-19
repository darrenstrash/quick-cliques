#ifndef _DJS_MISC_H_
#define _DJS_MISC_H_

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

#include <list>
#include <vector>
#include <string>
#include <stdio.h>

class Algorithm;

/*! \file Tools.h

    \brief see misc.c

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

#include <cmath>
//#define max(x,y) (x > y? x:y)
//#define min(x,y) (x < y? x:y)

int nodeComparator(int node1, int node2);

void printArray(int* array, int size);

void printArrayWithIndexArrows(int* array, int size, int index1, int index2, int index3);

void printArrayOfLinkedLists(std::vector<std::list<int>> const &listOfLists, int size);

void destroyCliqueResults(std::list<std::list<int>> &cliques);

std::vector<std::list<int>> readInGraphAdjList(int* n, int* m);

std::vector<std::list<int>> readInGraphAdjList(int &n, int &m, std::string const &fileName);
std::vector<std::list<int>> readInGraphAdjListEdgesPerLine(int &n, int &m, std::string const &fileName);

void runAndPrintStatsMatrix(long (*function)(char**,
                                             int),
                            const char* algName,
                            char** adjMatrix,
                            int n );

void RunAndPrintStats(Algorithm* pAlgorithm, std::list<std::list<int>> &cliques, bool const outputLatex);

void printListAbbv(std::list<int> const &linkedList, void (*printFunc)(int));

/*! \brief process a clique, which may include printing it in
           one of several formats and/or adding the 
           clique to a linked list.

    \param clique the clique to add to the linked list

*/

inline void processClique(std::list<int> const &clique)
{
    #ifdef PRINT_CLIQUES_TOMITA_STYLE
    printf("c ");
    #endif
}

void DescribeVertex(int const lineNumber, int *vertexSets, int *vertexLookup, int const size, int const vertex, int const beginX, int const beginD, int const beginP, int const beginR);

void DescribeSet(std::string const &setName, int const begin, int const end);

void DescribeState(int const lineNumber, int *vertexSets, int *vertexLookup, int const size, int const beginX, int const beginD, int const beginP, int const beginR);

void CheckConsistency(int const lineNumber, size_t const recursionNumber, int *vertexSets, int *vertexLookup, int const size);

void CheckReverseConsistency(int const lineNumber, size_t const recursionNumber, int *vertexSets, int *vertexLookup, int const size);

bool IsMaximalClique(std::list<int> const &clique, std::vector<std::vector<int>> const &adjacencyList);

namespace Tools
{
    void printList(std::list<int> const &linkedList, void (*printFunc)(int));
    void printInt(int integer);
    std::vector<int> ReadMetisOrdering(std::string const &filename);
    std::string GetTimeInSeconds(clock_t delta, bool brackets=true);
};

#endif

