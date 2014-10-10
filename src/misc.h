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

#include"LinkedList.h"

/*! \file misc.h

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

#define max(x,y) (x > y? x:y)
#define min(x,y) (x < y? x:y)

int nodeComparator(void* node1, void* node2);

void printArray(int* array, int size);

void printArrayOfLinkedLists(LinkedList** listOfLists, int size);

void printInt(void* integer);

void destroyCliqueResults(LinkedList* cliques);

LinkedList** readInGraphAdjList(int* n, int* m);

void runAndPrintStatsMatrix(long (*function)(char**,
                                             #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                             LinkedList*,
                                             #endif
                                             int),
                            const char* algName,
                            char** adjMatrix,
                            #ifdef RETURN_CLIQUES_ONE_BY_ONE
                            LinkedList* cliques,
                            #endif
                            int n );

void runAndPrintStatsListList( long (*function)(LinkedList**, 
                                                int**, 
                                                #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                                LinkedList*,
                                                #endif
                                                int*, int),
                               const char* algName,
                               LinkedList** adjListLinked,
                               int** adjListArray,
                               #ifdef RETURN_CLIQUES_ONE_BY_ONE
                               LinkedList* cliques,
                               #endif
                               int* degree,
                               int n );


inline void processClique(
                          #ifdef RETURN_CLIQUES_ONE_BY_ONE
                          LinkedList *cliques,
                          #endif
                          LinkedList *clique);
#endif

