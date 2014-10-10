#ifndef _DJS_TOMITA_ALGORITHM_H
#define _DJS_TOMITA_ALGORITHM_H

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

#include<assert.h>
#include<stdio.h>
#include<stdlib.h>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"

/*! \file tomita_algorithm.h

    \brief see tomita_algorithm.c

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

long listAllMaximalCliquesMatrix( char** adjacencyMatrix,
                                 #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                 LinkedList* cliques,
                                 #endif
                                 int    numVertices );


int findBestPivotNonNeighborsMatrix( int** pivotNonNeighbors, int* numNonNeighbors,
                                     char** adjacencyMatrix,
                                     int* vertexSets, int* vertexLookup, int size,
                                     int beginX, int beginP, int beginR);

void listAllMaximalCliquesMatrixRecursive( long* cliqueCount,
                                           #ifdef RETURN_CLIQUES_ONE_BY_ONE
                                           LinkedList* cliques,
                                           #endif
                                           LinkedList* partialClique, 
                                           char** adjacencyMatrix,
                                           int* vertexSets, int* vertexLookup, int size,
                                           int beginX, int beginP, int beginR );


#endif
