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

// local includes
#include "Tools.h"
#include "MemoryManager.h"

#include "Algorithm.h"

// system includes
#include <list>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstdlib>

/*! \file TomitaAlgorithm.h

    \brief see TomitaAlgorithm.cpp

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

class TomitaAlgorithm : public Algorithm
{
public:
    TomitaAlgorithm(char **ppAdjacencyMatrix, int const numVertices);
    virtual ~TomitaAlgorithm();

    virtual long Run(std::list<std::list<int>> &cliques);

    long listAllMaximalCliquesMatrix(char** adjacencyMatrix, int numVertices);


    int findBestPivotNonNeighborsMatrix(int** pivotNonNeighbors, int* numNonNeighbors,
                                     char** adjacencyMatrix,
                                     int* vertexSets, int* vertexLookup, int size,
                                     int beginX, int beginP, int beginR);

    void listAllMaximalCliquesMatrixRecursive(long* cliqueCount,
                                           std::list<int> &partialClique, 
                                           char** adjacencyMatrix,
                                           int* vertexSets, int* vertexLookup, int size,
                                           int beginX, int beginP, int beginR, long& stepsSinceLastReportedClique);

private:
    char **m_ppAdjacencyMatrix;
    int  m_iNumVertices;
};

#endif
