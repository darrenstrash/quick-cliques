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

#include "Tools.h"
#include "DegeneracyTools.h"
#include <list>
#include <vector>
#include "MemoryManager.h"

/*! \file compdegen.cpp

   \brief print the degeneracy of the input graph

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

using namespace std;

int main()
{

    int n; // number of vertices
    int m; // 2x number of edges

    vector<list<int>> adjacencyList = readInGraphAdjList(&n,&m);

    int maximumDegree(0);
    for (list<int> const &neighbors : adjacencyList) {
        if (maximumDegree < neighbors.size())
            maximumDegree = neighbors.size();
    }

    int d = computeDegeneracy(adjacencyList, n);

    adjacencyList.clear(); 

    fprintf(stderr, "Degeneracy is %d\n", d);
    fprintf(stderr, "MaxDegree  is %d\n", maximumDegree);
    return 0;
}
