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
#include "TomitaAlgorithm.h"
#include "AdjacencyListAlgorithm.h"
#include "HybridAlgorithm.h"
#include "DegeneracyAlgorithm.h"

// system includes
#include <list>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

/*! \file main.cpp

    \brief Main entry point for quick cliques software. This is where we parse the command line options, read the inputs, and decide which clique method to run.

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

bool isValidAlgorithm(string const &name)
{
    return (name == "tomita" || name == "adjlist" ||
            name == "hybrid" || name == "degeneracy");
}

int main(int argc, char** argv)
{
    int failureCode(0);

    if (argc==1 || !isValidAlgorithm(argv[1])) {
        cout << "usage: " << argv[0] << " <tomita|adjlist|hybrid|degeneracy>" << endl;
    }

    string const name(argv[1]);
    MaximalCliquesAlgorithm *pAlgorithm(nullptr);

    int n; // number of vertices
    int m; // 2x number of edges

    vector<list<int>> const adjacencyList = readInGraphAdjList(&n,&m);

    bool const bComputeAdjacencyMatrix(name == "tomita");

    char** adjacencyMatrix(nullptr);

    if (bComputeAdjacencyMatrix) {
        adjacencyMatrix = (char**)Calloc(n, sizeof(char*));

        for(int i=0; i<n; i++) {
            adjacencyMatrix[i] = (char*)Calloc(n, sizeof(char));
            for(int const neighbor : adjacencyList[i]) {
                adjacencyMatrix[i][neighbor] = 1; 
            }
        }
    }

    bool const bComputeAdjacencyArray(name == "adjlist");

    vector<vector<int>> adjacencyArray;

    if (bComputeAdjacencyArray) {
        adjacencyArray.resize(n);
        for (int i=0; i<n; i++) {
            adjacencyArray[i].resize(adjacencyList[i].size());
            int j = 0;
            for (int const neighbor : adjacencyList[i]) {
                adjacencyArray[i][j++] = neighbor;
            }
        }
    }

    if (name == "tomita") {
        pAlgorithm = new TomitaAlgorithm(adjacencyMatrix, n);
    } else if (name == "adjlist"){
        pAlgorithm = new AdjacencyListAlgorithm(adjacencyArray);
    } else if (name == "hybrid"){
        pAlgorithm = new HybridAlgorithm(adjacencyList);
    } else if (name == "degeneracy"){
        pAlgorithm = new DegeneracyAlgorithm(adjacencyList);
    } else {
        cout << "ERROR: unrecognized algorithm name " << name << endl;
    }

    // Run algorithm
    list<list<int>> cliques;
    RunAndPrintStats(pAlgorithm, cliques);

    cliques.clear();

    if (adjacencyMatrix != nullptr) {
        int i = 0;
        while(i<n) {
            Free(adjacencyMatrix[i]);
            i++;
        }
        Free(adjacencyMatrix); 
        adjacencyMatrix = nullptr;
    }

    delete pAlgorithm; pAlgorithm = nullptr;

    ////CommandLineOptions options = ParseCommandLineOptions(argc, argv);

    ////if (options.verify) {
    ////}

    return 0;
}
