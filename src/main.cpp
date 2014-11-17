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
#include "BronKerboschAlgorithm.h"
#include "TimeDelayAdjacencyListAlgorithm.h"
#include "TimeDelayMaxDegreeAlgorithm.h"
#include "TimeDelayDegeneracyAlgorithm.h"
#include "HybridAlgorithm.h"
#include "DegeneracyAlgorithm.h"
#include "FasterDegeneracyAlgorithm.h"

#include "AdjacencyListVertexSets.h"
#include "DegeneracyVertexSets.h"
#include "CacheEfficientDegeneracyVertexSets.h"
#include "IndependentSets.h"
#include "DegeneracyIndependentSets.h"
#include "MaximumCliqueAlgorithm.h"

// system includes
#include <map>
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
    return (name == "tomita" || name == "adjlist" || name == "generic-adjlist" || name == "timedelay-adjlist" || name == "timedelay-maxdegree" || 
            name == "hybrid" || name == "degeneracy" || name == "timedelay-degeneracy" || name == "faster-degeneracy" || name == "generic-degeneracy" || name == "cache-degeneracy" || name == "mis" || name == "degeneracy-mis");
}

void ProcessCommandLineArgs(int argc, char** argv, map<string,string> &mapCommandLineArgs)
{
    for (int i = 1; i < argc; ++i) {
        string const argument(argv[i]);
        size_t const positionOfEquals(argument.find_first_of('='));
        if (positionOfEquals != string::npos) {
            string const key  (argument.substr(0,positionOfEquals));
            string const value(argument.substr(positionOfEquals+1));
////            cout << "Parsed: " << key << "=" << value << endl;
            mapCommandLineArgs[key] = value;
        } else {
            mapCommandLineArgs[argument] = "";
        }
    }
}

int main(int argc, char** argv)
{
    int failureCode(0);


    map<string,string> mapCommandLineArgs;

    ProcessCommandLineArgs(argc, argv, mapCommandLineArgs);

    bool   const outputLatex(mapCommandLineArgs.find("--latex") != mapCommandLineArgs.end());
    string const inputFile((mapCommandLineArgs.find("--input-file") != mapCommandLineArgs.end()) ? mapCommandLineArgs["--input-file"] : "");
    string const algorithm((mapCommandLineArgs.find("--algorithm") != mapCommandLineArgs.end()) ? mapCommandLineArgs["--algorithm"] : "");

    if (inputFile.empty()) {
        cout << "ERROR: Missing input file " << endl;
        // ShowUsageMessage();
        // return 1; // TODO/DS
    }

    if (algorithm.empty()) {
        cout << "ERROR: Missing algorithm" << endl;
        // ShowUsageMessage();
        // return 1; // TODO/DS
    }

    if (argc <= 1 || !isValidAlgorithm(algorithm)) {
        cout << "usage: " << argv[0] << " <tomita|adjlist|hybrid|degeneracy|*> [--latex]" << endl;
    }

    string const name(algorithm);
    Algorithm *pAlgorithm(nullptr);

    int n; // number of vertices
    int m; // 2x number of edges

    vector<list<int>> adjacencyList;
    if (inputFile.find(".graph") != string::npos) {
        adjacencyList = readInGraphAdjListEdgesPerLine(n, m, inputFile);
    } else {
        adjacencyList = readInGraphAdjList(n, m, inputFile);
    }

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

    bool const bComputeAdjacencyArray(name == "adjlist" || name == "timedelay-adjlist" || name == "generic-adjlist" ||name == "timedelay-maxdegree" || name == "timedelay-degeneracy" || name == "faster-degeneracy" || name == "generic-degeneracy" || name == "cache-degeneracy" || name == "mis" || name == "degeneracy-mis");

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
        adjacencyList.clear(); // does this free up memory? probably some...
    }

    if (name == "tomita") {
        pAlgorithm = new TomitaAlgorithm(adjacencyMatrix, n);
    } else if (name == "adjlist") {
        pAlgorithm = new AdjacencyListAlgorithm(adjacencyArray);
    } else if (name == "generic-adjlist") {
        AdjacencyListVertexSets *pSets = new AdjacencyListVertexSets(adjacencyArray);
        pAlgorithm = new BronKerboschAlgorithm(pSets);
    } else if (name == "generic-degeneracy") {
        DegeneracyVertexSets *pSets = new DegeneracyVertexSets(adjacencyArray);
        pAlgorithm = new BronKerboschAlgorithm(pSets);
    } else if (name == "timedelay-adjlist") {
        pAlgorithm = new TimeDelayAdjacencyListAlgorithm(adjacencyArray);
    } else if (name == "timedelay-maxdegree") {
        pAlgorithm = new TimeDelayMaxDegreeAlgorithm(adjacencyArray);
    } else if (name == "hybrid") {
        pAlgorithm = new HybridAlgorithm(adjacencyList);
    } else if (name == "degeneracy") {
        pAlgorithm = new DegeneracyAlgorithm(adjacencyList);
    } else if (name == "faster-degeneracy") {
        pAlgorithm = new FasterDegeneracyAlgorithm(adjacencyArray);
    } else if (name == "cache-degeneracy") {
        CacheEfficientDegeneracyVertexSets *pSets = new CacheEfficientDegeneracyVertexSets(adjacencyArray);
        pAlgorithm = new BronKerboschAlgorithm(pSets);
    } else if (name == "mis") {
        IndependentSets *pSets = new IndependentSets(adjacencyArray);
        pAlgorithm = new BronKerboschAlgorithm(pSets);
    } else if (name == "degeneracy-mis") {
        DegeneracyIndependentSets *pSets = new DegeneracyIndependentSets(adjacencyArray);
        pAlgorithm = new BronKerboschAlgorithm(pSets);
    } else if (name == "timedelay-degeneracy") {
        pAlgorithm = new TimeDelayDegeneracyAlgorithm(adjacencyList);
    } else {
        cout << "ERROR: unrecognized algorithm name " << name << endl;
        return 1;
    }

    // Run algorithm
    list<list<int>> cliques;
    RunAndPrintStats(pAlgorithm, cliques, outputLatex);

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
