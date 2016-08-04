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
#include "CliqueTools.h"

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

    \copyright Copyright (c) 2011-2016 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

bool isValidAlgorithm(string const &name)
{
    return (name == "tomita" || name == "adjlist" || name == "hybrid" || name == "degeneracy");
}

void ProcessCommandLineArgs(int const argc, char** argv, map<string,string> &mapCommandLineArgs)
{
    for (int i = 1; i < argc; ++i) {
////        cout << "Processing argument " << i << endl;
        string const argument(argv[i]);
////        cout << "    Argument is " << argument << endl;
        size_t const positionOfEquals(argument.find_first_of("="));
////        cout << "    Position of = " << positionOfEquals << endl;
        if (positionOfEquals != string::npos) {
            string const key  (argument.substr(0,positionOfEquals));
            string const value(argument.substr(positionOfEquals+1));
////            cout << "    Parsed1: " << key << "=" << value << endl;
            mapCommandLineArgs[key] = value;
        } else {
////            cout << "    Parsed2: " << argument << endl;
            mapCommandLineArgs[argument] = "";
        }
    }
}

void PrintDebugWarning()
{
    cout << "\n\n\n\n\n" << flush;
    cout << "#########################################################################" << endl << flush;
    cout << "#                                                                       #" << endl << flush;
    cout << "#    WARNING: Debugging is turned on. Don't believe the run times...    #" << endl << flush;
    cout << "#                                                                       #" << endl << flush;
    cout << "#########################################################################" << endl << flush;
    cout << "\n\n\n\n\n" << flush;
}

void PrintHeader()
{
    cout << "NOTE: Quick Cliques v2.0beta." << endl;
}

string basename(string const &fileName)
{
    string sBaseName(fileName);

    size_t const lastSlash(sBaseName.find_last_of("/\\"));
    if (lastSlash != string::npos) {
        sBaseName = sBaseName.substr(lastSlash+1);
    }

    size_t const lastDot(sBaseName.find_last_of("."));
    if (lastDot != string::npos) {
        sBaseName = sBaseName.substr(0, lastDot);
    }

    return sBaseName;
}

int main(int argc, char** argv)
{
    int failureCode(0);

    map<string,string> mapCommandLineArgs;

    ProcessCommandLineArgs(argc, argv, mapCommandLineArgs);

    bool   const bQuiet(mapCommandLineArgs.find("--verbose") == mapCommandLineArgs.end());
    bool   const bOutputLatex(mapCommandLineArgs.find("--latex") != mapCommandLineArgs.end());
    bool   const bOutputTable(mapCommandLineArgs.find("--table") != mapCommandLineArgs.end());
    string const inputFile((mapCommandLineArgs.find("--input-file") != mapCommandLineArgs.end()) ? mapCommandLineArgs["--input-file"] : "");
    string const algorithm((mapCommandLineArgs.find("--algorithm") != mapCommandLineArgs.end()) ? mapCommandLineArgs["--algorithm"] : "");
    bool   const staging(mapCommandLineArgs.find("--staging") != mapCommandLineArgs.end());

    bool   const bTableMode(bOutputLatex || bOutputTable);

    if (!bTableMode) {
        PrintHeader();
#ifdef DEBUG_MESSAGE
        PrintDebugWarning();
#endif //DEBUG_MESSAGE
    }

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

    if (argc <= 1 || !isValidAlgorithm(algorithm) || inputFile.empty()) {
        cout << "USAGE: " << argv[0] << " --input-file=<filename> --algorithm=<tomita|adjlist|degeneracy|hybrid>" << endl;
        return 1;
    }

    string const name(algorithm);
    Algorithm *pAlgorithm(nullptr);

    int n; // number of vertices
    int m; // 2x number of edges

    vector<list<int>> adjacencyList;
    if (inputFile.find(".graph") != string::npos) {
        if (!bTableMode) cout << "Reading .graph file format. " << endl << flush;
        adjacencyList = readInGraphAdjListEdgesPerLine(n, m, inputFile);
    } else {
        if (!bTableMode) cout << "Reading .edges file format. " << endl << flush;
        adjacencyList = readInGraphAdjList(n, m, inputFile);
    }

    bool const bComputeAdjacencyMatrix(adjacencyList.size() < 20000);
    bool const bShouldComputeAdjacencyMatrix(name == "tomita");

    if (bShouldComputeAdjacencyMatrix && !bComputeAdjacencyMatrix) {
        cout << "ERROR!: unable to compute adjacencyMatrix, since the graph is too large: " << adjacencyList.size() << " vertices." << endl << flush;
        exit(1);
    }

    char** adjacencyMatrix(nullptr);

    vector<vector<char>> vAdjacencyMatrix;

    if (bComputeAdjacencyMatrix) {
        adjacencyMatrix = (char**)Calloc(n, sizeof(char*));
        vAdjacencyMatrix.resize(n);

        for(int i=0; i<n; i++) {
            adjacencyMatrix[i] = (char*)Calloc(n, sizeof(char));
            vAdjacencyMatrix[i].resize(n);
            for(int const neighbor : adjacencyList[i]) {
                adjacencyMatrix[i][neighbor]  = 1; 
                vAdjacencyMatrix[i][neighbor] = 1; 
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
        adjacencyList.clear(); // does this free up memory? probably some...
    }

    if (name == "tomita") {
        pAlgorithm = new TomitaAlgorithm(adjacencyMatrix, n);
    } else if (name == "adjlist") {
        pAlgorithm = new AdjacencyListAlgorithm(adjacencyArray);
    } else if (name == "degeneracy") {
        pAlgorithm = new DegeneracyAlgorithm(adjacencyList);
    } else if (name == "hybrid") {
        pAlgorithm = new HybridAlgorithm(adjacencyList);
    } else {
        cout << "ERROR: unrecognized algorithm name " << name << endl;
        return 1;
    }

#ifdef PRINT_CLIQUES_ONE_BY_ONE
    auto printClique = [](list<int> const &clique) {
        Tools::printList(clique, &Tools::printInt);
    };

    pAlgorithm->AddCallBack(printClique);
#endif //PRINT_CLIQUES_ONE_BY_ONE

    auto verifyMaximalCliqueArray = [&adjacencyArray](list<int> const &clique) {
        bool const isIS = CliqueTools::IsMaximalClique(adjacencyArray, clique, true /* verbose */);
        if (!isIS) {
            cout << "ERROR: Set " << (isIS ? "is" : "is not" ) << " a maximal clique!" << endl;
        }
    };

    auto verifyCliqueMatrix = [&vAdjacencyMatrix](list<int> const &clique) {
        bool const isIS = CliqueTools::IsClique(vAdjacencyMatrix, clique, true /* verbose */);
        if (!isIS) {
            cout << "ERROR: Set " << (isIS ? "is" : "is not" ) << " a clique!" << endl;
        }
    };

    auto printCliqueSize = [](list<int> const &clique) {
        cout << "Found clique of size " << clique.size() << endl << flush;
    };

////    pAlgorithm->AddCallBack(printCliqueSize);
////    pAlgorithm->AddCallBack(printClique);

    if (!bComputeAdjacencyMatrix) {
////        pAlgorithm->AddCallBack(verifyIndependentSetArray);
////        pAlgorithm->AddCallBack(verifyMaximalCliqueArray);
    } else {
////        pAlgorithm->AddCallBack(verifyCliqueMatrix);
////        pAlgorithm->AddCallBack(verifyIndependentSetMatrix);
    }

    // Run algorithm
    list<list<int>> cliques;

#ifdef RETURN_CLIQUES_ONE_BY_ONE
    auto storeCliqueInList = [&cliques](list<int> const &clique) {
        cliques.push_back(clique);
    };
    pAlgorithm->AddCallBack(storeCliqueInList);
#endif //RETURN_CLIQUES_ONE_BY_ONE

    pAlgorithm->SetQuiet(bQuiet);

    RunAndPrintStats(pAlgorithm, cliques, bTableMode);

////    cout << "Last clique has size: " << cliques.back().size() << endl << flush;

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

#ifdef DEBUG_MESSAGE
    PrintDebugWarning();
#endif

    ////CommandLineOptions options = ParseCommandLineOptions(argc, argv);

    ////if (options.verify) {
    ////}

    return 0;
}
