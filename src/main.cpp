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
#include "ArraySet.h"
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
#include "GraphTools.h"

#include "TesterMISS.h"
////#include "ComparisonFullMISS.h"
#include "ConnectedComponentMISS2.h"
#include "ConnectedComponentMISS.h"
#include "LightWeightReductionDominationMISR.h"
#include "LightWeightReductionDominationMISQ.h"
#include "LightWeightReductionSparseFullMISS.h"
#include "LightWeightReductionSparseStaticOrderMISS.h"
#include "LightWeightReductionSparseMISR.h"
#include "LightWeightReductionSparseMISQ.h"
#include "LightWeightReductionStaticOrderMISS.h"
#include "LightWeightReductionFullMISS.h"
#include "LightWeightReductionMISR.h"
#include "LightWeightReductionMISQ.h"
#include "LightWeightFullMISS.h"
#include "LightWeightStaticOrderMISS.h"
#include "LightWeightMISR.h"
#include "LightWeightMISQ.h"
#include "LightWeightFullMCS.h"
#include "LightWeightStaticOrderMCS.h"
#include "LightWeightMCR.h"
#include "LightWeightSparseMCQ.h"
#include "LightWeightMCQ.h"
#include "AdjacencyMatrixVertexSetsMax.h"
#include "AdjacencyMatrixVertexSets.h"
#include "AdjacencyListVertexSetsMax.h"
#include "ReverseDegeneracyVertexSets.h"
#include "AdjacencyListVertexSets.h"
#include "DegeneracyVertexSets.h"
#include "CacheEfficientDegeneracyVertexSets.h"
#include "IndependentSets.h"
#include "IndependentSetsReduction.h"
#include "ExperimentalReduction.h"
#include "DegeneracyIndependentSets.h"
#include "DegeneracyIndependentSets2.h"
#include "MaximumCliqueAlgorithm.h"
#include "MinimumCliqueAlgorithm.h"
#include "PartialMatchDegeneracyVertexSets.h"
#include "CliqueGraphAlgorithm.h"
#include "CliqueTools.h"
#include "Experiments.h"

#include "Staging.h"

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
    return (name == "tomita" || name == "generic-adjmatrix" || name == "generic-adjmatrix-max" || name == "mcq" || name == "mcr" || name == "static-order-mcs" || name == "mcs" || name == "misq" || name == "reduction-misq" || name == "reduction-misr" || name == "reduction-static-order-miss" || name == "reduction-miss" || name == "reduction-domination-misq" || name == "reduction-domination-misr" || name == "reduction-sparse-misq" || name == "reduction-sparse-misr" || name == "reduction-sparse-static-order-miss" || name == "reduction-sparse-miss" || name == "misr" || name == "static-order-miss" || name == "miss" || name == "adjlist" || name == "generic-adjlist" || name == "generic-adjlist-max" || name == "timedelay-adjlist" || name == "timedelay-maxdegree" || name == "hybrid" || name == "degeneracy" || name == "timedelay-degeneracy" || name == "faster-degeneracy" || name == "generic-degeneracy" || name == "cache-degeneracy" || name == "mis" || name == "degeneracy-mis" || name == "partial-match-degeneracy" || name == "reverse-degeneracy" || name == "degeneracy-min" || name == "degeneracy-mis-2" || name == "reduction-mis" || name == "experimental-mis" || name == "sparse-mcq" || name == "connected-component-miss" || name == "connected-component-miss2"/* || name == "comparison-miss"*/ || name == "tester-miss");
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

void PrintExperimentalWarning()
{
    cout << "WARNING: Quick Cliques v2.0beta. (Experimental)" << endl;
    cout << "WARNING: " << endl;
    cout << "WARNING: Proceed with caution: this software is currently in an experimental state." << endl;
    cout << "WARNING: This software may be slow, the algorithm may be unstable and the results may be incorrect." << endl;
    cout << "WARNING: If you care about this sort of thing, please consider using version 1.0 of this software, which has been more thoroughly tested." << endl;
}

void RunUnitTests()
{
    std::cout << "Running unit tests..." << std::endl;
    ArraySet::Test();
    GraphTools::TestMatchingCount();
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

void RunExperiment(string const &sInputFile, string const &sExperimentName, bool const bOutputLatex, bool const bPrintHeader, vector<vector<int>> const &adjacencyArray, vector<vector<char>> const &vAdjacencyMatrix, double const dTimeout)
{
    string const dataSetName(basename(sInputFile));
    Experiments experiments(dataSetName, dTimeout, bOutputLatex, bPrintHeader, adjacencyArray, vAdjacencyMatrix);

    if (sExperimentName=="kernel-size") {
        experiments.RunKernelSize();
    } else if (sExperimentName=="kernel-sparse-miss") {
        experiments.KernelizeAndRunReductionSparseMISS();
    } else if (sExperimentName=="kernel-component-sparse-miss") {
        experiments.KernelizeAndRunComponentWiseReductionSparseMISS();
    } else if (sExperimentName=="kernel-component-no-reduction-miss") {
        experiments.KernelizeAndRunComponentWiseMISS();
    } else if (sExperimentName=="exact-search") {
        experiments.RunExactSearch();
    } else if (sExperimentName=="standard-search") {
        experiments.RunStandardSearch();
    } else if (sExperimentName=="components-miss") {
        experiments.RunComponentsMISS();
    } else if (sExperimentName=="components-standard-search") {
        experiments.RunComponentsStandardSearch();
    } else if (sExperimentName=="components-forward-search") {
        experiments.RunComponentsForwardSearch();
    } else if (sExperimentName=="forward-search") {
        experiments.RunForwardSearch();
    } else {
        cout << "ERROR!: Invalid experiment name: " << sExperimentName << endl << flush;
    }
}

int main(int argc, char** argv)
{
    int failureCode(0);

    map<string,string> mapCommandLineArgs;

    ProcessCommandLineArgs(argc, argv, mapCommandLineArgs);

    bool   const bQuiet(mapCommandLineArgs.find("--verbose") == mapCommandLineArgs.end());
    bool   const bOutputLatex(mapCommandLineArgs.find("--latex") != mapCommandLineArgs.end());
    string const inputFile((mapCommandLineArgs.find("--input-file") != mapCommandLineArgs.end()) ? mapCommandLineArgs["--input-file"] : "");
    string const algorithm((mapCommandLineArgs.find("--algorithm") != mapCommandLineArgs.end()) ? mapCommandLineArgs["--algorithm"] : "");
    bool   const computeCliqueGraph(mapCommandLineArgs.find("--compute-clique-graph") != mapCommandLineArgs.end());
    bool   const staging(mapCommandLineArgs.find("--staging") != mapCommandLineArgs.end());
    bool   const findMaximumOnly(mapCommandLineArgs.find("--maximum") != mapCommandLineArgs.end());
    string const sExperimentName((mapCommandLineArgs.find("--experiment") != mapCommandLineArgs.end()) ? mapCommandLineArgs["--experiment"] : "");
    bool   const bPrintHeader(mapCommandLineArgs.find("--header") != mapCommandLineArgs.end());
    string const sTimeout(mapCommandLineArgs.find("--timeout") != mapCommandLineArgs.end() ? mapCommandLineArgs["--timeout"] : "");
    double dTimeout(0.0);
    if (!sTimeout.empty()) {
        try {
            dTimeout = stod(sTimeout);
        } catch(...) {
            cout << "ERROR!: Invalid --timeout argument, please enter valid double value." << endl << flush;
        }
    }
    bool   const bRunExperiment(!sExperimentName.empty());

    if (!bOutputLatex) {
        PrintExperimentalWarning();
#ifdef DEBUG_MESSAGE
        PrintDebugWarning();
#endif //DEBUG_MESSAGE
        RunUnitTests();
    }


    if (inputFile.empty()) {
        cout << "ERROR: Missing input file " << endl;
        // ShowUsageMessage();
        // return 1; // TODO/DS
    }

    if (!(staging || computeCliqueGraph || bRunExperiment) && algorithm.empty()) {
        cout << "ERROR: Missing algorithm" << endl;
        // ShowUsageMessage();
        // return 1; // TODO/DS
    }

    if (argc <= 1 || (!(computeCliqueGraph || staging || bRunExperiment) && !isValidAlgorithm(algorithm))) {
        cout << "usage: " << argv[0] << " --input-file=<filename> --algorithm=<algorithm-name> [--latex] [--experiment=<experiment-name>] [--header]" << endl;
    }

    string const name(algorithm);
    Algorithm *pAlgorithm(nullptr);

    int n; // number of vertices
    int m; // 2x number of edges

    vector<list<int>> adjacencyList;
    if (inputFile.find(".graph") != string::npos) {
        if (!bOutputLatex) cout << "Reading .graph file format. " << endl << flush;
        adjacencyList = readInGraphAdjListEdgesPerLine(n, m, inputFile);
    } else {
        if (!bOutputLatex) cout << "Reading .edges file format. " << endl << flush;
        adjacencyList = readInGraphAdjList(n, m, inputFile);
    }

    bool const bComputeAdjacencyMatrix(adjacencyList.size() < 20000);
    bool const bShouldComputeAdjacencyMatrix(name == "tomita" || name == "generic-adjmatrix" || name == "generic-adjmatrix-max" || name == "mcq" || name == "mcr" || name == "static-order-mcs" || name == "mcs" || name == "misq" || name == "reduction-misq" || name == "reduction-misr" || name == "reduction-static-order-miss" || name == "reduction-miss" || name == "misr" || name == "static-order-miss" || name == "miss" || name == "reduction-domination-misq" || name == "reduction-domination-misr" || name == "connected-component-miss" || name == "connected-component-miss2" /* || name == "comparison-miss"*/ || name == "tester-miss");

    bool const addDiagonals(bRunExperiment || name == "misq" || name == "misr" || name == "static-order-miss" || name == "miss" || name == "reduction-domination-misq" || name == "reduction-domination-misr" || name == "reduction-misq" || name == "reduction-misr" || name == "reduction-static-order-miss" || name == "reduction-miss" || name == "connected-component-miss" || name == "connected-component-miss2" || name == "comparison-miss" || name == "tester-miss");

    if (bShouldComputeAdjacencyMatrix && !bComputeAdjacencyMatrix) {
        cout << "ERROR!: unable to compute adjacencyMatrix, since the graph is too large: " << adjacencyList.size() << endl << flush;
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
            if (addDiagonals) {
                adjacencyMatrix[i][i]  = 1;
                vAdjacencyMatrix[i][i] = 1;
            }
        }
    }

    bool const bComputeAdjacencyArray(staging || computeCliqueGraph || bRunExperiment || name == "adjlist" || name == "timedelay-adjlist" || name == "generic-adjlist" || name == "generic-adjlist-max" ||name == "timedelay-maxdegree" || name == "timedelay-degeneracy" || name == "faster-degeneracy" || name == "generic-degeneracy" || name == "cache-degeneracy" || name == "mis" || name == "degeneracy-mis" || name == "partial-match-degeneracy" || name == "reverse-degeneracy" || name == "degeneracy-min" || name == "degeneracy-mis-2" || name == "reduction-mis" || name == "experimental-mis" || name == "reduction-misq" || name == "reduction-misr" || name == "reduction-static-order-miss" || name == "reduction-miss" || name == "reduction-sparse-misq" || name == "reduction-sparse-misr" || name == "reduction-sparse-static-order-miss" || name == "reduction-sparse-miss" || name == "reduction-domination-misq" || name == "reduction-domination-misr" || name == "sparse-mcq" || name == "connected-component-miss" || name == "connected-component-miss2"/* || name == "comparison-miss"*/ || name == "tester-miss");

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

    if (bRunExperiment) {
        RunExperiment(inputFile, sExperimentName, bOutputLatex, bPrintHeader, adjacencyArray, vAdjacencyMatrix, dTimeout);
        return 0;
    }

    VertexSets *pSets(nullptr);

    if (name == "tomita") {
        pAlgorithm = new TomitaAlgorithm(adjacencyMatrix, n);
    } else if (name == "generic-adjmatrix") {
        pSets = new AdjacencyMatrixVertexSets(vAdjacencyMatrix);
    } else if (name == "generic-adjmatrix-max") {
        pSets = new AdjacencyMatrixVertexSetsMax(vAdjacencyMatrix);
    } else if (name == "mcq") {
        pAlgorithm = new LightWeightMCQ(vAdjacencyMatrix);
    } else if (name == "sparse-mcq") {
        pAlgorithm = new LightWeightSparseMCQ(adjacencyArray);
    } else if (name == "mcr") {
        pAlgorithm = new LightWeightMCR(vAdjacencyMatrix);
    } else if (name == "static-order-mcs") {
        pAlgorithm = new LightWeightStaticOrderMCS(vAdjacencyMatrix);
    } else if (name == "mcs") {
        pAlgorithm = new LightWeightFullMCS(vAdjacencyMatrix);
    } else if (name == "misq") {
        pAlgorithm = new LightWeightMISQ(vAdjacencyMatrix);
    } else if (name == "misr") {
        pAlgorithm = new LightWeightMISR(vAdjacencyMatrix);
    } else if (name == "static-order-miss") {
        pAlgorithm = new LightWeightStaticOrderMISS(vAdjacencyMatrix);
    } else if (name == "miss") {
        pAlgorithm = new LightWeightFullMISS(vAdjacencyMatrix);
    } else if (name == "reduction-misq") {
        pAlgorithm = new LightWeightReductionMISQ(vAdjacencyMatrix, adjacencyArray);
    } else if (name == "reduction-misr") {
        pAlgorithm = new LightWeightReductionMISR(vAdjacencyMatrix, adjacencyArray);
    } else if (name == "reduction-static-order-miss") {
        pAlgorithm = new LightWeightReductionStaticOrderMISS(vAdjacencyMatrix, adjacencyArray);
    } else if (name == "reduction-miss") {
        pAlgorithm = new LightWeightReductionFullMISS(vAdjacencyMatrix, adjacencyArray);
////    } else if (name == "comparison-miss") {
////        pAlgorithm = new ComparisonFullMISS(vAdjacencyMatrix, adjacencyArray);
    } else if (name == "tester-miss") {
        pAlgorithm = new TesterMISS(vAdjacencyMatrix, adjacencyArray);
    } else if (name == "reduction-domination-misq") {
        pAlgorithm = new LightWeightReductionDominationMISQ(vAdjacencyMatrix, adjacencyArray);
    } else if (name == "reduction-domination-misr") {
        pAlgorithm = new LightWeightReductionDominationMISR(vAdjacencyMatrix, adjacencyArray);
    } else if (name == "reduction-sparse-misq") {
        pAlgorithm = new LightWeightReductionSparseMISQ(adjacencyArray);
    } else if (name == "reduction-sparse-misr") {
        pAlgorithm = new LightWeightReductionSparseMISR(adjacencyArray);
    } else if (name == "reduction-sparse-static-order-miss") {
        pAlgorithm = new LightWeightReductionSparseStaticOrderMISS(adjacencyArray);
    } else if (name == "reduction-sparse-miss") {
        pAlgorithm = new LightWeightReductionSparseFullMISS(adjacencyArray);
    } else if (name == "connected-component-miss") {
        pAlgorithm = new ConnectedComponentMISS(vAdjacencyMatrix, adjacencyArray);
    } else if (name == "connected-component-miss2") {
        pAlgorithm = new ConnectedComponentMISS2(vAdjacencyMatrix, adjacencyArray);
    } else if (name == "adjlist") {
        pAlgorithm = new AdjacencyListAlgorithm(adjacencyArray);
    } else if (name == "generic-adjlist") {
        pSets = new AdjacencyListVertexSets(adjacencyArray);
    } else if (name == "generic-adjlist-max") {
        pSets = new AdjacencyListVertexSetsMax(adjacencyArray);
        if (!findMaximumOnly) { cout << "Will only find maximum subgraph, algorithm does not support listing all subgraphs." << endl; }
        pAlgorithm = new MaximumCliqueAlgorithm(pSets);
    } else if (name == "generic-degeneracy") {
        pSets = new DegeneracyVertexSets(adjacencyArray);
    } else if (name == "reverse-degeneracy") {
        pSets = new ReverseDegeneracyVertexSets(adjacencyArray);
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
    } else if (name == "mis") {
        pSets = new IndependentSets(adjacencyArray);
    } else if (name == "reduction-mis") {
        pSets = new IndependentSetsReduction(adjacencyArray);
    } else if (name == "experimental-mis") {
        pSets = new ExperimentalReduction(adjacencyArray);
        if (!findMaximumOnly) { cout << "Will only find maximum subgraph, algorithm does not support listing all subgraphs." << endl; }
        pAlgorithm = new MaximumCliqueAlgorithm(pSets);
    } else if (name == "degeneracy-mis") {
        pSets = new DegeneracyIndependentSets(adjacencyArray);
    } else if (name == "degeneracy-mis-2") {
        pSets = new DegeneracyIndependentSets2(adjacencyArray);
    } else if (name == "degeneracy-min") {
        pSets = new DegeneracyIndependentSets(adjacencyArray);
        cout << "Will only find minimum subgraph, algorithm does not support listing all subgraphs, or finding maximum subgraph." << endl;
        pAlgorithm = new MinimumCliqueAlgorithm(pSets);
    } else if (name == "timedelay-degeneracy") {
        pAlgorithm = new TimeDelayDegeneracyAlgorithm(adjacencyList);
    } else if (name == "partial-match-degeneracy") {
        pSets = new PartialMatchDegeneracyVertexSets(adjacencyArray);
    } else if (!computeCliqueGraph && !staging && !bRunExperiment){
        cout << "ERROR: unrecognized algorithm name " << name << endl;
        return 1;
    }

    if (computeCliqueGraph) {
        pAlgorithm = new CliqueGraphAlgorithm(new DegeneracyVertexSets(adjacencyArray), adjacencyArray);
    }

    if (staging) {
        pAlgorithm = new Staging(adjacencyArray);
    }

    if (!pAlgorithm && findMaximumOnly) {
        pAlgorithm = new MaximumCliqueAlgorithm(pSets);
    } else if (!pAlgorithm && !findMaximumOnly) {
        pAlgorithm = new BronKerboschAlgorithm(pSets);
    }


#ifdef PRINT_CLIQUES_ONE_BY_ONE
    auto printClique = [](list<int> const &clique) {
        Tools::printList(clique, &Tools::printInt);
    };

    pAlgorithm->AddCallBack(printClique);
#endif //PRINT_CLIQUES_ONE_BY_ONE

    auto verifyMaximalCliqueArray = [&adjacencyArray](list<int> const &clique) {
////        bool const isIS = CliqueTools::IsIndependentSet(adjacencyArray, clique, true /* verbose */);
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

    auto verifyIndependentSetMatrix = [&vAdjacencyMatrix](list<int> const &clique) {
        bool const isIS = CliqueTools::IsIndependentSet(vAdjacencyMatrix, clique, true /* verbose */);
        if (!isIS) {
            cout << "ERROR: Set " << (isIS ? "is" : "is not" ) << " an independent set!" << endl;
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

    RunAndPrintStats(pAlgorithm, cliques, bOutputLatex);

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
