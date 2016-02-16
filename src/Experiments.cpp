#include "Experiments.h"

#include "Isolates4.h"
#include "SparseArraySet.h"
#include "TesterMISS.h"
#include "ForwardSearchMISS.h"
#include "LightWeightFullMISS.h"
#include "LightWeightReductionSparseFullMISS.h"
#include "GraphTools.h"
#include "Tools.h"
#include "OrderingTools.h"
#include "CliqueTools.h"

#include <list>
#include <cstdlib>

using namespace std;

Experiments::Experiments(string const sDataSetName, double const dTimeout, bool const bOutputLatex, bool const bPrintHeader, vector<vector<int>> const &adjacencyArray, vector<vector<char>> const &vAdjacencyMatrix)
 : m_sDataSetName(sDataSetName)
 , m_dTimeout(dTimeout)
 , m_bOutputLatex(bOutputLatex)
 , m_bPrintHeader(bPrintHeader)
 , m_AdjacencyArray(adjacencyArray)
 , m_AdjacencyMatrix(vAdjacencyMatrix)
{
}

void Experiments::RunKernelSize() const
{
    if (m_bPrintHeader) {
        if (m_bOutputLatex) {
            cout << "Graph Name & $n$ & $m$ & $t$ & $k$ & $c$ & $l$ \\\\ \\hline" << endl << flush;
        } else {
            cout << "Graph Name\tn\tm\tt\tk\tc\tl" << endl << flush;
        }
    }

    clock_t startTime(clock());

    size_t const numVertices(m_AdjacencyArray.size());
    size_t numEdges(0);
    for (vector<int> const &neighbors : m_AdjacencyArray) {
        numEdges+= neighbors.size();
    }
    numEdges >>=1;

    Isolates4<SparseArraySet> isolates(m_AdjacencyArray);
    vector<int> vIsolates;
    vector<int> vRemoved;
    vector<Reduction> vReductions;
    isolates.RemoveAllIsolates(0, vIsolates, vRemoved, vReductions, true /* consider all vertices for removal */);

    size_t const kernelSize(isolates.GetInGraph().Size());

    clock_t endTime(clock());

    vector<vector<int>> vComponents;
    GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());

    size_t const numComponents(vComponents.size());
    size_t largestComponentSize(0);
    for (vector<int> const &vComponent : vComponents) {
        largestComponentSize = max(vComponent.size(), largestComponentSize);
    }

    if (m_bOutputLatex) {
        cout << m_sDataSetName << " & " << numVertices << " & " << numEdges << " & " << Tools::GetTimeInSeconds(endTime-startTime) << "&" << kernelSize << " & " << numComponents << " & " << largestComponentSize << " \\\\ " << endl << flush;
    } else {
        cout << m_sDataSetName << "\t" << numVertices << "\t" << numEdges << "\t" << Tools::GetTimeInSeconds(endTime-startTime) << "\t" << kernelSize << "\t" << numComponents << "\t" << largestComponentSize << endl << flush;
    }
}

void Experiments::KernelizeAndRunReductionSparseMISS() const
{
    if (m_bPrintHeader) {
        if (m_bOutputLatex) {
            cout << "Graph Name & $n$ & $m$ & $t$ & $k$ & OPT \\\\ \\hline" << endl << flush;
        } else {
            cout << "Graph Name\tn\tm\tt\tk\tOPT" << endl << flush;
        }
    }

    clock_t startTime(clock());

    size_t const numVertices(m_AdjacencyArray.size());
    size_t numEdges(0);
    for (vector<int> const &neighbors : m_AdjacencyArray) {
        numEdges+= neighbors.size();
    }
    numEdges >>=1;

    Isolates4<SparseArraySet> isolates(m_AdjacencyArray);
    vector<int> vIsolates;
    vector<int> vRemoved;
    vector<Reduction> vReductions;
    isolates.RemoveAllIsolates(0, vIsolates, vRemoved, vReductions, true /* consider all vertices for removal */);

    size_t const kernelSize(isolates.GetInGraph().Size());

    vector<vector<int>> vAdjacencyArray;
    vAdjacencyArray.resize(isolates.GetInGraph().Size());

    map<int,int> vertexRemap;
    size_t uNewIndex = 0;

    for (int const vertex : isolates.GetInGraph()) {
        vertexRemap[vertex] = uNewIndex++;
    }

    for (pair<int,int> const &mapPair : vertexRemap) {
        int const oldVertex(mapPair.first);
        int const newVertex(mapPair.second);
        for (int const neighbor : isolates.Neighbors()[oldVertex]) {
            if (vertexRemap.find(neighbor) != vertexRemap.end()) {
////                cout << "vAdjacencyArray.size=" << vAdjacencyArray.size() << endl;
////                cout << "newVertex           =" << newVertex << endl; 
                vAdjacencyArray[newVertex].push_back(vertexRemap[neighbor]);
            }
        }
    }

    LightWeightReductionSparseFullMISS algorithm(vAdjacencyArray);
    algorithm.SetQuiet(true);

    list<list<int>> indsets;
    algorithm.Run(indsets);

    clock_t endTime(clock());

    int const solutionSize(indsets.back().size() + vReductions.size());

    if (m_bOutputLatex) {
        cout << m_sDataSetName << " & " << numVertices << " & " << numEdges << " & " << Tools::GetTimeInSeconds(endTime-startTime) << "&" << kernelSize << " & " << solutionSize << " \\\\ " << endl << flush;
    } else {
        cout << m_sDataSetName << "\t" << numVertices << "\t" << numEdges << "\t" << Tools::GetTimeInSeconds(endTime-startTime) << "\t" << kernelSize << "\t" << solutionSize << endl << flush;
    }
}

void Experiments::KernelizeAndRunComponentWiseMISS() const
{
    if (m_bPrintHeader) {
        if (m_bOutputLatex) {
            cout << "Graph Name & $n$ & $m$ & $t$ & $k$ & c & l & OPT \\\\ \\hline" << endl << flush;
        } else {
            cout << "Graph Name\tn\tm\tt\tk\tc\tl\tOPT" << endl << flush;
        }
    }

    clock_t startTime(clock());

    size_t const numVertices(m_AdjacencyArray.size());
    size_t numEdges(0);
    for (vector<int> const &neighbors : m_AdjacencyArray) {
        numEdges+= neighbors.size();
    }
    numEdges >>=1;

    Isolates4<SparseArraySet> isolates(m_AdjacencyArray);
    vector<int> vIsolates;
    vector<int> vRemoved;
    vector<Reduction> vReductions;
    isolates.RemoveAllIsolates(0, vIsolates, vRemoved, vReductions, true /* consider all vertices for removal */);

    size_t const kernelSize(isolates.GetInGraph().Size());

    vector<vector<int>> vComponents;
    GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());

    size_t const numComponents(vComponents.size());
    size_t largestComponentSize(0);
    for (vector<int> const &vComponent : vComponents) {
        largestComponentSize = max(vComponent.size(), largestComponentSize);
    }

    size_t solutionDelta(0);

    for (vector<int> const &vComponent : vComponents) {

        if (vComponent.size() > 20000) {
            cerr << "ERROR!: unable to compute adjacencyMatrix, since the graph is too large: " << vComponent.size() << endl << flush;
            exit(1);
        }

        vector<vector<char>> vAdjacencyMatrix;
        vAdjacencyMatrix.resize(vComponent.size());

        map<int,int> vertexRemap;
////        map<int,int> reverseMap;
        size_t uNewIndex = 0;

        for (int const vertex : vComponent) {
            vertexRemap[vertex] = uNewIndex++;
////            reverseMap[uNewIndex-1] = vertex;
        }

        for (pair<int,int> const &mapPair : vertexRemap) {
            int const oldVertex(mapPair.first);
            int const newVertex(mapPair.second);

            vAdjacencyMatrix[newVertex].resize(vComponent.size(), 0);
            for (int const neighbor : isolates.Neighbors()[oldVertex]) {
                if (vertexRemap.find(neighbor) != vertexRemap.end()) {
                    //cout << "newVertex           =" << newVertex << endl; 
                    vAdjacencyMatrix[newVertex][vertexRemap[neighbor]] = 1;
                }
            }
        }

        ////cout << "vAdjacencyArray.size=" << vAdjacencyArray.size() << endl;
        LightWeightFullMISS algorithm(vAdjacencyMatrix);
        algorithm.SetQuiet(true);

        list<list<int>> indsets;
        algorithm.Run(indsets);

#ifdef VERIFY
        cout << "Verifying independent set:" << endl;
        set<int> const indepset(indsets.back().begin(), indsets.back().end());
        for (int const vertex : indepset) {
            for (int const otherVertex : indepset) {
                if (vAdjacencyMatrix[vertex][otherVertex])
                    cout << "    ERROR: " << vertex << " and " << otherVertex << " are neighbors" << endl << flush;
            }
        }
#endif // VERIFY

        solutionDelta += indsets.back().size();
    }

    clock_t endTime(clock());

    int const solutionSize(solutionDelta + vReductions.size());

    if (m_bOutputLatex) {
        cout << m_sDataSetName << " & " << numVertices << " & " << numEdges << " & " << Tools::GetTimeInSeconds(endTime-startTime) << "&" << kernelSize << " & " << numComponents << " & " << largestComponentSize << " & " << solutionSize << " \\\\ " << endl << flush;
    } else {
        cout << m_sDataSetName << "\t" << numVertices << "\t" << numEdges << "\t" << Tools::GetTimeInSeconds(endTime-startTime) << "\t" << kernelSize << "\t" << numComponents << "\t" << largestComponentSize << "\t" << solutionSize << endl << flush;
    }

}

void Experiments::KernelizeAndRunComponentWiseReductionSparseMISS() const
{
    if (m_bPrintHeader) {
        if (m_bOutputLatex) {
            cout << "Graph Name & $n$ & $m$ & $t$ & $k$ & c & l & OPT \\\\ \\hline" << endl << flush;
        } else {
            cout << "Graph Name\tn\tm\tt\tk\tc\tl\tOPT" << endl << flush;
        }
    }

    clock_t startTime(clock());

    size_t const numVertices(m_AdjacencyArray.size());
    size_t numEdges(0);
    for (vector<int> const &neighbors : m_AdjacencyArray) {
        numEdges+= neighbors.size();
    }
    numEdges >>=1;

    Isolates4<SparseArraySet> isolates(m_AdjacencyArray);
    vector<int> vIsolates;
    vector<int> vRemoved;
    vector<Reduction> vReductions;
    isolates.RemoveAllIsolates(0, vIsolates, vRemoved, vReductions, true /* consider all vertices for removal */);

    size_t const kernelSize(isolates.GetInGraph().Size());

    vector<vector<int>> vComponents;
    GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());

    size_t const numComponents(vComponents.size());
    size_t largestComponentSize(0);
    for (vector<int> const &vComponent : vComponents) {
        largestComponentSize = max(vComponent.size(), largestComponentSize);
    }

    size_t solutionDelta(0);

    for (vector<int> const &vComponent : vComponents) {

        vector<vector<int>> vAdjacencyArray;
        vAdjacencyArray.resize(vComponent.size());

        map<int,int> vertexRemap;
////        map<int,int> reverseMap;
        size_t uNewIndex = 0;

        for (int const vertex : vComponent) {
            vertexRemap[vertex] = uNewIndex++;
////            reverseMap[uNewIndex-1] = vertex;
        }

        for (pair<int,int> const &mapPair : vertexRemap) {
            int const oldVertex(mapPair.first);
            int const newVertex(mapPair.second);
            for (int const neighbor : isolates.Neighbors()[oldVertex]) {
                if (vertexRemap.find(neighbor) != vertexRemap.end()) {
                    //cout << "newVertex           =" << newVertex << endl; 
                    vAdjacencyArray[newVertex].push_back(vertexRemap[neighbor]);
                }
            }
        }

        cout << "vAdjacencyArray.size=" << vAdjacencyArray.size() << endl;
        LightWeightReductionSparseFullMISS algorithm(vAdjacencyArray);

        list<list<int>> indsets;
        algorithm.Run(indsets);

        cout << "Verifying independent set:" << endl;
        set<int> const indepset(indsets.back().begin(), indsets.back().end());
        for (int const vertex : indepset) {
            for (int const neighbor : vAdjacencyArray[vertex]) {
                if (indepset.find(neighbor) != indepset.end()) {
                    cout << "    ERROR: " << vertex << " and " << neighbor << " are neighbors" << endl << flush;
                }
            }
        }

        solutionDelta += indsets.back().size();
    }

    clock_t endTime(clock());

    int const solutionSize(solutionDelta + vReductions.size());

    if (m_bOutputLatex) {
        cout << m_sDataSetName << " & " << numVertices << " & " << numEdges << " & " << Tools::GetTimeInSeconds(endTime-startTime) << "&" << kernelSize << " & " << numComponents << " & " << largestComponentSize << " & " << solutionSize << " \\\\ " << endl << flush;
    } else {
        cout << m_sDataSetName << "\t" << numVertices << "\t" << numEdges << "\t" << Tools::GetTimeInSeconds(endTime-startTime) << "\t" << kernelSize << "\t" << numComponents << "\t" << largestComponentSize << "\t" << solutionSize << endl << flush;
    }
}

void Experiments::RunStandardSearch() const
{
    list<list<int>> lCliques;
    TesterMISS algorithm(m_AdjacencyMatrix, m_AdjacencyArray);
    if (m_dTimeout > 0.0) {
        cout << "Using timeout: " << m_dTimeout << endl << flush;
        algorithm.SetTimeOutInSeconds(m_dTimeout);
    }
    algorithm.Run(lCliques);
    cout << "Found independent set of size " << lCliques.back().size() << endl << flush;
}

void Experiments::RunExactSearch() const
{
    if (m_bPrintHeader) {
        if (m_bOutputLatex) {
            cout << "Graph Name & $\\alpha(G)$ ${\\sc RMISS}$ & {\\sc MISS} \\\\ \\hline" << endl << flush;
        } else {
            cout << "Graph Name\ta(G)\tRMISS\tMISS" << endl << flush;
        }
    }

    cout << m_sDataSetName;
    RunComponentsStandardSearch();
////    if (m_bOutputLatex) cout << " & ";
////    else cout << "\t";
    RunComponentsMISS();
    if (m_bOutputLatex) cout << " \\\\" << endl;
    else cout << endl;
}

void Experiments::RunComponentsStandardSearch() const
{
    Isolates4<SparseArraySet> isolates(m_AdjacencyArray);
    vector<int> vIsolates;
    vector<int> vRemoved;
    vector<Reduction> vReductions;
    isolates.RemoveAllIsolates(0, vIsolates, vRemoved, vReductions, true /* consider all vertices for removal */);

    clock_t startTime(clock());
    vector<vector<int>> vComponents;
    GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());
////    cerr << "# vertices remaining in graph: " << isolates.GetInGraph().Size() << "/" << m_AdjacencyArray.size() << endl << flush;
////    cerr << "# connected components       : " << vComponents.size() << endl << flush;
////    cerr << "size of connected components : ";
////    cout << "[ ";
////    for (size_t index = 0; index < vComponents.size(); ++index) {
////        vector<int> const& vComponent(vComponents[index]);
////        cout << vComponent.size() << " ";
////    }
////    cout << "]" << endl;

    size_t totalCliqueSize(isolates.size() + isolates.GetFoldedVertexCount());

////    cout << "Initial independent set size: " << totalCliqueSize << endl << flush;

    for (vector<int> const &vComponent : vComponents) {
        if (vComponent.empty()) continue;
        size_t const previousTotalCliqueSize(totalCliqueSize);
        vector<vector<int>> componentAdjacencyList;
        set<int> componentSet(vComponent.begin(), vComponent.end());
        map<int,int> mapUnused;
        GraphTools::ComputeInducedSubgraphIsolates(isolates, componentSet, componentAdjacencyList, mapUnused);

////        cout << "component-size=" << componentAdjacencyList.size() << " " << Tools::GetTimeInSeconds(clock() - startTime) << endl << flush;

        vector<vector<char>> componentAdjacencyMatrix(componentAdjacencyList.size());
        for (size_t index = 0; index < componentAdjacencyList.size(); ++index) {
            componentAdjacencyMatrix[index].resize(componentAdjacencyList.size(), 0);
            componentAdjacencyMatrix[index][index] = 1;
            for (int const neighbor : componentAdjacencyList[index]) {
                componentAdjacencyMatrix[index][neighbor] = 1;
            }
        }

        clock_t newTime(clock());
        TesterMISS algorithm(componentAdjacencyMatrix, componentAdjacencyList);
        algorithm.SetQuiet(false); 
            ////        algorithm.SetOnlyVertex(vOrdering[splitPoint]);
        algorithm.SetTimeOutInSeconds(m_dTimeout);
        list<list<int>> cliques;
////        algorithm.SetQuiet(true);
        algorithm.Run(cliques);

////            if (algorithm)

        if (previousTotalCliqueSize + cliques.back().size() > totalCliqueSize) {
            totalCliqueSize = cliques.back().size() + previousTotalCliqueSize;
////            cout << "Found a better independent set: size=" << totalCliqueSize << "=" << cliques.back().size() << "+" << previousTotalCliqueSize << endl;
        }
    }

    clock_t const endTime(clock());

////    cout << "Found independent set of size " << totalCliqueSize << endl << flush;

    if (m_bOutputLatex) {
        cout << " & " << totalCliqueSize << " & " << Tools::GetTimeInSeconds(endTime - startTime, false /* no brackets */) << flush;
    } else {
        cout << "\t" << totalCliqueSize <<  "\t" << Tools::GetTimeInSeconds(endTime - startTime, false /* no brackets */) << flush;
    }
}

void Experiments::RunComponentsMISS() const
{
    Isolates4<SparseArraySet> isolates(m_AdjacencyArray);
    vector<int> vIsolates;
    vector<int> vRemoved;
    vector<Reduction> vReductions;
    isolates.RemoveAllIsolates(0, vIsolates, vRemoved, vReductions, true /* consider all vertices for removal */);

    clock_t startTime(clock());
    vector<vector<int>> vComponents;
    GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());
////    cerr << "# vertices remaining in graph: " << isolates.GetInGraph().Size() << "/" << m_AdjacencyArray.size() << endl << flush;
////    cerr << "# connected components       : " << vComponents.size() << endl << flush;
////    cerr << "size of connected components : ";
////    cout << "[ ";
////    for (size_t index = 0; index < vComponents.size(); ++index) {
////        vector<int> const& vComponent(vComponents[index]);
////        cout << vComponent.size() << " ";
////    }
////    cout << "]" << endl;

    size_t totalCliqueSize(isolates.size() + isolates.GetFoldedVertexCount());

////    cout << "Initial independent set size: " << totalCliqueSize << endl << flush;

    for (vector<int> const &vComponent : vComponents) {
        if (vComponent.empty()) continue;
        size_t const previousTotalCliqueSize(totalCliqueSize);
        vector<vector<int>> componentAdjacencyList;
        set<int> componentSet(vComponent.begin(), vComponent.end());
        map<int,int> mapUnused;
        GraphTools::ComputeInducedSubgraphIsolates(isolates, componentSet, componentAdjacencyList, mapUnused);

////        cout << "component-size=" << componentAdjacencyList.size() << " " << Tools::GetTimeInSeconds(clock() - startTime) << endl << flush;

        vector<vector<char>> componentAdjacencyMatrix(componentAdjacencyList.size());
        for (size_t index = 0; index < componentAdjacencyList.size(); ++index) {
            componentAdjacencyMatrix[index].resize(componentAdjacencyList.size(), 0);
            componentAdjacencyMatrix[index][index] = 1;
            for (int const neighbor : componentAdjacencyList[index]) {
                componentAdjacencyMatrix[index][neighbor] = 1;
            }
        }

        clock_t newTime(clock());
        LightWeightFullMISS algorithm(componentAdjacencyMatrix);
        algorithm.SetQuiet(false); 
            ////        algorithm.SetOnlyVertex(vOrdering[splitPoint]);
        algorithm.SetTimeOutInSeconds(m_dTimeout);
        list<list<int>> cliques;
////        algorithm.SetQuiet(true);
        algorithm.Run(cliques);

////            if (algorithm)

        if (previousTotalCliqueSize + cliques.back().size() > totalCliqueSize) {
            totalCliqueSize = cliques.back().size() + previousTotalCliqueSize;
////            cout << "Found a better independent set: size=" << totalCliqueSize << "=" << cliques.back().size() << "+" << previousTotalCliqueSize << endl;
        }
    }

    clock_t const endTime(clock());

////    cout << "Found independent set of size " << totalCliqueSize << endl << flush;

    if (m_bOutputLatex) {
        cout << " & " << Tools::GetTimeInSeconds(endTime - startTime, false /* no brackets */) << flush;
    } else {
        cout << "\t" << Tools::GetTimeInSeconds(endTime - startTime, false /* no brackets */) << flush;
    }
}

void Experiments::RunComponentsForwardSearch() const
{

    if (m_bPrintHeader) {
        if (m_bOutputLatex) {
            cout << "Graph Name & $i$ & time \\\\ \\hline" << endl << flush;
        } else {
            cout << "Graph Name\ti\ttime" << endl << flush;
        }
    }

    Isolates4<SparseArraySet> isolates(m_AdjacencyArray);
    vector<int> vIsolates;
    vector<int> vRemoved;
    vector<Reduction> vReductions;
    isolates.RemoveAllIsolates(0, vIsolates, vRemoved, vReductions, true /* consider all vertices for removal */);

    clock_t startTime(clock());
    vector<vector<int>> vComponents;
    GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());
////    cerr << "# vertices remaining in graph: " << isolates.GetInGraph().Size() << "/" << m_AdjacencyArray.size() << endl << flush;
////    cerr << "# connected components       : " << vComponents.size() << endl << flush;
////    cerr << "size of connected components : ";
////    cout << "[ ";
////    for (size_t index = 0; index < vComponents.size(); ++index) {
////        vector<int> const& vComponent(vComponents[index]);
////        cout << vComponent.size() << " ";
////    }
////    cout << "]" << endl;

    size_t totalCliqueSize(isolates.size() + isolates.GetFoldedVertexCount());

////    cout << "Initial independent set size: " << totalCliqueSize << endl << flush;

    for (vector<int> const &vComponent : vComponents) {
        if (vComponent.empty()) continue;
        size_t const previousTotalCliqueSize(totalCliqueSize);
        vector<vector<int>> componentAdjacencyList;
        set<int> componentSet(vComponent.begin(), vComponent.end());
        map<int,int> mapUnused;
        GraphTools::ComputeInducedSubgraphIsolates(isolates, componentSet, componentAdjacencyList, mapUnused);

////        cout << "component-size=" << componentAdjacencyList.size() << " " << Tools::GetTimeInSeconds(clock() - startTime) << endl << flush;

        vector<vector<char>> componentAdjacencyMatrix;
////        (componentAdjacencyList.size());
////        for (size_t index = 0; index < componentAdjacencyList.size(); ++index) {
////            componentAdjacencyMatrix[index].resize(componentAdjacencyList.size(), 0);
////            componentAdjacencyMatrix[index][index] = 1;
////            for (int const neighbor : componentAdjacencyList[index]) {
////                componentAdjacencyMatrix[index][neighbor] = 1;
////            }
////        }

        vector<int> vOrdering;
        vector<int> vColoring;
        size_t cliqueSize(0);
        OrderingTools::InitialOrderingMISR(componentAdjacencyList, vOrdering, vColoring, cliqueSize);

        if (previousTotalCliqueSize + cliqueSize > totalCliqueSize) {
            totalCliqueSize = previousTotalCliqueSize + cliqueSize;
////            cout << "After ordering, found new independent set of size: " << totalCliqueSize << endl << flush;
        }

        //    size_t pivot(0);
        //    for (size_t afterIndex = vColoring.size(); afterIndex > 0; afterIndex--) {
        //        if (vColoring[afterIndex-1] <= cliqueSize) { pivot = afterIndex-1; break; }
        //    }

        ////    set<int>     setVertices;
////        size_t const startVertex(min(vComponent.size()-1, static_cast<size_t>(100)));  ////static_cast<size_t>(vOrdering.size()*0.90)));
////        size_t const startVertex(vOrdering.size()-1);
        ////    size_t const startVertex(0);
////        set<int>     setVertices(vOrdering.begin() + startVertex, vOrdering.end());
        size_t minIndex(0), maxIndex(vOrdering.size());
////        for (size_t splitPoint = startVertex/*vOrdering.size()/2*/; splitPoint < vOrdering.size(); splitPoint++) {
        while (true) {
            ////    for (size_t splitPoint = startVertex/*vOrdering.size()/2*/; splitPoint > 0; splitPoint--) {
////            setVertices.insert(vOrdering[splitPoint]);
            size_t const currentIndex((minIndex + maxIndex)/2);
            cout << "Evaluating index: " << currentIndex << " in [" << minIndex << "," << maxIndex << "]" << endl << flush;
            set<int> const setVertices(vOrdering.begin(), vOrdering.begin() + currentIndex);
            vector<vector<int>> subgraphAdjacencyList;

            map<int,int> mapUnused2;
            GraphTools::ComputeInducedSubgraph(componentAdjacencyList, setVertices, subgraphAdjacencyList, mapUnused2);

////            cout << "Subgraph size=" << subgraphAdjacencyList.size() << " " << Tools::GetTimeInSeconds(clock() - startTime) << endl << flush;
////            cout << "subgraph=";
////            for (int const vertex : setVertices) {
////                cout << vertex << " ";
////            }
////            cout << endl;

            list<list<int>> cliques;
            list<int> realClique;
            size_t cliqueDelta(0);
            bool timedOut(false);
#if 0
            vector<vector<char>> subgraphAdjacencyMatrix; ////(subgraphAdjacencyList.size());
////            for (size_t index = 0; index < subgraphAdjacencyList.size(); ++index) {
////                subgraphAdjacencyMatrix[index].resize(subgraphAdjacencyList.size(), 0);
////                subgraphAdjacencyMatrix[index][index] = 1;
////                for (int const neighbor : subgraphAdjacencyList[index]) {
////                    subgraphAdjacencyMatrix[index][neighbor] = 1;
////                }
////            }

            clock_t newTime(clock());
            TesterMISS algorithm(subgraphAdjacencyMatrix, subgraphAdjacencyList);
////            TesterMISS algorithm(componentAdjacencyMatrix, componentAdjacencyList);
                algorithm.SetQuiet(true); 
            ////        algorithm.SetOnlyVertex(vOrdering[splitPoint]);
            algorithm.SetTimeOutInSeconds(m_dTimeout);
            algorithm.Run(cliques);
            timedOut = algorithm.GetTimedOut();

////            if (algorithm)

            if (previousTotalCliqueSize + cliques.back().size() + cliqueDelta > totalCliqueSize) {
                totalCliqueSize = cliques.back().size() + cliqueDelta + previousTotalCliqueSize;
                cout << "At " << setVertices.size() << "/" << vOrdering.size() << ", found a better independent set: size=" << totalCliqueSize << "=" << cliques.back().size() << "+" << cliqueDelta << "+" << previousTotalCliqueSize << endl;
            }
            realClique = std::move(cliques.back());
#else
            Isolates4<SparseArraySet> newIsolates(subgraphAdjacencyList);
            vector<int> vNewRemoved;
            vector<Reduction> vNewReductions;
            newIsolates.RemoveAllIsolates(0, vNewRemoved, vNewRemoved, vNewReductions, true);

////            size_t const cliqueDelta(newIsolates.size() + newIsolates.GetFoldedVertexCount());
            cliqueDelta = (newIsolates.size() + newIsolates.GetFoldedVertexCount());

            vector<vector<int>> vNewComponents;
            GraphTools::ComputeConnectedComponents(newIsolates, vNewComponents, subgraphAdjacencyList.size());
////            cerr << "# vertices remaining in graph: " << newIsolates.GetInGraph().Size() << "/" << subgraphAdjacencyList.size() << endl << flush;
////            cerr << "Vertices: ";
////            for (int const vertex : newIsolates.GetInGraph()) {
////                cout << vertex << " ";
////            }
////            cout << endl;
////            cerr << "# connected components       : " << vNewComponents.size() << endl << flush;
////            cerr << "size of connected components : ";
////            cout << "[ ";
////
////            for (size_t index = 0; index < vNewComponents.size(); ++index) {
////                vector<int> const& vNewComponent(vNewComponents[index]);
////                cout << vNewComponent.size() << " ";
////            }
////
////            cout << "]" << endl;
            for (vector<int> const &vNewComponent : vNewComponents) {
                if (vNewComponent.empty()) continue;
                cliques.clear();
                map<int,int> mapUnused3;
                vector<vector<int>> newComponentAdjacencyList;
                set<int> setVertices3(vNewComponent.begin(), vNewComponent.end());
                GraphTools::ComputeInducedSubgraphIsolates(newIsolates, setVertices3, newComponentAdjacencyList, mapUnused3);

                vector<vector<char>> newComponentAdjacencyMatrix;
////                (newComponentAdjacencyList.size());
////                for (size_t index = 0; index < newComponentAdjacencyList.size(); ++index) {
////                    newComponentAdjacencyMatrix[index].resize(newComponentAdjacencyList.size(), 0);
////                    newComponentAdjacencyMatrix[index][index] = 1;
////                    for (int const neighbor : newComponentAdjacencyList[index]) {
////                        newComponentAdjacencyMatrix[index][neighbor] = 1;
////                    }
////                }
                TesterMISS algorithm(newComponentAdjacencyMatrix, newComponentAdjacencyList);
                ////            TesterMISS algorithm(componentAdjacencyMatrix, componentAdjacencyList);
                algorithm.SetQuiet(true); 
                ////        algorithm.SetOnlyVertex(vOrdering[splitPoint]);
                algorithm.SetTimeOutInSeconds(m_dTimeout);
                algorithm.Run(cliques);
////                cout << "Clique size increase: " << realClique.size() << " to ";
                realClique.insert(realClique.end(), cliques.back().begin(), cliques.back().end());
////                cout << realClique.size() << endl << flush;
                ////cliqueDelta+=cliques.back().size();
                cliques.clear();
////                cliques.push_back(list<int>());
                timedOut = timedOut || algorithm.GetTimedOut();
            }
#endif // 0

////            if (componentAdjacencyMatrix.size() == 86) {
////                map<int,int> forwardMap;
////                for (pair<int,int> const &mapPair : mapUnused) {
////                    forwardMap[mapPair.second] = mapPair.first;
////                }
////                PrintSubgraphInEdgesFormat(m_AdjacencyList, forwardMap);
////            }

            if (previousTotalCliqueSize + realClique.size() + cliqueDelta > totalCliqueSize) {
                totalCliqueSize = realClique.size() + cliqueDelta + previousTotalCliqueSize;
                cout << "At " << setVertices.size() << "/" << vOrdering.size() << ", found a better independent set: size=" << totalCliqueSize << "=" << realClique.size() << "+" << cliqueDelta << "+" << previousTotalCliqueSize << endl;
            }

            if (minIndex == maxIndex) break;

            if (currentIndex == minIndex) {
                minIndex = maxIndex;
                if (timedOut) break;
            } else if (timedOut) {
                maxIndex = currentIndex;
            } else {
                minIndex = currentIndex;
            }

////            if (clock() - newTime > (0.1*CLOCKS_PER_SEC)) {
////                break;
////            }
        }
    }

    clock_t const endTime(clock());

    cout << m_sDataSetName;
    if (m_bOutputLatex) {
        cout << " & " << totalCliqueSize << " & " << Tools::GetTimeInSeconds(endTime - startTime, false /* no brackets */) << " \\\\" << endl << flush;
    } else {
        cout << "\t" << totalCliqueSize << "\t" << Tools::GetTimeInSeconds(endTime - startTime, false /* no brackets */) << endl << flush;
    }
}

void Experiments::RunForwardSearch() const
{
    if (m_bPrintHeader) {
        if (m_bOutputLatex) {
            cout << "Graph Name & $i$ & time \\\\ \\hline" << endl << flush;
        } else {
            cout << "Graph Name\ti\ttime" << endl << flush;
        }
    }

    clock_t startTime(clock());
    Isolates4<SparseArraySet> isolates(m_AdjacencyArray);
    vector<int> vIsolates;
    vector<int> vRemoved;
    vector<Reduction> vReductions;
    isolates.RemoveAllIsolates(0, vIsolates, vRemoved, vReductions, true /* consider all vertices for removal */);

    map<int,int> mapUnused;
    vector<vector<int>> subgraphAdjacencyList;
    set<int> setVertices(isolates.GetInGraph().begin(), isolates.GetInGraph().end());
    GraphTools::ComputeInducedSubgraphIsolates(isolates, setVertices, subgraphAdjacencyList, mapUnused);

    vector<vector<char>> subgraphAdjacencyMatrix;

    auto printCliqueSize = [](list<int> const &clique) {
        cout << "Found clique of size " << clique.size() << endl << flush;
    };

////    auto printClique = [](list<int> const &clique) {
////        cout << "Clique: ";
////        for (int const vertex : clique) {
////            cout << vertex << " ";
////        }
////        cout << endl;
////    };

    list<list<int>> cliques;
    ForwardSearchMISS algorithm(subgraphAdjacencyMatrix, subgraphAdjacencyList);
    algorithm.SetQuiet(true); 
////    algorithm.SetOnlyVertex(vOrdering[splitPoint]);
    algorithm.SetTimeOutInSeconds(m_dTimeout);
    algorithm.AddCallBack(printCliqueSize);
    algorithm.Run(cliques);

    clock_t const endTime(clock());

    size_t totalCliqueSize = isolates.size() + isolates.GetFoldedVertexCount() + cliques.back().size();

    cout << m_sDataSetName;
    if (m_bOutputLatex) {
        cout << " & " << totalCliqueSize << " & " << Tools::GetTimeInSeconds(endTime - startTime, false /* no brackets */) << " \\\\" << endl << flush;
    } else {
        cout << "\t" << totalCliqueSize << "\t" << Tools::GetTimeInSeconds(endTime - startTime, false /* no brackets */) << endl << flush;
    }
}

void Experiments::ComputeCriticalIndependentSet() const
{
    if (m_bPrintHeader) {
        if (m_bOutputLatex) {
            cout << "Graph Name & $n$ & $m$ & $t$ & $k$ & c & l \\\\ \\hline" << endl << flush;
        } else {
            cout << "Graph Name\tn\tm\tt\tk\tc\tl" << endl << flush;
        }
    }


    size_t const numVertices(m_AdjacencyArray.size());
    size_t numEdges(0);
    for (vector<int> const &neighbors : m_AdjacencyArray) {
        numEdges+= neighbors.size();
    }
    numEdges >>=1;

    Isolates4<SparseArraySet> isolates(m_AdjacencyArray);

    clock_t startTime(clock());

    set<int> const criticalSet(CliqueTools::ComputeCriticalIndependentSet(m_AdjacencyArray));
    cout << "Critical set (" << criticalSet.size() << " elements):" << endl;

#if 0
    size_t const kernelSize(isolates.GetInGraph().Size());

    vector<vector<int>> vComponents;
    GraphTools::ComputeConnectedComponents(isolates, vComponents, m_AdjacencyArray.size());

    size_t const numComponents(vComponents.size());
    size_t largestComponentSize(0);
    for (vector<int> const &vComponent : vComponents) {
        largestComponentSize = max(vComponent.size(), largestComponentSize);
    }

    if (m_bOutputLatex) {
        cout << m_sDataSetName << " & " << numVertices << " & " << numEdges << " & " << Tools::GetTimeInSeconds(endTime-startTime) << "&" << kernelSize << " & " << numComponents << " & " << largestComponentSize << " \\\\ " << endl << flush;
    } else {
        cout << m_sDataSetName << "\t" << numVertices << "\t" << numEdges << "\t" << Tools::GetTimeInSeconds(endTime-startTime) << "\t" << kernelSize << "\t" << numComponents << "\t" << largestComponentSize << endl << flush;
    }
#endif
}

void Experiments::ComputeCriticalIndependentSetKernel() const
{
    if (m_bPrintHeader) {
        if (m_bOutputLatex) {
            cout << "Graph Name & $n$ & $m$ & $t$ & $k$ \\\\ \\hline" << endl << flush;
        } else {
            cout << "Graph Name\tn\tm\tt\tk" << endl << flush;
        }
    }


    size_t const numVertices(m_AdjacencyArray.size());
    size_t numEdges(0);
    for (vector<int> const &neighbors : m_AdjacencyArray) {
        numEdges+= neighbors.size();
    }
    numEdges >>=1;

    clock_t startTime(clock());

    set<int> const remainingVertices(CliqueTools::IterativelyRemoveCriticalIndependentSets(m_AdjacencyArray));
////    cout << "Remaining graph (" << remainingVertices.size() << " elements):" << endl;

    clock_t endTime(clock());

    if (m_bOutputLatex) {
        cout << m_sDataSetName << " & " << numVertices << " & " << numEdges << " & " << Tools::GetTimeInSeconds(endTime-startTime) << "&" << remainingVertices.size()/2 << " \\\\ " << endl << flush;
    } else {
        cout << m_sDataSetName << "\t" << numVertices << "\t" << numEdges << "\t" << Tools::GetTimeInSeconds(endTime-startTime) << "\t" << remainingVertices.size()/2 << endl << flush;
    }
}

void Experiments::ComputeMaximumCriticalIndependentSetKernel() const
{
    if (m_bPrintHeader) {
        if (m_bOutputLatex) {
            cout << "Graph Name & $n$ & $m$ & $t$ & $|I|$ & $k$ \\\\ \\hline" << endl << flush;
        } else {
            cout << "Graph Name\tn\tm\tt\t|I|\tk" << endl << flush;
        }
    }


    size_t const numVertices(m_AdjacencyArray.size());
    size_t numEdges(0);
    for (vector<int> const &neighbors : m_AdjacencyArray) {
        numEdges+= neighbors.size();
    }
    numEdges >>=1;

    clock_t startTime(clock());

    set<int> independentVertices;
    set<int> const remainingVertices(CliqueTools::IterativelyRemoveMaximumCriticalIndependentSets(m_AdjacencyArray, independentVertices));
////    cout << "Remaining graph (" << remainingVertices.size() << " elements):" << endl;

    clock_t endTime(clock());

    if (m_bOutputLatex) {
        cout << m_sDataSetName << " & " << numVertices << " & " << numEdges << " & " << Tools::GetTimeInSeconds(endTime-startTime) << " & " << independentVertices.size() << " & " << remainingVertices.size()/2 << " \\\\ " << endl << flush;
    } else {
        cout << m_sDataSetName << "\t" << numVertices << "\t" << numEdges << "\t" << Tools::GetTimeInSeconds(endTime-startTime) << "\t" << independentVertices.size() << "\t" << remainingVertices.size()/2 << endl << flush;
    }
}

