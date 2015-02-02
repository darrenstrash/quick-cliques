#include "LightWeightReductionSparseMISQ.h"
#include "OrderingTools.h"
#include "GraphTools.h"

#include <cmath>
#include <iostream>

using namespace std;

LightWeightReductionSparseMISQ::LightWeightReductionSparseMISQ(vector<vector<char>> const &vAdjacencyMatrix, vector<vector<int>> const &vAdjacencyArray)
: Algorithm("MISQ")
, m_AdjacencyMatrix(vAdjacencyMatrix)
, m_AdjacencyArray(vAdjacencyArray)
, m_uMaximumCliqueSize(0)
, stackP(vAdjacencyMatrix.size())
, stackColors(vAdjacencyMatrix.size())
, stackOrder(vAdjacencyMatrix.size())
, nodeCount(0)
, isolates(vAdjacencyArray)
, coloringStrategy(isolates.Neighbors())
////, m_bInvert(0)
{
}

////void LightWeightReductionSparseMISQ::SetInvert(bool const invert)
////{
////    m_bInvert = invert;
////}

void LightWeightReductionSparseMISQ::InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors)
{
    P = std::move(GraphTools::OrderVerticesByDegree(isolates.GetInGraph(), isolates.Neighbors(), true /* non-decreasing*/));

    size_t maxDegree(0);
#ifdef SPARSE
    for (SparseArraySet const &neighborSet : isolates.Neighbors()) {
#else
    for (ArraySet const &neighborSet : isolates.Neighbors()) {
#endif //SPARSE
        maxDegree = max(maxDegree, neighborSet.Size());
    }

    vColors.reserve(P.size());
    vColors.clear();
    for (int degree = 1; degree <= maxDegree; degree++ ) {
        vColors.push_back(degree);
    }

    vColors.resize(P.size(), maxDegree + 1);

    vVertexOrder = P;
}

long LightWeightReductionSparseMISQ::Run(list<std::list<int>> &cliques)
{
    R.reserve(m_AdjacencyMatrix.size());

    // don't reserve for 0-th vectors, they get std::move'd by InitialOrdering
    for (int index = 1; index < stackP.size(); ++index) {
        stackP.reserve(m_AdjacencyMatrix.size());
        stackColors.reserve(m_AdjacencyMatrix.size());
        stackOrder.reserve(m_AdjacencyMatrix.size());
    }

    // Initial coloring should be 1 to maxDegree, then the color the rest maxDegree+1.
    vector<int> &P(stackP[0]);
    vector<int> &vColors(stackColors[0]);
    vector<int> &vVertexOrder(stackOrder[0]);

    vector<int> &vCliqueVertices(R);
    vector<int> vOtherVerticesNotUsed;
    vector<pair<int,int>> vAddedEdgesNotUsed;
////    isolates.RemoveAllIsolates(0 /* unused*/, vCliqueVertices, vOtherVerticesNotUsed, vAddedEdgesNotUsed, true /* evaluate all vertices */);

    cliques.push_back(list<int>());

    InitializeOrder(P, vVertexOrder, vColors);

    if (R.size() < m_uMaximumCliqueSize) {
        cliques.back().clear();
        cliques.back().insert(cliques.back().end(), P.begin(), P.begin() + m_uMaximumCliqueSize);
        ExecuteCallBacks(cliques.back());
    }

#ifdef REMOVE_INITIAL_ISOLATES
    isolates.RemoveAllIsolates(0 /* unused*/, vCliqueVertices, vOtherVerticesNotUsed, vAddedEdgesNotUsed, true /* evaluate all vertices */);

    cout << "Initial R has " << R.size() << " elements" << endl;

    size_t uNewIndex(0);
    for (size_t index = 0; index < P.size(); ++index) {
        int const vertex(P[index]);
        if (isolates.GetInGraph().Contains(vertex)) {
            P[uNewIndex] = P[index];
            vColors[uNewIndex] = vColors[index];
            uNewIndex++;
        }
    }
    P.resize(uNewIndex);
    vColors.resize(uNewIndex);

    uNewIndex = 0;
    for (size_t index = 0; index < vVertexOrder.size(); ++index) {
        int const vertex(vVertexOrder[index]);
        if (isolates.GetInGraph().Contains(vertex)) {
            vVertexOrder[uNewIndex++] = vVertexOrder[index];
        }
    }
    vVertexOrder.resize(uNewIndex);

#endif //REMOVE_INITIAL_ISOLATES

    if (R.size() > m_uMaximumCliqueSize) {
        cliques.back().clear();
        cliques.back().insert(cliques.back().end(), R.begin(), R.end());
        ExecuteCallBacks(cliques.back());
    }

    RunRecursive(P, vVertexOrder, cliques, vColors);
    return cliques.size();
}

void LightWeightReductionSparseMISQ::Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors)
{
    coloringStrategy.Color(isolates.Neighbors(), vVertexOrder/* evaluation order */, vVerticesToReorder /* color order */, vColors);
}

void LightWeightReductionSparseMISQ::GetNewOrder(vector<int> &vNewVertexOrder, vector<int> &vVertexOrder, vector<int> const &P, int const chosenVertex, vector<int> &vCliqueVertices, vector<int> &vRemoved)
{
    vNewVertexOrder.resize(P.size());
    size_t uNewIndex(0);
    for (int const candidate : P) {
////    if (!m_bInvert && m_AdjacencyMatrix[chosenVertex][candidate]) vNewVertexOrder[uNewIndex++] = candidate;
        if (isolates.GetInGraph().Contains(candidate)) vNewVertexOrder[uNewIndex++] = candidate;
    }
    vNewVertexOrder.resize(uNewIndex);
}

void LightWeightReductionSparseMISQ::RunRecursive(vector<int> &P, vector<int> &vVertexOrder, list<list<int>> &cliques, vector<int> &vColors)
{
    int const recursionNode(nodeCount);
////    cout << "Node " << recursionNode << ":" << endl;
    nodeCount++;
    vector<int> &vNewP(stackP[R.size() + 1]);
    vector<int> &vNewColors(stackColors[R.size() + 1]);
    vector<int> vCliqueVerticesToReplace;
    vector<int> vRemovedVerticesToReplace;
    vector<pair<int,int>> vAddedEdgesUnused;

    while (!P.empty()) {
        int const largestColor(vColors.back());
        if (R.size() + largestColor <= m_uMaximumCliqueSize) {
            isolates.ReplaceAllRemoved(vCliqueVerticesToReplace);
            isolates.ReplaceAllRemoved(vRemovedVerticesToReplace);

            for (int const cliqueVertex : vCliqueVerticesToReplace) {
                R.pop_back();
            }
            return;
        }

        vColors.pop_back();
        int const nextVertex(P.back()); P.pop_back();
////        R.push_back(nextVertex); // gets pushed back with other clique vertices...
////        cout << "Adding " << nextVertex << " to clique" << endl;

        vector<int> &vNewVertexOrder(stackOrder[R.size()]);
        vector<int> vCliqueVertices; vCliqueVertices.push_back(nextVertex);
        vector<int> vRemoved;
        isolates.RemoveVertexAndNeighbors(nextVertex, vRemoved);
////        isolates.RemoveAllIsolates(0/*unused*/, vCliqueVertices, vRemoved, vAddedEdgesUnused /* unused */, false /* only consider updated vertices */);
        GetNewOrder(vNewVertexOrder, vVertexOrder, P, nextVertex, vCliqueVertices, vRemoved);
        R.insert(R.end(), vCliqueVertices.begin(), vCliqueVertices.end());

        if (!vNewVertexOrder.empty()) {
            vNewP.resize(vNewVertexOrder.size());
            vNewColors.resize(vNewVertexOrder.size());
            Color(vNewVertexOrder/* evaluation order */, vNewP /* color order */, vNewColors);
            RunRecursive(vNewP, vNewVertexOrder, cliques, vNewColors);
        } else if (R.size() > m_uMaximumCliqueSize) {
            cliques.back().clear();
            cliques.back().insert(cliques.back().end(), R.begin(), R.end());
            ExecuteCallBacks(cliques.back());
            m_uMaximumCliqueSize = R.size();
        }

        PostProcessOrder(vVertexOrder, nextVertex);

        // return R back to the state it was at the beginning of the loop.
        // remove reduction clique vertices from R
        for (int const vertex : vCliqueVertices)
            R.pop_back();

////        // one more for first vertex added to clique
////        R.pop_back();

        // put back removed vertices.
////        cout << "Putting back into graph: " << endl;
////        for (int const cliqueVertex : vCliqueVertices) {
////            cout << cliqueVertex << " ";
////        }
////        cout << " | ";
////        for (int const otherVertex : vRemoved) {
////            cout << otherVertex << " ";
////        }
////        cout << endl;
        isolates.ReplaceAllRemoved(vCliqueVertices);
        isolates.ReplaceAllRemoved(vRemoved);

        isolates.RemoveVertex(nextVertex);
        vRemovedVerticesToReplace.push_back(nextVertex);

        // remove vertices
        vector<int> vTempCliqueVertices;
        vector<int> vTempRemovedVertices;
        isolates.RemoveAllIsolates(0 /*unused*/,vTempCliqueVertices, vTempRemovedVertices, vAddedEdgesUnused, false /* only consider changed vertices */);

        R.insert(R.end(), vTempCliqueVertices.begin(), vTempCliqueVertices.end());

        vCliqueVerticesToReplace.insert(vCliqueVerticesToReplace.end(), vTempCliqueVertices.begin(), vTempCliqueVertices.end());
        vRemovedVerticesToReplace.insert(vRemovedVerticesToReplace.end(), vTempRemovedVertices.begin(), vTempRemovedVertices.end());

////        if (P.size() != isolates.GetInGraph().Size()) {
////            cout << "Node " << recursionNode << ": Consistancy Error" << endl;
////            cout << "P: " << endl;
////            for (int const vertexInP : P) {
////                cout << vertexInP << "(" << isolates.GetInGraph().Contains(vertexInP) << ") ";
////            }
////            cout << endl;
////        }

        if (P.size() != isolates.GetInGraph().Size()) { // if false?
////            cout << __LINE__ << ": This should not be triggered!" << endl;
            size_t uNewIndex(0);
            // pull vertices out of P and vColors
            for (size_t index = 0; index < P.size(); ++index) {
                if (isolates.GetInGraph().Contains(P[index])) {
                    P[uNewIndex] = P[index];
                    vColors[uNewIndex] = vColors[index];
                    uNewIndex++;
                }
            }
            P.resize(uNewIndex);
            vColors.resize(uNewIndex);

            // pull vertices out of vVertexOrder
            uNewIndex = 0;
            for (size_t index = 0; index < vVertexOrder.size(); ++index) {
                if (isolates.GetInGraph().Contains(vVertexOrder[index])) {
                    vVertexOrder[uNewIndex++] = vVertexOrder[index];
                }
            }
            vVertexOrder.resize(uNewIndex);
////            P.resize(uNewIndex);
////            vColors.resize(uNewIndex);
////            Color(vVertexOrder/* evaluation order */, P /* color order */, vColors);

            if (R.size() > m_uMaximumCliqueSize) {
                cliques.back().clear();
                cliques.back().insert(cliques.back().end(), R.begin(), R.end());
                ExecuteCallBacks(cliques.back());
                m_uMaximumCliqueSize = R.size();
            }

        }
    }


    isolates.ReplaceAllRemoved(vCliqueVerticesToReplace);
    isolates.ReplaceAllRemoved(vRemovedVerticesToReplace);

    for (int const cliqueVertex : vCliqueVerticesToReplace) {
        R.pop_back();
    }

    vNewColors.clear();
    vNewP.clear();
}


LightWeightReductionSparseMISQ::~LightWeightReductionSparseMISQ()
{
    cerr << "Largest Clique     : " << m_uMaximumCliqueSize << endl;
    cerr << "Node    Count      : " << nodeCount << endl;
}
