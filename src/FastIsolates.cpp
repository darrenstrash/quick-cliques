#include "FastIsolates.h"
#include "ArraySet.h"
#include "SparseArraySet.h"

#include <vector>
#include <set>
#include <iostream>
#include <ctime>
#include <cassert>
#include <climits>

using namespace std;

#define NEW_ISOLATE_CHECKS
////#define DOMINATION_CHECKS
#define NON_PATH

template <typename NeighborSet>
FastIsolates<NeighborSet>::FastIsolates(vector<vector<int>> const &adjacencyArray)
 : m_AdjacencyArray(adjacencyArray)
 , neighbors(adjacencyArray.size())
 , m_vbInGraph(adjacencyArray.size(), true)
 , remaining(adjacencyArray.size())
 , vMarkedVertices(adjacencyArray.size(), false)
#ifdef TIMERS
 , timer(0)
 , removeTimer(0)
 , replaceTimer(0)
 , sortDuringNextTimer(0)
 , removeOneDuringNextTimer(0)
 , removeDuringNextTimer(0)
 , replaceDuringNextTimer(0)
 #endif // TIMERS
 , m_bConnectedComponentMode(false)
 , foldedVertexCount(0)
 , m_bAllowVertexFolds(true)
 , m_uReductionCount(0)
 , m_uRemainingGraphSize(adjacencyArray.size())
{
    for (size_t u=0; u < adjacencyArray.size(); ++u) {
        remaining.Insert(u);
#ifdef SPARSE
        neighbors[u].Resize(m_AdjacencyArray[u].size());
#else
        neighbors[u].InitializeFromAdjacencyArray(m_AdjacencyArray, u);
#endif // SPARSE
        for (int const vertex : m_AdjacencyArray[u]) {
            neighbors[u].Insert(vertex);
        }
    }
}

template <typename NeighborSet>
FastIsolates<NeighborSet>::~FastIsolates()
{

#ifdef TIMERS
    cout << "Total time spent computing next vertex: " << (timer/(double)CLOCKS_PER_SEC) << endl;
    cout << "    Sort      : " << (sortDuringNextTimer/(double)CLOCKS_PER_SEC) << endl;
    cout << "    RemoveOne : " << (removeOneDuringNextTimer/(double)CLOCKS_PER_SEC) << endl;
    cout << "    RemoveRest: " << (removeDuringNextTimer/(double)CLOCKS_PER_SEC) << endl;
    cout << "    ReplaceAll: " << (replaceDuringNextTimer/(double)CLOCKS_PER_SEC) << endl;
   
    cout << "Total time spent applying reductions  : " << (removeTimer/(double)CLOCKS_PER_SEC) << endl;
    cout << "Total time spent undoing  reductions  : " << (replaceTimer/(double)CLOCKS_PER_SEC) << endl;
#endif // TIMERS
}

template <typename NeighborSet>
bool FastIsolates<NeighborSet>::RemoveIsolatedClique(int const vertex, vector<int> &vIsolateVertices, vector<int> &vOtherRemovedVertices, vector<Reduction> &vReductions)
{
////    bool const debug(vertex==40653);
////    if (vertex == 31) {
////        cout << "Removing isolated clique with vertex " << vertex << endl << flush;
////        if (!inGraph.Contains(31)) cout << "Vertex 31 is not in the graph!" << endl << flush;
////    }

    size_t neighborCount(0);
    for (int const neighbor : neighbors[vertex]) {
        if (neighbors[neighbor].Size() < neighbors[vertex].Size()) {
            return false;
        }
    }

    bool superSet(true);

////    bool debug(vertex == 113);
////
////    if (debug) {
////        cout << "vertex" << vertex << " has " << neighbors[vertex].size() << " neighbors." << endl << flush;
////    }

#if 1
////    bool const debug(true); //vertex==12);
////    if (debug) {
////        cout << "vertex " << vertex << " has neighbors: ";
////    }
////    if (debug) {
////        for (int const neighbor : neighbors[vertex]) {
////            cout << neighbor << " " << flush;
////        }
////        cout << endl << flush;
////        cout << "(";
////
////        for (int const neighbor : m_AdjacencyArray[vertex]) {
////            cout << neighbor << " " << flush;
////        }
////        cout << ")" << endl << flush;
////    }
////
    for (int const neighbor : neighbors[vertex]) {
////        if (debug) {
////            cout << "Considering neighbor " <<  neighbor  << flush;
////        }
        neighborCount++;

        for (int const nNeighbor : neighbors[neighbor]) {
            vMarkedVertices[nNeighbor] = true;
        }
        vMarkedVertices[neighbor] = true;

        for (int const neighbor2 : neighbors[vertex]) {
////            if (debug && superSet && !vMarkedVertices[neighbor2]) {
////                cout << "(missing neighbor " << neighbor2 << ") ";
////            }
            superSet = superSet && vMarkedVertices[neighbor2];
        }

        //superSet = superSet && (neighborCount == neighbors[vertex].size()); // TODO/DS : put back?

        for (int const nNeighbor : neighbors[neighbor]) {
            vMarkedVertices[nNeighbor] = false;
        }
        vMarkedVertices[neighbor] = false;

        if (!superSet) {
////            if (debug) {
////                cout << "[x] " << endl;
////            }
            return false;
        }

////        if (debug) cout << endl << flush;

        ////            cout << "Iteration(" << neighborCount << ": Intersection=";
        ////            for (int const vertex : intersection) {
        ////                cout << vertex << " ";
        ////            }
        ////            cout << endl << flush;
    }
#else

    for (int const neighbor : neighbors[vertex]) {
////        cout << "Evaluating neighbor: " << neighbor << endl << flush;
        neighborCount++;

        for (int const nNeighbor : neighbors[neighbor]) {
////            cout << "Marking nNeighbor: " << nNeighbor << endl << flush;
            vMarkedVertices[nNeighbor] = true;
        }
////        cout << "Marking neighbor: " << neighbor << endl << flush;
        vMarkedVertices[neighbor] = true;

        for (int const neighbor2 : neighbors[vertex]) {
////            cout << "Checking neighbor: " << neighbor2 << endl << flush;
            superSet = superSet && vMarkedVertices[neighbor2];
        }

        for (int const nNeighbor : neighbors[neighbor]) {
////            cout << "Unmarking nNeighbor: " << nNeighbor << endl << flush;
            vMarkedVertices[nNeighbor] = false;
        }
////        cout << "Unmarking neighbor: " << neighbor << endl << flush;
        vMarkedVertices[neighbor] = false;

        if (!superSet) return false;

        ////            cout << "Iteration(" << neighborCount << ": Intersection=";
        ////            for (int const vertex : intersection) {
        ////                cout << vertex << " ";
        ////            }
        ////            cout << endl << flush;
    }
#endif

////    cout << "Done evaluating neighbors" << endl << flush;

    ////        cout << "Iteration(" << neighborCount << ": Intersection=";
    ////        for (int const vertex : intersection) {
    ////            cout << vertex << " ";
    ////        }
    ////        cout << endl << flush;

    if (superSet) {
////        if (debug) {
////            cout << "Removing Clique: [" << vertex;
////            for (int const neighbor : neighbors[vertex]) {
////                cout << " " << neighbor;
////            }
////            cout << "]" << endl;
////        }
////        for (int const neighbor : neighbors[vertex]) {
////            cout << neighbor << " ";
////        }
////        cout << endl;

#ifdef SLOW
        m_uReductionCount++;
        Reduction reduction(ISOLATED_VERTEX);
        reduction.SetVertex(vertex);
#endif // SLOW

////        cout << "Removing isolated vertex " << vertex << " with neighbors ";
        for (int const neighbor : neighbors[vertex]) {
////            if (!inGraph.Contains(neighbor)) {
////                cout << "Trying to remove non-existant neighbor " << neighbor << " from " << vertex << " neighbor list!" << endl << flush;
////            }
////            cout << __LINE__ << ": calling remove" << endl << flush;
            m_vbInGraph[neighbor] = false;
            m_uRemainingGraphSize--;
            remaining.Remove(neighbor);
#ifdef SLOW
            reduction.AddNeighbor(neighbor);
            reduction.AddRemovedEdge(vertex,   neighbor);
            reduction.AddRemovedEdge(neighbor, vertex);
#endif // SLOW
////            cout << neighbor << " ";
        }
////        cout << endl << flush;

////        if (!inGraph.Contains(vertex)) {
////            cout << "Trying to remove non-existant vertex " << vertex << " from graph!" << endl << flush;
////        }

////        cout << __LINE__ << ": calling remove" << endl << flush;
        m_vbInGraph[vertex] = false;
        m_uRemainingGraphSize--;
        vIsolateVertices.push_back(vertex);
////        if (vertex == 0) {
////            cout << __LINE__ << ": Removing " << vertex << " from graph." << endl << flush;
////        }
        vOtherRemovedVertices.insert(vOtherRemovedVertices.end(), neighbors[vertex].begin(), neighbors[vertex].end());

        for (int const neighbor : neighbors[vertex]) {
            ////                cout << "   Expunging neighbor " << neighbor << " with " << neighbors[neighbor].size() << " neighbors" << endl << flush;
            for (int const nNeighbor : neighbors[neighbor]) {
                if (m_vbInGraph[nNeighbor]) {
                    remaining.Insert(nNeighbor);
                }

                if (nNeighbor != vertex) {
                    neighbors[nNeighbor].Remove(neighbor);
#ifdef SLOW
                    reduction.AddRemovedEdge(nNeighbor, neighbor);
                    reduction.AddRemovedEdge(neighbor, nNeighbor);
#endif
                }
            }
            neighbors[neighbor].Clear();
        }
        neighbors[vertex].Clear();

////    if (vertex == 21952) {
////        cout << "vertex" << vertex << " is being removed." << endl << flush;
////    }

#ifdef SLOW
        vReductions.emplace_back(std::move(reduction));
#endif // SLOW
        
        return true;
    }
    return false;
}

#ifdef NEW_ISOLATE_CHECKS
template <typename NeighborSet>
bool FastIsolates<NeighborSet>::FoldVertex(int const vertex, vector<int> &vIsolateVertices, vector<int> &vOtherRemovedVertices, vector<Reduction> &vReductions)
{
#ifdef NON_PATH
    if (neighbors[vertex].Size() != 2) return false;
    if (neighbors[neighbors[vertex][0]].Contains(neighbors[vertex][1])) return false; // neighbors can't be adjacent.

    foldedVertexCount++;

////    cout << "Folding vertex " << vertex << ":" << endl << flush;

    int const vertex1(neighbors[vertex][0]);
    int const vertex2(neighbors[vertex][1]);

#ifdef SLOW
    Reduction reduction(FOLDED_VERTEX);
    reduction.SetVertex(vertex);
    reduction.AddNeighbor(vertex1);
    reduction.AddNeighbor(vertex2);
#else
    m_uReductionCount++;
#endif // SLOW


////    bool const debug(vertex == 1480);
////    if (debug) {
////        cout << "Folding: " << vertex << ", " << vertex1 << ", " << vertex2 << endl;
////        cout << "BEFORE:" << endl;
////        cout << vertex << ":";
////        for (int const neighbor : neighbors[vertex]) {
////            cout << " " << neighbor;
////        }
////        cout << endl;
////
////        cout << vertex1 << ":";
////        for (int const neighbor : neighbors[vertex1]) {
////            cout << " " << neighbor;
////        }
////        cout << endl;
////
////        cout << vertex2 << ":";
////        for (int const neighbor : neighbors[vertex2]) {
////            cout << " " << neighbor;
////        }
////        cout << endl;
////    }

    neighbors[vertex].Clear();
    neighbors[vertex].Resize(neighbors[vertex1].Size() + neighbors[vertex2].Size());
    neighbors[vertex1].Remove(vertex);
    neighbors[vertex2].Remove(vertex);

#ifdef SLOW
    reduction.AddRemovedEdge(vertex, vertex1);
    reduction.AddRemovedEdge(vertex1, vertex);
    reduction.AddRemovedEdge(vertex, vertex2);
    reduction.AddRemovedEdge(vertex2, vertex);
#endif // SLOW

    for (int const neighbor1 : neighbors[vertex1]) {
        if (neighbor1 == vertex) continue;
        neighbors[neighbor1].Remove(vertex1);
#ifdef SLOW
        reduction.AddRemovedEdge(neighbor1, vertex1);
        reduction.AddRemovedEdge(vertex1, neighbor1);
#endif // SLOW
////        neighbors[neighbor1].Insert(vertex);
        neighbors[vertex].Insert(neighbor1);
        remaining.Insert(neighbor1);
    }
    neighbors[vertex1].Clear();

////    cout << "vertex " << vertex << " contains self?=" << (neighbors[vertex].Contains(vertex)) << endl;

    for (int const neighbor2 : neighbors[vertex2]) {
        if (neighbor2 == vertex) continue;
        neighbors[neighbor2].Remove(vertex2);
#ifdef SLOW
        reduction.AddRemovedEdge(neighbor2, vertex2);
        reduction.AddRemovedEdge(vertex2, neighbor2);
#endif // SLOW
////        neighbors[neighbor2].Insert(vertex);
        neighbors[vertex].Insert(neighbor2);
        remaining.Insert(neighbor2);
    }

////    cout << "vertex " << vertex << " contains self?=" << (neighbors[vertex].Contains(vertex)) << endl;

    neighbors[vertex2].Clear();
    for (int const neighbor : neighbors[vertex]) {
////        if (neighbor == vertex) {
////            cout << "Something is very wrong..." << endl;
////        }
        neighbors[neighbor].Insert(vertex);
    }

////    cout << "vertex " << vertex << " contains self?=" << (neighbors[vertex].Contains(vertex)) << endl;

    remaining.Insert(vertex);

#ifdef SLOW
    vReductions.emplace_back(std::move(reduction));
#endif // SLOW

    // which one goes in independent set? Need to update...
    vOtherRemovedVertices.push_back(vertex1);
    vOtherRemovedVertices.push_back(vertex2);

    remaining.Remove(vertex1);
    remaining.Remove(vertex2);
    m_vbInGraph[vertex1] = false;
    m_vbInGraph[vertex2] = false;
    m_uRemainingGraphSize--;
    m_uRemainingGraphSize--;

////    if (debug) {
////        cout << "AFTER:" << endl;
////        cout << vertex << ":";
////        for (int const neighbor : neighbors[vertex]) {
////            cout << " " << neighbor;
////        }
////        cout << endl;
////
////        cout << vertex1 << ":";
////        for (int const neighbor : neighbors[vertex1]) {
////            cout << " " << neighbor;
////        }
////        cout << endl;
////
////        cout << vertex2 << ":";
////        for (int const neighbor : neighbors[vertex2]) {
////            cout << " " << neighbor;
////        }
////        cout << endl;
////    }
////    cout << "Done folding..." << endl << flush;
    return true;
#else
    return false;
#endif 
}

#endif // NEW_ISOLATE_CHECKS

template <typename NeighborSet>
void FastIsolates<NeighborSet>::RemoveAllIsolates(int const independentSetSize, vector<int> &vIsolateVertices, vector<int> &vOtherRemovedVertices, vector<Reduction> &vReductions, bool bConsiderAllVertices)
{
////    if (find (vIsolateVertices.begin(), vIsolateVertices.end(), 31) != vIsolateVertices.end())
////        cout << "Calling RemoveAllIsolates with 31 in the isolate set!" << endl;
#ifdef TIMERS
    clock_t startClock = clock();
#endif // TIMERS
////    remaining = inGraph; // TODO/DS : We can optimize this by knowing which vertex (and neighbors where removed last.
////    if (vOtherRemovedVertices.empty()) {

        // TODO/DS: Put this in; it saves us from having to consider obvious non-candidates. Only works if we establish
        // the invariant that the graph contains no vertices that can be reduced.
////        if (bConsiderAllVertices) { ////true) { //bConsiderAllVertices) {
////        if (true) { //bConsiderAllVertices) {
////            remaining.Clear();
////            for (int const vertex : inGraph) {
////                remaining.Insert(vertex);
////            }
////        }
////    } else {
////        remaining.clear();
////        for (int const removedVertex : vOtherRemovedVertices) {
////            remaining.insert(neighbors[removedVertex].begin(), neighbors[removedVertex].end());
////        }
////    }
////    cout << "Removing all isolates." << endl << flush;
    int iterations(0);
    while (!remaining.Empty()) {
////        size_t const numRemoved(vIsolateVertices.size() + vOtherRemovedVertices.size());
////        if ((vIsolateVertices.size() + vOtherRemovedVertices.size()) %10000 == 0)
////        cout << "Progress: Removed: " << vIsolateVertices.size() << " isolates, and " << vOtherRemovedVertices.size() << " others" << endl << flush;
        int const vertex = *(remaining.begin());
        remaining.Remove(vertex);

    // can prune out high-degree vertices too. // but this is slow right now.

////        if (inGraph.size() - neighbors[vertex].size() < independentSetSize) {
////            remaining.insert(neighbors[vertex].begin(), neighbors[vertex].end());
////            RemoveVertex(vertex);
////            vOtherRemovedVertices.push_back(vertex);
////            continue;
////        }

////        cout << "Attempting to remove vertex " << vertex << endl << flush;

        bool reduction = RemoveIsolatedClique(vertex, vIsolateVertices, vOtherRemovedVertices, vReductions);
#ifdef NEW_ISOLATE_CHECKS
#if 1 ////def FOLD_VERTEX
        if (!reduction) {
////            reduction = RemoveIsolatedPath(vertex, vIsolateVertices, vOtherRemovedVertices, vAddedEdges);
            reduction = FoldVertex(vertex, vIsolateVertices, vOtherRemovedVertices, vReductions);
        }
#endif // 0
#ifdef DOMINATION_CHECKS
        if (!reduction) {
            reduction = RemoveDominatedVertex(vertex, vIsolateVertices, vOtherRemovedVertices, vReductions);
        }
#endif //DOMINATION_CHECKS

#endif // NEW_ISOLATE_CHECKS

////    if (find (vIsolateVertices.begin(), vIsolateVertices.end(), 31) != vIsolateVertices.end())
////        cout << "31 was added to the isolate set!" << endl;
////        if (!inGraph.Contains(31)) cout << "And it's not in the graph..." << endl << flush;

        iterations++;

////        size_t const numNewRemoved(vIsolateVertices.size() + vOtherRemovedVertices.size());
////        if (numNewRemoved != numRemoved) {
////            cout << "Progress: Removed: " << vIsolateVertices.size() << " isolates, and " << vOtherRemovedVertices.size() << " others, in " << iterations << " iterations." << endl << flush;
////            iterations = 0;
////            cout << "Remaining graph has " << inGraph.Size() << " vertices." << endl << flush;
////        }
    }

////    cout << "Removed " << isolates.size() - isolateSize << " isolates." << endl << flush;
////    cout << "Removed: " << vIsolateVertices.size() << " isolates, and " << vOtherRemovedVertices.size() << " others, in " << iterations << " iterations." << endl << flush;

#ifdef TIMERS
    clock_t endClock = clock();
    removeTimer += (endClock - startClock);
#endif // TIMERS
}

////template class FastIsolates<ArraySet>;
template class FastIsolates<SparseArraySet>;

