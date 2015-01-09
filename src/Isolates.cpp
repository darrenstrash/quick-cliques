#include "Isolates.h"

#include <vector>
#include <set>
#include <iostream>
#include <ctime>
#include <cassert>

using namespace std;

Isolates::Isolates(vector<vector<int>> &adjacencyArray)
 : m_AdjacencyArray(adjacencyArray)
 , neighbors(adjacencyArray.size())
 , inGraph(adjacencyArray.size())
 , isolates()
 , remaining(adjacencyArray.size())
 , vMarkedVertices(adjacencyArray.size(), false)
 , vvRemovedVertices()
 , m_AlternativeVertices()
 , timer(0)
 , removeTimer(0)
 , replaceTimer(0)
{
    for (size_t u=0; u < adjacencyArray.size(); ++u) {
        remaining.Insert(u);
        inGraph.Insert(u);
        neighbors[u].insert(m_AdjacencyArray[u].begin(), m_AdjacencyArray[u].end());
    }
}

Isolates::~Isolates()
{
    cout << "Total time spent computing next vertex: " << (timer/(double)CLOCKS_PER_SEC) << endl;
    cout << "Total time spent applying reductions  : " << (removeTimer/(double)CLOCKS_PER_SEC) << endl;
    cout << "Total time spent undoing  reductions  : " << (replaceTimer/(double)CLOCKS_PER_SEC) << endl;
}

void Isolates::RemoveEdges(vector<pair<int,int>> const &vEdges)
{
    for (pair<int,int> const &edge : vEdges) {
        neighbors[edge.first].erase(edge.second);
        neighbors[edge.second].erase(edge.first);

        m_AdjacencyArray[edge.first].pop_back();
        m_AdjacencyArray[edge.second].pop_back();
    }
}

bool Isolates::RemoveIsolatedClique(int const vertex, vector<int> &vIsolateVertices, vector<int> &vOtherRemovedVertices)
{
    size_t neighborCount(0);
    bool superSet(true);
    for (int const neighbor : neighbors[vertex]) {
        if (neighbors[neighbor].size() < neighbors[vertex].size()) {
            superSet = false; break;
        }
    }

////    bool debug(vertex == 113);
////
////    if (debug) {
////        cout << "vertex" << vertex << " has " << neighbors[vertex].size() << " neighbors." << endl << flush;
////    }

    if (!superSet) return false;

#if 1
////    bool const debug(true); //vertex==12);
////    if (debug) {
////        cout << "vertex " << vertex << " has neighbors: ";
////    }
////    if (debug) {
////        for (int const neighbor : neighbors[vertex]) {
////            cout << neighbor << " " << flush;
////        }
////////        cout << endl << flush;
////        cout << "(";
////
////        for (int const neighbor : m_AdjacencyArray[vertex]) {
////            cout << neighbor << " " << flush;
////        }
////        cout << ")" << endl << flush;
////    }

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
////        cout << "Removing Clique: [" << vertex;
////        for (int const neighbor : neighbors[vertex]) {
////            cout << " " << neighbor;
////        }
////        cout << "]" << endl;
////        for (int const neighbor : neighbors[vertex]) {
////            cout << neighbor << " ";
////        }
////        cout << endl;

        for (int const neighbor : neighbors[vertex]) {
            inGraph.Remove(neighbor);
            remaining.Remove(neighbor);
        }
        inGraph.Remove(vertex);
        isolates.insert(vertex);
        vIsolateVertices.push_back(vertex);
////        if (vertex == 0) {
////            cout << __LINE__ << ": Removing " << vertex << " from graph." << endl << flush;
////        }
        vOtherRemovedVertices.insert(vOtherRemovedVertices.end(), neighbors[vertex].begin(), neighbors[vertex].end());

        for (int const neighbor : neighbors[vertex]) {
            ////                cout << "   Expunging neighbor " << neighbor << " with " << neighbors[neighbor].size() << " neighbors" << endl << flush;
            for (int const nNeighbor : neighbors[neighbor]) {
                if (inGraph.Contains(nNeighbor)) {
                    remaining.Insert(nNeighbor);
                }

                if (nNeighbor != vertex) {
                    neighbors[nNeighbor].erase(neighbor);
                }
            }
            neighbors[neighbor].clear();
        }
        neighbors[vertex].clear();

////    if (vertex == 21952) {
////        cout << "vertex" << vertex << " is being removed." << endl << flush;
////    }
        
        return true;
    }
    return false;
}

// TODO/DS: need to remember added edge, so we can remove it later.
// TODO/DS: not currently working, proceeding without it.
bool Isolates::RemoveIsolatedPath(int const vertex, vector<int> &vIsolateVertices, vector<int> &vOtherRemovedVertices, vector<pair<int,int>> &vAddedEdges)
{
    if (neighbors[vertex].size() != 2) return false;

    for (int const neighbor : neighbors[vertex]) {
        if (neighbors[neighbor].size() == 2) {
            isolates.insert(vertex);

            int endVertex1(-1);
            int endVertex2(-1);

            // remove from other neighbor's list...
            for (int const nNeighbor : neighbors[neighbor]) {
                if (nNeighbor != vertex) {
                    endVertex1 = nNeighbor;
                    neighbors[nNeighbor].erase(neighbor);
////                    cout << __LINE__ << ": Removing edge " << nNeighbor << "," << neighbor << endl;
                    remaining.Insert(nNeighbor);
                    break;
                }
            }

            for (int const otherNeighbor : neighbors[vertex]) {
                if (otherNeighbor != neighbor) {
                    endVertex2 = otherNeighbor;
                    neighbors[otherNeighbor].erase(vertex);
////                    cout << __LINE__ << ": Removing edge " << otherNeighbor << "," << vertex << endl;
                    remaining.Insert(otherNeighbor);

                    if (neighbors[otherNeighbor].find(endVertex1) == neighbors[otherNeighbor].end()) {

                        m_AdjacencyArray[otherNeighbor].push_back(endVertex1);
                        m_AdjacencyArray[endVertex1].push_back(otherNeighbor);

                        neighbors[otherNeighbor].insert(endVertex1);
                        neighbors[endVertex1].insert(otherNeighbor);
                        vAddedEdges.push_back(make_pair(endVertex1, otherNeighbor));
                        break;
                    }
                }
            }

            neighbors[vertex].clear();
            neighbors[neighbor].clear();

            vIsolateVertices.push_back(vertex);
////        if (vertex == 0) {
////            cout << __LINE__ << ": Removing " << vertex << " from graph." << endl << flush;
////        }
            vOtherRemovedVertices.push_back(neighbor);
////        if (neighbor == 0) {
////            cout << __LINE__ << ": Removing " << neighbor << " from graph." << endl << flush;
////        }

            remaining.Remove(vertex);
            remaining.Remove(neighbor);

            m_AlternativeVertices[vertex]   = neighbor;
            m_AlternativeVertices[neighbor] = vertex;
            ////cout << "Vertex " << vertex << " has alternative " << neighbor << endl;

            inGraph.Remove(vertex);
            inGraph.Remove(neighbor);
////            cout << "Removing Path: [" << vertex << " " << neighbor << "]" << " (" << endVertex1 << " " << endVertex2 << ")" << endl << flush;

            return true;
        }
    }
    return false;
}

void Isolates::RemoveVertex(int const vertex)
{
////    cout << __LINE__ << ": Removing vertex " << vertex << endl << flush;
    inGraph.Remove(vertex);
    remaining.Remove(vertex);
    isolates.erase(vertex);

    for (int const neighbor : neighbors[vertex]) {
        neighbors[neighbor].erase(vertex);
    }

    neighbors[vertex].clear();
}


void Isolates::RemoveVertexAndNeighbors(int const vertex, vector<int> &vRemoved)
{
////    cout << __LINE__ << ": Removing vertex " << vertex << endl << flush;
    inGraph.Remove(vertex);
    remaining.Remove(vertex);
    isolates.insert(vertex);
    vRemoved.push_back(vertex);
////    if (vertex == 0) {
////        cout << __LINE__ << ": Removing " << vertex << " from graph." << endl << flush;
////    }
    inGraph.Remove(vertex);

    for (int const neighbor : neighbors[vertex]) {
////        cout << __LINE__ << ":     Removing neighbor " << neighbor << endl << flush;
        inGraph.Remove(neighbor);
        remaining.Remove(neighbor);
        vRemoved.push_back(neighbor);
////        if (neighbor == 0) {
////            cout << __LINE__ << ": Removing " << neighbor << " from graph." << endl << flush;
////        }

        for (int const nNeighbor : neighbors[neighbor]) {
////            cout << __LINE__ << ":         Removing from neighbor's neighbor "<< nNeighbor << endl << flush;
            if (nNeighbor == vertex) continue;
            neighbors[nNeighbor].erase(neighbor);
////            cout << __LINE__ << ":         Done removing" << endl << flush;
            if (inGraph.Contains(nNeighbor))
                remaining.Insert(nNeighbor);
        }
////        cout << __LINE__ << ":     Done Removing neighbor" << neighbor << endl << flush;
        neighbors[neighbor].clear();
////        cout << __LINE__ << ": Cleared neighbors: " << neighbor << endl << flush;
    }

////    cout << __LINE__ << ": Done Removing vertex " << vertex << endl << flush;
    neighbors[vertex].clear();
////    cout << __LINE__ << ": Cleared neighbors: " << vertex << endl << flush;
}

void Isolates::RemoveAllIsolates(int const independentSetSize, vector<int> &vIsolateVertices, vector<int> &vOtherRemovedVertices, vector<pair<int,int>> &vAddedEdges)
{
    clock_t startClock = clock();
////    remaining = inGraph; // TODO/DS : We can optimize this by knowing which vertex (and neighbors where removed last.
////    if (vOtherRemovedVertices.empty()) {
        remaining.Clear();
        for (int const vertex : inGraph) {
            remaining.Insert(vertex);
        }
////    } else {
////        remaining.clear();
////        for (int const removedVertex : vOtherRemovedVertices) {
////            remaining.insert(neighbors[removedVertex].begin(), neighbors[removedVertex].end());
////        }
////    }
    //cout << "Removing all isolates." << endl << flush;
    int isolateSize(isolates.size());
    while (!remaining.Empty()) {
////            if (remaining.size() %10000 == 0)
        //cout << "Remaining: " << remaining.size() << ", Removed: " << removed.size() << endl << flush;
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

        bool reduction = RemoveIsolatedClique(vertex, vIsolateVertices, vOtherRemovedVertices);
////        if (!reduction) {
////            reduction = RemoveIsolatedPath(vertex, vIsolateVertices, vOtherRemovedVertices, vAddedEdges);
////        }
    }

////    cout << "Removed " << isolates.size() - isolateSize << " isolates." << endl << flush;

    clock_t endClock = clock();

    removeTimer += (endClock - startClock);
}

void Isolates::ReplaceAllRemoved(vector<int> const &vRemoved)
{
    clock_t startClock = clock();
    for (int const removedVertex : vRemoved) {
        inGraph.Insert(removedVertex);
        isolates.erase(removedVertex);
        m_AlternativeVertices.erase(removedVertex);
    }
////    cout << "Replacing all removed vertices." << endl << flush;
    for (int const removedVertex : vRemoved) {
////        if (remaining.size() %10000 == 0)
////            cout << "Remaining: " << remaining.size() << ", Removed: " << removed.size() << endl << flush;
        for (int const neighbor : m_AdjacencyArray[removedVertex]) {
            if (!inGraph.Contains(neighbor)) continue;
            neighbors[removedVertex].insert(neighbor);
            neighbors[neighbor].insert(removedVertex);
        }
    }
////    cout << "Done replacing vertices..." << endl << flush;

    clock_t endClock = clock();
    replaceTimer += (endClock - startClock);
}

int Isolates::NextVertexToRemove(std::vector<int> &vVertices)
{
    clock_t startClock = clock();

    remaining.Clear();
    for (int const vertex : inGraph) {
        remaining.Insert(vertex); // only consider neighbors of removed vertices... makes it faster.
    }

    vector<int> &vVerticesOrderedByDegree(vVertices);

    auto sortByDegree = [this](int const &leftVertex, int const &rightVertex) { return neighbors[leftVertex].size() > neighbors[rightVertex].size(); };
    sort(vVerticesOrderedByDegree.begin(), vVerticesOrderedByDegree.end(), sortByDegree);

    if (vVerticesOrderedByDegree.empty()) return -1;

    int vertexWithMaxReductions(vVerticesOrderedByDegree[0]); // default to max degree vertex
    int maxIsolates(-1);
    int maxRemoved(-1);
    int maxIndex(-1);
    size_t index(0); // default to index of maximum degree vertex

////    cout << "Testing remaining " << inGraph.size() << " vertices, to maximize isolate removal" << endl << flush;
    for (int const vertex : vVerticesOrderedByDegree) {
////        cout << "Starting loop..." << endl << flush;
////        set<int> tempRemaining(remaining); // TODO/DS: put back?

        // TODO/DS: Remove when vertex selection becomes faster, this is a compromise between fast vertex selection and good vertex selection.
        if (index >= vVerticesOrderedByDegree.size()/10) break;

        remaining.Clear();

        inGraph.SaveState();

        set<int> tempIsolates(isolates);

        vector<int> vRemoved;
////        cout << "Try removing vertex " << vertex << endl << flush;
#ifdef OPTIMIZE_BRANCHING
        RemoveVertex(vertex); vRemoved.push_back(vertex);
#else
        RemoveVertexAndNeighbors(vertex, vRemoved);
#endif // OPTIMIZE_BRANCHING
////        cout << "And the new isolates..." << vertex << endl << flush;
        vector<pair<int,int>> vAddedEdges;

        clock_t startremove(clock());
        RemoveAllIsolates(0, vRemoved, vRemoved, vAddedEdges);
        removeTimer -= (clock() - startremove);

        int const newNumRemoved(m_AdjacencyArray.size() - inGraph.Size());

        inGraph.RestoreState();
        int const oldNumRemoved(m_AdjacencyArray.size() - inGraph.Size());
        int const deltaRemoved(newNumRemoved - oldNumRemoved);
#if 0 //def OPTIMIZE_ISOLATES
        if (static_cast<int>(isolates.size() - tempIsolates.size()) > maxIsolates) {
            maxIsolates = static_cast<int>(isolates.size() - tempIsolates.size());
            vertexWithMaxReductions = vertex;
        }
#else // optimize removed vertices
        if (deltaRemoved > maxRemoved) {
            maxRemoved = deltaRemoved;
            vertexWithMaxReductions = vertex;
            maxIndex = index;
        }
#endif

////        remaining = std::move(tempRemaining); // TODO/DS: put back?

        isolates  = std::move(tempIsolates);

        // need to do this after restoring variable inGraph
        RemoveEdges(vAddedEdges); // TODO/DS: Put back...
        clock_t startreplace(clock());
        ReplaceAllRemoved(vRemoved);
        replaceTimer -= (clock() - startreplace);

        index++;

////        cout << "Done with loop" << endl;
    }


    clock_t endClock = clock();

    timer += (endClock - startClock);

    // inefficient, should be a better way to elmininate it before getting here.
    for (int &vertex : vVertices) {
        if (vertex == vertexWithMaxReductions) {
            vertex = vVertices.back();
            vVertices.pop_back();
            return vertexWithMaxReductions;
        }
    }

    assert(0);
    return vertexWithMaxReductions;
}

int Isolates::NextVertexToRemove()
{
    vector<int> vVertices;
    vVertices.insert(vVertices.end(), inGraph.begin(), inGraph.end());
    return NextVertexToRemove(vVertices);
}

int Isolates::GetAlternativeVertex(int const vertex) const
{
    map<int,int>::const_iterator cit(m_AlternativeVertices.find(vertex));
    if (cit != m_AlternativeVertices.end()) {
        return cit->second;
    }
    return -1;
}

