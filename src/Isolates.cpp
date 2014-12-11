#include "Isolates.h"

#include <vector>
#include <set>
#include <iostream>

using namespace std;

Isolates::Isolates(vector<vector<int>> &adjacencyArray)
 : m_AdjacencyArray(adjacencyArray)
 , neighbors(adjacencyArray.size())
 , inGraph()
 , isolates()
 , removed()
 , remaining()
 , vMarkedVertices(adjacencyArray.size(), false)
 , vvRemovedVertices()
 , m_AlternativeVertices()
{
    for (size_t u=0; u < adjacencyArray.size(); ++u) {
        remaining.insert(u);
        inGraph.insert(u);
        neighbors[u].insert(m_AdjacencyArray[u].begin(), m_AdjacencyArray[u].end());
    }
}

Isolates::~Isolates()
{
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
        
        removed.insert(neighbors[vertex].begin(), neighbors[vertex].end());
        removed.insert(vertex);
        isolates.insert(vertex);
        vIsolateVertices.push_back(vertex);
        inGraph.erase(vertex);
////        if (vertex == 0) {
////            cout << __LINE__ << ": Removing " << vertex << " from graph." << endl << flush;
////        }
        vOtherRemovedVertices.insert(vOtherRemovedVertices.end(), neighbors[vertex].begin(), neighbors[vertex].end());

        for (int const neighbor : neighbors[vertex]) {
            remaining.erase(neighbor);
            inGraph.erase(neighbor);
            ////                cout << "   Expunging neighbor " << neighbor << " with " << neighbors[neighbor].size() << " neighbors" << endl << flush;
            for (int const nNeighbor : neighbors[neighbor]) {
                if (removed.find(nNeighbor) == removed.end()) {
                    remaining.insert(nNeighbor);
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
bool Isolates::RemoveIsolatedPath(int const vertex, vector<int> &vIsolateVertices, vector<int> &vOtherRemovedVertices, vector<pair<int,int>> &vAddedEdges)
{
    if (neighbors[vertex].size() != 2) return false;

    for (int const neighbor : neighbors[vertex]) {
        if (neighbors[neighbor].size() == 2) {
            removed.insert(vertex);
            removed.insert(neighbor);
            isolates.insert(vertex);

            int endVertex1(-1);
            int endVertex2(-1);

            // remove from other neighbor's list...
            for (int const nNeighbor : neighbors[neighbor]) {
                if (nNeighbor != vertex) {
                    endVertex1 = nNeighbor;
                    neighbors[nNeighbor].erase(neighbor);
////                    cout << __LINE__ << ": Removing edge " << nNeighbor << "," << neighbor << endl;
                    remaining.insert(nNeighbor);
                    break;
                }
            }

            for (int const otherNeighbor : neighbors[vertex]) {
                if (otherNeighbor != neighbor) {
                    endVertex2 = otherNeighbor;
                    neighbors[otherNeighbor].erase(vertex);
////                    cout << __LINE__ << ": Removing edge " << otherNeighbor << "," << vertex << endl;
                    remaining.insert(otherNeighbor);

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

            remaining.erase(vertex);
            remaining.erase(neighbor);

            m_AlternativeVertices[vertex]   = neighbor;
            m_AlternativeVertices[neighbor] = vertex;
            ////cout << "Vertex " << vertex << " has alternative " << neighbor << endl;

            inGraph.erase(vertex);
            inGraph.erase(neighbor);
////            cout << "Removing Path: [" << vertex << " " << neighbor << "]" << " (" << endVertex1 << " " << endVertex2 << ")" << endl << flush;

            return true;
        }
    }
    return false;
}

void Isolates::RemoveVertex(int const vertex)
{
////    cout << __LINE__ << ": Removing vertex " << vertex << endl << flush;
    removed.insert(vertex);
    remaining.erase(vertex);
    inGraph.erase(vertex);
    isolates.erase(vertex);

    for (int const neighbor : neighbors[vertex]) {
        neighbors[neighbor].erase(vertex);
    }

    neighbors[vertex].clear();
}


void Isolates::RemoveVertexAndNeighbors(int const vertex, vector<int> &vRemoved)
{
////    cout << __LINE__ << ": Removing vertex " << vertex << endl << flush;
    removed.insert(vertex);
    remaining.erase(vertex);
    isolates.insert(vertex);
    vRemoved.push_back(vertex);
    inGraph.erase(vertex);
////    if (vertex == 0) {
////        cout << __LINE__ << ": Removing " << vertex << " from graph." << endl << flush;
////    }
    inGraph.erase(vertex);

    for (int const neighbor : neighbors[vertex]) {
////        cout << __LINE__ << ":     Removing neighbor " << neighbor << endl << flush;
        removed.insert(neighbor);
        remaining.erase(neighbor);
        vRemoved.push_back(neighbor);
        inGraph.erase(neighbor);
////        if (neighbor == 0) {
////            cout << __LINE__ << ": Removing " << neighbor << " from graph." << endl << flush;
////        }
        inGraph.erase(neighbor);

        for (int const nNeighbor : neighbors[neighbor]) {
////            cout << __LINE__ << ":         Removing from neighbor's neighbor "<< nNeighbor << endl << flush;
            if (nNeighbor == vertex) continue;
            neighbors[nNeighbor].erase(neighbor);
////            cout << __LINE__ << ":         Done removing" << endl << flush;
            if (removed.find(nNeighbor) == removed.end())
                remaining.insert(nNeighbor);
        }
////        cout << __LINE__ << ":     Done Removing neighbor" << neighbor << endl << flush;
        neighbors[neighbor].clear();
////        cout << __LINE__ << ": Cleared neighbors: " << neighbor << endl << flush;
    }

////    cout << __LINE__ << ": Done Removing vertex " << vertex << endl << flush;
    neighbors[vertex].clear();
////    cout << __LINE__ << ": Cleared neighbors: " << vertex << endl << flush;
}

void Isolates::RemoveAllIsolates(vector<int> &vIsolateVertices, vector<int> &vOtherRemovedVertices, vector<pair<int,int>> &vAddedEdges)
{
    remaining = inGraph; // TODO/DS : Put back?
    //cout << "Removing all isolates." << endl << flush;
    int isolateSize(isolates.size());
    while (!remaining.empty()) {
////            if (remaining.size() %10000 == 0)
        //cout << "Remaining: " << remaining.size() << ", Removed: " << removed.size() << endl << flush;
        set<int>::iterator sit = remaining.begin();
        int const vertex = *sit;
        remaining.erase(sit);

////        cout << "Attempting to remove vertex " << vertex << endl << flush;

        bool reduction = RemoveIsolatedClique(vertex, vIsolateVertices, vOtherRemovedVertices);
////        if (!reduction) {
////            reduction = RemoveIsolatedPath(vertex, vIsolateVertices, vOtherRemovedVertices, vAddedEdges);
////        }
    }

////    cout << "Removed " << isolates.size() - isolateSize << " isolates." << endl << flush;
}

void Isolates::ReplaceAllRemoved(vector<int> const &vRemoved)
{
    inGraph.insert(vRemoved.begin(), vRemoved.end());
    for (int const removedVertex : vRemoved) {
        isolates.erase(removedVertex);
        removed.erase(removedVertex);
        m_AlternativeVertices.erase(removedVertex);
    }
////    cout << "Replacing all removed vertices." << endl << flush;
    for (int const removedVertex : vRemoved) {
////        if (remaining.size() %10000 == 0)
////            cout << "Remaining: " << remaining.size() << ", Removed: " << removed.size() << endl << flush;
        for (int const neighbor : m_AdjacencyArray[removedVertex]) {
            if (inGraph.find(neighbor) == inGraph.end()) continue;
            neighbors[removedVertex].insert(neighbor);
            neighbors[neighbor].insert(removedVertex);
        }
    }
////    cout << "Done replacing vertices..." << endl << flush;
}

int Isolates::NextVertexToRemove()
{
#if 1
    int vertexWithMaxReductions(-1);

    remaining = inGraph;

    set<int> tempInGraph(inGraph);

    int maxIsolates(-1);
////    cout << "Testing remaining " << inGraph.size() << " vertices, to maximize isolate removal" << endl << flush;
    for (int const vertex : tempInGraph) {
////        cout << "Starting loop..." << endl << flush;
        set<int> tempRemaining(remaining);
        set<int> tempRemoved(removed);
        set<int> tempIsolates(isolates);

        vector<int> vRemoved;
////        cout << "Try removing vertex " << vertex << endl << flush;
        RemoveVertexAndNeighbors(vertex, vRemoved);
////        cout << "And the new isolates..." << vertex << endl << flush;
        vector<pair<int,int>> vAddedEdges;
        RemoveAllIsolates(vRemoved, vRemoved, vAddedEdges);
        if (static_cast<int>(isolates.size() - tempIsolates.size()) > maxIsolates) {
            maxIsolates = static_cast<int>(isolates.size() - tempIsolates.size());
            vertexWithMaxReductions = vertex;
        }

        remaining = std::move(tempRemaining);
        removed   = std::move(tempRemoved);
        isolates  = std::move(tempIsolates);
        inGraph   = tempInGraph;

        // need to do this after restoring variable inGraph
        RemoveEdges(vAddedEdges); // TODO/DS: Put back...
        ReplaceAllRemoved(vRemoved);

////        cout << "Done with loop" << endl;
    }

    // if no vertex was found, then return max degree vertex
    if (vertexWithMaxReductions != -1) {
        cout << "Removing vertex " << vertexWithMaxReductions << " maximizes isolate removal (" << maxIsolates << " isolates)." << endl;
        return vertexWithMaxReductions;
    }
#endif

    int maxDegreeVertex(-1);
    int maxDegree(INT_MIN);
    for (int i = 0; i < neighbors.size(); ++i) {
        if (static_cast<int>(neighbors[i].size()) > maxDegree && removed.find(i) == removed.end()) {
            maxDegree = static_cast<int>(neighbors[i].size());
            maxDegreeVertex = i;
        }
    }

    return maxDegreeVertex;
}

int Isolates::GetAlternativeVertex(int const vertex) const
{
    map<int,int>::const_iterator cit(m_AlternativeVertices.find(vertex));
    if (cit != m_AlternativeVertices.end()) {
        return cit->second;
    }
    return -1;
}

