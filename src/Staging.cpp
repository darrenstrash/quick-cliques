// local includes
#include "Staging.h"
#include "CliqueTools.h"

// system includes
#include <vector>
#include <list>
#include <set>
#include <iostream>
#include <algorithm>
#include <climits>

using namespace std;

Staging::Staging(std::vector<std::vector<int>> &adjacencyList)
 : Algorithm("staging")
 , m_AdjacencyList(adjacencyList)
{
}

Staging::~Staging()
{
}

// SLOW remove isolated cliques:
#if 0
    cout << "Removing isolated cliques" << endl << flush;

    set<int> isolates;
    set<int> removed;
    set<int> remaining;

    vector<set<int>> neighbors(m_AdjacencyList.size());
    for (size_t u = 0; u < m_AdjacencyList.size(); ++u) {
        neighbors[u].insert(m_AdjacencyList[u].begin(), m_AdjacencyList[u].end());
        remaining.insert(u);
    }

////    for (vector<int> &neighbors : m_AdjacencyList) {
////        sort(neighbors.begin(), neighbors.end());
////    }

////    vector<bool> markedVertices(m_AdjacencyList.size(), false);

    while (!remaining.empty()) {
        cout << "Remaining: " << remaining.size() << ", Removed: " << removed.size() << endl << flush;
        set<int>::iterator sit = remaining.begin();
        int const vertex = *sit;
        remaining.erase(sit);

        vector<int> intersection(neighbors[vertex].begin(), neighbors[vertex].end());
        int neighborCount(0);
        vector<int> neighborhood;
        for (int const neighbor : neighbors[vertex]) {
            neighborhood.push_back(neighbor);

            set<int> hub(neighbors[neighbor]);
            hub.insert(neighbor);

////            cout << "Iteration(" << neighborCount << ": Intersection=";
////            for (int const vertex : intersection) {
////                cout << vertex << " ";
////            }
////            cout << endl << flush;
            neighborCount++;
            vector<int> intersectionTemp(intersection.size(),-1);
            vector<int>::iterator it = set_intersection(hub.begin(), hub.end(), intersection.begin(), intersection.end(), intersectionTemp.begin());
            intersectionTemp.resize(it - intersectionTemp.begin());
            intersection = std::move(intersectionTemp);
        }

////        cout << "Iteration(" << neighborCount << ": Intersection=";
////        for (int const vertex : intersection) {
////            cout << vertex << " ";
////        }
////        cout << endl << flush;

        vector<int> commonNeighborhood(intersection.size(), -1);
        vector<int>::iterator it = set_intersection(neighborhood.begin(), neighborhood.end(), intersection.begin(), intersection.end(), commonNeighborhood.begin());

        if ((it - commonNeighborhood.begin()) == neighborhood.size()) {
            cout << "Removing " << vertex << " and its " << neighborCount << " remaining neighbors." << endl << flush;
            isolates.insert(vertex);
            removed.insert(neighbors[vertex].begin(), neighbors[vertex].end());
            removed.insert(vertex);
            for (int const neighbor : neighbors[vertex]) {
                remaining.erase(neighbor);
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
        }
    }

    cout << "Method 1: Removed " << removed.size() << " vertices." << endl << flush;

#endif

#if 0
    isolates.clear();
    removed.clear();

    for (int i = 0; i < m_AdjacencyList.size(); ++i) {
        if (removed.find(i) != removed.end()) continue;
////        cout << "Evaluating " << i << endl << flush;
////        bool isIsolate(true);
        vector<int> intersection(m_AdjacencyList[i]);
        int neighborCount(0);
        vector<int> neighborhood;
        for (int const neighbor : m_AdjacencyList[i]) {
////            cout << "    Neighbor?: " << neighbor << endl << flush;
            if (removed.find(neighbor) != removed.end()) continue;
////            cout << "    Neighbor: " << neighbor << endl << flush;
            neighborhood.push_back(neighbor);

            vector<int> hub(m_AdjacencyList[neighbor]);
            hub.push_back(neighbor);
            sort(hub.begin(), hub.end());

////            cout << "Iteration(" << neighborCount << ": Intersection=";
////            for (int const vertex : intersection) {
////                cout << vertex << " ";
////            }
////            cout << endl << flush;
            neighborCount++;
            vector<int> intersectionTemp(intersection.size(),-1);
            vector<int>::iterator it = set_intersection(hub.begin(), hub.end(), intersection.begin(), intersection.end(), intersectionTemp.begin());
            intersectionTemp.resize(it - intersectionTemp.begin());
            intersection = std::move(intersectionTemp);
        }

////        cout << "Iteration(" << neighborCount << ": Intersection=";
////        for (int const vertex : intersection) {
////            cout << vertex << " ";
////        }
////        cout << endl << flush;

        vector<int> commonNeighborhood(intersection.size(), -1);

        vector<int>::iterator it = set_intersection(neighborhood.begin(), neighborhood.end(), intersection.begin(), intersection.end(), commonNeighborhood.begin());

        if ((it - commonNeighborhood.begin()) == neighborhood.size()) {
////            cout << "Removing " << i << " and its " << neighborCount << " remaining neighbors." << endl << flush;
            isolates.insert(i);
            removed.insert(m_AdjacencyList[i].begin(), m_AdjacencyList[i].end());
            removed.insert(i);
            i = -1; // restart search from beginning, TODO/DS: make more efficient
        }
    }

    cout << "Removed " << removed.size() << " vertices." << endl << flush;
#endif //0


bool RemoveIsolatedClique(int const vertex, vector<set<int>> &neighbors, set<int> &remaining, set<int> &isolates, set<int> &removed, vector<bool> &vMarkedVertices)
{
    size_t neighborCount(0);
    bool superSet(true);
    for (int const neighbor : neighbors[vertex]) {
        if (neighbors[neighbor].size() < neighbors[vertex].size()) {
            superSet = false; break;
        }
    }

    if (vertex == 21952) {
        cout << "vertex" << vertex << " has " << neighbors[vertex].size() << " neighbors." << endl << flush;
    }

    if (!superSet) return false;

    for (int const neighbor : neighbors[vertex]) {
        neighborCount++;

        for (int const nNeighbor : neighbors[neighbor]) {
            vMarkedVertices[nNeighbor] = true;
        }
        vMarkedVertices[neighbor] = true;

        for (int const neighbor2 : neighbors[vertex]) {
            superSet = superSet && vMarkedVertices[neighbor2];
        }

        for (int const nNeighbor : neighbors[neighbor]) {
            vMarkedVertices[nNeighbor] = false;
        }
        vMarkedVertices[neighbor] = false;

        if (!superSet) return false;

        ////            cout << "Iteration(" << neighborCount << ": Intersection=";
        ////            for (int const vertex : intersection) {
        ////                cout << vertex << " ";
        ////            }
        ////            cout << endl << flush;
    }

    ////        cout << "Iteration(" << neighborCount << ": Intersection=";
    ////        for (int const vertex : intersection) {
    ////            cout << vertex << " ";
    ////        }
    ////        cout << endl << flush;

    if (superSet) {
////        cout << "Removing " << vertex << " and its " << neighborCount << " remaining neighbors." << endl << flush;
        removed.insert(neighbors[vertex].begin(), neighbors[vertex].end());
        removed.insert(vertex);
        isolates.insert(vertex);
        for (int const neighbor : neighbors[vertex]) {
            remaining.erase(neighbor);
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

    if (vertex == 21952) {
        cout << "vertex" << vertex << " is being removed." << endl << flush;
    }
        
        return true;
    }
    return false;
}

bool RemoveIsolatedPath(int const vertex, vector<set<int>> &neighbors, set<int> &remaining, set<int> &isolates, set<int> &removed, vector<bool> &vMarkedVertices)
{
    if (neighbors[vertex].size() != 2) return false;

    for (int const neighbor : neighbors[vertex]) {
        if (neighbors[neighbor].size() == 2) {
////            cout << "Removing path vertices!" << endl << flush;
            removed.insert(vertex);
            removed.insert(neighbor);
            isolates.insert(vertex);

            int endVertex(-1);

            // remove from other neighbor's list...
            for (int const nNeighbor : neighbors[neighbor]) {
                if (nNeighbor != vertex) {
                    endVertex = nNeighbor;
                    neighbors[nNeighbor].erase(neighbor);
                    remaining.insert(nNeighbor);
                    break;
                }
            }

            for (int const otherNeighbor : neighbors[vertex]) {
                if (otherNeighbor != neighbor) {
                    neighbors[otherNeighbor].erase(neighbor);
                    remaining.insert(otherNeighbor);

                    neighbors[otherNeighbor].insert(endVertex);
                    neighbors[endVertex].insert(otherNeighbor);
                    break;
                }
            }

            neighbors[vertex].clear();
            neighbors[neighbor].clear();

            remaining.erase(vertex);
            remaining.erase(neighbor);

            return true;
        }
    }
    return false;
}

bool RemoveVertexAndNeighbors(int const vertex, vector<set<int>> &neighbors, set<int> &remaining, set<int> &isolates, set<int> &removedVertices)
{
//    cout << __LINE__ << ": Removing vertex " << vertex << endl << flush;
    removedVertices.insert(vertex);
    remaining.erase(vertex);
    isolates.insert(vertex);

    for (int const neighbor : neighbors[vertex]) {
//        cout << __LINE__ << ":     Removing neighbor " << neighbor << endl << flush;
        removedVertices.insert(neighbor);
        remaining.erase(neighbor);
        for (int const nNeighbor : neighbors[neighbor]) {
//            cout << __LINE__ << ":         Removing from neighbor's neighbor "<< nNeighbor << endl << flush;
            if (nNeighbor == vertex) continue;
            neighbors[nNeighbor].erase(neighbor);
//            cout << __LINE__ << ":         Done removing" << endl << flush;
            if (removedVertices.find(nNeighbor) == removedVertices.end())
                remaining.insert(nNeighbor);
        }
//        cout << __LINE__ << ":     Done Removing neighbor" << neighbor << endl << flush;
        neighbors[neighbor].clear();
//        cout << __LINE__ << ": Cleared neighbors: " << neighbor << endl << flush;
    }

//    cout << __LINE__ << ": Done Removing vertex " << vertex << endl << flush;
    neighbors[vertex].clear();
//    cout << __LINE__ << ": Cleared neighbors: " << vertex << endl << flush;
    return true;
}


bool RemoveMinDegreeVertex(vector<set<int>> &neighbors, set<int> &remaining, set<int> &isolates, set<int> &removed, vector<bool> &vMarkedVertices)
{
    if (removed.size() == neighbors.size()) return false;

    int vertex(-1);
    int minDegree(INT_MAX);
    for (int i = 0; i < neighbors.size(); ++i) {
        if (neighbors[i].size() < minDegree && removed.find(i) == removed.end()) {
            minDegree = neighbors[i].size();
            vertex = i;
        }
    }

    //cout << "Removing vertex " << vertex << " with degree " << minDegree << endl << flush;

    RemoveVertexAndNeighbors(vertex, neighbors, remaining, isolates, removed);

    return true;
}

bool RemoveMaxDegreeVertex(vector<set<int>> &neighbors, set<int> &remaining, set<int> &isolates, set<int> &removed, vector<bool> &vMarkedVertices)
{
    if (removed.size() == neighbors.size()) return false;

    int vertex(-1);
    int maxDegree(INT_MIN);
    for (int i = 0; i < neighbors.size(); ++i) {
        if (static_cast<int>(neighbors[i].size()) > maxDegree && removed.find(i) == removed.end()) {
            maxDegree = static_cast<int>(neighbors[i].size());
            vertex = i;
        }
    }

    cout << "Removing vertex " << vertex << " with degree " << maxDegree << endl << flush;

    RemoveVertexAndNeighbors(vertex, neighbors, remaining, isolates, removed);

    return true;
}

#if 0
// compute maximum common neighbors.
    size_t maxCommonNeighbors(0);

    for (vector<int> &neighbors : m_AdjacencyList) {
        sort(neighbors.begin(), neighbors.end());
    }

    vector<int> vCommon(m_AdjacencyList.size(), -1);

    for (size_t u = 0; u < m_AdjacencyList.size(); ++u) {
        for (size_t v = u+1; v < m_AdjacencyList.size(); ++v) {
            size_t const commonNeighbors = set_difference(m_AdjacencyList[u].begin(), m_AdjacencyList[u].end(), m_AdjacencyList[v].begin(), m_AdjacencyList[v].end(), vCommon.begin()) - vCommon.begin();
            if (commonNeighbors > maxCommonNeighbors)
                maxCommonNeighbors = commonNeighbors;
        }
    }

    cout << "Maximum common neighbors: " << maxCommonNeighbors << endl << flush;
#endif

bool RemoveAllIsolates(vector<set<int>> &neighbors, set<int> &remaining, set<int> &isolates, set<int> &removed, vector<bool> &vMarkedVertices)
{
    //cout << "Removing all isolates." << endl << flush;
    int isolateSize(isolates.size());
    while (!remaining.empty()) {
        ////            if (remaining.size() %10000 == 0)
        //cout << "Remaining: " << remaining.size() << ", Removed: " << removed.size() << endl << flush;
        set<int>::iterator sit = remaining.begin();
        int const vertex = *sit;
        remaining.erase(sit);

        vector<int> intersection(neighbors[vertex].begin(), neighbors[vertex].end());

        bool reduction = RemoveIsolatedClique(vertex, neighbors, remaining, isolates, removed, vMarkedVertices);
        if (!reduction) {
            reduction = RemoveIsolatedPath(vertex, neighbors, remaining, isolates, removed, vMarkedVertices);
        }
    }

    cout << "Removed " << isolates.size() - isolateSize << " isolates." << endl << flush;
    return true;
}


template <typename C> bool ReplaceRemovedVertices(vector<vector<int>> const &adjacencyArray, vector<set<int>> &neighbors, set<int> const &inGraph, C const &removedVertices)
{
    //cout << "Replacing all removed vertices." << endl << flush;
    for (typename C::value_type const removedVertex : removedVertices) {
        ////            if (remaining.size() %10000 == 0)
////        cout << "Remaining: " << remaining.size() << ", Removed: " << removed.size() << endl << flush;
        for (int const neighbor : adjacencyArray[removedVertex]) {
            if (inGraph.find(neighbor) == inGraph.end()) continue;
            neighbors[removedVertex].insert(neighbor);
            neighbors[neighbor].insert(removedVertex);
        }
    }
    //cout << "Replaced vertices..." << endl << flush;
    return true;
}



void Staging::Run()
{
    cout << "Running Staging..." << endl << flush;

    //CliqueTools::FindMaximalIndependentSetInCliqueGraph(m_AdjacencyList);


    cout << "Applying Reductions..." << endl << flush;

    set<int> isolates;
    set<int> removed;
    set<int> remaining;
    set<int> inGraph;

    vector<set<int>> neighbors(m_AdjacencyList.size());
    for (size_t u = 0; u < m_AdjacencyList.size(); ++u) {
        neighbors[u].insert(m_AdjacencyList[u].begin(), m_AdjacencyList[u].end());
        remaining.insert(u);
        inGraph.insert(u);
    }

    vector<bool> vMarkedVertices(m_AdjacencyList.size(), false);

    cout << "Removing isolates..." << endl;
    RemoveAllIsolates(neighbors, remaining, isolates, removed, vMarkedVertices);
    for (int const vertex : removed) {
        inGraph.erase(vertex);
    }
    cout << "# vertices remaining in graph: " << inGraph.size() << "/" << m_AdjacencyList.size() << endl << flush;

    // as of now, this loop helps find a really small Independent set, this is to help limit the recursion depth.
    // need to expand this into a branch-and-bound, pick the vertices that constrain the search the most, to keep
    // the search at a reasonable depth. And try to "peel off" as many vertices as possible through reductions.
    while (removed.size() != m_AdjacencyList.size()) {
#if 1
        // find vertex whose removal maximizes the number of vertices moved to IS by reduction.
        int vertexWithMaxReductions(-1);

        remaining = inGraph;

        int maxIsolates(-1);
        cout << "Testing remaining " << inGraph.size() << " vertices, to maximize isolate removal" << endl << flush;
        for (int const vertex : inGraph) {
                set<int> tempRemaining(remaining);
                set<int> tempRemoved(removed);
                set<int> tempIsolates(isolates);

                RemoveVertexAndNeighbors(vertex, neighbors, tempRemaining, tempIsolates, tempRemoved);
                RemoveAllIsolates(neighbors, tempRemaining, tempIsolates, tempRemoved, vMarkedVertices);
                if (static_cast<int>(tempIsolates.size() - isolates.size()) > maxIsolates) {
                    maxIsolates = static_cast<int>(tempIsolates.size() - isolates.size());
                    vertexWithMaxReductions = vertex;
                }

                vector<int> newlyRemovedVertices(tempRemoved.size() - removed.size(), -1);

                vector<int>::iterator it = set_difference(tempRemoved.begin(), tempRemoved.end(), removed.begin(), removed.end(), newlyRemovedVertices.begin());

                newlyRemovedVertices.resize(it - newlyRemovedVertices.begin());

                ReplaceRemovedVertices(m_AdjacencyList, neighbors, inGraph, newlyRemovedVertices);
        }

        cout << "Removing vertex " << vertexWithMaxReductions << " maximizes isolate removal (" << maxIsolates << " isolates)." << endl << flush;
        //cout << "Removed " << removed.size() << "/" << m_AdjacencyList.size() << " vertices" << endl << flush;

        // don't actually remove anything, just run like before. // TODO/DS actually remove vertex!
        set<int> tempRemoved(removed);
        if (vertexWithMaxReductions == -1) {
            //RemoveAllIsolates    (neighbors, remaining, isolates, tempRemoved, vMarkedVertices); // TODO: remove, does nothing
            RemoveMaxDegreeVertex(neighbors, remaining, isolates, tempRemoved, vMarkedVertices); // TODO: remove max degree instead?
        } else {
            RemoveVertexAndNeighbors(vertexWithMaxReductions, neighbors, remaining, isolates, tempRemoved);
            RemoveAllIsolates(neighbors, remaining, isolates, tempRemoved, vMarkedVertices);
        }

        cout << "# vertices remaining in graph: " << inGraph.size() << "/" << m_AdjacencyList.size() << endl << flush;

        vector<int> newlyRemovedVertices(tempRemoved.size() - removed.size(), -1);
        vector<int>::iterator it = set_difference(tempRemoved.begin(), tempRemoved.end(), removed.begin(), removed.end(), newlyRemovedVertices.begin()); // TODO: This is an expensive diff, just grab the added vertices before we get to this point.
        newlyRemovedVertices.resize(it - newlyRemovedVertices.begin());

        for (int const removedVertex : newlyRemovedVertices) {
            //cout << "Removing " << removedVertex << " from graph " << endl;
            inGraph.erase(removedVertex);
        }
//        inGraph.erase(tempRemoved.begin(), tempRemoved.end());
        removed.insert(newlyRemovedVertices.begin(), newlyRemovedVertices.end());
#else
        RemoveAllIsolates(neighbors, remaining, isolates, removed, vMarkedVertices);
        RemoveMinDegreeVertex(neighbors, remaining, isolates, removed, vMarkedVertices);
#endif

    }

    cout << "Removed " << removed.size() << " vertices." << endl << flush;
    cout << "Found independent set of size: " << isolates.size() << endl << flush;
}
