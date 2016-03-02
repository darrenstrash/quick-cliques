
#include "MatchingTools.h"
#include "GraphTools.h"
#include "BiDoubleGraph.h"

#include <set>
#include <vector>
#include <list>

#define UNMATCHED_VERTEX -1

using namespace std;

set<int> MatchingTools::ComputeLeftMIS(vector<vector<int>> const &biDoubleGraph)
{
    // each edge in matching is indexed by first vertex, and contains the second vertex
    // initially there are no edges in matching (thus invalid second vertex in matching).
    vector<int> matching(biDoubleGraph.size(), UNMATCHED_VERTEX);

    // residual graph is represented the same as edge in matching.
    // we can iterate over neighbors (excluding neighbor in matching
    // to find residual paths.
////    vector<int> residualGraph(matching);

    cout << "Computing BiDoubleMIS..." << endl << flush;

    auto inLeftSide = [&biDoubleGraph] (int const vertex) {
        return (vertex < biDoubleGraph.size()/2);
    };

    // finds a path from start to finish with excess capacity.
    auto findPath = [&inLeftSide, &biDoubleGraph] (vector<int> const &matching) {
        vector<bool> inStack(matching.size(), false);
        vector<bool> evaluated(matching.size(), false);
        list<int> stack;

        // depth first search, starting from imaginary start vertex
        // i.e., search starts from left-side.
        for (size_t index = 0; index < matching.size()/2; ++index) {
////        for (size_t index = matching.size()/2; index > 0; --index) {
            // only insert vertices without edges in matching, otherwise
            // imaginary first vertex has residual capacity 0 to that vertex.
            if (matching[index] != UNMATCHED_VERTEX) continue;
            stack.push_back(index);
            inStack[index] = true;
        }

        vector<int> vPreviousVertexOnPath(biDoubleGraph.size(), UNMATCHED_VERTEX);
        int endVertex(UNMATCHED_VERTEX);

        bool foundPath(false);
        while (!stack.empty() && !foundPath) {
            int const vertex = stack.back(); stack.pop_back();
            evaluated[vertex] = true;
            inStack[vertex] = false;
            for (int const neighbor : biDoubleGraph[vertex]) {
                // evaluate neighbor if the edge to that neighbor has residual capacity.
                if (evaluated[neighbor] || inStack[neighbor]) continue;

                // forward edge with residual capacity
                if (inLeftSide(vertex) && matching[vertex] != neighbor) {
                    vPreviousVertexOnPath[neighbor] = vertex;
                    stack.push_back(neighbor);
                    inStack[neighbor] = true;

                    if (!inLeftSide(neighbor) && matching[neighbor] == UNMATCHED_VERTEX) { //found path
                        foundPath = true;
                        endVertex = neighbor;
                        break;
                    }
                }
                // backward edge that we can "undo" by pushing flow back...
                else if (inLeftSide(neighbor) && matching[neighbor] == vertex) {
////                else if (!inLeftSide(vertex) && matching[neighbor] == vertex) {
                    vPreviousVertexOnPath[neighbor] = vertex;
                    stack.push_back(neighbor);
                    inStack[neighbor] = true;
                }
            }
        }

        vector<int> vPath;
        if (endVertex == UNMATCHED_VERTEX) return vPath;
        vPath.push_back(endVertex);
        while (vPreviousVertexOnPath[endVertex] != UNMATCHED_VERTEX) {
            vPath.push_back(vPreviousVertexOnPath[endVertex]);
            endVertex = vPreviousVertexOnPath[endVertex];
        }

        std::reverse(vPath.begin(), vPath.end());

////        cout << "Path through residual graph: ";
////        for (int const vertex : vPath) {
////            cout << vertex << " ";
////        }
////        cout << endl;
        return vPath;
    };

    vector<int> path;
    path = findPath(matching);
    size_t iterations(0);
    while (!path.empty()) {
        iterations++;
////        cout << "Found path of length " << path.size() << ", iteration: " << iterations << endl << flush;
        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const backwardEdge(vertex1 > vertex2);
            if (backwardEdge) {
                matching[vertex1] = UNMATCHED_VERTEX;
                matching[vertex2] = UNMATCHED_VERTEX;
////                cout << "Remove backward edge " << vertex1 << "->" << vertex2 << " from matching" << endl << flush;
            }
////            if (matching[vertex1] != UNMATCHED_VERTEX) {
////                matching[matching[vertex1]] = UNMATCHED_VERTEX;
////            }
////            if (matching[vertex2] != UNMATCHED_VERTEX) {
////                matching[matching[vertex2]] = UNMATCHED_VERTEX;
////            }
////            matching[vertex1] = UNMATCHED_VERTEX;
////            matching[vertex2] = UNMATCHED_VERTEX;
        }

        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const forwardEdge(vertex1 < vertex2);
            if (forwardEdge) {
                matching[vertex1] = vertex2;
                matching[vertex2] = vertex1;

////                cout << "Add    forward  edge " << vertex1 << "->" << vertex2 << " to   matching" << endl << flush;
            }
        }

        path = findPath(matching);
    }

    // From each unmatched vertex, calculate the alternating paths.
    // -- Paths that alternate matched edges and non-matched edges.
    // Mark all vertices that are connected by an alternating path.

#ifdef VERIFY
    bool bVerified(true);
    for (int vertex = 0; vertex < matching.size(); ++vertex) {
        if (matching[vertex] != UNMATCHED_VERTEX) {
            if (matching[matching[vertex]] != vertex) {
                cout << "ERROR! mismatch: " << vertex << " -> " << matching[vertex] << " -> " << matching[matching[vertex]] << endl << flush;
                bVerified = false;
                break;
            }
        }
    }

    cout << "Verification of matching " << (bVerified ? "passed!" : "failed!") << endl << flush;
#endif // VERIFY

    enum LastEdge {NO_LAST_EDGE, MATCHED_LAST_EDGE, UNMATCHED_LAST_EDGE, BOTH_LAST_EDGE, NULL_LAST_EDGE};

    // TODO/DS: Make more efficient. Stop search when after marking a vertex in each matchings.
    // i.e., number of vertices marked = number of matchings.
    vector<LastEdge> marked(matching.size(), NO_LAST_EDGE);

    // first, mark unmatched vertices, and their neighbors.
    vector<int> matchedVertices;

    vector<int> previousVertex(matching.size(), -1);


    // TODO/DS: Check that lambda works 
    auto isMatchedEdge = [&matching] (int const vertex, int const neighbor) {
        return matching[neighbor] == vertex; 
    };

    int coverVertexCount(0);
    //Mark all unmatched vertices on left-hand side
    //mark their matched neighbors.
    for (int vertex = 0; vertex < matching.size()/2; ++vertex) {
        if (matching[vertex] == UNMATCHED_VERTEX) {
////            cout << "Mark " << vertex << " as unmatched!" << endl << flush;
            marked[vertex] = LastEdge::NULL_LAST_EDGE;
            coverVertexCount++;
            for (int const neighbor : biDoubleGraph[vertex]) {
                if (marked[neighbor] == LastEdge::NO_LAST_EDGE)
                    coverVertexCount++;
                marked[neighbor] = LastEdge::UNMATCHED_LAST_EDGE;
                previousVertex[neighbor] = vertex;
////                if (isMatchedEdge(vertex, neighbor)) {
////                    cout << "ERROR! (" << __LINE__ << ")" << endl << flush;
////                }
            }
        } else {
           matchedVertices.push_back(vertex); 
        }
    }

    cout << "Matching       has " << matchedVertices.size() << " vertices" << endl << flush;
    cout << "Cover count     is " << coverVertexCount                      << endl << flush;


    int iteration(0);
    while (iteration < matching.size()) {
        int const savedCoverCount(coverVertexCount);
        //cout << "Iteration: " << iteration << "/" << matching.size() << endl << flush;
        iteration++;
////        for (int const matchedVertex : matchedVertices) {
        for (int matchedVertex = 0; matchedVertex < matching.size(); matchedVertex++) {
            if (marked[matchedVertex] == NO_LAST_EDGE) continue;

            LastEdge const edgeIntoVertex(marked[matchedVertex]);

            // If the last edge was matched, follow unmatched edges
            if (edgeIntoVertex == LastEdge::MATCHED_LAST_EDGE) {
////                if (!inLeftSide(matchedVertex))
////                    cout << "ERROR! (" << __LINE__ << ") Vertex " << matchedVertex << " has   matched last edge " << endl << flush;
                for (int const neighbor : biDoubleGraph[matchedVertex]) {
                    if (!isMatchedEdge(matchedVertex, neighbor)) {
                        if (marked[neighbor] == LastEdge::NO_LAST_EDGE) {
                            marked[neighbor] = LastEdge::UNMATCHED_LAST_EDGE;
                            coverVertexCount++;
                        }
                    }
////                    if (marked[neighbor] == LastEdge::MATCHED_LAST_EDGE)
////                        cout << "ERROR! Vertex " << neighbor << " on " << (inLeftSide(neighbor)? "left": "right") << " has   matched last edge." << endl << flush;
                }

            // If last edge was unmatched, follow matched edges.
            } else if (edgeIntoVertex == LastEdge::UNMATCHED_LAST_EDGE) {
////                if (inLeftSide(matchedVertex))
////                    cout << "ERROR! (" << __LINE__ << ") Vertex " << matchedVertex << " has unmatched last edge " << endl << flush;
                for (int const neighbor : biDoubleGraph[matchedVertex]) {
                    if (isMatchedEdge(matchedVertex, neighbor)) {
                        if (marked[neighbor] == NO_LAST_EDGE) {
                            marked[neighbor] = LastEdge::MATCHED_LAST_EDGE;
                            coverVertexCount++;
                        }
////                        else if (marked[neighbor] == LastEdge::UNMATCHED_LAST_EDGE)
////                            cout << "ERROR! Vertex " << neighbor << " on " << (inLeftSide(neighbor)? "left": "right") << " has unmatched last edge." << endl << flush;
                    }
                }
            }
        }

        if (savedCoverCount != coverVertexCount) {
            cout << "Cover increased to " << coverVertexCount << endl << flush;
        }
        // TODO/DS: doesn't work, I'm not sure why...
////        else {
////            break;
////        }
    }

    int matchVertices(0);
    for (int const vertex : matching) {
        if (vertex != UNMATCHED_VERTEX) {
            matchVertices++;
        }
    }

    matchVertices >>= 1;

    vector<bool> mvc(matching.size(), false);
    size_t numMVCVertices(0);
    set<int> misToReturn;
    size_t numTotalMISVertices(0);

    // put vertices in mvc
    for (int vertex = 0; vertex < matching.size(); ++vertex) {
        if ((inLeftSide(vertex) && (marked[vertex] == LastEdge::NO_LAST_EDGE)) ||
            (!inLeftSide(vertex) && (marked[vertex] != LastEdge::NO_LAST_EDGE))) {
#ifdef VERIFY
            if (matching[vertex] == UNMATCHED_VERTEX) {
                cout << "ERROR! vertex " << vertex << " is not in matching!" << endl << flush;
            }
#endif //VERIFY
            mvc[vertex] = true;
            numMVCVertices++;

            if (inLeftSide(vertex) && marked[vertex] == LastEdge::UNMATCHED_LAST_EDGE) {
                cout << "ERROR!" << endl << flush;
            }
            if (!inLeftSide(vertex) && marked[vertex] == LastEdge::MATCHED_LAST_EDGE) {
                cout << "ERROR!" << endl << flush;
            }
        }
        else {
            numTotalMISVertices++;
            if (inLeftSide(vertex)) {
                misToReturn.insert(vertex);
            }

            // TODO/DS: I don't know if this belongs or not...
////            else {
////                misToReturn.insert(vertex - matching.size()/2);
////            }
        }
    }

#ifdef VERIFY
    for (int const vertex : misToReturn) {
        for (int const neighbor : biDoubleGraph[vertex]) {
            if (misToReturn.find(neighbor) != misToReturn.end()) {
                cout << "ERROR! Failed independent set validation." << endl << flush;
            }
        }
    }
#endif // VERIFY


    cout << "MVC             Vertices: " << numMVCVertices << endl << flush;
    cout << "MIS             Vertices: " << misToReturn.size() << endl << flush;
    cout << "Cover and Match Vertices: " << coverVertexCount << "/" << matchVertices << endl << flush;

    if (numMVCVertices != matchVertices) {
        cout << "ERROR! MVC != size of matching" << endl << flush;
    }

    // first mark unmatched vertices
    // mark them as having the last edge matched, because the only option
    // for the next edge is an unmatched edge.
////    for (int vertex = 0; vertex < matching.size(); ++vertex) {
////        if (matching[vertex] == UNMATCHED_VERTEX) {
////            marked[vertex] = LastEdge::MATCHED_LAST_EDGE;
////        }
////    }
////    for (int vertex = 0; vertex < matching.size(); ++vertex) {

////    for (int const vertex : matching) {
////        if (inLeftSide(vertex) && vertex == UNMATCHED_VERTEX) {
////            misToReturn.push_back(vertex);
////        }
////    }

    return misToReturn;
}

void MatchingTools::ComputeAlternatingPaths(BiDoubleGraph const &biDouble, vector<int> const &vMatching, vector<MatchingTools::LastEdge> &vOnAlternatingPath)
{
    // first, mark unmatched vertices, and their neighbors.
    vector<int> matchedVertices;

    vector<int> previousVertex(vMatching.size(), -1);

    auto isMatchedEdge = [&vMatching] (int const vertex, int const neighbor) {
        return vMatching[neighbor] == vertex; 
    };

    int coverVertexCount(0);
    //Mark all unmatched vertices on left-hand side
    //mark their matched neighbors.
    for (int vertex = 0; vertex < vMatching.size()/2; ++vertex) {
        if (vMatching[vertex] == UNMATCHED_VERTEX) {
////            cout << "Mark " << vertex << " as unmatched!" << endl << flush;
            vOnAlternatingPath[vertex] = MatchingTools::LastEdge::NULL_LAST_EDGE;
            coverVertexCount++;
            for (int const neighbor : biDouble.Neighbors(vertex)) {
                if (vOnAlternatingPath[neighbor] == MatchingTools::LastEdge::NO_LAST_EDGE)
                    coverVertexCount++;
                vOnAlternatingPath[neighbor] = MatchingTools::LastEdge::UNMATCHED_LAST_EDGE;
                previousVertex[neighbor] = vertex;
////                if (isMatchedEdge(vertex, neighbor)) {
////                    cout << "ERROR! (" << __LINE__ << ")" << endl << flush;
////                }
            }
        } else {
           matchedVertices.push_back(vertex); 
        }
    }

    cout << "Matching       has " << matchedVertices.size() << " vertices" << endl << flush;
    cout << "Cover count     is " << coverVertexCount                      << endl << flush;


    int iteration(0);
    while (iteration < vMatching.size()) {
        int const savedCoverCount(coverVertexCount);
        //cout << "Iteration: " << iteration << "/" << vMatching.size() << endl << flush;
        iteration++;
////        for (int const matchedVertex : matchedVertices) {
        for (int matchedVertex = 0; matchedVertex < vMatching.size(); matchedVertex++) {
            if (vOnAlternatingPath[matchedVertex] == NO_LAST_EDGE) continue;

            LastEdge const edgeIntoVertex(vOnAlternatingPath[matchedVertex]);

            // If the last edge was matched, follow unmatched edges
            if (edgeIntoVertex == MatchingTools::LastEdge::MATCHED_LAST_EDGE) {
////                if (!inLeftSide(matchedVertex))
////                    cout << "ERROR! (" << __LINE__ << ") Vertex " << matchedVertex << " has   matched last edge " << endl << flush;
                for (int const neighbor : biDouble.Neighbors(matchedVertex)) {
                    if (!isMatchedEdge(matchedVertex, neighbor)) {
                        if (vOnAlternatingPath[neighbor] == MatchingTools::LastEdge::NO_LAST_EDGE) {
                            vOnAlternatingPath[neighbor] = MatchingTools::LastEdge::UNMATCHED_LAST_EDGE;
                            coverVertexCount++;
                        }
                    }
////                    if (vOnAlternatingPath[neighbor] == MatchingTools::LastEdge::MATCHED_LAST_EDGE)
////                        cout << "ERROR! Vertex " << neighbor << " on " << (inLeftSide(neighbor)? "left": "right") << " has   matched last edge." << endl << flush;
                }

            // If last edge was unmatched, follow matched edges.
            } else if (edgeIntoVertex == MatchingTools::LastEdge::UNMATCHED_LAST_EDGE) {
////                if (inLeftSide(matchedVertex))
////                    cout << "ERROR! (" << __LINE__ << ") Vertex " << matchedVertex << " has unmatched last edge " << endl << flush;
                for (int const neighbor : biDouble.Neighbors(matchedVertex)) {
                    if (isMatchedEdge(matchedVertex, neighbor)) {
                        if (vOnAlternatingPath[neighbor] == NO_LAST_EDGE) {
                            vOnAlternatingPath[neighbor] = MatchingTools::LastEdge::MATCHED_LAST_EDGE;
                            coverVertexCount++;
                        }
////                        else if (vOnAlternatingPath[neighbor] == MatchingTools::LastEdge::UNMATCHED_LAST_EDGE)
////                            cout << "ERROR! Vertex " << neighbor << " on " << (inLeftSide(neighbor)? "left": "right") << " has unmatched last edge." << endl << flush;
                    }
                }
            }
        }

        if (savedCoverCount != coverVertexCount) {
            cout << "Cover increased to " << coverVertexCount << endl << flush;
        }
        // TODO/DS: doesn't work, I'm not sure why...
////        else {
////            break;
////        }
    }

    cout << "Cover           Vertices: " << coverVertexCount << endl << flush;
}

void MatchingTools::ComputeAlternatingPathsOptimized(BiDoubleGraph const &biDouble, vector<int> const &vMatching, vector<MatchingTools::LastEdge> &vOnAlternatingPath)
{
    auto isMatchedEdge = [&vMatching] (int const vertex, int const neighbor) {
        return vMatching[neighbor] == vertex; 
    };

    vector<bool> vEvaluated(vMatching.size(), false);
    vector<int> vStack;
    vStack.reserve(vMatching.size());

    // First mark all unmatched vertices on left-hand side and their neighbors.
    for (int vertex = 0; vertex < vMatching.size()/2; ++vertex) {
        if (vMatching[vertex] != UNMATCHED_VERTEX) continue;
        vEvaluated[vertex] = true;
        vOnAlternatingPath[vertex] = MatchingTools::LastEdge::NULL_LAST_EDGE;
        for (int const neighbor : biDouble.Neighbors(vertex)) {
            vOnAlternatingPath[neighbor] = MatchingTools::LastEdge::UNMATCHED_LAST_EDGE;
            if (!vEvaluated[neighbor]) vStack.push_back(neighbor);
            vEvaluated[neighbor] = true;
        }
    }

    while (!vStack.empty()) {
        int const matchedVertex(vStack.back()); vStack.pop_back();

        LastEdge const edgeIntoVertex(vOnAlternatingPath[matchedVertex]);

        // If the last edge was matched, follow unmatched edges
        if (edgeIntoVertex == MatchingTools::LastEdge::MATCHED_LAST_EDGE) {
            for (int const neighbor : biDouble.Neighbors(matchedVertex)) {
                if (!isMatchedEdge(matchedVertex, neighbor)) {
                    if (vOnAlternatingPath[neighbor] == MatchingTools::LastEdge::NO_LAST_EDGE) {
                        vOnAlternatingPath[neighbor] = MatchingTools::LastEdge::UNMATCHED_LAST_EDGE;

                        if (!vEvaluated[neighbor]) {
                            vStack.push_back(neighbor);
                            vEvaluated[neighbor] = true;
                        }
                    }
                }
            }

        // If last edge was unmatched, follow matched edge.
        } else if (edgeIntoVertex == MatchingTools::LastEdge::UNMATCHED_LAST_EDGE) {
            int const neighbor(vMatching[matchedVertex]);
            if (neighbor != UNMATCHED_VERTEX) {
                if (vOnAlternatingPath[neighbor] == NO_LAST_EDGE) {
                    vOnAlternatingPath[neighbor] = MatchingTools::LastEdge::MATCHED_LAST_EDGE;

                    if (!vEvaluated[neighbor]) {
                        vStack.push_back(neighbor);
                        vEvaluated[neighbor] = true;
                    }
                }
            }
        }
    }
}

set<int> MatchingTools::GetLeftVerticesOnAlternatingPaths(BiDoubleGraph const &biDouble, vector<int> const &vMatching, vector<MatchingTools::LastEdge> &vOnAlternatingPath)
{
    set<int> misToReturn;
    for (int vertex = 0; vertex < vMatching.size()/2; ++vertex) {
        if (vOnAlternatingPath[vertex] != MatchingTools::LastEdge::NO_LAST_EDGE) misToReturn.insert(vertex);
    }

    return misToReturn;
}

// An attempt to create a more efficient matching
set<int> MatchingTools::ComputeCriticalSet(vector<vector<int>> const &adjacencyList)
{
    // each edge in matching is indexed by first vertex, and contains the second vertex
    // initially there are no edges in matching (thus invalid second vertex in matching).
    vector<int> matching(adjacencyList.size()*2, UNMATCHED_VERTEX);

    cout << "Computing critical set..." << endl << flush;

    BiDoubleGraph biDouble(adjacencyList);
    biDouble.ComputeMaximumMatching(matching);

    // TODO/DS: Make more efficient. Stop search when after marking a vertex in each matchings.
    // i.e., number of vertices marked = number of matchings.
    vector<LastEdge> vOnAlternatingPath(matching.size(), MatchingTools::LastEdge::NO_LAST_EDGE);

    ComputeAlternatingPathsOptimized(biDouble, matching, vOnAlternatingPath);

    set<int> misToReturn;

    misToReturn = std::move(GetLeftVerticesOnAlternatingPaths(biDouble, matching, vOnAlternatingPath));

    int matchVertices(0);
    for (int const vertex : matching) {
        if (vertex != UNMATCHED_VERTEX) {
            matchVertices++;
        }
    }

    matchVertices >>= 1;

    cout << "MIS             Vertices: " << misToReturn.size() << endl << flush;
    cout << "Match           Vertices: " << matchVertices << endl << flush;

    return misToReturn;
}


set<int> MatchingTools::ComputeBiDoubleMIS(vector<vector<int>> const &biDoubleGraph, vector<bool> const &vInGraph, set<int> const &setInGraph)
{
    // each edge in matching is indexed by first vertex, and contains the second vertex
    // initially there are no edges in matching (thus invalid second vertex in matching).
    vector<int> matching(biDoubleGraph.size(), UNMATCHED_VERTEX);

    // residual graph is represented the same as edge in matching.
    // we can iterate over neighbors (excluding neighbor in matching
    // to find residual paths.
////    vector<int> residualGraph(matching);

////    cout << "Computing BiDoubleMIS..." << endl << flush;

    auto inLeftSide = [&biDoubleGraph] (int const vertex) {
        return (vertex < biDoubleGraph.size()/2);
    };

    auto inGraph = [&vInGraph] (int const vertex) {
        return vInGraph[vertex];
    };

    // finds a path from start to finish with excess capacity.
    auto findPath = [&inLeftSide, &biDoubleGraph, &inGraph, &setInGraph] (vector<int> const &matching) {
        vector<bool> inStack(matching.size(), false);
        vector<bool> evaluated(matching.size(), false);
        list<int> stack;

        // depth first search, starting from imaginary start vertex
        // i.e., search starts from left-side.
        for (int const vertex : setInGraph) {
////        for (size_t index = matching.size()/2; index > 0; --index) {
            // only insert vertices without edges in matching, otherwise
            // imaginary first vertex has residual capacity 0 to that vertex.
            if (matching[vertex] != UNMATCHED_VERTEX || !inLeftSide(vertex)) continue;
            stack.push_back(vertex);
            inStack[vertex] = true;
        }

        vector<int> vPreviousVertexOnPath(biDoubleGraph.size(), UNMATCHED_VERTEX);
        int endVertex(UNMATCHED_VERTEX);

        bool foundPath(false);
        while (!stack.empty() && !foundPath) {
            int const vertex = stack.back(); stack.pop_back();
            evaluated[vertex] = true;
            inStack[vertex] = false;
            for (int const neighbor : biDoubleGraph[vertex]) {
                // evaluate neighbor if the edge to that neighbor has residual capacity.
                if (!inGraph(neighbor) || evaluated[neighbor] || inStack[neighbor]) continue;

                // forward edge with residual capacity
                if (inLeftSide(vertex) && matching[vertex] != neighbor) {
                    vPreviousVertexOnPath[neighbor] = vertex;
                    stack.push_back(neighbor);
                    inStack[neighbor] = true;

                    if (!inLeftSide(neighbor) && matching[neighbor] == UNMATCHED_VERTEX) { //found path
                        foundPath = true;
                        endVertex = neighbor;
                        break;
                    }
                }
                // backward edge that we can "undo" by pushing flow back...
                else if (inLeftSide(neighbor) && matching[neighbor] == vertex) {
////                else if (!inLeftSide(vertex) && matching[neighbor] == vertex) {
                    vPreviousVertexOnPath[neighbor] = vertex;
                    stack.push_back(neighbor);
                    inStack[neighbor] = true;
                }
            }
        }

        vector<int> vPath;
        if (endVertex == UNMATCHED_VERTEX) return vPath;
        vPath.push_back(endVertex);
        while (vPreviousVertexOnPath[endVertex] != UNMATCHED_VERTEX) {
            vPath.push_back(vPreviousVertexOnPath[endVertex]);
            endVertex = vPreviousVertexOnPath[endVertex];
        }

        std::reverse(vPath.begin(), vPath.end());

////        cout << "Path through residual graph: ";
////        for (int const vertex : vPath) {
////            cout << vertex << " ";
////        }
////        cout << endl;
        return vPath;
    };

    vector<int> path;
    path = findPath(matching);
    size_t iterations(0);
    while (!path.empty()) {
        iterations++;
////        cout << "Found path of length " << path.size() << ", iteration: " << iterations << endl << flush;
        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const backwardEdge(vertex1 > vertex2);
            if (backwardEdge) {
                matching[vertex1] = UNMATCHED_VERTEX;
                matching[vertex2] = UNMATCHED_VERTEX;
////                cout << "Remove backward edge " << vertex1 << "->" << vertex2 << " from matching" << endl << flush;
            }
////            if (matching[vertex1] != UNMATCHED_VERTEX) {
////                matching[matching[vertex1]] = UNMATCHED_VERTEX;
////            }
////            if (matching[vertex2] != UNMATCHED_VERTEX) {
////                matching[matching[vertex2]] = UNMATCHED_VERTEX;
////            }
////            matching[vertex1] = UNMATCHED_VERTEX;
////            matching[vertex2] = UNMATCHED_VERTEX;
        }

        for (size_t index = 1; index < path.size(); ++index) {
            int const vertex1(path[index-1]);
            int const vertex2(path[index]);
            bool const forwardEdge(vertex1 < vertex2);
            if (forwardEdge) {
                matching[vertex1] = vertex2;
                matching[vertex2] = vertex1;

////                cout << "Add    forward  edge " << vertex1 << "->" << vertex2 << " to   matching" << endl << flush;
            }
        }

        path = findPath(matching);
    }

    // From each unmatched vertex, calculate the alternating paths.
    // -- Paths that alternate matched edges and non-matched edges.
    // Mark all vertices that are connected by an alternating path.

#ifdef VERIFY
    bool bVerified(true);
    for (int const vertex : setInGraph) {
        if (matching[vertex] != UNMATCHED_VERTEX) {
            if (matching[matching[vertex]] != vertex) {
                cout << "ERROR! mismatch: " << vertex << " -> " << matching[vertex] << " -> " << matching[matching[vertex]] << endl << flush;
                bVerified = false;
                break;
            }
        }
    }

    cout << "Verification of matching " << (bVerified ? "passed!" : "failed!") << endl << flush;
#endif // VERIFY

    enum LastEdge {NO_LAST_EDGE, MATCHED_LAST_EDGE, UNMATCHED_LAST_EDGE, BOTH_LAST_EDGE, NULL_LAST_EDGE};

    // TODO/DS: Make more efficient. Stop search when after marking a vertex in each matchings.
    // i.e., number of vertices marked = number of matchings.
    vector<LastEdge> marked(matching.size(), NO_LAST_EDGE);

    // first, mark unmatched vertices, and their neighbors.
    vector<int> matchedVertices;

    vector<int> previousVertex(matching.size(), -1);

    auto isMatchedEdge = [&matching] (int const vertex, int const neighbor) {
        return matching[neighbor] == vertex; 
    };

    int coverVertexCount(0);
    //Mark all unmatched vertices on left-hand side
    //mark their matched neighbors.
    for (int const vertex : setInGraph) {
        if (!inLeftSide(vertex)) continue;
        if (matching[vertex] == UNMATCHED_VERTEX) {
////            cout << "Mark " << vertex << " as unmatched!" << endl << flush;
            marked[vertex] = LastEdge::NULL_LAST_EDGE;
            coverVertexCount++;
            for (int const neighbor : biDoubleGraph[vertex]) {
                if (!inGraph(neighbor)) continue;
                if (marked[neighbor] == LastEdge::NO_LAST_EDGE)
                    coverVertexCount++;
                marked[neighbor] = LastEdge::UNMATCHED_LAST_EDGE;
                previousVertex[neighbor] = vertex;
////                if (isMatchedEdge(vertex, neighbor)) {
////                    cout << "ERROR! (" << __LINE__ << ")" << endl << flush;
////                }
            }
        } else {
           matchedVertices.push_back(vertex); 
        }
    }

////    cout << "Matching       has " << matchedVertices.size() << " vertices" << endl << flush;
////    cout << "Cover count     is " << coverVertexCount                      << endl << flush;


    int iteration(0);
    while (iteration < setInGraph.size()) {
        int const savedCoverCount(coverVertexCount);
        //cout << "Iteration: " << iteration << "/" << matching.size() << endl << flush;
        iteration++;
////        for (int const matchedVertex : matchedVertices) {
        for (int const matchedVertex : setInGraph) {
            if (marked[matchedVertex] == NO_LAST_EDGE) continue;

            LastEdge const edgeIntoVertex(marked[matchedVertex]);

            // If the last edge was matched, follow unmatched edges
            if (edgeIntoVertex == LastEdge::MATCHED_LAST_EDGE) {
////                if (!inLeftSide(matchedVertex))
////                    cout << "ERROR! (" << __LINE__ << ") Vertex " << matchedVertex << " has   matched last edge " << endl << flush;
                for (int const neighbor : biDoubleGraph[matchedVertex]) {
                    if (!inGraph(neighbor)) continue;
                    if (!isMatchedEdge(matchedVertex, neighbor)) {
                        if (marked[neighbor] == LastEdge::NO_LAST_EDGE) {
                            marked[neighbor] = LastEdge::UNMATCHED_LAST_EDGE;
                            coverVertexCount++;
                        }
                    }
////                    if (marked[neighbor] == LastEdge::MATCHED_LAST_EDGE)
////                        cout << "ERROR! Vertex " << neighbor << " on " << (inLeftSide(neighbor)? "left": "right") << " has   matched last edge." << endl << flush;
                }

            // If last edge was unmatched, follow matched edges.
            } else if (edgeIntoVertex == LastEdge::UNMATCHED_LAST_EDGE) {
////                if (inLeftSide(matchedVertex))
////                    cout << "ERROR! (" << __LINE__ << ") Vertex " << matchedVertex << " has unmatched last edge " << endl << flush;
                for (int const neighbor : biDoubleGraph[matchedVertex]) {
                    if (!inGraph(neighbor)) continue;
                    if (isMatchedEdge(matchedVertex, neighbor)) {
                        if (marked[neighbor] == NO_LAST_EDGE) {
                            marked[neighbor] = LastEdge::MATCHED_LAST_EDGE;
                            coverVertexCount++;
                        }
////                        else if (marked[neighbor] == LastEdge::UNMATCHED_LAST_EDGE)
////                            cout << "ERROR! Vertex " << neighbor << " on " << (inLeftSide(neighbor)? "left": "right") << " has unmatched last edge." << endl << flush;
                    }
                }
            }
        }

////        if (savedCoverCount != coverVertexCount) {
////            cout << "Cover increased to " << coverVertexCount << endl << flush;
////        }
        // TODO/DS: doesn't work, I'm not sure why...
////        else {
////            break;
////        }
    }

    int matchVertices(0);
    for (int const vertex : setInGraph) {
        if (matching[vertex] != UNMATCHED_VERTEX) {
            matchVertices++;
        }
    }

    matchVertices >>= 1;

    vector<bool> mvc(matching.size(), false);
    set<int> misToReturn;

    size_t numMVCVertices(0);
    size_t numTotalMISVertices(0);

    // put vertices in mvc
    for (int const vertex : setInGraph) {
        if ((inLeftSide(vertex) && (marked[vertex] == LastEdge::NO_LAST_EDGE)) ||
                (!inLeftSide(vertex) && (marked[vertex] != LastEdge::NO_LAST_EDGE))) {
#ifdef VERIFY
            if (matching[vertex] == UNMATCHED_VERTEX) {
                cout << "ERROR! vertex " << vertex << " is not in matching!" << endl << flush;
            }
#endif //VERIFY
            mvc[vertex] = true;
            numMVCVertices++;

            if (inLeftSide(vertex) && marked[vertex] == LastEdge::UNMATCHED_LAST_EDGE) {
                cout << "ERROR!" << endl << flush;
            }
            if (!inLeftSide(vertex) && marked[vertex] == LastEdge::MATCHED_LAST_EDGE) {
                cout << "ERROR!" << endl << flush;
            }
        }
        else {
            numTotalMISVertices++;
            if (inLeftSide(vertex)) {
                misToReturn.insert(vertex);
            }

            // TODO/DS: I don't know if this belongs or not...
            else {
                misToReturn.insert(vertex);
            }
        }
    }

#ifdef VERIFY
    for (int vertex=0; vertex < mvc.size(); ++vertex) {
        if (!inGraph(vertex)) continue;
        for (int const neighbor : biDoubleGraph[vertex]) {
            if (!inGraph(neighbor)) continue;
            if (!mvc[neighbor] && !mvc[vertex]) {
                cout << "ERROR! MVC is not a vertex cover!, edge (" << vertex << "," << neighbor << ") is not covered." << endl << flush;
            }
        }
    }

    for (int const vertex : misToReturn) {
        if (!inGraph(vertex)) {
            cout << "ERROR! vertex " << vertex << " in MIS is not in graph!" << endl;
        }
        for (int const neighbor : biDoubleGraph[vertex]) {
            if (misToReturn.find(neighbor) != misToReturn.end()) {
                cout << "ERROR! Failed independent set validation." << endl << flush;
            }
        }
    }
#endif // VERIFY


////    cout << "Graph Vertices: " << setInGraph.size() << endl << flush;          
////    cout << "Match Vertices: " << matchVertices << endl << flush;          
////    cout << "MVC   Vertices: " << numMVCVertices << endl << flush;
////    cout << "MIS   Vertices: " << misToReturn.size() << endl << flush;

    if (numMVCVertices + numTotalMISVertices != setInGraph.size()) {
        cout << "ERROR! MVC + MIS != Graph" << endl << flush;
    }

    if (numMVCVertices != matchVertices) {
        cout << "ERROR! MVC != size of matching" << endl << flush;
    }

    return misToReturn;
}

set<int> MatchingTools::ComputeLeftMIS(vector<vector<int>> const &biDoubleGraph, vector<bool> const &vInGraph, set<int> const &setInGraph)
{
    set<int> mis(MatchingTools::ComputeBiDoubleMIS(biDoubleGraph, vInGraph, setInGraph));
    for (set<int>::iterator it = mis.begin(); it != mis.end(); ++it) {
        if (*it >= biDoubleGraph.size()/2) {
            mis.erase(it, mis.end());
            break;
        }
    }

    return mis;
}

void MatchingTools::ComputeAlternatingPathsOptimized(BiDoubleGraph const &biDouble, vector<int> const &vMatching, vector<MatchingTools::LastEdge> &vOnAlternatingPath, vector<bool> const &vInGraph, set<int> const &setInGraph)
{
    auto isMatchedEdge = [&vMatching] (int const vertex, int const neighbor) {
        return vMatching[neighbor] == vertex; 
    };

    auto inGraph = [&vInGraph] (int const vertex) {
        return vInGraph[vertex];
    };

    vector<bool> vEvaluated(vMatching.size(), false);
    vector<int> vStack;
    vStack.reserve(vMatching.size());

    // First mark all unmatched vertices on left-hand side and their neighbors.
////    for (int vertex = 0; vertex < vMatching.size()/2; ++vertex) {
    for (int const vertex : setInGraph) {
        if (vMatching[vertex] != UNMATCHED_VERTEX || !biDouble.InLeftSide(vertex)) continue;
        vEvaluated[vertex] = true;
        vOnAlternatingPath[vertex] = MatchingTools::LastEdge::NULL_LAST_EDGE;
        for (int const neighbor : biDouble.Neighbors(vertex)) {
            if (!inGraph(neighbor)) continue;
            vOnAlternatingPath[neighbor] = MatchingTools::LastEdge::UNMATCHED_LAST_EDGE;
            if (!vEvaluated[neighbor]) vStack.push_back(neighbor);
            vEvaluated[neighbor] = true;
        }
    }

    while (!vStack.empty()) {
        int const matchedVertex(vStack.back()); vStack.pop_back();

        LastEdge const edgeIntoVertex(vOnAlternatingPath[matchedVertex]);

        // If the last edge was matched, follow unmatched edges
        if (edgeIntoVertex == MatchingTools::LastEdge::MATCHED_LAST_EDGE) {
            for (int const neighbor : biDouble.Neighbors(matchedVertex)) {
                if (!inGraph(neighbor)) continue;
                if (!isMatchedEdge(matchedVertex, neighbor)) {
                    if (vOnAlternatingPath[neighbor] == MatchingTools::LastEdge::NO_LAST_EDGE) {
                        vOnAlternatingPath[neighbor] = MatchingTools::LastEdge::UNMATCHED_LAST_EDGE;

                        if (!vEvaluated[neighbor]) {
                            vStack.push_back(neighbor);
                            vEvaluated[neighbor] = true;
                        }
                    }
                }
            }

        // If last edge was unmatched, follow matched edge.
        } else if (edgeIntoVertex == MatchingTools::LastEdge::UNMATCHED_LAST_EDGE) {
            int const neighbor(vMatching[matchedVertex]);
            if (neighbor != UNMATCHED_VERTEX) {
                if (!inGraph(neighbor)) {cout << "ERROR! line: " << __LINE__ << endl << flush;}
                if (vOnAlternatingPath[neighbor] == NO_LAST_EDGE) {
                    vOnAlternatingPath[neighbor] = MatchingTools::LastEdge::MATCHED_LAST_EDGE;

                    if (!vEvaluated[neighbor]) {
                        vStack.push_back(neighbor);
                        vEvaluated[neighbor] = true;
                    }
                }
            }
        }
    }
}

set<int> MatchingTools::GetVerticesOnAlternatingPaths(BiDoubleGraph const &biDouble, vector<int> const &vMatching, vector<MatchingTools::LastEdge> &vOnAlternatingPath, set<int> const &setInGraph)
{
    set<int> misToReturn;
////    for (int vertex = 0; vertex < vMatching.size()/2; ++vertex) {
////        if (vOnAlternatingPath[vertex] != MatchingTools::LastEdge::NO_LAST_EDGE) misToReturn.insert(vertex);
    for (int const vertex : setInGraph) {
        if ((biDouble.InLeftSide(vertex) && (vOnAlternatingPath[vertex] == MatchingTools::LastEdge::NO_LAST_EDGE)) ||
                (!biDouble.InLeftSide(vertex) && (vOnAlternatingPath[vertex] != MatchingTools::LastEdge::NO_LAST_EDGE))) {
        } else {
            misToReturn.insert(vertex);
        }
    }

    return misToReturn;
}

set<int> MatchingTools::ComputeBiDoubleMISOptimized(BiDoubleGraph &biDouble, vector<bool> const &vInGraph, set<int> const &setInGraph)
{
    // each edge in matching is indexed by first vertex, and contains the second vertex
    // initially there are no edges in matching (thus invalid second vertex in matching).
    vector<int> matching(biDouble.Size(), UNMATCHED_VERTEX);

    // residual graph is represented the same as edge in matching.
    // we can iterate over neighbors (excluding neighbor in matching
    // to find residual paths.

////    cout << "Computing BiDoubleMIS..." << endl << flush;

    auto inLeftSide = [&biDouble] (int const vertex) {
        return biDouble.InLeftSide(vertex);
    };

    auto inGraph = [&vInGraph] (int const vertex) {
        return vInGraph[vertex];
    };

    biDouble.ComputeMaximumMatching(matching, vInGraph, setInGraph);

#ifdef VERIFY
    bool bVerified(true);
    for (int const vertex : setInGraph) {
        if (matching[vertex] != UNMATCHED_VERTEX) {
            if (matching[matching[vertex]] != vertex) {
                cout << "ERROR! mismatch: " << vertex << " -> " << matching[vertex] << " -> " << matching[matching[vertex]] << endl << flush;
                bVerified = false;
                break;
            }
        }
    }

    cout << "Verification of matching " << (bVerified ? "passed!" : "failed!") << endl << flush;
#endif // VERIFY

    vector<MatchingTools::LastEdge> marked(matching.size(), MatchingTools::LastEdge::NO_LAST_EDGE);

    // TODO/DS combine last two steps into one. -> more efficient.
    MatchingTools::ComputeAlternatingPathsOptimized(biDouble, matching, marked, vInGraph, setInGraph);

    set<int> const misToReturn(std::move(MatchingTools::GetVerticesOnAlternatingPaths(biDouble, matching, marked, setInGraph)));

    return misToReturn;
}

set<int> MatchingTools::ComputeBiDoubleMISOptimizedWithMatching(BiDoubleGraph &biDouble, vector<int> &matching, vector<bool> const &vInGraph, set<int> const &setInGraph)
{
    // each edge in matching is indexed by first vertex, and contains the second vertex
    // initially there are no edges in matching (thus invalid second vertex in matching).
////    vector<int> matching(biDouble.Size(), UNMATCHED_VERTEX);

    // residual graph is represented the same as edge in matching.
    // we can iterate over neighbors (excluding neighbor in matching
    // to find residual paths.

////    cout << "Computing BiDoubleMIS..." << endl << flush;

    auto inLeftSide = [&biDouble] (int const vertex) {
        return biDouble.InLeftSide(vertex);
    };

    auto inGraph = [&vInGraph] (int const vertex) {
        return vInGraph[vertex];
    };

    biDouble.ComputeMaximumMatching(matching, vInGraph, setInGraph);

#ifdef VERIFY
    bool bVerified(true);
    for (int const vertex : setInGraph) {
        if (matching[vertex] != UNMATCHED_VERTEX) {
            if (matching[matching[vertex]] != vertex) {
                cout << "ERROR! mismatch: " << vertex << " -> " << matching[vertex] << " -> " << matching[matching[vertex]] << endl << flush;
                bVerified = false;
                break;
            }
        }
    }

    cout << "Verification of matching " << (bVerified ? "passed!" : "failed!") << endl << flush;
#endif // VERIFY

    vector<MatchingTools::LastEdge> marked(matching.size(), MatchingTools::LastEdge::NO_LAST_EDGE);

    // TODO/DS combine last two steps into one. -> more efficient.
    MatchingTools::ComputeAlternatingPathsOptimized(biDouble, matching, marked, vInGraph, setInGraph);

    set<int> const misToReturn(std::move(MatchingTools::GetVerticesOnAlternatingPaths(biDouble, matching, marked, setInGraph)));

    return misToReturn;
}

set<int> MatchingTools::ComputeLeftMISOptimizedWithMatching(BiDoubleGraph &biDouble, vector<int> &matching, vector<bool> const &vInGraph, set<int> const &setInGraph)
{
    set<int> mis(MatchingTools::ComputeBiDoubleMISOptimizedWithMatching(biDouble, matching, vInGraph, setInGraph));
    for (set<int>::iterator it = mis.begin(); it != mis.end(); ++it) {
        if (*it >= biDouble.Size()/2) {
            mis.erase(it, mis.end());
            break;
        }
    }

    return mis;
}
