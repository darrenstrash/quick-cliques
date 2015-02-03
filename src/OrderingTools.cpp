#include "OrderingTools.h"
#include "DegeneracyTools.h"
#include "GraphTools.h"
#include "CliqueColoringStrategy.h"
#include "SparseIndependentSetColoringStrategy.h"

#include <vector>
#include <iostream>

using namespace std;

void OrderingTools::InitialOrderingMCQ(vector<vector<char>> const &adjacencyMatrix, vector<int> &vOrderedVertices, vector<int> &vColoring)
{
    size_t maxDegree(0);
    {
        vector<int> vDegree(adjacencyMatrix.size(), 0);
        for (size_t u = 0; u < adjacencyMatrix.size(); ++u) {
            for (size_t v = 0; v < adjacencyMatrix.size(); ++v) {
                if (adjacencyMatrix[u][v]) vDegree[u]++;
            }
            maxDegree = max(maxDegree, static_cast<size_t>(vDegree[u]));
        }
        vOrderedVertices = std::move(GraphTools::OrderVerticesByDegree(adjacencyMatrix, vDegree, false /* non-increasing order */));
    }

    vColoring.reserve(adjacencyMatrix.size());
    vColoring.clear();
    for (int degree = 1; degree <= maxDegree; degree++ ) {
        vColoring.push_back(degree);
    }

    vColoring.resize(adjacencyMatrix.size(), maxDegree + 1);
}

void OrderingTools::InitialOrderingMISQ(vector<vector<char>> const &adjacencyMatrix, vector<int> &vOrderedVertices, vector<int> &vColoring)
{
    size_t maxDegree(0);
    {
        vector<int> vDegree(adjacencyMatrix.size(), 0);
        for (size_t u = 0; u < adjacencyMatrix.size(); ++u) {
            for (size_t v = 0; v < adjacencyMatrix.size(); ++v) {
                if (!adjacencyMatrix[u][v]) vDegree[u]++;
            }
            maxDegree = max(maxDegree, static_cast<size_t>(vDegree[u]));
        }
        vOrderedVertices = std::move(GraphTools::OrderVerticesByDegree(adjacencyMatrix, vDegree, false /* non-increasing order */));
    }

    vColoring.reserve(adjacencyMatrix.size());
    vColoring.clear();
    for (int degree = 1; degree <= maxDegree; degree++ ) {
        vColoring.push_back(degree);
    }

    vColoring.resize(adjacencyMatrix.size(), maxDegree + 1);
}

void OrderingTools::InitialOrderingMCR(vector<vector<char>> const &adjacencyMatrix, vector<int> &vOrderedVertices, vector<int> &vColoring, size_t &cliqueSize)
{
    vOrderedVertices.resize(adjacencyMatrix.size(), -1);
    vColoring.resize(adjacencyMatrix.size(), -1);

    // create an adjacencyArray, much faster for degeneracy ordering.
    vector<vector<int>> adjacencyArray(adjacencyMatrix.size());
    for (int vertex = 0; vertex < adjacencyMatrix.size(); ++vertex) {
        for (int otherVertex = 0; otherVertex < adjacencyMatrix.size(); ++otherVertex) {
            if (adjacencyMatrix[vertex][otherVertex]) {
                adjacencyArray[vertex].push_back(otherVertex);
            }
        }
    }

    size_t const size(adjacencyArray.size());

    // array of lists of vertices, indexed by degree
    vector<list<int>> verticesByDegree(size);

    // array of lists of vertices, indexed by degree
    vector<list<int>::iterator> vertexLocator(size);

    vector<int> degree(size);

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    size_t maxDegree(0);
    for(size_t i = 0; i < size; i++) {
        degree[i] = adjacencyArray[i].size();
        verticesByDegree[degree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[degree[i]].begin();

        maxDegree = max(maxDegree, adjacencyArray[i].size());
    }

    int currentDegree = 0;
    int numVerticesRemoved = 0;

    vector<NeighborListArray> vOrderingArray(size);

    while (numVerticesRemoved < size) {
        if (!verticesByDegree[currentDegree].empty()) {

            int vertex(-1);
            if (verticesByDegree[currentDegree].size() > 1) {

                // if remaining graph is regular.
                if ((size - numVerticesRemoved) == verticesByDegree[currentDegree].size()) {
                    cout << "Remaining graph of " << verticesByDegree[currentDegree].size() << " vertices is regular with degree " << currentDegree << endl;
                    cout << "Vertices: ";
                    for (int const lastVertex : verticesByDegree[currentDegree]) {
////                        if (vColoring.size() - index < 10)
                            cout << lastVertex << " ";
                    }
                    cout << endl;


                    // if regular, and degree is # vertices - 1, then it's a clique.
                    if (static_cast<int>(verticesByDegree[currentDegree].size()) == currentDegree + 1) {
                        cliqueSize = currentDegree + 1;
                    }

                    vector<int> remainingVertices(verticesByDegree[currentDegree].begin(), verticesByDegree[currentDegree].end());
                    vector<int> remainingColors(remainingVertices.size(), 0);
                    CliqueColoringStrategy coloringStrategy(adjacencyMatrix);
                    coloringStrategy.Color(adjacencyMatrix, remainingVertices /* evaluation order */, remainingVertices /* color order */, remainingColors);
                    //copy initial ordering to output arrays

////                    int maxColor(0);
////                    size_t index(0);
                    for (size_t index = 0; index < remainingVertices.size(); ++index) {
                        vOrderedVertices[index] = remainingVertices[index];
                        vColoring[index] = remainingColors[index];
////                        maxColor = max(maxColor, vColoring[index]);
                    }

                    // simpler
#if 1
                    for (size_t index = remainingColors.size(); index < vColoring.size(); ++index) {
                        vColoring[index] = min(vColoring[index-1] + 1, static_cast<int>(maxDegree + 1));
                    }

                    cout << "All vertices: ";
                    for (size_t index = 0; index < vOrderedVertices.size(); ++index) {
////                        if (vColoring.size() - index < 10)
                            cout << vOrderedVertices[index] << " ";
                    }
                    cout << endl;
                    cout << "All colors: ";
                    for (size_t index = 0; index < vColoring.size(); ++index) {
////                        if (vColoring.size() - index < 10)
                            cout << vColoring[index] << " ";
                    }
                    cout << endl;
#else
                    
                    int const lastIndexWithSmallerColor(min(remainingVertices.size() + maxDegree - maxColor, size-1));
                    int currentColor = maxColor + 1;
                    while (index <= lastIndexWithSmallerColor) {
                        vColoring[index] = currentColor;
                        currentColor++;
                        index++;
                    }
                    while (index < size) {
                        vColoring[index] = maxDegree + 1;
                        index++;
                    }
#endif //0

                    return;
                } else {
                    // break ties by neighborhood-degree
                    size_t minNeighborhoodDegree(ULONG_MAX);
                    int    chosenVertex=verticesByDegree[currentDegree].front();
                    for (int const candidate : verticesByDegree[currentDegree]) {
                        size_t neighborhoodDegree(0);
                        for (int const neighbor : adjacencyArray[candidate]) {
                            if (degree[neighbor] != -1) {
                                neighborhoodDegree += degree[neighbor];
                            }
                        }

                        if (neighborhoodDegree < minNeighborhoodDegree || (neighborhoodDegree == minNeighborhoodDegree && candidate > chosenVertex)) {
                            minNeighborhoodDegree = neighborhoodDegree;
                            chosenVertex = candidate;
                        }
                    }
                    vertex = chosenVertex;
                    verticesByDegree[currentDegree].erase(vertexLocator[vertex]);
                    cout << "Choosing multi vertex: " << vertex << " with degree " << degree[vertex] << " and neighborhood " << minNeighborhoodDegree << endl << flush;
                }
            } else {
                vertex = verticesByDegree[currentDegree].front();
                verticesByDegree[currentDegree].pop_front();
                cout << "Choosing solo  vertex: " << vertex << " with degree " << degree[vertex] << endl << flush;
            }

            vOrderedVertices[vOrderedVertices.size() - numVerticesRemoved - 1] = vertex;

            degree[vertex] = -1;

            for (int const neighbor : adjacencyArray[vertex]) {
                if (degree[neighbor] > -1)
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    degree[neighbor]--;
                    verticesByDegree[degree[neighbor]].push_front(neighbor);
                    vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                }
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else {
            currentDegree++;
        }
    }
}

void OrderingTools::InitialOrderingMISR(vector<vector<int>> const &adjacencyArray, vector<int> &vOrderedVertices, vector<int> &vColoring, size_t &cliqueSize)
{
    vOrderedVertices.resize(adjacencyArray.size(), -1);
    vColoring.resize(adjacencyArray.size(), -1);

    vector<bool> vMarkedVertices(adjacencyArray.size());

    size_t const size(adjacencyArray.size());

    // array of lists of vertices, indexed by degree
    vector<list<int>> verticesByDegree(size);

    // array of lists of vertices, indexed by degree
    vector<list<int>::iterator> vertexLocator(size);

    vector<int> degree(size);

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    size_t maxDegree(0);
    size_t maxCoDegree(0);
    for(size_t i = 0; i < size; i++) {
        degree[i] = adjacencyArray[i].size();
        verticesByDegree[degree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[degree[i]].begin();

        maxDegree = max(maxDegree, adjacencyArray[i].size());
        maxCoDegree = max(maxCoDegree, size - adjacencyArray[i].size() - 1);
    }

    // perform degeneracy ordering in complement graph.

    int currentDegree = maxDegree;
    int numVerticesRemoved = 0;

    vector<NeighborListArray> vOrderingArray(size);

    while (numVerticesRemoved < size) {
        if (!verticesByDegree[currentDegree].empty()) {

            int vertex(-1);
            if (verticesByDegree[currentDegree].size() > 1) {

                // if remaining graph is regular.
                if (currentDegree == 0) { ////(size - numVerticesRemoved) == verticesByDegree[currentDegree].size()) {
                    cout << "Remaining graph of " << verticesByDegree[currentDegree].size() << " vertices is regular with degree " << (size - numVerticesRemoved - currentDegree - 1) << endl;
                    cout << "Vertices: ";
                    for (int const lastVertex : verticesByDegree[currentDegree]) {
////                        if (vColoring.size() - index < 10)
                            cout << lastVertex << " ";
                    }
                    cout << endl;

                    // if regular, and degree is # vertices - 1, then it's a clique.
                    if (currentDegree == 0) { /////static_cast<int>(verticesByDegree[currentDegree].size()) == (size - (currentDegree + 1))) { // TODO/DS put in the proper check for independent set...
                        cliqueSize = verticesByDegree[currentDegree].size();
                    }

                    vector<int> remainingVertices(verticesByDegree[currentDegree].begin(), verticesByDegree[currentDegree].end());
                    vector<int> remainingColors(remainingVertices.size(), 0);
////                    CliqueColoringStrategy coloringStrategy(adjacencyMatrix);
                    SparseIndependentSetColoringStrategy coloringStrategy(adjacencyArray); // TODO/DS: Validate.
                    coloringStrategy.Color(adjacencyArray, remainingVertices /* evaluation order */, remainingVertices /* color order */, remainingColors);
                    //copy initial ordering to output arrays

////                    int maxColor(0);
////                    size_t index(0);
                    for (size_t index = 0; index < remainingVertices.size(); ++index) {
                        vOrderedVertices[index] = remainingVertices[index];
                        vColoring[index] = remainingColors[index];
////                        maxColor = max(maxColor, vColoring[index]);
                    }

#if 1
                    // simpler
                    for (size_t index = remainingColors.size(); index < vColoring.size(); ++index) {
                        vColoring[index] = min(vColoring[index-1] + 1, static_cast<int>(maxCoDegree + 1));
                    }

                    cout << "All vertices: ";
                    for (size_t index = 0; index < vOrderedVertices.size(); ++index) {
////                        if (vColoring.size() - index < 10)
                            cout << vOrderedVertices[index] << " ";
                    }
                    cout << endl;
                    cout << "All colors: ";
                    for (size_t index = 0; index < vColoring.size(); ++index) {
////                        if (vColoring.size() - index < 10)
                            cout << vColoring[index] << " ";
                    }
                    cout << endl;
#else
                    int const lastIndexWithSmallerColor(min(remainingVertices.size() + maxDegree - maxColor, size-1));
                    int currentColor = maxColor + 1;
                    while (index <= lastIndexWithSmallerColor) {
                        vColoring[index] = currentColor;
                        currentColor++;
                        index++;
                    }
                    while (index < size) {
                        vColoring[index] = maxDegree + 1;
                        index++;
                    }
#endif //0

                    return;
                } else {
                    // break ties by neighborhood-degree
                    size_t minNeighborhoodDegree(ULONG_MAX);
////                    cout << "currentDegree= " << currentDegree << endl << flush;
                    int    chosenVertex=verticesByDegree[currentDegree].front();
                    for (int const candidate : verticesByDegree[currentDegree]) {
                        size_t coNeighborhoodDegree(0);
                        for (int const neighbor : adjacencyArray[candidate]) {
                            vMarkedVertices[neighbor] = true;
                        }

                        // iterate over all non-neighbors, and add up their non neighbors
                        for (int nonNeighbor = 0; nonNeighbor < adjacencyArray.size(); ++nonNeighbor) {
                            if (degree[nonNeighbor] != -1 && !vMarkedVertices[nonNeighbor]) {
                                coNeighborhoodDegree += (size - numVerticesRemoved - degree[nonNeighbor] - 1);
                            }
                        }

                        for (int const neighbor : adjacencyArray[candidate]) {
                            vMarkedVertices[neighbor] = false;
                        }

                        if (coNeighborhoodDegree < minNeighborhoodDegree || (coNeighborhoodDegree == minNeighborhoodDegree && candidate > chosenVertex)) {
                            minNeighborhoodDegree = coNeighborhoodDegree;
                            chosenVertex = candidate;
                        }
                    }
                    vertex = chosenVertex;
                    verticesByDegree[currentDegree].erase(vertexLocator[vertex]);
                    cout << "Choosing multi vertex: " << vertex << " with degree " << (size - numVerticesRemoved - degree[vertex] - 1) << " and neighborhood " << minNeighborhoodDegree << endl << flush;
                }
            } else {
                vertex = verticesByDegree[currentDegree].front();
                verticesByDegree[currentDegree].pop_front();
                cout << "Choosing solo  vertex: " << vertex << " with degree " << (size - numVerticesRemoved - degree[vertex] - 1) << endl << flush;
            }

            vOrderedVertices[vOrderedVertices.size() - numVerticesRemoved - 1] = vertex;

            degree[vertex] = -1;

            for (int const neighbor : adjacencyArray[vertex]) {
                if (degree[neighbor] > -1)
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    degree[neighbor]--;
                    verticesByDegree[degree[neighbor]].push_front(neighbor);
                    vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                }
            }

            numVerticesRemoved++;
            ////currentDegree = numVertices-1; // degrees can't grow...
        }
        else {
            currentDegree--;
        }
    }
}

void OrderingTools::InitialOrderingMISR(vector<vector<char>> const &adjacencyMatrix, vector<int> &vOrderedVertices, vector<int> &vColoring, size_t &cliqueSize)
{
#if 0
    // create an adjacencyArray, much faster for degeneracy ordering.

    vector<vector<int>> adjacencyArray(adjacencyMatrix.size());
    for (int vertex = 0; vertex < adjacencyMatrix.size(); ++vertex) {
        for (int otherVertex = 0; otherVertex < adjacencyMatrix.size(); ++otherVertex) {
            if (vertex == otherVertex) continue;
            if (!adjacencyMatrix[vertex][otherVertex]) {
                adjacencyArray[vertex].push_back(otherVertex);
            }
        }
    }

    OrderingTools::InitialOrderingMISR(adjacencyArray, vOrderedVertices, vColoring, cliqueSize);
#else
    vector<vector<char>> coAdjacencyMatrix(adjacencyMatrix.size());
    for (int u = 0; u < adjacencyMatrix.size(); ++u) {
        coAdjacencyMatrix[u].resize(adjacencyMatrix.size());
        for (int v = 0; v < adjacencyMatrix.size(); ++v) {
            if (u==v) continue;
            coAdjacencyMatrix[u][v] = !adjacencyMatrix[u][v];
        }
    }

    OrderingTools::InitialOrderingMCR(coAdjacencyMatrix, vOrderedVertices, vColoring, cliqueSize);
#endif //0
}
