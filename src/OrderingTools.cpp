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

void OrderingTools::InitialOrderingMCQ(vector<vector<int>> const &adjacencyArray, vector<int> &vOrderedVertices, vector<int> &vColoring)
{
    size_t maxDegree(0);
    for (vector<int> const &vNeighbor : adjacencyArray) {
        maxDegree = max(maxDegree, vNeighbor.size());
    }

    vOrderedVertices = std::move(GraphTools::OrderVerticesByDegree(adjacencyArray, false /* non-increasing order */));
    vColoring.reserve(adjacencyArray.size());
    vColoring.clear();
    for (int degree = 1; degree <= maxDegree; degree++ ) {
        vColoring.push_back(degree);
    }

    vColoring.resize(adjacencyArray.size(), maxDegree + 1);
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
        vertexLocator[i] = verticesByDegree[degree[i]].insert(verticesByDegree[degree[i]].end(), i);

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
////                    cout << "Remaining graph of " << verticesByDegree[currentDegree].size() << " vertices is regular with degree " << currentDegree << endl;
////                    cout << "Vertices: ";
////                    for (int const lastVertex : verticesByDegree[currentDegree]) {
////////                        if (vColoring.size() - index < 10)
////                            cout << lastVertex << " ";
////                    }
////                    cout << endl;


                    // if regular, and degree is # vertices - 1, then it's a clique.
                    if (static_cast<int>(verticesByDegree[currentDegree].size()) == currentDegree + 1) {
                        cliqueSize = currentDegree + 1;
                    }

                    vector<int> remainingVertices(verticesByDegree[currentDegree].begin(), verticesByDegree[currentDegree].end());
                    vector<int> remainingColors(remainingVertices.size(), 0);
                    CliqueColoringStrategy coloringStrategy(adjacencyMatrix);
                    coloringStrategy.Color(adjacencyMatrix, remainingVertices /* evaluation order */, remainingVertices /* color order */, remainingColors);
                    //copy initial ordering to output arrays

                    int maxColor(0);
                    size_t index(0);
                    for (/*size_t index = 0*/; index < remainingVertices.size(); ++index) {
                        vOrderedVertices[index] = remainingVertices[index];
                        vColoring[index] = remainingColors[index];
                        maxColor = max(maxColor, vColoring[index]);
                    }

                    // simpler
#if 0
                    for (size_t index = remainingColors.size(); index < vColoring.size(); ++index) {
                        vColoring[index] = min(vColoring[index-1] + 1, static_cast<int>(maxDegree + 1));
                    }

////                    cout << "All vertices: ";
////                    for (size_t index = 0; index < vOrderedVertices.size(); ++index) {
////////                        if (vColoring.size() - index < 10)
////                            cout << vOrderedVertices[index] << " ";
////                    }
////                    cout << endl;
////                    cout << "All colors: ";
////                    for (size_t index = 0; index < vColoring.size(); ++index) {
////////                        if (vColoring.size() - index < 10)
////                            cout << vColoring[index] << " ";
////                    }
////                    cout << endl;
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

                        if (neighborhoodDegree < minNeighborhoodDegree) { //// || (neighborhoodDegree == minNeighborhoodDegree && candidate < chosenVertex)) {
                            minNeighborhoodDegree = neighborhoodDegree;
                            chosenVertex = candidate;
                        }
                    }
                    vertex = chosenVertex;
                    verticesByDegree[currentDegree].erase(vertexLocator[vertex]);
////                    cout << "Choosing multi vertex: " << vertex << " with degree " << degree[vertex] << " and neighborhood " << minNeighborhoodDegree << endl << flush;
                }
            } else {
                vertex = verticesByDegree[currentDegree].front();
                verticesByDegree[currentDegree].pop_front();
////                cout << "Choosing solo  vertex: " << vertex << " with degree " << degree[vertex] << endl << flush;
            }

            vOrderedVertices[vOrderedVertices.size() - numVerticesRemoved - 1] = vertex;

            degree[vertex] = -1;

            for (int const neighbor : adjacencyArray[vertex]) {
                if (degree[neighbor] > -1)
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    degree[neighbor]--;
                    vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].insert(verticesByDegree[degree[neighbor]].end(), neighbor);
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

#if 0
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
#else
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

    vector<int> coDegree(size);

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    size_t maxDegree(0);
    size_t maxCoDegree(0);
    for(size_t i = 0; i < size; i++) {
        coDegree[i] = size - adjacencyArray[i].size() - 1;
////        vertexLocator[i] = verticesByDegree[coDegree[i]].insert(verticesByDegree[coDegree[i]].end(), i);
        verticesByDegree[coDegree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[coDegree[i]].begin();

        maxDegree = max(maxDegree, adjacencyArray[i].size());
        maxCoDegree = max(maxCoDegree, static_cast<size_t>(coDegree[i]));
    }

    // perform degeneracy ordering in complement graph.

    int currentDegree = 0;
    int numVerticesRemoved = 0;

    vector<NeighborListArray> vOrderingArray(size);

    while (numVerticesRemoved < size) {
////        cout << "Ordered " << numVerticesRemoved << "/" << size << " vertices" << endl;
        if (!verticesByDegree[currentDegree].empty()) {

            int vertex(-1);
            if (verticesByDegree[currentDegree].size() > 1) {

                // if remaining graph is regular.
                if ((size - numVerticesRemoved) == verticesByDegree[currentDegree].size()) {
////                    cout << "Remaining graph of " << verticesByDegree[currentDegree].size() << " vertices is regular with degree " << currentDegree << endl;
////                    cout << "Vertices: ";
////                    for (int const lastVertex : verticesByDegree[currentDegree]) {
////////                        if (vColoring.size() - index < 10)
////                            cout << lastVertex << " ";
////                    }
////                    cout << endl;

                    // if regular, and degree is # vertices - 1, then it's a clique.
                    if (static_cast<int>(verticesByDegree[currentDegree].size()) == (currentDegree + 1)) { // TODO/DS put in the proper check for independent set...
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

////                    cout << "All vertices: ";
////                    for (size_t index = 0; index < vOrderedVertices.size(); ++index) {
////////                        if (vColoring.size() - index < 10)
////                            cout << vOrderedVertices[index] << " ";
////                    }
////                    cout << endl;
////                    cout << "All colors: ";
////                    for (size_t index = 0; index < vColoring.size(); ++index) {
////////                        if (vColoring.size() - index < 10)
////                            cout << vColoring[index] << " ";
////                    }
////                    cout << endl;
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
                            if (coDegree[nonNeighbor] != -1 && !vMarkedVertices[nonNeighbor]) {
                                coNeighborhoodDegree += (coDegree[nonNeighbor]);
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
////                    cout << "Choosing multi vertex: " << vertex << " with degree " << coDegree[vertex]  << " and neighborhood " << minNeighborhoodDegree << endl << flush;
                }
            } else {
                vertex = verticesByDegree[currentDegree].front();
                verticesByDegree[currentDegree].pop_front();
////                cout << "Choosing solo  vertex: " << vertex << " with degree " << coDegree[vertex] << endl << flush;
            }

            vOrderedVertices[vOrderedVertices.size() - numVerticesRemoved - 1] = vertex;

            coDegree[vertex] = -1;

#if 0
            for (int const neighbor : adjacencyArray[vertex]) {
                if (coDegree[neighbor] > -1)
                {
                    verticesByDegree[coDegree[neighbor]].erase(vertexLocator[neighbor]);
                    coDegree[neighbor]--;
                    verticesByDegree[coDegree[neighbor]].push_front(neighbor);
                    vertexLocator[neighbor] = verticesByDegree[coDegree[neighbor]].begin();
                }
            }
#else
            // need to iterate over non-neighbors
            for (int const neighbor : adjacencyArray[vertex]) {
                vMarkedVertices[neighbor] = true;
            }

            // iterate over all non-neighbors, subtract from their co-Degree
            for (int nonNeighbor = 0; nonNeighbor < adjacencyArray.size(); ++nonNeighbor) {
                if (coDegree[nonNeighbor] != -1 && !vMarkedVertices[nonNeighbor]) {
                    verticesByDegree[coDegree[nonNeighbor]].erase(vertexLocator[nonNeighbor]);
                    coDegree[nonNeighbor]--;
////                    vertexLocator[nonNeighbor] = verticesByDegree[coDegree[nonNeighbor]].insert(verticesByDegree[coDegree[nonNeighbor]].end(), nonNeighbor);
                    verticesByDegree[coDegree[nonNeighbor]].push_front(nonNeighbor);
                    vertexLocator[nonNeighbor] = verticesByDegree[coDegree[nonNeighbor]].begin();
                }
            }

            for (int const neighbor : adjacencyArray[vertex]) {
                vMarkedVertices[neighbor] = false;
            }

#endif // 0
            numVerticesRemoved++;
            currentDegree--; // degrees can't grow...
            if (currentDegree == -1) currentDegree = 0;
       } else {
            currentDegree++;
        }
    }
}
#endif //0

void OrderingTools::InitialOrderingMISR(vector<vector<int>> const &adjacencyArray, Isolates4<SparseArraySet> const &isolates, vector<int> &vOrderedVertices, vector<int> &vColoring, size_t &cliqueSize)
{
    size_t const numVertices(isolates.GetInGraph().Size());
    vOrderedVertices.resize(numVertices, -1);
    vColoring.resize(numVertices, -1);

    size_t const size(adjacencyArray.size());
    vector<bool> vMarkedVertices(size);

    // array of lists of vertices, indexed by degree
    vector<list<int>> verticesByDegree(size);

    // array of lists of vertices, indexed by degree
    vector<list<int>::iterator> vertexLocator(size);

    vector<int> coDegree(size);

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    size_t maxDegree(0);
    size_t maxCoDegree(0);
    for(int const i : isolates.GetInGraph()) {
        coDegree[i] = size - isolates.Neighbors()[i].Size() - 1;
////        vertexLocator[i] = verticesByDegree[coDegree[i]].insert(verticesByDegree[coDegree[i]].end(), i);
        verticesByDegree[coDegree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[coDegree[i]].begin();

        maxDegree = max(maxDegree, isolates.Neighbors()[i].Size());
        maxCoDegree = max(maxCoDegree, static_cast<size_t>(coDegree[i]));
    }

    // perform degeneracy ordering in complement graph.

    int currentDegree = 0;
    int numVerticesRemoved = 0;

    vector<NeighborListArray> vOrderingArray(size);

    while (numVerticesRemoved < numVertices) {
////        cout << "Ordered " << numVerticesRemoved << "/" << size << " vertices" << endl;
        if (!verticesByDegree[currentDegree].empty()) {

            int vertex(-1);
            if (verticesByDegree[currentDegree].size() > 1) {

                // if remaining graph is regular.
                if ((numVertices - numVerticesRemoved) == verticesByDegree[currentDegree].size()) {
////                    cout << "Remaining graph of " << verticesByDegree[currentDegree].size() << " vertices is regular with degree " << currentDegree << endl;
////                    cout << "Vertices: ";
////                    for (int const lastVertex : verticesByDegree[currentDegree]) {
////////                        if (vColoring.size() - index < 10)
////                            cout << lastVertex << " ";
////                    }
////                    cout << endl;

                    // if regular, and degree is # vertices - 1, then it's a clique.
                    if (static_cast<int>(verticesByDegree[currentDegree].size()) == (currentDegree + 1)) { // TODO/DS put in the proper check for independent set...
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

////                    cout << "All vertices: ";
////                    for (size_t index = 0; index < vOrderedVertices.size(); ++index) {
////////                        if (vColoring.size() - index < 10)
////                            cout << vOrderedVertices[index] << " ";
////                    }
////                    cout << endl;
////                    cout << "All colors: ";
////                    for (size_t index = 0; index < vColoring.size(); ++index) {
////////                        if (vColoring.size() - index < 10)
////                            cout << vColoring[index] << " ";
////                    }
////                    cout << endl;
#else
                    int const lastIndexWithSmallerColor(min(remainingVertices.size() + maxDegree - maxColor, numVertices-1));
                    int currentColor = maxColor + 1;
                    while (index <= lastIndexWithSmallerColor) {
                        vColoring[index] = currentColor;
                        currentColor++;
                        index++;
                    }
                    while (index < numVertices) {
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
                        for (int const neighbor : isolates.Neighbors()[candidate]) {
                            vMarkedVertices[neighbor] = true;
                        }

                        // iterate over all non-neighbors, and add up their non neighbors
                        for (int nonNeighbor = 0; nonNeighbor < numVertices; ++nonNeighbor) {
                            if (coDegree[nonNeighbor] != -1 && !vMarkedVertices[nonNeighbor]) {
                                coNeighborhoodDegree += (coDegree[nonNeighbor]);
                            }
                        }

                        for (int const neighbor : isolates.Neighbors()[candidate]) {
                            vMarkedVertices[neighbor] = false;
                        }

                        if (coNeighborhoodDegree < minNeighborhoodDegree || (coNeighborhoodDegree == minNeighborhoodDegree && candidate > chosenVertex)) {
                            minNeighborhoodDegree = coNeighborhoodDegree;
                            chosenVertex = candidate;
                        }
                    }
                    vertex = chosenVertex;
                    verticesByDegree[currentDegree].erase(vertexLocator[vertex]);
////                    cout << "Choosing multi vertex: " << vertex << " with degree " << coDegree[vertex]  << " and neighborhood " << minNeighborhoodDegree << endl << flush;
                }
            } else {
                vertex = verticesByDegree[currentDegree].front();
                verticesByDegree[currentDegree].pop_front();
////                cout << "Choosing solo  vertex: " << vertex << " with degree " << coDegree[vertex] << endl << flush;
            }

            vOrderedVertices[vOrderedVertices.size() - numVerticesRemoved - 1] = vertex;

            coDegree[vertex] = -1;

#if 0
            for (int const neighbor : adjacencyArray[vertex]) {
                if (coDegree[neighbor] > -1)
                {
                    verticesByDegree[coDegree[neighbor]].erase(vertexLocator[neighbor]);
                    coDegree[neighbor]--;
                    verticesByDegree[coDegree[neighbor]].push_front(neighbor);
                    vertexLocator[neighbor] = verticesByDegree[coDegree[neighbor]].begin();
                }
            }
#else
            // need to iterate over non-neighbors
            for (int const neighbor : isolates.Neighbors()[vertex]) {
                vMarkedVertices[neighbor] = true;
            }

            // iterate over all non-neighbors, subtract from their co-Degree
            for (int nonNeighbor = 0; nonNeighbor < numVertices; ++nonNeighbor) {
                if (coDegree[nonNeighbor] != -1 && !vMarkedVertices[nonNeighbor]) {
                    verticesByDegree[coDegree[nonNeighbor]].erase(vertexLocator[nonNeighbor]);
                    coDegree[nonNeighbor]--;
////                    vertexLocator[nonNeighbor] = verticesByDegree[coDegree[nonNeighbor]].insert(verticesByDegree[coDegree[nonNeighbor]].end(), nonNeighbor);
                    verticesByDegree[coDegree[nonNeighbor]].push_front(nonNeighbor);
                    vertexLocator[nonNeighbor] = verticesByDegree[coDegree[nonNeighbor]].begin();
                }
            }

            for (int const neighbor : isolates.Neighbors()[vertex]) {
                vMarkedVertices[neighbor] = false;
            }

#endif // 0
            numVerticesRemoved++;
            currentDegree--; // degrees can't grow...
            if (currentDegree == -1) currentDegree = 0;
       } else {
            currentDegree++;
        }
    }
}

void OrderingTools::InitialOrderingMISR(vector<vector<char>> const &adjacencyMatrix, vector<int> &vOrderedVertices, vector<int> &vColoring, size_t &cliqueSize)
{
#if 1
    // create an adjacencyArray, much faster for degeneracy ordering.

    vector<vector<int>> adjacencyArray(adjacencyMatrix.size());
    for (int vertex = 0; vertex < adjacencyMatrix.size(); ++vertex) {
        for (int otherVertex = 0; otherVertex < adjacencyMatrix.size(); ++otherVertex) {
            if (vertex == otherVertex) continue;
            if (adjacencyMatrix[vertex][otherVertex]) {
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


void OrderingTools::InitialOrderingConnectedComponent(string const &filename, vector<vector<char>> const &adjacencyMatrix, vector<int> &vOrderedVertices, vector<int> &vColoring)
{

////#if 1 ////def FACEBOOK
////    cout << "ERROR: Only works for yeast.955.kernel subgraph!" << endl;
////    vector<int> const vVertexOrder = Tools::ReadMetisOrdering("/Users/strash/graphs/kernels/biogrid-yeast.955.kernel.graph.iperm");
////#else
////    cout << "ERROR: Only works for 1et.1024.kernel subgraph!" << endl;
////    vector<int> const vVertexOrder = Tools::ReadMetisOrdering("/Users/strash/graphs/kernels/1et.1024.component.graph.iperm");
////#endif //FACEBOOK

    vector<int> const vVertexOrder = Tools::ReadMetisOrdering(filename);

    vOrderedVertices.resize(vVertexOrder.size(), 0);

    for (size_t index = 0; index < vVertexOrder.size(); ++index) {
////        vOrderedVertices[vOrderedVertices.size() - vVertexOrder[index] - 1] = index;
        vOrderedVertices[vVertexOrder[index]] = index;
    }

////    std::reverse(vOrderedVertices.begin(), vOrderedVertices.end());

    cout << "Snippet of ordering: ";
    for (size_t index = 0; index < vOrderedVertices.size(); ++index) {
#ifdef SNIPPET
        if (index == 10 && vOrderedVertices.size() > 21) {
            cout << "...";
            index = vOrderedVertices.size() - 10;
        }
#endif // SNIPPET

        cout << vOrderedVertices[index] << " ";
    }
    cout << endl;

#if 0
    size_t maxDegree(0);
    {
        vector<int> vDegree(adjacencyMatrix.size(), 0);
        for (size_t u = 0; u < adjacencyMatrix.size(); ++u) {
            for (size_t v = 0; v < adjacencyMatrix.size(); ++v) {
                if (adjacencyMatrix[u][v]) vDegree[u]++;
            }
            maxDegree = max(maxDegree, static_cast<size_t>(vDegree[u]));
        }
    }

    vColoring.reserve(adjacencyMatrix.size());
    vColoring.clear();
    for (int degree = 1; degree <= maxDegree; degree++ ) {
        vColoring.push_back(degree);
    }

    vColoring.resize(adjacencyMatrix.size(), maxDegree + 1);
#endif // 0
}
