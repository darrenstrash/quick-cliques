#include "OrderingTools.h"
#include "DegeneracyTools.h"
#include "GraphTools.h"

#include <vector>

using namespace std;

vector<int> OrderingTools:: InitialOrderingMCQ(vector<vector<char>> const &adjacencyMatrix, vector<int> const &degree)
{
    return std::move(GraphTools::OrderVerticesByDegree(adjacencyMatrix, degree, false /* non-increasing order */));
}

vector<int> OrderingTools:: InitialOrderingMCR(vector<vector<char>> const &adjacencyMatrix)
{
    vector<int> vOrderedVertices(adjacencyMatrix.size(), -1);

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

    for(size_t i = 0; i < size; i++) {
        degree[i] = adjacencyArray[i].size();
        verticesByDegree[degree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[degree[i]].begin();
    }

    int currentDegree = 0;
    int numVerticesRemoved = 0;

    vector<NeighborListArray> vOrderingArray(size);

    while (numVerticesRemoved < size) {
        if (!verticesByDegree[currentDegree].empty()) {
            int const vertex = verticesByDegree[currentDegree].front();
            verticesByDegree[currentDegree].pop_front();

            vOrderingArray[vertex].vertex = vertex;
            vOrderingArray[vertex].orderNumber = numVerticesRemoved;
            vOrderedVertices[vOrderedVertices.size() - numVerticesRemoved - 1] = vertex;

            degree[vertex] = -1;

            // will swap later neighbors to end of neighbor list
            vector<int> &neighborList = adjacencyArray[vertex];

            int splitPoint(neighborList.size());
            for(int i=0; i < splitPoint; ++i) {
                int const neighbor(neighborList[i]);
                // if later neighbor, swap to end of neighborList (there are few of these)
                if(degree[neighbor]!=-1) {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);

                    neighborList[i] = neighborList[--splitPoint];
                    neighborList[splitPoint] = neighbor;
                    i--;

                    degree[neighbor]--;

                    if (degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
                // earlier neighbor, do nothing.
            }

            // create space for later neighbors to ordering
            vOrderingArray[vertex].laterDegree = neighborList.size() - splitPoint;
            vOrderingArray[vertex].later.resize(neighborList.size() - splitPoint);

            // create space for earlier neighbors to ordering
            vOrderingArray[vertex].earlierDegree = splitPoint;
            vOrderingArray[vertex].earlier.resize(splitPoint);

            // fill in earlier and later neighbors
            for (int i = 0; i < splitPoint; ++i) {
////            cout << "earlier: " << vOrderingArray[vertex].earlier.size() << endl;
////            cout << "split  : " << splitPoint << endl;
////            cout << "earlier[" << i << "]" << endl;
                vOrderingArray[vertex].earlier[i] = neighborList[i];
            }

            for (int i = splitPoint; i < neighborList.size(); ++i) {
////            cout << "later  : " << vOrderingArray[vertex].later.size() << endl;
////            cout << "split  : " << splitPoint << endl;
////            cout << "later [" << i - splitPoint << "]" << endl;
                vOrderingArray[vertex].later[i-splitPoint] = neighborList[i];
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else {
            currentDegree++;
        }
    }

    return vOrderedVertices;
}
