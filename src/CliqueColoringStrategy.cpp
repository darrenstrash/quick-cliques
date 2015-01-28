#include "CliqueColoringStrategy.h"
#include "DegeneracyTools.h"

#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

CliqueColoringStrategy::CliqueColoringStrategy(vector<vector<char>> const &adjacencyMatrix) : ColoringStrategy(), m_ColorToVertices()////, m_Colors(adjacencyList.size(), -1) , m_VertexOrder()
{
    m_ColorToVertices.resize(adjacencyMatrix.size());
    for (vector<int> vVertices : m_ColorToVertices) {
        vVertices.reserve(adjacencyMatrix.size());
    }
////    m_VertexOrder = std::move(GetVerticesInDegeneracyOrder(adjacencyList));
}

void CliqueColoringStrategy::Color(vector<vector<char>> const &adjacencyMatrix, vector<int> &vVerticesToReorder, vector<int> &vColors)
{
    if (vVerticesToReorder.empty()) return;

#if 0
    cout << "Coloring (in ): ";
////    for (int const vertex : vVerticesToReorder) {
    for (int index = 0; index < vVerticesToReorder.size(); index++) {
        int const vertex(vVerticesToReorder[index]);
        cout << vertex << "(" << vColors[index] << ") ";
    }
    cout << endl;
#endif // 0

    size_t maxColor(0);

    size_t const numVerticesToReorder(vVerticesToReorder.size());

    for (int const vertex : vVerticesToReorder) {
        size_t color = 0;
        for (vector<int> const &verticesWithColor : m_ColorToVertices) {
            bool hasNeighborWithColor(false);
            if (verticesWithColor.empty()) break;
            for (int const coloredVertex : verticesWithColor) {
                // can be made more efficient?
                hasNeighborWithColor = adjacencyMatrix[vertex][coloredVertex];
                if (hasNeighborWithColor) {
                    color++;
                    break;
                }
            }

            if (!hasNeighborWithColor) {
                break;
            }
        }

        m_ColorToVertices[color].push_back(vertex);
        maxColor = max(maxColor, color);
    }

////    cout << "maxColor=" << maxColor << ", numVertices=" << vVerticesToReorder.size() << endl;

    int currentIndex(0);
    int currentColor(0);
    for (int currentColor = 0; currentColor <= maxColor; ++currentColor) {
        for (int const vertex : m_ColorToVertices[currentColor]) {
            vVerticesToReorder[currentIndex] = vertex;
            vColors[currentIndex] = currentColor+1;
            currentIndex++;
        }
        m_ColorToVertices[currentColor].clear();
    }

#if 0
    cout << "Coloring (out): ";
    for (int index = 0; index < vVerticesToReorder.size(); index++) {
        int const vertex(vVerticesToReorder[index]);
        cout << vertex << "(" << vColors[index] << ") ";
    }
    cout << endl;
#endif // 0

// verify that it is a valid coloring.
#if DEBUG ////def DEBUG
    vector<int> vColor(adjacencyList.size(), -1);
    for (size_t index = 0; index < vVerticesToReorder.size(); ++index) {
        vColor[vVerticesToReorder[index]] = vColors[index];
    }
    for (int const vertex : vVerticesToReorder) {
        for (int const neighbor : adjacencyList[vertex]) {
            if (vColor[vertex] == vColor[neighbor]) {
                cout << "Consistency Error: vertex " << vertex << " has the same color as neighbor " << neighbor << ", color=" << vColor[vertex] << endl << flush;
            }
        }
    }
#endif // DEBUG
}
