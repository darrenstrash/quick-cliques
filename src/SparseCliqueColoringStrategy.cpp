#include "SparseCliqueColoringStrategy.h"
#include "DegeneracyTools.h"

#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

SparseCliqueColoringStrategy::SparseCliqueColoringStrategy(vector<vector<int>> const &adjacencyList) : ColoringStrategy(), m_vlColorToVertices(), m_vVertexToColor(adjacencyList.size(), -1), m_vNeighborHasColor()
{
////    m_VertexOrder = std::move(GetVerticesInDegeneracyOrder(adjacencyList));
}

void SparseCliqueColoringStrategy::Color(vector<vector<int>> const &adjacencyList, vector<int> const &vVertexOrder, vector<int> &vVerticesToReorder, vector<int> &vColors)
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

    m_vlColorToVertices.resize(vVerticesToReorder.size());
    m_vNeighborHasColor.resize(vVerticesToReorder.size(), false);

    size_t const numVerticesToReorder(vVerticesToReorder.size());

    for (int const vertex : vVertexOrder) {
        size_t color = 0;
        for (int const neighbor : adjacencyList[vertex]) {
            if (m_vVertexToColor[neighbor] != -1) {
                m_vNeighborHasColor[m_vVertexToColor[neighbor]] = true;
            }
        }

        for (size_t color = 0; color < m_vNeighborHasColor.size(); ++color) {
            if (!m_vNeighborHasColor[color]) {
                m_vlColorToVertices[color].push_back(vertex);
                m_vVertexToColor[vertex] = color;
                maxColor = max(maxColor, color);
                break;
            }
        }

        for (int const neighbor : adjacencyList[vertex]) {
            if (m_vVertexToColor[neighbor] != -1) {
                m_vNeighborHasColor[m_vVertexToColor[neighbor]] = false;
            }
        }
    }

////    cout << "maxColor=" << maxColor << ", numVertices=" << vVerticesToReorder.size() << endl;

    int currentIndex(0);
    int currentColor(0);
    for (int currentColor = 0; currentColor <= maxColor; ++currentColor) {
        for (int const vertex : m_vlColorToVertices[currentColor]) {
            vVerticesToReorder[currentIndex] = vertex;
            vColors[currentIndex] = currentColor + 1;
            currentIndex++;

            m_vVertexToColor[vertex] = -1; // reset all vertex colorings to be invalid
        }
        m_vlColorToVertices[currentColor].clear();
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

void SparseCliqueColoringStrategy::Color(vector<vector<int>> const &adjacencyList, vector<int> &vVerticesToReorder, vector<int> &vColors)
{
    Color(adjacencyList, vVerticesToReorder, vVerticesToReorder, vColors);
}
