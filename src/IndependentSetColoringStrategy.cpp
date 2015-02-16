#include "IndependentSetColoringStrategy.h"
#include "DegeneracyTools.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

IndependentSetColoringStrategy::IndependentSetColoringStrategy(vector<vector<char>> const &adjacencyMatrix) : ColoringStrategy(), m_AdjacencyMatrix(adjacencyMatrix), m_vvVerticesWithColor()////, m_Colors(adjacencyList.size(), -1) , m_VertexOrder()
{
    m_vvVerticesWithColor.resize(adjacencyMatrix.size());
    for (vector<int> vVertices : m_vvVerticesWithColor) {
        vVertices.reserve(adjacencyMatrix.size());
    }
////    m_VertexOrder = std::move(GetVerticesInDegeneracyOrder(adjacencyList));
}

int IndependentSetColoringStrategy::ColorWithoutReorder(vector<vector<char>> const &adjacencyMatrix, vector<int> const &vVertexOrder, vector<int> &vVerticesToReorder, vector<int> &vColors)
{
    if (vVerticesToReorder.empty()) return 0;

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

    vVerticesToReorder = vVertexOrder;
    vColors.clear();
    vColors.reserve(vVerticesToReorder.size());

    size_t const numVerticesToReorder(vVerticesToReorder.size());

    for (int const vertex : vVertexOrder) {
        size_t color = 0;
        for (vector<int> const &verticesWithColor : m_vvVerticesWithColor) {
            bool hasNeighborWithColor(false);
            if (verticesWithColor.empty()) break;
            for (int const coloredVertex : verticesWithColor) {
                // can be made more efficient?
                hasNeighborWithColor = !adjacencyMatrix[vertex][coloredVertex];
                if (hasNeighborWithColor) {
                    color++;
                    break;
                }
            }

            if (!hasNeighborWithColor) {
                break;
            }
        }

        vColors.push_back(color+1);
        m_vvVerticesWithColor[color].push_back(vertex);
        maxColor = max(maxColor, color);
    }


////    cout << "maxColor=" << maxColor << ", numVertices=" << vVerticesToReorder.size() << endl;

    for (int currentColor = 0; currentColor <= maxColor; ++currentColor) {
        m_vvVerticesWithColor[currentColor].clear();
    }

    return maxColor + 1;

#if 0
    cout << "Coloring (out): ";
    for (int index = 0; index < vVerticesToReorder.size(); index++) {
        int const vertex(vVerticesToReorder[index]);
        cout << vertex << "(" << vColors[index] << ") ";
    }
    cout << endl;
#endif // 0

// verify that it is a valid coloring.
#ifdef DEBUG
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


void IndependentSetColoringStrategy::Color(vector<vector<char>> const &adjacencyMatrix, vector<int> const &vVertexOrder, vector<int> &vVerticesToReorder, vector<int> &vColors)
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

    for (int const vertex : vVertexOrder) {
        size_t color = 0;
        for (vector<int> const &verticesWithColor : m_vvVerticesWithColor) {
            bool hasNeighborWithColor(false);
            if (verticesWithColor.empty()) break;
            for (int const coloredVertex : verticesWithColor) {
                // can be made more efficient?
                hasNeighborWithColor = !adjacencyMatrix[vertex][coloredVertex];
                if (hasNeighborWithColor) {
                    color++;
                    break;
                }
            }

            if (!hasNeighborWithColor) {
                break;
            }
        }

        m_vvVerticesWithColor[color].push_back(vertex);
        maxColor = max(maxColor, color);
    }

////    cout << "maxColor=" << maxColor << ", numVertices=" << vVerticesToReorder.size() << endl;

    int currentIndex(0);
    int currentColor(0);
    for (int currentColor = 0; currentColor <= maxColor; ++currentColor) {
        for (int const vertex : m_vvVerticesWithColor[currentColor]) {
            vVerticesToReorder[currentIndex] = vertex;
            vColors[currentIndex] = currentColor+1;
            currentIndex++;
        }
        m_vvVerticesWithColor[currentColor].clear();
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
#ifdef DEBUG
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

void IndependentSetColoringStrategy::Recolor(vector<vector<char>> const &adjacencyMatrix, vector<int> const &vVertexOrder, vector<int> &vVerticesToReorder, vector<int> &vColors, int const currentBestCliqueSize, int const currentCliqueSize)
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

    int iBestCliqueDelta(currentBestCliqueSize - currentCliqueSize);

    size_t const numVerticesToReorder(vVerticesToReorder.size());

    for (int const vertex : vVertexOrder) {
        size_t color = 0;
        for (vector<int> const &verticesWithColor : m_vvVerticesWithColor) {
            bool hasNeighborWithColor(false);
            if (verticesWithColor.empty()) break;
            for (int const coloredVertex : verticesWithColor) {
                // can be made more efficient?
                hasNeighborWithColor = !adjacencyMatrix[vertex][coloredVertex];
                if (hasNeighborWithColor) {
                    color++;
                    break;
                }
            }

            if (!hasNeighborWithColor) {
                break;
            }
        }

        m_vvVerticesWithColor[color].push_back(vertex);
        maxColor = max(maxColor, color);
        if (color+1 > iBestCliqueDelta && /*m_vvVerticesWithColor[color].size() == 1*/ color == maxColor) {
            Repair(vertex, color, iBestCliqueDelta);
            if (m_vvVerticesWithColor[maxColor].empty())
                maxColor--;
        }
    }

////    cout << "maxColor=" << maxColor << ", numVertices=" << vVerticesToReorder.size() << endl;

    int currentIndex(0);
    int currentColor(0);
    for (int currentColor = 0; currentColor <= maxColor; ++currentColor) {
        for (int const vertex : m_vvVerticesWithColor[currentColor]) {
            vVerticesToReorder[currentIndex] = vertex;
            vColors[currentIndex] = currentColor+1;
            currentIndex++;
        }
        m_vvVerticesWithColor[currentColor].clear();
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
#ifdef DEBUG
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

bool IndependentSetColoringStrategy::HasConflict(int const vertex, vector<int> const &vVerticesWithColor)
{
    if (vVerticesWithColor.empty()) return false;
    for (int const coloredVertex : vVerticesWithColor) {
        // can be made more efficient?
        if (!m_AdjacencyMatrix[vertex][coloredVertex]) {
            return true;
        }
    }

    return false;
}

int IndependentSetColoringStrategy::GetConflictingVertex(int const vertex, vector<int> const &vVerticesWithColor)
{
    int conflictingVertex(-1);
    int count(0);
    for (int const candidateVertex : vVerticesWithColor) {
        if (!m_AdjacencyMatrix[vertex][candidateVertex]) {
            conflictingVertex = candidateVertex;
            count++;
            if (count > 1) return -1;
        }
    }
    return conflictingVertex;
}

bool IndependentSetColoringStrategy::Repair(int const vertex, int const color, int const iBestCliqueDelta)
{
    for (int newColor = 0; newColor <= iBestCliqueDelta-1; newColor++) {
        int const conflictingVertex(GetConflictingVertex(vertex, m_vvVerticesWithColor[newColor]));
        if (conflictingVertex < 0) continue;
        for (int nextColor = newColor+1; nextColor <= iBestCliqueDelta; nextColor++) {
            if (HasConflict(conflictingVertex, m_vvVerticesWithColor[nextColor])) continue;
            m_vvVerticesWithColor[color].erase(find(m_vvVerticesWithColor[color].begin(), m_vvVerticesWithColor[color].end(), vertex));
            m_vvVerticesWithColor[newColor].erase(find(m_vvVerticesWithColor[newColor].begin(), m_vvVerticesWithColor[newColor].end(), conflictingVertex));
            m_vvVerticesWithColor[newColor].push_back(vertex);
            m_vvVerticesWithColor[nextColor].push_back(conflictingVertex);
////            cout << "Repairing vertices " << vertex << " and " << conflictingVertex << endl;
            return true;
        }
    }
    return false;
}
