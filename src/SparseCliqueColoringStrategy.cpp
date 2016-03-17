#include "SparseCliqueColoringStrategy.h"
#include "DegeneracyTools.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

SparseCliqueColoringStrategy::SparseCliqueColoringStrategy(vector<vector<int>> const &adjacencyArray)
 : ColoringStrategy()
 , m_AdjacencyArray(adjacencyArray)
 , m_vlColorToVertices()
 , m_vVertexToColor(adjacencyArray.size(), -1)
 , m_vNeighborHasColor()
 , m_vvVerticesWithColor(adjacencyArray.size())
 , m_vNeighborColorCount(adjacencyArray.size(), 0)
 , m_vbNeighbors(adjacencyArray.size(), false)
 , m_vbConflictNeighbors(adjacencyArray.size(), false)
{
////    m_VertexOrder = std::move(GetVerticesInDegeneracyOrder(adjacencyArray));
}

// TODO/DS (17.03.2016): Is the sparse independent set version better?
//                       Using lists isn't very efficient. Use m_vvVerticesWithColor instead
void SparseCliqueColoringStrategy::Color(vector<vector<int>> const &adjacencyArray, vector<int> const &vVertexOrder, vector<int> &vVerticesToReorder, vector<int> &vColors)
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
        for (int const neighbor : adjacencyArray[vertex]) {
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

        for (int const neighbor : adjacencyArray[vertex]) {
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
    vector<int> vColor(adjacencyArray.size(), -1);
    for (size_t index = 0; index < vVerticesToReorder.size(); ++index) {
        vColor[vVerticesToReorder[index]] = vColors[index];
    }
    for (int const vertex : vVerticesToReorder) {
        for (int const neighbor : adjacencyArray[vertex]) {
            if (vColor[vertex] == vColor[neighbor]) {
                cout << "Consistency Error: vertex " << vertex << " has the same color as neighbor " << neighbor << ", color=" << vColor[vertex] << endl << flush;
            }
        }
    }
#endif // DEBUG
}

void SparseCliqueColoringStrategy::Color(vector<vector<int>> const &adjacencyArray, vector<int> &vVerticesToReorder, vector<int> &vColors)
{
    Color(adjacencyArray, vVerticesToReorder, vVerticesToReorder, vColors);
}


// TODO/DS: Convert to Sparse
void SparseCliqueColoringStrategy::Recolor(vector<vector<int>> const &adjacencyArray, vector<int> const &vVertexOrder, vector<int> &vVerticesToReorder, vector<int> &vColors, int const currentBestCliqueSize, int const currentCliqueSize)
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

    int maxColor(-1);

    int iBestCliqueDelta(currentBestCliqueSize - currentCliqueSize);

    size_t const numVerticesToReorder(vVerticesToReorder.size());

    for (int const vertex : vVertexOrder) {
#if 1
        int uSmallestFreeColor = maxColor + 1;
        // first count the number of neighbors with a given color
        for (int const neighbor : m_AdjacencyArray[vertex]) {
            m_vbNeighbors[neighbor] = true;
            if (m_vVertexToColor[neighbor] != -1) {
                m_vNeighborColorCount[m_vVertexToColor[neighbor]]++;
            }
        }

        // compare color counts to total number of vertices with the color
        // if there is a difference, then there exists a non-neighbor with
        // that color; otherwise the color is free. Pick the smallest such
        // free color
        for (int const neighbor : m_AdjacencyArray[vertex]) {
            int const currentColor(m_vVertexToColor[neighbor]);
////            if (debug && neighbor==28) {
////                cout << " neighbor 28 has color " << currentColor << ", there are " << m_vvVerticesWithColor[currentColor].size() << " vertices with that color, and " << m_vNeighborColorCount[currentColor] << " neighbors with that color" << endl;
////            }
            if (currentColor != -1 && static_cast<int>(m_vvVerticesWithColor[currentColor].size()) == m_vNeighborColorCount[currentColor]) {
                uSmallestFreeColor = min(currentColor, static_cast<int>(uSmallestFreeColor));
            }
        }

////        // put color counts back to 0.
////        for (int const neighbor : m_AdjacencyArray[vertex]) {
////            if (m_vVertexToColor[neighbor] != -1) {
////                m_vNeighborColorCount[m_vVertexToColor[neighbor]] = 0;
////            }
////        }

#else
        int uSmallestFreeColor = 0;
        for (vector<int> const &verticesWithColor : m_vvVerticesWithColor) {
            bool hasNeighborWithColor(false);
            if (verticesWithColor.empty()) break;
            for (int const coloredVertex : verticesWithColor) {
                // can be made more efficient?
                hasNeighborWithColor = (find(adjacencyArray[vertex].begin(), adjacencyArray[vertex].end(), coloredVertex) == adjacencyArray[vertex].end());
                if (hasNeighborWithColor) {
                    uSmallestFreeColor++;
                    break;
                }
            }

            if (!hasNeighborWithColor) {
                break;
            }
        }
#endif // 1
////        cout << "vertex " << vertex << " gets initial color " << uSmallestFreeColor << endl;

        m_vvVerticesWithColor[uSmallestFreeColor].push_back(vertex);
        m_vVertexToColor[vertex] = uSmallestFreeColor;
        maxColor = max(maxColor, uSmallestFreeColor);
        if (uSmallestFreeColor +1 > iBestCliqueDelta && /*m_vvVerticesWithColor[color].size() == 1*/ uSmallestFreeColor == maxColor) {
            Repair(vertex, uSmallestFreeColor, iBestCliqueDelta);
            if (m_vvVerticesWithColor[maxColor].empty())
                maxColor--;
        }

        // put color counts back to 0. Needs to come after repair, repair uses the counts.
        for (int const neighbor : m_AdjacencyArray[vertex]) {
            m_vbNeighbors[neighbor] = false;
            if (m_vVertexToColor[neighbor] != -1) {
                m_vNeighborColorCount[m_vVertexToColor[neighbor]] = 0;
            }
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

            m_vVertexToColor[vertex] = -1;
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
    vector<int> vColor(adjacencyArray.size(), -1);
    for (size_t index = 0; index < vVerticesToReorder.size(); ++index) {
        vColor[vVerticesToReorder[index]] = vColors[index];
    }
    for (int const vertex : vVerticesToReorder) {
        for (int const neighbor : adjacencyArray[vertex]) {
            if (vColor[vertex] == vColor[neighbor]) {
                cout << "Consistency Error: vertex " << vertex << " has the same color as neighbor " << neighbor << ", color=" << vColor[vertex] << endl << flush;
            }
        }
    }
#endif // DEBUG
}

bool SparseCliqueColoringStrategy::HasConflict(int const vertex, vector<int> const &vVerticesWithColor)
{
    if (vVerticesWithColor.empty()) return false;
    for (int const coloredVertex : vVerticesWithColor) {
        // can be made more efficient?
        if (m_vbConflictNeighbors[coloredVertex]) {
            return true;
        }
    }

    return false;
}

int SparseCliqueColoringStrategy::GetConflictingVertex(int const vertex, vector<int> const &vVerticesWithColor)
{
    int conflictingVertex(-1);
    int count(0);
    for (int const candidateVertex : vVerticesWithColor) {
        if (m_vbNeighbors[candidateVertex]) { ////find(m_AdjacencyArray[vertex].begin(), m_AdjacencyArray[vertex].end(), candidateVertex) == m_AdjacencyArray[vertex].end()) {
            conflictingVertex = candidateVertex;
            count++;
            if (count > 1) return -1;
        }
    }
    return conflictingVertex;
}

bool SparseCliqueColoringStrategy::Repair(int const vertex, int const color, int const iBestCliqueDelta)
{
    for (int newColor = 0; newColor <= iBestCliqueDelta-1; newColor++) {
        int const conflictingVertex(GetConflictingVertex(vertex, m_vvVerticesWithColor[newColor]));
        if (conflictingVertex < 0) continue;

        for (int const neighbor : m_AdjacencyArray[conflictingVertex]) {
            m_vbConflictNeighbors[neighbor] = true;
        }

        // TODO/DS: put conflicting neighbors into array for quick testing.
        for (int nextColor = newColor+1; nextColor <= iBestCliqueDelta; nextColor++) {
            if (HasConflict(conflictingVertex, m_vvVerticesWithColor[nextColor])) continue;
            m_vvVerticesWithColor[color].erase(find(m_vvVerticesWithColor[color].begin(), m_vvVerticesWithColor[color].end(), vertex));
            m_vvVerticesWithColor[newColor].erase(find(m_vvVerticesWithColor[newColor].begin(), m_vvVerticesWithColor[newColor].end(), conflictingVertex));
            m_vvVerticesWithColor[newColor].push_back(vertex);
            m_vvVerticesWithColor[nextColor].push_back(conflictingVertex);

            m_vVertexToColor[vertex] = newColor;
            m_vVertexToColor[conflictingVertex] = nextColor;
////            cout << "Repairing vertices " << vertex << " and " << conflictingVertex << endl;

            for (int const neighbor : m_AdjacencyArray[conflictingVertex]) {
                m_vbConflictNeighbors[neighbor] = false;
            }
            return true;
        }

        for (int const neighbor : m_AdjacencyArray[conflictingVertex]) {
            m_vbConflictNeighbors[neighbor] = false;
        }
    }
    return false;
}
