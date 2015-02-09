#ifndef COLORING_STRATEGY_H
#define COLORING_STRATEGY_H

#include <vector>

class ColoringStrategy
{
public:
    ColoringStrategy() {}
    virtual void Color(std::vector<std::vector<int>>  const &adjacencyList, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors) {};
    virtual void Color(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors) {};

    virtual void Color(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors) {};
    virtual void Recolor(std::vector<std::vector<char>> const &adjacencyMatrix, std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors) {};
    virtual void Recolor(std::vector<std::vector<int>> const &adjacencyArray, std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors) {};
////    virtual void Recolor() = 0;
////    virtual void RemoveVertex(int const vertex) = 0;
////    virtual void PeekAtNextVertexAndColor(int &vertex, int &color) = 0;
};

#endif //COLORING_STRATEGY_H
