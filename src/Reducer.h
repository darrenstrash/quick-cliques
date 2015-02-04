#ifndef REDUCER_H
#define REDUCER_H

#include "Isolates2.h"
#include "ArraySet.h"

#include <vector>

class Reducer
{
public:
    Reducer(std::vector<std::vector<int>> const &);
    ~Reducer();
    void InitialReduce(std::vector<int> &vCliqueVertices);
    void Reduce(std::vector<int> const &vAlreadyConsideredVertices, std::vector<int> &vCliqueVertices, std::vector<int> &vOtherRemovedVertices);
    void RemoveVertex(int const vertex);
    void RemoveVertexAndNeighbors(int const vertex, std::vector<int> &vOtherRemovedVertices);
    void ReplaceVertices(std::vector<int> const &vVertices);
    void ReplaceVertex(int const vertex);
    bool InRemainingGraph(int const vertex) const;

    size_t RemainingGraphSize() const { return isolates.GetInGraph().Size(); }

private:
    std::vector<std::vector<int>> const &m_AdjacencyArray;
    Isolates2<ArraySet> isolates;
    std::vector<bool> vMarkedVertices;
    ArraySet remaining;
};
#endif // REDUCER_H
