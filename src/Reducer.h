#ifndef REDUCER_H
#define REDUCER_H

#include "Isolates2.h"
#include "ArraySet.h"

#include <vector>

class IsolateReducer
{
public:
    IsolateReducer(std::vector<std::vector<int>> const &);
    virtual ~IsolateReducer();
    void InitialReduce(std::vector<int> &vCliqueVertices);
    virtual void Reduce(std::vector<int> &vAlreadyConsideredVertices, std::vector<int> &vCliqueVertices, std::vector<int> &vOtherRemovedVertices);
    void RemoveVertex(int const vertex);
    void RemoveVertexAndNeighbors(int const vertex, std::vector<int> &vOtherRemovedVertices);
    void ReplaceVertices(std::vector<int> const &vVertices);
    void ReplaceVertex(int const vertex);
    bool InRemainingGraph(int const vertex) const;

    size_t RemainingGraphSize() const { return isolates.GetInGraph().Size(); }

protected:
    std::vector<std::vector<int>> const &m_AdjacencyArray;
    std::vector<bool> vMarkedVertices;
    Isolates2<ArraySet> isolates;
};

class IsolateDominationReducer : public IsolateReducer
{
public:
    IsolateDominationReducer(std::vector<std::vector<int>> const &);
    virtual ~IsolateDominationReducer();
    virtual void Reduce(std::vector<int> &vAlreadyConsideredVertices, std::vector<int> &vCliqueVertices, std::vector<int> &vOtherRemovedVertices);

protected:
    ArraySet remaining;
};

#endif // REDUCER_H
