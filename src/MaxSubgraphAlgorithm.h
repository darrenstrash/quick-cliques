#ifndef MAX_SUBGRAPH_ALGORITHM_H
#define MAX_SUBGRAPH_ALGORITHM_H

#include "Algorithm.h"
#include "IndependentSetColoringStrategy.h"

#include <vector>
#include <list>
#include <ctime>

////#define PREPRUNE
////#define REMOVE_ISOLATES_BEFORE_ONLY
////#define ALWAYS_REMOVE_ISOLATES_AFTER
////#define NO_ISOLATES_P_LEFT_10

class MaxSubgraphAlgorithm : public Algorithm
{
public:
    MaxSubgraphAlgorithm(std::string const &name);
    virtual ~MaxSubgraphAlgorithm();

    virtual long Run(std::list<std::list<int>> &cliques);
////    virtual long Run(vector<int> const &startingVertices, std::list<std::list<int>> &cliques);

    virtual void Color(std::vector<int> const &vVertexOrder, std::vector<int> &vVerticesToReorder, std::vector<int> &vColors) = 0;

    virtual void InitializeOrder(std::vector<int> &P, std::vector<int> &vVertexOrder, std::vector<int> &vColors) = 0;
    virtual void GetNewOrder(std::vector<int> &vNewVertexOrder, std::vector<int> &vVertexOrder, std::vector<int> const &P, int const chosenVertex) = 0;
    virtual void ProcessOrderAfterRecursion(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors, int const chosenVertex) = 0;

    virtual void ProcessOrderBeforeReturn(std::vector<int> &vVertexOrder, std::vector<int> &P, std::vector<int> &vColors) = 0;

    virtual void RunRecursive(std::vector<int> &P, std::vector<int> &vVertexOrder, std::list<std::list<int>> &cliques, std::vector<int> &vColors);

    virtual void SetQuiet(bool const quiet) { m_bQuiet = quiet; }

    virtual void PrintState() const;

    virtual void SetNodeCount(size_t const count) { nodeCount = count; }
    virtual size_t GetNodeCount() { return nodeCount; }

    void SetR(std::vector<int> const &newR) { R = newR; }
    void SetMaximumCliqueSize(size_t const newCliqueSize) { m_uMaximumCliqueSize = newCliqueSize; }

    void SetOnlyVertex(int const vertex) { m_iOnlyVertex = vertex; }
protected:
    size_t m_uMaximumCliqueSize;
    std::vector<int> R;
    std::vector<std::vector<int>> stackP;
    std::vector<std::vector<int>> stackColors;
    std::vector<std::vector<int>> stackOrder;
    size_t nodeCount;
    int depth;
    clock_t startTime;
    clock_t timeToLargestClique;
    bool    m_bQuiet;
    std::vector<bool> stackEvaluatedHalfVertices;
////    bool m_bInvert;
    int m_iOnlyVertex;
};
#endif // MAX_SUBGRAPH_ALGORITHM_H
