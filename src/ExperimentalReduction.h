#ifndef EXPERIMENTAL_REDUCTION_H
#define EXPERIMENTAL_REDUCTION_H

// local includes
#include "Set.h"
#include "SetsXPR.h"
#include "ArraySetsXPR.h"
#include "VertexSets.h"
#include "Isolates.h"

// system includes
#include <vector>
#include <cstring> // memcopy
#include <iostream>
#include <utility>

#define REMOVE_ISOLATES
////#define DO_PIVOT

class ExperimentalReduction : public VertexSets
{
public:
    ExperimentalReduction(std::vector<std::vector<int>> &adjacencyList);
    virtual ~ExperimentalReduction();

    ExperimentalReduction           (ExperimentalReduction const &sets) = delete;
    ExperimentalReduction& operator=(ExperimentalReduction const &sets) = delete;

    void AddToAdjacencyList     (std::vector<std::pair<int,int>> const &vEdges);
    void RemoveFromAdjacencyList(std::vector<std::pair<int,int>> const &vEdges);

    virtual void MoveFromPToR(int const vertexInP) __attribute__((always_inline));
    virtual void MoveFromPToR(std::list<int> &partialClique, int const vertexInP);
    virtual void MoveFromRToX(int const vertexInP) __attribute__((always_inline));
    virtual void MoveFromRToX(std::list<int> &partialClique, int const vertexInP);

    virtual void ReturnVerticesToP(std::vector<int> const &vVertices) __attribute__((always_inline));

    std::vector<int> ChoosePivot() const __attribute__((always_inline));
    bool InP(int const vertex) const __attribute__((always_inline));

    bool PIsEmpty() const __attribute__((always_inline));
    bool XAndPAreEmpty() const __attribute__((always_inline));

    virtual size_t SizeOfX() const { return m_Sets.SizeOfX(); }
    virtual size_t SizeOfP() const { return m_Sets.SizeOfP(); }

    virtual size_t GetGraphSize() const { return m_AdjacencyList.size(); }

    void Initialize();

    virtual void PrintSummary(int const line) const;

    void SaveState() __attribute__((always_inline));
    void RestoreState(std::list<int> &partialClique) __attribute__((always_inline));

    void UndoReductions(std::list<int> &partialClique) __attribute__((always_inline));
    void ApplyReductions(int const vertex, std::list<int> &partialClique) __attribute__((always_inline));

    virtual bool GetNextTopLevelPartition();

    virtual void GetTopLevelPartialClique(std::list<int> &partialClique) const
    {
        for (int const cliqueVertex : isolates.GetIsolates()) {
            partialClique.push_back(cliqueVertex);
        }
    }

    virtual int GetNextVertexToEvaluate(std::vector<int> &vVertices)
    {
#ifdef REMOVE_ISOLATES
        int const vertexToRemove(isolates.NextVertexToRemove(vVertices));
        // inefficient, should be a better way to elmininate it before getting here.
        for (int &vertex : vVertices) {
            if (vertex == vertexToRemove) {
                vertex = vVertices.back();
                vVertices.pop_back();
                return vertexToRemove;
            }
        }
#endif
        return VertexSets::GetNextVertexToEvaluate(vVertices);
    }

    std::string CheckP()
    {
        std::string errorString;
        std::set<int> const &inGraph(isolates.GetInGraph());
        for (int const vertex : m_Sets.GetP()) {
            if (inGraph.find(vertex) == inGraph.end()) {
                errorString += "Mismatch: " + std::to_string(vertex) + " in P, but not graph. ";
            }
        }

        if (m_Sets.SizeOfP() != inGraph.size()) {
            errorString += "Mismatch: Size of P=" + std::to_string(m_Sets.SizeOfP()) + ", Graph Size=" + std::to_string(inGraph.size());
            for (int const vertex : inGraph) {
                if (!InP(vertex)) {
                    errorString += "Mismatch: " + std::to_string(vertex) + " in graph, but not P. ";
                }
            }
        }

        if (!errorString.empty()) return errorString;
        return "Ok";
    }

protected: // methods

    void StoreGraphChanges();
    void RetrieveGraphChanges();

protected: // members
    std::vector<std::vector<int>> &m_AdjacencyList;
    std::vector<int> vertexSets;
    std::vector<int> vertexLookup;
    std::vector<int> degree;
    std::vector<SetDelineator> m_lDelineators;
    Isolates isolates;
    std::vector<std::vector<int>> vvCliqueVertices;
    std::vector<std::vector<int>> vvOtherRemovedVertices;
    std::vector<std::vector<std::pair<int,int>>> vvAddedEdges;
    ArraySetsXPR m_Sets;
    std::vector<int> m_vCliqueVertices;
};

inline void ExperimentalReduction::ReturnVerticesToP(std::vector<int> const &vVertices)
{
#ifndef DO_PIVOT
    for (int const vertex : vVertices) {
        m_Sets.MoveFromXToP(vertex);
    }

////    std::cout << "Returning vertices to P" << std::endl;
////    PrintSummary(__LINE__);

#ifdef REMOVE_ISOLATES
    isolates.ReplaceAllRemoved(vVertices);
#endif
#endif // DO_PIVOT
////    PrintSummary(__LINE__);
////    std::cout << __LINE__ << ": CheckP = " << CheckP() << std::endl;

////    std::cout << "EndCheckPoint: ";
////    PrintSummary(__LINE__);
}

inline void ExperimentalReduction::MoveFromPToR(int const vertex)
{
    std::cout << "FATAL ERROR: This function should not be called" << std::endl;
    exit(1);
}

// DONE, need to verify
inline void ExperimentalReduction::MoveFromPToR(std::list<int> &partialClique, int const vertex)
{
////    std::cout << "Moving " << vertex << " from P to R " << std::endl;
////    PrintSummary(__LINE__);

    ApplyReductions(vertex, partialClique);

////    bool const debug(vertex == 37 || vertex == 36);

    SaveState();

    for (int const cliqueVertex : m_vCliqueVertices) {
        m_Sets.MoveFromPToR(cliqueVertex);
////        partialClique.push_back(cliqueVertex);
////        PrintSummary(__LINE__);
        for (int const neighbor : m_AdjacencyList[cliqueVertex]) {
////            std::cout << "Evaluating neighbor " << neighbor << std::endl;
            if (m_Sets.InX(neighbor)) {
////                std::cout << "    Neighbor is in X, Removing..." << std::endl;
                m_Sets.RemoveFromX(neighbor);
            } else if (m_Sets.InP(neighbor)) {
////                std::cout << "    Neighbor is in P, Removing..." << std::endl;
                m_Sets.RemoveFromP(neighbor);
            }
        }
    }

////    PrintSummary(__LINE__);
////    std::cout << __LINE__ << ": CheckP = " << CheckP() << std::endl;
}

inline void ExperimentalReduction::MoveFromRToX(int const vertex)
{
    std::cout << "FATAL ERROR: This function should not be called" << std::endl;
    exit(1);
}

// DONE: need to verify
inline void ExperimentalReduction::MoveFromRToX(std::list<int> &partialClique, int const vertex)
{
    RestoreState(partialClique);

    // after Restoring State, the vertex is back in P
#ifndef DO_PIVOT
    m_Sets.MoveFromPToX(vertex);
    isolates.RemoveVertex(vertex);
#endif // DO_PIVOT

////    std::cout << "Moving " << vertex << " from R to X " << std::endl;
////    PrintSummary(__LINE__);
////    std::cout << __LINE__ << ": CheckP = " << CheckP() << std::endl;
}

/*! \brief Computes the vertex v in P union X that has the most neighbors in P,
  and places P \ {neighborhood of v} in an array. These are the 
  vertices to consider adding to the partial clique during the current
  recursive call of the algorithm.

  \param pivotNeighbors  An intially empty vector, which will contain the set 
  P \ {neighborhood of v} when this function completes.
 */

#ifndef DO_PIVOT
#define NOT_DONE
#endif // DO_PIVOT

// TODO/DS: Choose pivot to maximize number of non-neighbors.
inline std::vector<int> ExperimentalReduction::ChoosePivot() const
{

////    std::cout << "Checkpoint: ";
////    PrintSummary(__LINE__);

#ifdef NOT_DONE
    std::vector<int> vVerticesInP;
    for (int const vertex : m_Sets.GetP()) {
        vVerticesInP.push_back(vertex);
    }

    return vVerticesInP;
#else
    int pivot = -1;
    int maxIntersectionSize = -1;

    // loop through all vertices in P union X
    for (int const vertex : m_Sets.GetP()) {
        int neighborCount = 0;
        int numNeighbors = 0;

        // count the number of neighbors vertex has in P.
        // only count them if the degree of the vertex
        // is greater than the the best count so far.
        // count the number of neighbors vertex has in P.
        for (int const neighbor : m_AdjacencyList[vertex]) {
            if (m_Sets.InP(neighbor))
                neighborCount++;
        }

        // if vertex has more neighbors in P, then update the pivot
        if (static_cast<int>(m_Sets.SizeOfP() - neighborCount) > maxIntersectionSize) {
            maxIntersectionSize = m_Sets.SizeOfP() - neighborCount;
            pivot = vertex;
        }
    }

    std::vector<int> pivotNeighbors;

    if (pivot == -1) {
        return pivotNeighbors;
    }

    std::set<int> verticesToConsider;

    // mark neighbors of pivot that are in P.
    for (int const neighbor : m_AdjacencyList[pivot]) {
        // if the neighbor is in P, put it in pivot nonNeighbors
        if (InP(neighbor)) {
            verticesToConsider.insert(neighbor);
        }
    }

    pivotNeighbors.insert(pivotNeighbors.end(), verticesToConsider.begin(), verticesToConsider.end());

    pivotNeighbors.push_back(pivot); // pivot is in P

#ifdef DEBUG
    std::cout << " - : ";
    for (int const neighbor: pivotNeighbors)
        std::cout << neighbor << " ";
    std::cout << std::endl;
#endif // DEBUG

    return pivotNeighbors;
#endif // NOT_DONE
}

inline bool ExperimentalReduction::InP(int const vertex) const
{
    return m_Sets.InP(vertex);
}

inline bool ExperimentalReduction::PIsEmpty() const
{
    return m_Sets.PIsEmpty();
}

inline bool ExperimentalReduction::XAndPAreEmpty() const
{
    return m_Sets.XIsEmpty() && m_Sets.PIsEmpty();
}

inline void ExperimentalReduction::ApplyReductions(int const vertex, std::list<int> &partialClique)
{
    std::vector<int> &vCliqueVertices(vvCliqueVertices.back());
    std::vector<int> &vOtherRemoved(vvOtherRemovedVertices.back());
    std::vector<std::pair<int,int>> &vAddedEdges(vvAddedEdges.back());

    isolates.RemoveVertexAndNeighbors(vertex, vOtherRemoved);
    vCliqueVertices.push_back(vertex);

#ifdef REMOVE_ISOLATES
    isolates.RemoveAllIsolates(0, vCliqueVertices, vOtherRemoved, vAddedEdges);
////    std::cout << "Removed " << vCliqueVertices.size() + vOtherRemoved.size() << " vertices in reduction" << std::endl;
////    std::cout << "    Vertices: ";

////    for (int const cliqueVertex : vCliqueVertices) {
////        std::cout << cliqueVertex << " ";
////    }
////
////    for (int const otherVertex : vOtherRemoved) {
////        std::cout << otherVertex << " ";
////    }
////    std::cout << std::endl;

    partialClique.insert(partialClique.end(), vCliqueVertices.begin(), vCliqueVertices.end());
    AddToAdjacencyList(vAddedEdges);
#else
    partialClique.push_back(vertex);
#endif

}

inline void ExperimentalReduction::UndoReductions(std::list<int> &partialClique)
{
    RetrieveGraphChanges();

    // restore graph to contain everything in P.
    std::vector<int>                &vCliqueVertices(vvCliqueVertices.back());
    std::vector<int>                &vOtherRemoved  (vvOtherRemovedVertices.back());
    std::vector<std::pair<int,int>> &vEdges         (vvAddedEdges.back());

    RemoveFromAdjacencyList(vEdges); // TODO/DS: this probably isn't right

////    std::cout << "Returning " << vCliqueVertices.size() + vOtherRemoved.size() << " vertices to isolate graph" << std::endl;
    isolates.ReplaceAllRemoved(vCliqueVertices);
    for (int const vertex : vCliqueVertices) {
////        m_Sets.MoveFromRToP(vertex);
        partialClique.pop_back();
    }

    isolates.ReplaceAllRemoved(vOtherRemoved);
////    for (int const vertex : vOtherRemoved) {
////        m_Sets.MoveFromRToP(vertex);
////    }

    vCliqueVertices.clear();
    vOtherRemoved.clear();
    vEdges.clear();
}

inline void ExperimentalReduction::SaveState()
{
    m_Sets.SaveState();

////    std::cout << "Saving P : ";
////    for (int const vertexInP : P.GetStdSet()) {
////        std::cout << vertexInP << " ";
////    }
////    std::cout << std::endl;

    StoreGraphChanges();
}

inline void ExperimentalReduction::RestoreState(std::list<int> &partialClique)
{
    m_Sets.RestoreState();

////    std::cout << "Restoring P to : ";
////    for (int const vertexInP : P.GetStdSet()) {
////        std::cout << vertexInP << " ";
////    }
////    std::cout << std::endl;

    UndoReductions(partialClique);
}

#endif //EXPERIMENTAL_REDUCTION_H
