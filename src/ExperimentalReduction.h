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
    void RestoreState() __attribute__((always_inline));

    void UndoReductions(std::list<int> &partialClique) __attribute__((always_inline));
    void ApplyReductions(int const vertex, std::list<int> &partialClique) __attribute__((always_inline));

    virtual bool GetNextTopLevelPartition();

    virtual void GetTopLevelPartialClique(std::list<int> &partialClique) const
    {
////        for (int const cliqueVertex : isolates.GetIsolates()) {
////            partialClique.push_back(cliqueVertex);
////        }
    }

    virtual int GetNextVertexToEvaluate()
    {
////#ifdef REMOVE_ISOLATES
////        return isolates.NextVertexToRemove();
////#else
        return -1;
////#endif
    }

////    std::string CheckP()
////    {
////        std::string errorString;
////        std::set<int> const &inGraph(isolates.GetInGraph());
////        for (int const vertex : P) {
////            if (inGraph.find(vertex) == inGraph.end()) {
////                errorString += "Mismatch: " + std::to_string(vertex) + " in P, but not graph. ";
////            }
////        }
////
////        if (P.Size() != inGraph.size()) {
////            errorString += "Mismatch: Size of P=" + std::to_string(SizeOfP()) + ", Graph Size=" + std::to_string(inGraph.size());
////            for (int const vertex : inGraph) {
////                if (!InP(vertex)) {
////                    errorString += "Mismatch: " + std::to_string(vertex) + " in graph, but not P. ";
////                }
////            }
////        }
////
////        if (!errorString.empty()) return errorString;
////        return "Ok";
////    }

protected: // methods

    void StoreGraphChanges();
    void RetrieveGraphChanges();

protected: // members
    std::vector<std::vector<int>> &m_AdjacencyList;
    std::vector<int> vertexSets;
    std::vector<int> vertexLookup;
    std::vector<int> degree;
    std::vector<SetDelineator> m_lDelineators;
////    Isolates isolates;
    std::vector<std::vector<int>> vvCliqueVertices;
    std::vector<std::vector<int>> vvOtherRemovedVertices;
    std::vector<std::vector<std::pair<int,int>>> vvAddedEdges;
    ArraySetsXPR m_Sets;
    std::vector<int> m_vCliqueVertices;
};

inline void ExperimentalReduction::ReturnVerticesToP(std::vector<int> const &vVertices)
{
#if 1 // TODO/DS: PUT BACK!
    for (int const vertex : vVertices) {
        m_Sets.MoveFromXToP(vertex);
    }

////    std::cout << "Returning vertices to P" << std::endl;
////    PrintSummary(__LINE__);

#ifdef REMOVE_ISOLATES
////    isolates.ReplaceAllRemoved(vVertices);
#endif
#endif
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

    partialClique.push_back(vertex);
////    ApplyReductions(vertex, partialClique);

////    bool const debug(vertex == 37 || vertex == 36);

    SaveState();

////    for (int const cliqueVertex : m_vCliqueVertices) {
        int const cliqueVertex(vertex);
        m_Sets.MoveFromPToR(cliqueVertex);
////        PrintSummary(__LINE__);
        for (int const neighbor : m_AdjacencyList[cliqueVertex]) {
////            if (debug) std::cout << "Evaluating neighbor " << neighbor << std::endl;
            if (m_Sets.InX(neighbor)) {
////                if (debug) std::cout << "    Neighbor is in X, Removing..." << std::endl;
                m_Sets.RemoveFromX(neighbor);
            } else if (m_Sets.InP(neighbor)) {
////                if (debug) std::cout << "    Neighbor is in P, Removing..." << std::endl;
                m_Sets.RemoveFromP(neighbor);
            }
        }
////    }

////    PrintSummary(__LINE__);
}

inline void ExperimentalReduction::MoveFromRToX(int const vertex)
{
    std::cout << "FATAL ERROR: This function should not be called" << std::endl;
    exit(1);
}

// DONE: need to verify
inline void ExperimentalReduction::MoveFromRToX(std::list<int> &partialClique, int const vertex)
{
    partialClique.pop_back();

    RestoreState();

    // after Restoring State, the vertex is back in P
    m_Sets.MoveFromPToX(vertex);

////    std::cout << "Moving " << vertex << " from R to X " << std::endl;
////    PrintSummary(__LINE__);
}

/*! \brief Computes the vertex v in P union X that has the most neighbors in P,
  and places P \ {neighborhood of v} in an array. These are the 
  vertices to consider adding to the partial clique during the current
  recursive call of the algorithm.

  \param pivotNeighbors  An intially empty vector, which will contain the set 
  P \ {neighborhood of v} when this function completes.
 */

////#define NOT_DONE

// TODO/DS: Choose pivot to maximize number of non-neighbors.
inline std::vector<int> ExperimentalReduction::ChoosePivot() const
{

#ifdef NOT_DONE
    std::vector<int> vVerticesInP;
    for (int const vertex : m_Sets.GetP()) {
        vVerticesInP.push_back(vertex);
    }

#ifdef DEBUG
    std::cout << " P : ";
    for (int const neighbor: vVerticesInP)
        std::cout << neighbor << " ";
    std::cout << std::endl;
#endif //DEBUG

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

////inline void ExperimentalReduction::ApplyReductions(int const vertex, std::list<int> &partialClique)
////{
////    std::vector<int> vCliqueVertices;
////    std::vector<int> vOtherRemoved;
////    std::vector<std::pair<int,int>> vAddedEdges;
////
////    isolates.RemoveVertexAndNeighbors(vertex, vOtherRemoved);
////    vCliqueVertices.push_back(vertex);
////
////#ifdef REMOVE_ISOLATES
////    isolates.RemoveAllIsolates(0, vCliqueVertices, vOtherRemoved, vAddedEdges);
////////    std::cout << "Removed " << vCliqueVertices.size() + vOtherRemoved.size() << " vertices in reduction" << std::endl;
////
////    partialClique.insert(partialClique.end(), vCliqueVertices.begin(), vCliqueVertices.end());
////    AddToAdjacencyList(vAddedEdges);
////#else
////    partialClique.push_back(vertex);
////#endif
////}

////inline void ExperimentalReduction::UndoReductions(std::list<int> &partialClique)
////{
////    // restore graph to contain everything in P.
////    std::vector<int>           vCliqueVertices;
////    std::vector<int>           vOtherRemoved;
////    std::vector<std::pair<int,int>> vEdges;
////    RetrieveGraphChanges(vCliqueVertices, vOtherRemoved, vEdges);
////
////    RemoveFromAdjacencyList(vEdges); // TODO/DS: this probably isn't right
////
////    isolates.ReplaceAllRemoved(vCliqueVertices);
////    for (int const vertex : vCliqueVertices) {
////        R.MoveTo(vertex, P);
////    }
////
////    isolates.ReplaceAllRemoved(vOtherRemoved);
////    for (int const vertex : vOtherRemoved) {
////        P.Insert(vertex);
////    }
////}

inline void ExperimentalReduction::SaveState()
{
    m_Sets.SaveState();

////    std::cout << "Saving P : ";
////    for (int const vertexInP : P.GetStdSet()) {
////        std::cout << vertexInP << " ";
////    }
////    std::cout << std::endl;

////    StoreGraphChanges();
}

inline void ExperimentalReduction::RestoreState()
{
    m_Sets.RestoreState();

////    std::cout << "Restoring P to : ";
////    for (int const vertexInP : P.GetStdSet()) {
////        std::cout << vertexInP << " ";
////    }
////    std::cout << std::endl;

////    UndoReductions();
}

#endif //EXPERIMENTAL_REDUCTION_H
