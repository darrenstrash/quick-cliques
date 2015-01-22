#ifndef EXPERIMENTAL_REDUCTION_H
#define EXPERIMENTAL_REDUCTION_H

// local includes
#include "Set.h"
#include "SetsXPR.h"
#include "ArraySetsXPR.h"
#include "VertexSets.h"
#include "Isolates.h"
#include "Isolates2.h"

// system includes
#include <vector>
#include <set>
#include <cstring> // memcopy
#include <iostream>
#include <utility>

#define PERSISTENT_REMOVE_ISOLATES
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
    virtual void MoveFromRToX(std::list<int> &partialClique, std::vector<int> &vVerticesToEvaluate, int const vertexInP);

    virtual void ReturnVerticesToP(std::vector<int> const &vVertices) __attribute__((always_inline));
    virtual void ReturnVerticesToP(std::list<int> &partialClique, std::vector<int> const &vVertices) __attribute__((always_inline));

    std::vector<int> ChoosePivotNonConst() __attribute__((always_inline));
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
////            if (cliqueVertex == 40653) { std::cout << "vertex 40653 is in top-level clique." << std::endl << std::flush; }
            partialClique.push_back(cliqueVertex);
        }
    }

    virtual int GetNextVertexToEvaluate(std::vector<int> &vVertices)
    {
#ifdef REMOVE_ISOLATES
        return isolates.NextVertexToRemove(vVertices);
#else
        return VertexSets::GetNextVertexToEvaluate(vVertices);
#endif
    }

    std::string CheckP()
    {
        std::string errorString;
        ArraySet const &inGraph(isolates.GetInGraph());
        for (int const vertex : m_Sets.GetP()) {
            if (!inGraph.Contains(vertex)) {
                errorString += "Mismatch: " + std::to_string(vertex) + " in P, but not graph. ";
            }
        }

        if (m_Sets.SizeOfP() != inGraph.Size()) {
            errorString += "Mismatch: Size of P=" + std::to_string(m_Sets.SizeOfP()) + ", Graph Size=" + std::to_string(inGraph.Size());
            for (int const vertex : inGraph) {
                if (!InP(vertex)) {
                    errorString += "Mismatch: " + std::to_string(vertex) + " in graph, but not P. ";
                }
            }
        }

        if (!errorString.empty()) return errorString;
        return "Ok";
    }

////virtual void RemoveDominatedVerticesFromVector(std::vector<int> &vVerticesInP);

virtual void RemoveDominatedVertices(std::vector<int> &dominatedVertices)
{
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
    Isolates2 isolates;
    std::vector<std::vector<int>> vvCliqueVertices;
    std::vector<std::vector<int>> vvOtherRemovedVertices;
    std::vector<std::vector<std::pair<int,int>>> vvAddedEdges;
    std::vector<std::vector<int>> vvPersistentCliqueVertices;
    std::vector<std::vector<int>> vvPersistentOtherRemovedVertices;
    std::vector<std::vector<std::pair<int,int>>> vvPersistentAddedEdges;
    ArraySetsXPR m_Sets;
    std::vector<int> m_vCliqueVertices;
};

inline void ExperimentalReduction::ReturnVerticesToP(std::vector<int> const &vVertices)
{
    std::cout << "FATAL ERROR: This function should not be called" << std::endl;
    exit(1);
}

inline void ExperimentalReduction::ReturnVerticesToP(std::list<int> &partialClique, std::vector<int> const &vVertices)
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

#ifdef PERSISTENT_REMOVE_ISOLATES

    // remove clique vertices from partial clique.
    for (int const cliqueVertices : vvPersistentCliqueVertices.back()) {
        partialClique.pop_back();
    }

    // put all removed vertices back into P.
    // TODO/DS: not needed, because the higher level recursive call has correct P...?

    // put all removed vertices back into graph.
    isolates.ReplaceAllRemoved(vvPersistentCliqueVertices.back());
    isolates.ReplaceAllRemoved(vvPersistentOtherRemovedVertices.back());

    for (int const vertex : vvPersistentCliqueVertices.back()) {
        m_Sets.MoveFromRToP(vertex);
    }

    for (int const vertex : vvPersistentOtherRemovedVertices.back()) {
        m_Sets.MoveFromRToP(vertex);
    }


    ////    AddToAdjacencyList(vAddedEdges); // would need to add back to isolates graph too.

    // remove sets for this recursive call;
    vvPersistentCliqueVertices.pop_back();
    vvPersistentOtherRemovedVertices.pop_back();
    vvPersistentAddedEdges.pop_back();
#endif //PERSISTENT_REMOVE_ISOLATES

////    PrintSummary(__LINE__);
////    std::cout << __LINE__ << ": CheckP = " << CheckP() << std::endl;

////    std::cout << "EndCheckPoint: ";
////    std::set<int> P;
////    P.insert(m_Sets.GetP().begin(), m_Sets.GetP().end());
////    for (int const vertex : P) {
////        std::cout << vertex << " ";
////    }
////    std::cout << std::endl << std::flush;
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

////    std::cout << __LINE__ << ": Checking P..." << std::endl;
////    CheckP();


    ApplyReductions(vertex, partialClique);

////    std::cout << __LINE__ << ": Checking P..." << std::endl;
////    CheckP();
////    if (!InP(31)) std::cout << "vertex 31 is not in P..." << std::endl << std::flush;
////    ArraySet const &inGraph(isolates.GetInGraph());
////    if (!inGraph.Contains(31)) std::cout << "vertex 31 is not in Graph..." << std::endl << std::flush;

////    bool const debug(vertex == 37 || vertex == 36);

    SaveState();


////    std::cout << "Clique Vertices : ";
    for (int const cliqueVertex : m_vCliqueVertices) {
////        std::cout << cliqueVertex << " " << std::flush;
////        if (!InP(cliqueVertex)) {
////            std::cout << "Moving clique vertex " << cliqueVertex << ", when it is not in P!" << std::endl << std::flush; 
////            CheckP();

        m_Sets.MoveFromPToR(cliqueVertex);
////        partialClique.push_back(cliqueVertex);
////        PrintSummary(__LINE__);
        for (int const neighbor : m_AdjacencyList[cliqueVertex]) { // TODO/DS: validate this change
//        for (int const neighbor : m_vOtherVertices) { // TODO/DS: validate this change
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
////    std::cout << std::endl << std::flush;

////    PrintSummary(__LINE__);
////    std::cout << __LINE__ << ": CheckP = " << CheckP() << std::endl;
}

inline void ExperimentalReduction::MoveFromRToX(int const vertex)
{
    std::cout << "FATAL ERROR: This function should not be called" << std::endl;
    exit(1);
}

// DONE: need to verify
inline void ExperimentalReduction::MoveFromRToX(std::list<int> &partialClique, std::vector<int> &vVerticesToEvaluate, int const vertex)
{
    RestoreState(partialClique);    // after Restoring State, the vertex is back in P

#ifndef DO_PIVOT
    // Remove from P, and from subgraph
    m_Sets.MoveFromPToX(vertex);
    isolates.RemoveVertex(vertex);

#ifdef PERSISTENT_REMOVE_ISOLATES
    // apply reductions to new subgraph, add clique vertices to partialClique
    // save them so we can undo them when returning vertices to P.
    std::vector<int> vCliqueVertices;
    std::vector<int> vOtherRemoved;
    std::vector<std::pair<int,int>> &vAddedEdges(vvPersistentAddedEdges.back());
    isolates.RemoveAllIsolates(0, vCliqueVertices, vOtherRemoved, vAddedEdges, false /* consider only precomputed vertices */);

    // add clique vertices to partial clique
    partialClique.insert(partialClique.end(), vCliqueVertices.begin(), vCliqueVertices.end());

    // remove them so the vertices aren't considered in future recursive calls
    vVerticesToEvaluate.clear();
    vVerticesToEvaluate.insert(vVerticesToEvaluate.end(), isolates.GetInGraph().begin(), isolates.GetInGraph().end());

    // remove vertices from P
    for (int const cliqueVertex : vCliqueVertices) {
        m_Sets.RemoveFromP(cliqueVertex);
    }
    for (int const otherVertex : vOtherRemoved) {
        m_Sets.RemoveFromP(otherVertex);
    }

    // save removed vertices, so we can restore them later.
    vvPersistentCliqueVertices.back().insert(vvPersistentCliqueVertices.back().end(), vCliqueVertices.begin(), vCliqueVertices.end());
    vvPersistentOtherRemovedVertices.back().insert(vvPersistentOtherRemovedVertices.back().end(), vOtherRemoved.begin(), vOtherRemoved.end());
    vvPersistentAddedEdges.back().insert(vvPersistentAddedEdges.back().end(), vAddedEdges.begin(), vAddedEdges.end());
#endif // PERSISTENT_REMOVE_ISOLATES
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
inline std::vector<int> ExperimentalReduction::ChoosePivotNonConst()
{

#ifdef PERSISTENT_REMOVE_ISOLATES
    vvPersistentCliqueVertices.push_back(std::vector<int>());
    vvPersistentOtherRemovedVertices.push_back(std::vector<int>());
    vvPersistentAddedEdges.push_back(std::vector<std::pair<int,int>>());
#endif //PERSISTENT_REMOVE_ISOLATES

////    std::cout << "Checkpoint: ";
////    std::set<int> P;
////    P.insert(m_Sets.GetP().begin(), m_Sets.GetP().end());
////    for (int const vertex : P) {
////        std::cout << vertex << " ";
////    }
////    std::cout << std::endl << std::flush;
    ////PrintSummary(__LINE__);

#ifdef NOT_DONE
    std::vector<int> vVerticesInP;
    for (int const vertex : m_Sets.GetP()) {
        vVerticesInP.push_back(vertex);
    }

    return vVerticesInP;
#else
#ifdef MIN_DEGREE_PIVOT
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
#else
    int const pivot = isolates.NextVertexToRemove();
#endif // MIN_DEGREE_PIVOT

    std::vector<int> pivotNeighbors;

    if (pivot == -1) {
        return pivotNeighbors;
    }

    pivotNeighbors.push_back(pivot); // pivot is in P
    // mark neighbors of pivot that are in P.
    for (int const neighbor : m_AdjacencyList[pivot]) {
        // if the neighbor is in P, put it in pivot nonNeighbors
        if (InP(neighbor)) {
            pivotNeighbors.push_back(neighbor);
        }
    }
    
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

    int removed(vOtherRemoved.size() + 1);

#ifdef REMOVE_ISOLATES
    isolates.RemoveAllIsolates(0, vCliqueVertices, vOtherRemoved, vAddedEdges, false /* use precomputed list of vertices to reduce */);
////    std::cout << "Removed " << vCliqueVertices.size() + vOtherRemoved.size() - removed << "/" << m_Sets.SizeOfP() << " vertices in reduction" << std::endl;

////    std::set<int> removedVertices;
////    removedVertices.insert(vCliqueVertices.begin(), vCliqueVertices.end());
////    removedVertices.insert(vOtherRemoved.begin(), vOtherRemoved.end());
////    std::cout << "All removed vertices: ";
////    for (int const removedVertex : removedVertices) {
////        std::cout << removedVertex << " ";
////    }
////    std::cout << std::endl << std::flush;

////    std::cout << "    Clique Vertices: ";
////
////    for (int const cliqueVertex : vCliqueVertices) {
////        std::cout << cliqueVertex << " ";
////    }
////
////    for (int const otherVertex : vOtherRemoved) {
////        std::cout << otherVertex << " ";
////    }
////    std::cout << std::endl;

    partialClique.insert(partialClique.end(), vCliqueVertices.begin(), vCliqueVertices.end());
////    AddToAdjacencyList(vAddedEdges);
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
