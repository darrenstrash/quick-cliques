#include "Reducer.h"
#include "Isolates2.h"
#include "ArraySet.h"

#include <utility>

using namespace std;

IsolateReducer::IsolateReducer(std::vector<std::vector<int>> const &adjacencyArray)
: m_AdjacencyArray(adjacencyArray)
, vMarkedVertices(adjacencyArray.size(), false)
, isolates(adjacencyArray)
{
}

IsolateReducer::~IsolateReducer()
{
}

void IsolateReducer::InitialReduce(vector<int> &vCliqueVertices)
{
    vector<int> vRemovedUnused;
    vector<pair<int,int>> vAddedEdgesUnused;
    isolates.RemoveAllIsolates(0, vCliqueVertices, vRemovedUnused, vAddedEdgesUnused, true /* consider all vertices */);
}

void IsolateReducer::Reduce(vector<int> &vAlreadyConsideredVertices, vector<int> &vCliqueVertices, vector<int> &vOtherRemovedVertices)
{
    vector<pair<int,int>> vAddedEdgesUnused;
////    int const graphSize(isolates.GetInGraph().Size());
    isolates.RemoveAllIsolates(0, vCliqueVertices, vOtherRemovedVertices, vAddedEdgesUnused, false /* consider only changed vertices */);

////    cout << "Reduced " << (graphSize - isolates.GetInGraph().Size()) << "/" << graphSize << " vertices through isolate removal" << endl;

    for (int const cliqueVertex : vCliqueVertices) {
        vMarkedVertices[cliqueVertex] = true;
    }

    size_t uNewIndex(0);
    for (size_t index = 0; index < vAlreadyConsideredVertices.size(); ++index) {
        int const vertex(vAlreadyConsideredVertices[index]);
        bool keep(true);
        for (int const neighbor : m_AdjacencyArray[vertex]) {
            if (vMarkedVertices[neighbor]) { keep = false; break; }
        }

        if (keep) {
            vAlreadyConsideredVertices[uNewIndex++] = vAlreadyConsideredVertices[index];
        }
    }

    vAlreadyConsideredVertices.resize(uNewIndex);

    for (int const cliqueVertex : vCliqueVertices) {
        vMarkedVertices[cliqueVertex] = false;
    }
}

void IsolateReducer::RemoveVertex(int const vertex)
{
    isolates.RemoveVertex(vertex);
}

void IsolateReducer::RemoveVertexAndNeighbors(int const vertex, std::vector<int> &vOtherRemovedVertices)
{
    isolates.RemoveVertexAndNeighbors(vertex, vOtherRemovedVertices);
}

void IsolateReducer::ReplaceVertices(vector<int> const &vVertices)
{
    isolates.ReplaceAllRemoved(vVertices);
}

bool IsolateReducer::InRemainingGraph(int const vertex) const
{
    return isolates.GetInGraph().Contains(vertex);
}

IsolateDominationReducer::IsolateDominationReducer(vector<vector<int>> const &adjacencyArray)
: IsolateReducer(adjacencyArray)
, remaining(adjacencyArray.size())
{
}

IsolateDominationReducer::~IsolateDominationReducer()
{
}

void IsolateDominationReducer::Reduce(vector<int> &vAlreadyConsideredVertices, vector<int> &vCliqueVertices, vector<int> &vOtherRemovedVertices)
{
    IsolateReducer::Reduce(vAlreadyConsideredVertices, vCliqueVertices, vOtherRemovedVertices);

#if 0
    return;
#endif
    remaining = isolates.GetInGraph();
    if (vAlreadyConsideredVertices.empty()) return;

    // TODO/DS: remove dominated vertices...
    // for independent sets, $u$ dominates $v$ if $u$ is not a neighbor of $v$ and
    // $u$ has at least the same set of non-neighbors that $v$ does in the remaining graph.
    // thus if $u$ and $v$ are neighbors or $u$ has a neighbor that $v$ does not, then $u$ does
    // not dominate $v$.
    // check for dominatation. p in P can't have neighbors that x in X does not
    for (int const vertex : isolates.GetInGraph()) {
        remaining.Insert(vertex);
    }

////    cout << "Removing dominated vertices: ";

    while (!remaining.Empty()) {
        int const p = *(remaining.begin());
        remaining.Remove(p);
        if (!isolates.GetInGraph().Contains(p)) continue;
////        vMarkedVertices[p] = ;
        for (int const neighborP : isolates.Neighbors()[p]) {
            vMarkedVertices[neighborP] = true;
        }

        for (int const x : vAlreadyConsideredVertices) {
            bool dominates(true);
            for (int const neighborX : m_AdjacencyArray[x]) {
                if (InRemainingGraph(neighborX) && !vMarkedVertices[neighborX]) {
                    dominates = false;
                    break;
                }
            }

            for (int const neighborP : isolates.Neighbors()[p]) {
                if (dominates) remaining.Insert(neighborP); // need to re-evaluate neighbors in case they are now dominated.
                vMarkedVertices[neighborP] = false;
            }

            if (dominates) {
////                cout << p << "(dominated by " << x << ") ";
                RemoveVertex(p);
////                IsolateReducer::Reduce(vAlreadyConsideredVertices, vCliqueVertices, vOtherRemovedVertices);
                vOtherRemovedVertices.push_back(p);
                break;
            } 
        }
    }

////    cout << endl;


////    if (newGraphSize != isolates.GetInGraph().Size()) {
////        cout << "Reduced " << (newGraphSize - isolates.GetInGraph().Size()) << "/" << newGraphSize << " vertices through dominated removal" << endl;
////    }
}
