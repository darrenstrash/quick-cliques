#include "IndependentSetsReduction.h"

#include <vector>
#include <iostream>
#include <cstring> // memcpy

using namespace std;

IndependentSetsReduction::IndependentSetsReduction(vector<vector<int>> &adjacencyList)
: VertexSets("adjlist")
, beginX(0)
, beginP(0)
, beginR(adjacencyList.size())
, m_AdjacencyList(adjacencyList)
, vertexSets()
, vertexLookup()
, degree()
, m_lDelineators()
, isolates(adjacencyList)
, vvCliqueVertices()
, vvOtherRemovedVertices()
, vvAddedEdges()
{
}

IndependentSetsReduction::~IndependentSetsReduction()
{
}

void IndependentSetsReduction::Initialize()
{
    m_lDelineators.reserve(10);

    vertexSets  .resize(m_AdjacencyList.size(), 0);
    vertexLookup.resize(m_AdjacencyList.size(), 0);
    degree      .resize(m_AdjacencyList.size(), 0);

    for (size_t i = 0; i < m_AdjacencyList.size(); ++i) {
        degree[i] = m_AdjacencyList[i].size();
    }

    // indices indicating where each set P, X, R starts in 
    // vertexSets; initially, P contains all the vertices
    for (size_t i = 0; i < m_AdjacencyList.size(); ++i) {
        vertexSets[i] = i;
        vertexLookup[i] = i;
    }
}

void IndependentSetsReduction::PrintSummary(int const line) const
{
    cout << line << ": X[" << beginX << ":" << beginP << "), P[" << beginP << ":" << beginR << ")" << endl; 
}

bool IndependentSetsReduction::GetNextTopLevelPartition()
{
    if (m_bDoneWithTopLevelPartitions) return false;
    beginX = 0;
    beginP = 0;
    beginR = m_AdjacencyList.size();

    std::vector<int> vCliqueVertices;
    std::vector<int> vOtherRemoved;
    std::vector<std::pair<int,int>> vAddedEdges;
    isolates.RemoveAllIsolates(vCliqueVertices, vOtherRemoved, vAddedEdges);
    for (int const cliqueVertex : vCliqueVertices) {
        int const vertexLocation = vertexLookup[cliqueVertex];
        beginR--;
        vertexSets[vertexLocation] = vertexSets[beginR];
        vertexLookup[vertexSets[beginR]] = vertexLocation;
        vertexSets[beginR] = cliqueVertex;
        vertexLookup[cliqueVertex] = beginR;
    }

    for (int const otherVertex : vOtherRemoved) {
        int const vertexLocation = vertexLookup[otherVertex];
        beginR--;
        vertexSets[vertexLocation] = vertexSets[beginR];
        vertexLookup[vertexSets[beginR]] = vertexLocation;
        vertexSets[beginR] = otherVertex;
        vertexLookup[otherVertex] = beginR;
    }

    bool const returnValue(!m_bDoneWithTopLevelPartitions);
    m_bDoneWithTopLevelPartitions = true;

    return returnValue;
}

void IndependentSetsReduction::StoreGraphChanges(vector<int> &vCliqueVertices, vector<int> &vOtherRemoved, vector<pair<int,int>> &vAddedEdges)
{
    vvCliqueVertices      .emplace_back(std::move(vCliqueVertices));
    vvOtherRemovedVertices.emplace_back(std::move(vOtherRemoved));
    vvAddedEdges          .emplace_back(std::move(vAddedEdges));
}

void IndependentSetsReduction::RetrieveGraphChanges(vector<int> &vCliqueVertices, vector<int> &vOtherRemoved, vector<pair<int,int>> &vAddedEdges)
{
    vCliqueVertices = std::move(vvCliqueVertices.back());       vvCliqueVertices.pop_back();
    vOtherRemoved   = std::move(vvOtherRemovedVertices.back()); vvOtherRemovedVertices.pop_back();
    vAddedEdges     = std::move(vvAddedEdges.back());           vvAddedEdges.pop_back();
}

void IndependentSetsReduction::AddToAdjacencyList(vector<pair<int,int>> const &vEdges)
{
    for (pair<int,int> const &edge : vEdges) {
        m_AdjacencyList[edge.first].push_back(edge.second);
        m_AdjacencyList[edge.second].push_back(edge.first);
    }
}

void IndependentSetsReduction::RemoveFromAdjacencyList(vector<pair<int,int>> const &vEdges)
{
    for (pair<int,int> const &edge : vEdges) {
        vector<int>::iterator it = find(m_AdjacencyList[edge.first].begin(), m_AdjacencyList[edge.first].end(), edge.second);
        if (it == m_AdjacencyList[edge.first].end()) {
            cout << __LINE__ << ": found unexpected end of " << edge.first << "'s neighbor list" << endl;
        }
        m_AdjacencyList[edge.first].erase(it);
        it = find(m_AdjacencyList[edge.second].begin(), m_AdjacencyList[edge.second].end(), edge.first);
        if (it == m_AdjacencyList[edge.second].end()) {
            cout << __LINE__ << ": found unexpected end of " << edge.second << "'s neighbor list" << endl;
        }
        m_AdjacencyList[edge.second].erase(it);
    }
}
