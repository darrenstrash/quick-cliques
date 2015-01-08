#include "ExperimentalReduction.h"

#include <vector>
#include <iostream>
#include <cstring> // memcpy

using namespace std;

ExperimentalReduction::ExperimentalReduction(vector<vector<int>> &adjacencyList)
: VertexSets("adjlist")
, m_AdjacencyList(adjacencyList)
, vertexSets()
, vertexLookup()
, degree()
, m_lDelineators()
, isolates(adjacencyList)
, vvCliqueVertices(1)
, vvOtherRemovedVertices(1)
, vvAddedEdges(1)
, m_Sets(m_AdjacencyList.size())
, m_vCliqueVertices(vvCliqueVertices.back())
{
}

ExperimentalReduction::~ExperimentalReduction()
{
}

void ExperimentalReduction::Initialize()
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

void ExperimentalReduction::PrintSummary(int const line) const
{
    m_Sets.PrintSummary();
}

bool ExperimentalReduction::GetNextTopLevelPartition()
{
    if (m_bDoneWithTopLevelPartitions) return false;

    for (int vertex = 0; vertex < m_AdjacencyList.size(); ++vertex) {
        m_Sets.InsertIntoP(vertex);
    }

    std::vector<int> vCliqueVertices;
    std::vector<int> vOtherRemoved;
    std::vector<std::pair<int,int>> vAddedEdges;
    isolates.RemoveAllIsolates(0, vCliqueVertices, vOtherRemoved, vAddedEdges);
    for (int const cliqueVertex : vCliqueVertices) {
        m_Sets.RemoveFromP(cliqueVertex);
    }

    for (int const otherVertex : vOtherRemoved) {
        m_Sets.RemoveFromP(otherVertex);
    }

    bool const returnValue(!m_bDoneWithTopLevelPartitions);
    m_bDoneWithTopLevelPartitions = true;

    return returnValue;
}

void ExperimentalReduction::StoreGraphChanges()
{
    m_vCliqueVertices = vvCliqueVertices.back();

    vvCliqueVertices      .push_back(vector<int>());
    vvOtherRemovedVertices.push_back(vector<int>());
    vvAddedEdges          .push_back(vector<pair<int,int>>());
}

void ExperimentalReduction::RetrieveGraphChanges()
{
    vvCliqueVertices.pop_back();
    vvOtherRemovedVertices.pop_back();
    vvAddedEdges.pop_back();

    m_vCliqueVertices = vvCliqueVertices.back();
}

void ExperimentalReduction::AddToAdjacencyList(vector<pair<int,int>> const &vEdges)
{
    for (pair<int,int> const &edge : vEdges) {
        m_AdjacencyList[edge.first].push_back(edge.second);
        m_AdjacencyList[edge.second].push_back(edge.first);
    }
}

void ExperimentalReduction::RemoveFromAdjacencyList(vector<pair<int,int>> const &vEdges)
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
