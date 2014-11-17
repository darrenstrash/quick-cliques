#include "IndependentSets.h"

#include <vector>
#include <iostream>
#include <cstring> // memcpy

using namespace std;

IndependentSets::IndependentSets(vector<vector<int>> const &adjacencyList)
: VertexSets("adjlist")
, beginX(0)
, beginP(0)
, beginR(adjacencyList.size())
, m_AdjacencyList(adjacencyList)
, vertexSets()
, vertexLookup()
, degree()
, m_lDelineators()
{
}

IndependentSets::~IndependentSets()
{
}

void IndependentSets::Initialize()
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

void IndependentSets::PrintSummary(int const line) const
{
    cout << line << ": X[size=" << beginP-beginX << "], P[size=" << beginR-beginP << "]" << endl; 
}

bool IndependentSets::GetNextTopLevelPartition()
{
    beginX = 0;
    beginP = 0;
    beginR = m_AdjacencyList.size();

    bool const returnValue(!m_bDoneWithTopLevelPartitions);
    m_bDoneWithTopLevelPartitions = true;

    return returnValue;
}
