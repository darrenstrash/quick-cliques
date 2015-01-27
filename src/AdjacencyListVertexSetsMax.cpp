#include "AdjacencyListVertexSetsMax.h"
#include "CliqueColoringStrategy.h"
#include "DegeneracyTools.h"
#include "GraphTools.h"

#include <vector>
#include <iostream>
#include <cstring> // memcpy

using namespace std;

AdjacencyListVertexSetsMax::AdjacencyListVertexSetsMax(vector<vector<int>> &adjacencyList)
: VertexSets("adjlist")
, beginX(0)
, beginP(0)
, beginR(adjacencyList.size())
, m_AdjacencyList(adjacencyList)
, vertexSets()
, vertexLookup()
, degree()
, m_lDelineators()
, stackP()
, stackColor()
, coloringStrategy(adjacencyList)
{
}

AdjacencyListVertexSetsMax::~AdjacencyListVertexSetsMax()
{
}

void AdjacencyListVertexSetsMax::Initialize()
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

    // initial ordering of vertices is reverse degeneracy ordering
////    stackP.emplace_back(std::move(GetVerticesInDegeneracyOrder(m_AdjacencyList)));
////    std::reverse(stackP.back().begin(), stackP.back().end());
    stackP.emplace_back(std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyList, false /* descending*/)));
    stackColor.emplace_back(vector<int>(m_AdjacencyList.size(), 0));

    coloringStrategy.Color(m_AdjacencyList, stackP.back(), stackColor.back());
}

void AdjacencyListVertexSetsMax::PrintSummary(int const line) const
{
    cout << line << ": X[size=" << beginP-beginX << "], P[size=" << beginR-beginP << "]" << endl; 
}

bool AdjacencyListVertexSetsMax::GetNextTopLevelPartition()
{
    beginX = 0;
    beginP = 0;
    beginR = m_AdjacencyList.size();

    bool const returnValue(!m_bDoneWithTopLevelPartitions);
    m_bDoneWithTopLevelPartitions = true;

    return returnValue;
}
