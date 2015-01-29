#include "AdjacencyMatrixVertexSetsMax.h"
#include "CliqueColoringStrategy.h"
#include "DegeneracyTools.h"
#include "GraphTools.h"

#include <vector>
#include <iostream>
#include <cstring> // memcpy

using namespace std;

AdjacencyMatrixVertexSetsMax::AdjacencyMatrixVertexSetsMax(vector<vector<char>> const &adjacencyMatrix)
: VertexSets("adjmatrix-max")
, beginX(0)
, beginP(0)
, beginR(adjacencyMatrix.size())
, m_AdjacencyMatrix(adjacencyMatrix)
, vertexSets()
, vertexLookup()
, degree()
, m_lDelineators()
, stackP()
, stackColor()
, coloringStrategy(adjacencyMatrix)
{
}

AdjacencyMatrixVertexSetsMax::~AdjacencyMatrixVertexSetsMax()
{
}

void AdjacencyMatrixVertexSetsMax::Initialize()
{
    m_lDelineators.reserve(10);

    vertexSets  .resize(m_AdjacencyMatrix.size(), 0);
    vertexLookup.resize(m_AdjacencyMatrix.size(), 0);
    degree      .resize(m_AdjacencyMatrix.size(), 0);

    for (size_t i = 0; i < m_AdjacencyMatrix.size(); ++i) {
        for (size_t j = 0; i < m_AdjacencyMatrix[i].size(); ++i) {
            degree[i]++;
        }
    }

    // indices indicating where each set P, X, R starts in 
    // vertexSets; initially, P contains all the vertices
    for (size_t i = 0; i < m_AdjacencyMatrix.size(); ++i) {
        vertexSets[i] = i;
        vertexLookup[i] = i;
    }

    // initial ordering of vertices is reverse degeneracy ordering
////    stackP.emplace_back(std::move(GetVerticesInDegeneracyOrder(m_AdjacencyList)));
////    std::reverse(stackP.back().begin(), stackP.back().end());
    stackP.emplace_back(std::move(GraphTools::OrderVerticesByDegree(m_AdjacencyMatrix, degree, false /* descending*/)));
    stackColor.emplace_back(vector<int>(m_AdjacencyMatrix.size(), 0));

    coloringStrategy.Color(m_AdjacencyMatrix, stackP.back() /* evaluation order */, stackP.back() /* color order */, stackColor.back());
}

void AdjacencyMatrixVertexSetsMax::PrintSummary(int const line) const
{
    cout << line << ": X[size=" << beginP-beginX << "], P[size=" << beginR-beginP << "]" << endl; 
}

bool AdjacencyMatrixVertexSetsMax::GetNextTopLevelPartition()
{
    beginX = 0;
    beginP = 0;
    beginR = m_AdjacencyMatrix.size();

    bool const returnValue(!m_bDoneWithTopLevelPartitions);
    m_bDoneWithTopLevelPartitions = true;

    return returnValue;
}
