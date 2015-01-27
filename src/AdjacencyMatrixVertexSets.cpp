#include "AdjacencyMatrixVertexSets.h"

#include <vector>
#include <iostream>
#include <cstring> // memcpy

using namespace std;

AdjacencyMatrixVertexSets::AdjacencyMatrixVertexSets(vector<vector<char>> const &adjacencyMatrix)
: VertexSets("adjmatrix")
, beginX(0)
, beginP(0)
, beginR(adjacencyMatrix.size())
, m_AdjacencyMatrix(adjacencyMatrix)
, vertexSets()
, vertexLookup()
, degree()
, m_lDelineators()
{
}

AdjacencyMatrixVertexSets::~AdjacencyMatrixVertexSets()
{
}

void AdjacencyMatrixVertexSets::Initialize()
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
}

void AdjacencyMatrixVertexSets::PrintSummary(int const line) const
{
    cout << line << ": X[size=" << beginP-beginX << "], P[size=" << beginR-beginP << "]" << endl; 
}

bool AdjacencyMatrixVertexSets::GetNextTopLevelPartition()
{
    beginX = 0;
    beginP = 0;
    beginR = m_AdjacencyMatrix.size();

    bool const returnValue(!m_bDoneWithTopLevelPartitions);
    m_bDoneWithTopLevelPartitions = true;

    return returnValue;
}
