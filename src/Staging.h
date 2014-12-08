#ifndef QC_STAGING_H
#define QC_STAGING_H

// local includes
#include "Algorithm.h"

// system includes
#include <vector>
#include <list>

class Staging : public Algorithm
{
public:
    Staging(std::vector<std::vector<int>> &adjacencyList);
    ~Staging();

    virtual long Run(std::list<std::list<int>> &cliques) { Run(); return -1; }
    virtual void Run();

private:
    std::vector<std::vector<int>> &m_AdjacencyList;

};
#endif  //QC_STAGING_H
