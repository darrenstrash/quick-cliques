#ifndef MAXIMAL_CLIQUES_ALGORITHM_H
#define MAXIMAL_CLIQUES_ALGORITHM_H

// system includes
#include <list>
#include <string>

class MaximalCliquesAlgorithm
{
public:
    MaximalCliquesAlgorithm(std::string const &name);
    virtual ~MaximalCliquesAlgorithm();

    virtual long Run(std::list<std::list<int>> &cliques) = 0;

    std::string GetName() const;

private:
    std::string m_sName;
};

#endif //MAXIMAL_CLIQUES_ALGORITHM_H
