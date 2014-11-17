#ifndef ALGORITHM_H
#define ALGORITHM_H

// system includes
#include <list>
#include <string>

class Algorithm
{
public:
    Algorithm(std::string const &name);
    virtual ~Algorithm();

    virtual long Run(std::list<std::list<int>> &cliques) = 0;

    std::string GetName() const;

private:
    std::string m_sName;
};

#endif //ALGORITHM_H
