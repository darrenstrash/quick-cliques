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
    virtual void Run() {}

    std::string GetName() const;

    void SetQuiet(bool const quiet);
    bool GetQuiet() const;

private:
    std::string m_sName;
    bool m_bQuiet;
};

#endif //ALGORITHM_H
