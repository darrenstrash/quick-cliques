#ifndef ALGORITHM_H
#define ALGORITHM_H

// system includes
#include <vector>
#include <list>
#include <string>
#include <functional>

class Algorithm
{
public:
    Algorithm(std::string const &name);
    virtual ~Algorithm();

    virtual long Run(std::list<std::list<int>> &cliques) = 0;
    virtual void Run() {}

    void SetName(std::string const &name);
    std::string GetName() const;

    void AddCallBack(std::function<void(std::list<int> const&)> callback);

    void ExecuteCallBacks(std::list<int> const &vertexSet) const;

    void SetQuiet(bool const quiet);
    bool GetQuiet() const;

private:
    std::string m_sName;
    bool m_bQuiet;
    std::vector<std::function<void(std::list<int> const&)>> m_vCallBacks;
};

#endif //ALGORITHM_H
