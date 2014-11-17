#include "Algorithm.h"

// system includes
#include <string>

using namespace std;

Algorithm::Algorithm(std::string const &name)
 : m_sName(name)
{
}

Algorithm::~Algorithm()
{
}

string Algorithm::GetName() const
{
    return m_sName;
}
