#include "MaximalCliquesAlgorithm.h"

// system includes
#include <string>

using namespace std;

MaximalCliquesAlgorithm::MaximalCliquesAlgorithm(std::string const &name)
 : m_sName(name)
{
}

MaximalCliquesAlgorithm::~MaximalCliquesAlgorithm()
{
}

string MaximalCliquesAlgorithm::GetName() const
{
    return m_sName;
}
