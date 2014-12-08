#include "Algorithm.h"

// system includes
#include <string>

using namespace std;

Algorithm::Algorithm(std::string const &name)
 : m_sName(name)
 , m_bQuiet(false)
{
}

Algorithm::~Algorithm()
{
}

string Algorithm::GetName() const
{
    return m_sName;
}

void Algorithm::SetQuiet(bool const quiet)
{
    m_bQuiet = quiet;
}

bool Algorithm::GetQuiet() const
{
    return m_bQuiet;
}
