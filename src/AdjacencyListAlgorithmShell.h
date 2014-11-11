#ifndef _DJS_ADJLIST_ALGORITHM_SHELL_H_
#define _DJS_ADJLIST_ALGORITHM_SHELL_H_
/* 
    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by 
    the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version. 
 
    This program is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
    GNU General Public License for more details. 
 
    You should have received a copy of the GNU General Public License 
    along with this program.  If not, see <http://www.gnu.org/licenses/> 
*/

/*! \file AdjacencyListAlgorithmShell.h

    \brief see AdjacencyListAlgorithmShell.cpp

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

// local includes
#include "MaximalCliquesAlgorithm.h"
#include "AdjacencyListVertexSets.h"
#include "Tools.h"

// system includes
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <list>
#include <vector>

class AdjacencyListAlgorithmShell : public MaximalCliquesAlgorithm
{
public:
    AdjacencyListAlgorithmShell(AdjacencyListVertexSets &sets);
    virtual ~AdjacencyListAlgorithmShell();

    virtual long Run(std::list<std::list<int>> &cliques);

    virtual void RunRecursive(long &cliqueCount, std::list<std::list<int>> &cliques, std::list<int> &partialClique);

    AdjacencyListAlgorithmShell           (AdjacencyListAlgorithmShell const &) = delete;
    AdjacencyListAlgorithmShell& operator=(AdjacencyListAlgorithmShell const &) = delete;

private:
    AdjacencyListVertexSets &m_Sets;
};

#endif
