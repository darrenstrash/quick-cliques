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

#include <climits>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "Tools.h"
#include <list>
#include <vector>
#include "MemoryManager.h"
#include "DegeneracyTools.h"

using namespace std;

/*! \file DegeneracyTools.cpp

    \brief A collection of functions to compute degeneracy and degeneracy order.

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

/*! \brief

    \param list an input graph, represented as an array of linked lists of integers

    \param size the number of vertices in the graph

    \return the degeneracy of the input graph.
*/

int computeDegeneracy(vector<list<int>> const &adjList, int size)
{
    int i = 0;

    int degeneracy = 0;
    
    // array of lists of vertices, indexed by degree
    vector<list<int>> verticesByDegree(size);

    // array of lists of vertices, indexed by degree
    vector<std::list<int>::iterator> vertexLocator(size);

    vector<int> degree(size);

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    for(i=0; i<size; i++)
    {
        degree[i] = adjList[i].size();
        verticesByDegree[degree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[degree[i]].begin();
    }
    
    int currentDegree = 0;

    int numVerticesRemoved = 0;

    while(numVerticesRemoved < size)
    {
        if(!verticesByDegree[currentDegree].empty())
        {
            degeneracy = max(degeneracy,currentDegree);
            
            int const vertex = verticesByDegree[currentDegree].front();

            verticesByDegree[currentDegree].erase(vertexLocator[vertex]);

            degree[vertex] = -1;

            list<int> const &neighborList = adjList[vertex];

            for(int const neighbor : neighborList)
            {
                if(degree[neighbor]!=-1)
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);

                    degree[neighbor]--;

                    if(degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else
        {
            currentDegree++;
        }
    }

    verticesByDegree.clear();

    return degeneracy;
}

/*! \brief

    \param list an input graph, represented as an array of linked lists of integers

    \param size the number of vertices in the graph

    \return an array of NeighborLists representing a degeneracy ordering of the vertices.

    \see NeighborList
*/

NeighborList** computeDegeneracyOrderList(vector<list<int>> const &adjList, int size)
{

#ifdef DEBUG
    printf("degeneracy is %d\n", computeDegeneracy(list, size));
#endif

    NeighborList** ordering = (NeighborList**)Calloc(size, sizeof(NeighborList*));

    int i = 0;

    int degeneracy = 0;
    
    // array of lists of vertices, indexed by degree
    vector<list<int>> verticesByDegree(size);

    // array of lists of vertices, indexed by degree
    vector<list<int>::iterator> vertexLocator(size);

    vector<int> degree(size);

    for(i = 0; i < size; i++)
    {
        ordering[i] = (NeighborList*)Malloc(sizeof(NeighborList));
    }

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    for(i=0; i<size; i++)
    {
        degree[i] = adjList[i].size();
        //printf("degree[%d] = %d\n", i, degree[i]);
        verticesByDegree[degree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[degree[i]].begin();
    }
    
    int currentDegree = 0;

    int numVerticesRemoved = 0;

    while(numVerticesRemoved < size)
    {
        if(!verticesByDegree[currentDegree].empty())
        {
            degeneracy = max(degeneracy,currentDegree);
            
            int const vertex = verticesByDegree[currentDegree].front();
            verticesByDegree[currentDegree].erase(vertexLocator[vertex]);

            ordering[vertex]->vertex = vertex;
            ordering[vertex]->orderNumber = numVerticesRemoved;

            degree[vertex] = -1;

            list<int> const &neighborList = adjList[vertex];

            for(int const neighbor : neighborList)
            {
                //printf("Neighbor: %d\n", neighbor);

                if(degree[neighbor]!=-1)
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    ordering[vertex]->later.push_back(neighbor);

                    degree[neighbor]--;

                    if(degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
                else
                {
                    ordering[vertex]->earlier.push_back(neighbor);
                }
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else
        {
            currentDegree++;
        }
    }

    return ordering;
}

/*! \brief

    \param list an input graph, represented as an array of linked lists of integers

    \param size the number of vertices in the graph

    \return an array of NeighborListArrays representing a degeneracy ordering of the vertices.

    \see NeighborListArray
*/

NeighborListArray** computeDegeneracyOrderArray(vector<list<int>> const &adjList, int size)
{

    vector<NeighborList> vOrdering(size);

    int i = 0;

    int degeneracy = 0;
    
    // array of lists of vertices, indexed by degree
    vector<list<int>> verticesByDegree(size);

    // array of lists of vertices, indexed by degree
    vector<list<int>::iterator> vertexLocator(size);

    vector<int> degree(size);

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    for(i=0; i<size; i++)
    {
        degree[i] = adjList[i].size();
        verticesByDegree[degree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[degree[i]].begin();
    }

    int currentDegree = 0;

    int numVerticesRemoved = 0;

    while(numVerticesRemoved < size)
    {
        if(!verticesByDegree[currentDegree].empty())
        {
            degeneracy = max(degeneracy,currentDegree);
            
            int const vertex = verticesByDegree[currentDegree].front();
            verticesByDegree[currentDegree].pop_front();

            vOrdering[vertex].vertex = vertex;
            vOrdering[vertex].orderNumber = numVerticesRemoved;

            degree[vertex] = -1;

            list<int> const &neighborList = adjList[vertex];

            for(int const neighbor : neighborList)
            {
                if(degree[neighbor]!=-1)
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    (vOrdering[vertex].later).push_back(neighbor);

                    degree[neighbor]--;

                    if(degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
                else
                {
                    vOrdering[vertex].earlier.push_back(neighbor);
                }
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else
        {
            currentDegree++;
        }

    }

    NeighborListArray** orderingArray = (NeighborListArray**)Calloc(size, sizeof(NeighborListArray*));

    for(i = 0; i<size;i++)
    {
        orderingArray[i] = (NeighborListArray*)Malloc(sizeof(NeighborListArray));
        orderingArray[i]->vertex = vOrdering[i].vertex;
        orderingArray[i]->orderNumber = vOrdering[i].orderNumber;

        orderingArray[i]->laterDegree = vOrdering[i].later.size();
        orderingArray[i]->later.resize(orderingArray[i]->laterDegree);

        int j=0;
        for(int const laterNeighbor : vOrdering[i].later)
        {
            orderingArray[i]->later[j++] = laterNeighbor;
        }

        orderingArray[i]->earlierDegree = vOrdering[i].earlier.size();
        orderingArray[i]->earlier.resize(orderingArray[i]->earlierDegree);

        j=0;
        for (int const earlierNeighbor : vOrdering[i].earlier)
        {
            orderingArray[i]->earlier[j++] = earlierNeighbor;
        }
    }

    return orderingArray;
}

vector<NeighborListArray> computeDegeneracyOrderArray(vector<vector<int>> &adjArray, int size)
{
    int i = 0;

    // array of lists of vertices, indexed by degree
    vector<list<int>> verticesByDegree(size);

    // array of lists of vertices, indexed by degree
    vector<list<int>::iterator> vertexLocator(size);

    vector<int> degree(size);

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    for(i=0; i<size; i++)
    {
        degree[i] = adjArray[i].size();
        verticesByDegree[degree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[degree[i]].begin();
    }

    int currentDegree = 0;

    int numVerticesRemoved = 0;

    vector<NeighborListArray> vOrderingArray(size);

    while (numVerticesRemoved < size) {
        if (!verticesByDegree[currentDegree].empty()) {
            int const vertex = verticesByDegree[currentDegree].front();
            verticesByDegree[currentDegree].pop_front();

            vOrderingArray[vertex].vertex = vertex;
            vOrderingArray[vertex].orderNumber = numVerticesRemoved;

            degree[vertex] = -1;

            // will swap later neighbors to end of neighbor list
            vector<int> &neighborList = adjArray[vertex];

            int splitPoint(neighborList.size());
            for(int i=0; i < splitPoint; ++i) {
                int const neighbor(neighborList[i]);
                // if later neighbor, swap to end of neighborList (there are few of these)
                if(degree[neighbor]!=-1) {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);

                    neighborList[i] = neighborList[--splitPoint];
                    neighborList[splitPoint] = neighbor;
                    i--;

                    degree[neighbor]--;

                    if (degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
                // earlier neighbor, do nothing.
            }

            // create space for later neighbors to ordering
            vOrderingArray[vertex].laterDegree = neighborList.size() - splitPoint;
            vOrderingArray[vertex].later.resize(neighborList.size() - splitPoint);

            // create space for earlier neighbors to ordering
            vOrderingArray[vertex].earlierDegree = splitPoint;
            vOrderingArray[vertex].earlier.resize(splitPoint);

            // fill in earlier and later neighbors
            for (int i = 0; i < splitPoint; ++i) {
////            cout << "earlier: " << vOrderingArray[vertex].earlier.size() << endl;
////            cout << "split  : " << splitPoint << endl;
////            cout << "earlier[" << i << "]" << endl;
                vOrderingArray[vertex].earlier[i] = neighborList[i];
            }

            for (int i = splitPoint; i < neighborList.size(); ++i) {
////            cout << "later  : " << vOrderingArray[vertex].later.size() << endl;
////            cout << "split  : " << splitPoint << endl;
////            cout << "later [" << i - splitPoint << "]" << endl;
                vOrderingArray[vertex].later[i-splitPoint] = neighborList[i];
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else {
            currentDegree++;
        }
    }

    return vOrderingArray;
}
