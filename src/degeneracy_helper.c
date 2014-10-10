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

#include<limits.h>
#include<assert.h>
#include<stdio.h>
#include<stdlib.h>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include"degeneracy_helper.h"

/*! \file degeneracy_helper.c

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

int computeDegeneracy(LinkedList** list, int size)
{
    int i = 0;

    int degeneracy = 0;
    
    // array of lists of vertices, indexed by degree
    LinkedList** verticesByDegree = (LinkedList**) Calloc(size, sizeof(LinkedList*));

    // array of lists of vertices, indexed by degree
    Link** vertexLocator = (Link**) Calloc(size, sizeof(Link*));

    int* degree = (int*) Calloc(size, sizeof(int));

    for(i = 0; i < size; i++)
    {
        verticesByDegree[i] = createLinkedList();
    }

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    for(i=0; i<size; i++)
    {
        degree[i] = length(list[i]);
        vertexLocator[i] = addFirst(verticesByDegree[degree[i]], ((void*) i));
    }
    
    int currentDegree = 0;

    int numVerticesRemoved = 0;

    while(numVerticesRemoved < size)
    {
        if(!isEmpty(verticesByDegree[currentDegree]))
        {
            degeneracy = max(degeneracy,currentDegree);
            
            int vertex = (int) getFirst(verticesByDegree[currentDegree]);

            delete(vertexLocator[vertex]);

            degree[vertex] = -1;

            LinkedList* neighborList = list[vertex];

            Link* neighborLink = neighborList->head->next;

            while(!isTail(neighborLink))
            {
                int neighbor = (int) neighborLink->data;

                if(degree[neighbor]!=-1)
                {
                    delete(vertexLocator[neighbor]);

                    degree[neighbor]--;

                    if(degree[neighbor] != -1)
                    {
                        vertexLocator[neighbor] = 
                            addFirst(verticesByDegree[degree[neighbor]], 
                                     (void*) neighbor);
                    }
                }

                neighborLink = neighborLink->next;
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else
        {
            currentDegree++;
        }

    }

    for(i = 0; i<size;i++)
    {
        destroyLinkedList(verticesByDegree[i]);
    }

    Free(vertexLocator);
    Free(verticesByDegree);
    Free(degree);

    return degeneracy;
}

/*! \brief

    \param list an input graph, represented as an array of linked lists of integers

    \param size the number of vertices in the graph

    \return an array of NeighborLists representing a degeneracy ordering of the vertices.

    \see NeighborList
*/

NeighborList** computeDegeneracyOrderList(LinkedList** list, int size)
{

#ifdef DEBUG
    printf("degeneracy is %d\n", computeDegeneracy(list, size));
#endif

    NeighborList** ordering = Calloc(size, sizeof(NeighborList*));

    int i = 0;

    int degeneracy = 0;
    
    // array of lists of vertices, indexed by degree
    LinkedList** verticesByDegree = (LinkedList**) Calloc(size, sizeof(LinkedList*));

    // array of lists of vertices, indexed by degree
    Link** vertexLocator = (Link**) Calloc(size, sizeof(Link*));

    int* degree = (int*) Calloc(size, sizeof(int));

    for(i = 0; i < size; i++)
    {
        verticesByDegree[i] = createLinkedList();
        ordering[i] = Malloc(sizeof(NeighborList));
        ordering[i]->earlier = createLinkedList();
        ordering[i]->later = createLinkedList();
    }

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    for(i=0; i<size; i++)
    {
        degree[i] = length(list[i]);
        //printf("degree[%d] = %d\n", i, degree[i]);
        vertexLocator[i] = addFirst(verticesByDegree[degree[i]], ((void*) i));
    }
    
    int currentDegree = 0;

    int numVerticesRemoved = 0;

    while(numVerticesRemoved < size)
    {
        if(!isEmpty(verticesByDegree[currentDegree]))
        {
            degeneracy = max(degeneracy,currentDegree);
            
            int vertex = (int) getFirst(verticesByDegree[currentDegree]);

            delete(vertexLocator[vertex]);

            ordering[vertex]->vertex = vertex;
            ordering[vertex]->orderNumber = numVerticesRemoved;

            degree[vertex] = -1;

            LinkedList* neighborList = list[vertex];

            Link* neighborLink = neighborList->head->next;

            while(!isTail(neighborLink))
            {
                int neighbor = (int) neighborLink->data;
                //printf("Neighbor: %d\n", neighbor);

                if(degree[neighbor]!=-1)
                {
                    delete(vertexLocator[neighbor]);
                    addLast(ordering[vertex]->later, (void*)neighbor);

                    degree[neighbor]--;

                    if(degree[neighbor] != -1)
                    {
                        vertexLocator[neighbor] = 
                            addFirst(verticesByDegree[degree[neighbor]], 
                                     (void*) neighbor);
                    }
                }
                else
                {
                    addLast(ordering[vertex]->earlier, (void*) neighbor);
                }

                neighborLink = neighborLink->next;
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else
        {
            currentDegree++;
        }

    }

    for(i = 0; i<size;i++)
    {
        destroyLinkedList(verticesByDegree[i]);
    }

    Free(vertexLocator);
    Free(verticesByDegree);
    Free(degree);

    return ordering;
}

/*! \brief

    \param list an input graph, represented as an array of linked lists of integers

    \param size the number of vertices in the graph

    \return an array of NeighborListArrays representing a degeneracy ordering of the vertices.

    \see NeighborListArray
*/

NeighborListArray** computeDegeneracyOrderArray(LinkedList** list, int size)
{

    NeighborList** ordering = Calloc(size, sizeof(NeighborList*));

    int i = 0;

    int degeneracy = 0;
    
    // array of lists of vertices, indexed by degree
    LinkedList** verticesByDegree = (LinkedList**) Calloc(size, sizeof(LinkedList*));

    // array of lists of vertices, indexed by degree
    Link** vertexLocator = (Link**) Calloc(size, sizeof(Link*));

    int* degree = (int*) Calloc(size, sizeof(int));

    for(i = 0; i < size; i++)
    {
        verticesByDegree[i] = createLinkedList();
        ordering[i] = Malloc(sizeof(NeighborList));
        ordering[i]->earlier = createLinkedList();
        ordering[i]->later = createLinkedList();
    }

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    for(i=0; i<size; i++)
    {
        degree[i] = length(list[i]);
        vertexLocator[i] = addFirst(verticesByDegree[degree[i]], ((void*) i));
    }
    
    int currentDegree = 0;

    int numVerticesRemoved = 0;

    while(numVerticesRemoved < size)
    {
        if(!isEmpty(verticesByDegree[currentDegree]))
        {
            degeneracy = max(degeneracy,currentDegree);
            
            int vertex = (int) getFirst(verticesByDegree[currentDegree]);

            delete(vertexLocator[vertex]);

            ordering[vertex]->vertex = vertex;
            ordering[vertex]->orderNumber = numVerticesRemoved;

            degree[vertex] = -1;

            LinkedList* neighborList = list[vertex];

            Link* neighborLink = neighborList->head->next;

            while(!isTail(neighborLink))
            {
                int neighbor = (int) neighborLink->data;

                if(degree[neighbor]!=-1)
                {
                    delete(vertexLocator[neighbor]);
                    addLast(ordering[vertex]->later, (void*)neighbor);

                    degree[neighbor]--;

                    if(degree[neighbor] != -1)
                    {
                        vertexLocator[neighbor] = 
                            addFirst(verticesByDegree[degree[neighbor]], 
                                     (void*) neighbor);
                    }
                }
                else
                {
                    addLast(ordering[vertex]->earlier, (void*) neighbor);
                }

                neighborLink = neighborLink->next;
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else
        {
            currentDegree++;
        }

    }

    NeighborListArray** orderingArray = Calloc(size, sizeof(NeighborListArray*));

    for(i = 0; i<size;i++)
    {
        orderingArray[i] = Malloc(sizeof(NeighborListArray));
        orderingArray[i]->vertex = ordering[i]->vertex;
        orderingArray[i]->orderNumber = ordering[i]->orderNumber;

        orderingArray[i]->laterDegree = length(ordering[i]->later);
        orderingArray[i]->later = Calloc(orderingArray[i]->laterDegree, sizeof(int));

        int j=0;
        Link* curr = ordering[i]->later->head->next;
        while(!isTail(curr))
        {
            orderingArray[i]->later[j++] = (int)(curr->data);
            curr = curr->next;
        }

        orderingArray[i]->earlierDegree = length(ordering[i]->earlier);
        orderingArray[i]->earlier = Calloc(orderingArray[i]->earlierDegree, sizeof(int));

        j=0;
        curr = ordering[i]->earlier->head->next;
        while(!isTail(curr))
        {
            orderingArray[i]->earlier[j++] = (int)(curr->data);
            curr = curr->next;
        }
    }

    for(i = 0; i<size;i++)
    {

        Free(ordering[i]);
        destroyLinkedList(verticesByDegree[i]);
    }

    Free(ordering);

    Free(vertexLocator);
    Free(verticesByDegree);
    Free(degree);

    return orderingArray;
}
